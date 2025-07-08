import multiprocessing as mp
from multiprocessing import Queue, Process, Manager
import random
from collections import defaultdict
import time
import numpy as np
from itertools import islice
from queue import Empty  # Ajout de l'import correct pour Queue.Empty

class ParallelBLAST:
    def __init__(self, word_size=3, match_score=2, mismatch_score=-1, gap_score=-2, threshold=5):
        self.word_size = word_size
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        self.threshold = threshold
        self.nucleotides = ['A', 'T', 'C', 'G']
        
        total_cores = mp.cpu_count()
        # Optimisation finale du nombre de processus
        if total_cores <= 4:
            self.n_indexers = 1
            self.n_seed_finders = 1
            self.n_extenders = max(2, total_cores - 2)
        else:
            self.n_indexers = 2
            self.n_seed_finders = 2
            self.n_extenders = max(4, total_cores - 4)



    def parallel_index_worker(self, sequence_chunks, word_size, result_queue):
        """Worker optimisé pour la création parallèle de l'index"""
        local_index = {}  # Utiliser un dict au lieu de defaultdict pour la performance
        
        for chunk, chunk_offset in sequence_chunks:
            chunk_offset = int(chunk_offset)
            for i in range(len(chunk) - word_size + 1):
                word = chunk[i:i + word_size]
                if word not in local_index:
                    local_index[word] = []
                local_index[word].append(i + chunk_offset)
                    
        result_queue.put(local_index)



    def parallel_create_index(self, sequence, chunk_size=5000):
        """Création parallèle optimisée de l'index"""
        chunks = []
        offset = 0
        
        # Utiliser numpy pour un découpage plus efficace
        sequence_array = np.array(list(sequence))
        chunk_indices = np.arange(0, len(sequence), chunk_size)
        
        for start in chunk_indices:
            end = min(start + chunk_size, len(sequence))
            chunk = ''.join(sequence_array[start:end])
            # Convertir explicitement l'offset en int
            chunks.append((chunk, int(offset)))
            offset += len(chunk)
        
        # Utiliser un Manager pour partager l'index entre les processus
        manager = Manager()
        shared_index = manager.dict()
        processes = []
        
        chunk_sets = np.array_split(chunks, self.n_indexers)
        
        for chunk_set in chunk_sets:
            if len(chunk_set) > 0:
                p = Process(target=self._index_worker_optimized, 
                        args=(chunk_set, self.word_size, shared_index))
                processes.append(p)
                p.start()
        
        for p in processes:
            p.join()
        
        # Convertir l'index partagé en defaultdict
        final_index = defaultdict(list)
        for word, positions in shared_index.items():
            final_index[word].extend(positions)
        
        return final_index

    def _index_worker_optimized(self, chunk_set, word_size, shared_index):
        """Worker optimisé pour l'indexation"""
        local_index = {}
        
        for chunk, chunk_offset in chunk_set:
            # Assurer que chunk_offset est un int
            chunk_offset = int(chunk_offset)
            words = [chunk[i:i + word_size] for i in range(len(chunk) - word_size + 1)]
            for i, word in enumerate(words):
                if word not in local_index:
                    local_index[word] = []
                local_index[word].append(i + chunk_offset)
        
        # Mettre à jour l'index partagé en une seule fois
        for word, positions in local_index.items():
            if word not in shared_index:
                shared_index[word] = positions
            else:
                shared_index[word].extend(positions)


    def seed_finder_worker(self, query_chunks, subject_index, result_queue):
        """Worker pour la recherche parallèle des seeds"""
        local_seeds = []
        
        for chunk, chunk_offset in query_chunks:
            for i in range(len(chunk) - self.word_size + 1):
                word = chunk[i:i + self.word_size]
                if word in subject_index:
                    for j in subject_index[word]:
                        local_seeds.append((i + int(chunk_offset), j))  # Conversion en int
        
        result_queue.put(local_seeds)

    def parallel_find_seeds(self, query, subject_index, chunk_size=2000):
        """Recherche de seeds optimisée"""
        # Utiliser numpy pour une manipulation plus efficace des données
        query_array = np.array(list(query))
        n_chunks = max(1, len(query) // chunk_size)
        chunks = np.array_split(query_array, n_chunks)
        
        # Utiliser une Queue avec une capacité maximale pour éviter la surcharge mémoire
        result_queue = Queue(maxsize=self.n_seed_finders * 2)
        processes = []
        
        for i, chunk in enumerate(chunks):
            chunk_str = ''.join(chunk)
            p = Process(target=self._seed_finder_worker_optimized,
                    args=(chunk_str, i * len(chunk), subject_index, result_queue))
            processes.append(p)
            p.start()
        
        # Collecter les résultats de manière plus efficace
        seeds = []
        completed = 0
        while completed < len(processes):
            try:
                batch = result_queue.get(timeout=1)
                seeds.extend(batch)
                completed += 1
            except Empty:
                continue
        
        for p in processes:
            p.join()
        
        return seeds
    
    def _seed_finder_worker_optimized(self, chunk, offset, subject_index, result_queue):
        """Worker optimisé pour la recherche de seeds"""
        local_seeds = []
        word_size = self.word_size
        
        # Utiliser une compréhension de liste pour générer tous les mots possibles
        words = [chunk[i:i + word_size] for i in range(len(chunk) - word_size + 1)]
        
        # Traiter les mots par lots pour réduire les accès à la queue
        for i, word in enumerate(words):
            if word in subject_index:
                local_seeds.extend((i + offset, j) for j in subject_index[word])
        
        # Envoyer les résultats en une seule fois
        if local_seeds:
            result_queue.put(local_seeds)



    def extend_alignment_worker(self, query, subject, seeds_queue, results_queue, done_queue):
        """Worker optimisé pour l'extension des alignements"""
        while True:
            try:
                # Vérifier d'abord si le travail est terminé
                if not seeds_queue.empty():
                    try:
                        seeds = seeds_queue.get(timeout=0.1)
                    except Empty:
                        continue
                else:
                    # Vérifier si nous devons terminer
                    if not done_queue.empty():
                        break
                    time.sleep(0.1)
                    continue

                local_alignments = []
                for seed in seeds:
                    try:
                        alignment = self.extend_alignment(query, subject, seed)
                        if alignment:
                            local_alignments.append(alignment)
                    except Exception:
                        continue

                if local_alignments:
                    results_queue.put(local_alignments)

            except Exception:
                # Si une erreur se produit, vérifier si nous devons terminer
                if not done_queue.empty():
                    break
                continue

        # Assurer que nous ne perdons pas les derniers résultats
        if local_alignments:
            results_queue.put(local_alignments)



    def parallel_extend_alignments(self, query, subject, seeds, chunk_size=None):
        """Extension parallèle avec gestion améliorée des queues"""
        if chunk_size is None:
            chunk_size = max(100, len(seeds) // (self.n_extenders * 4))
        
        seeds_queue = Queue()
        results_queue = Queue()
        done_queue = Queue()
        
        # Distribuer les seeds dans la queue
        for i in range(0, len(seeds), chunk_size):
            chunk = seeds[i:i + chunk_size]
            if chunk:  # Vérifier que le chunk n'est pas vide
                seeds_queue.put(chunk)
        
        processes = []
        for _ in range(self.n_extenders):
            p = Process(target=self.extend_alignment_worker,
                    args=(query, subject, seeds_queue, results_queue, done_queue))
            p.start()
            processes.append(p)
        
        # Signaler la fin du travail aux workers
        for _ in range(self.n_extenders):
            done_queue.put(True)
        
        # Collecter les résultats
        alignments = []
        completed_workers = 0
        
        while completed_workers < self.n_extenders:
            try:
                result = results_queue.get(timeout=1.0)
                if result:
                    alignments.extend(result)
            except Empty:
                completed_workers += 1
        
        # Nettoyage
        for p in processes:
            p.join(timeout=1.0)
            if p.is_alive():
                p.terminate()
        
        return alignments


    def extend_alignment_worker_dynamic(self, query, subject, seeds_queue, results_queue, 
                                    done_queue, status_queue, worker_id):
        """Worker optimisé avec reporting de charge"""
        processed_count = 0
        batch_size = 10
        
        while True:
            try:
                try:
                    seeds = seeds_queue.get(timeout=0.1)
                except Empty:
                    if not done_queue.empty():
                        status_queue.put((worker_id, processed_count))
                        break
                    continue
                
                local_alignments = []
                for seed in seeds:
                    try:
                        alignment = self.extend_alignment(query, subject, seed)
                        if alignment:
                            local_alignments.append(alignment)
                        processed_count += 1
                        
                        if processed_count % batch_size == 0:
                            status_queue.put((worker_id, processed_count))
                            
                    except Exception as e:
                        print(f"Error processing seed in worker {worker_id}: {str(e)}")
                        continue
                
                if local_alignments:
                    results_queue.put(local_alignments)
                    
            except Exception as e:
                print(f"Error in worker {worker_id}: {str(e)}")
                continue


    def extend_alignment(self, query, subject, seed):
        """Étend un seed pour créer un alignement local"""
        query_pos, subject_pos = seed
        
        # Extension à gauche
        left_query_pos = query_pos
        left_subject_pos = subject_pos
        left_score = 0
        best_left_score = 0
        best_left_endpoints = (query_pos, subject_pos)
        
        while left_query_pos > 0 and left_subject_pos > 0:
            left_query_pos -= 1
            left_subject_pos -= 1
            if query[left_query_pos] == subject[left_subject_pos]:
                left_score += self.match_score
            else:
                left_score += self.mismatch_score
            
            if left_score > best_left_score:
                best_left_score = left_score
                best_left_endpoints = (left_query_pos, left_subject_pos)
            elif left_score < best_left_score - self.threshold:
                break
        
        # Extension à droite
        right_query_pos = query_pos + self.word_size
        right_subject_pos = subject_pos + self.word_size
        right_score = 0
        best_right_score = 0
        best_right_endpoints = (right_query_pos, right_subject_pos)
        
        while right_query_pos < len(query) and right_subject_pos < len(subject):
            if query[right_query_pos] == subject[right_subject_pos]:
                right_score += self.match_score
            else:
                right_score += self.mismatch_score
            
            if right_score > best_right_score:
                best_right_score = right_score
                best_right_endpoints = (right_query_pos + 1, right_subject_pos + 1)
            elif right_score < best_right_score - self.threshold:
                break
                
            right_query_pos += 1
            right_subject_pos += 1
        
        total_score = best_left_score + best_right_score
        
        # Ne retourner l'alignement que s'il dépasse le seuil
        if total_score >= self.threshold:
            return {
                'score': total_score,
                'query_start': best_left_endpoints[0],
                'query_end': best_right_endpoints[0],
                'subject_start': best_left_endpoints[1],
                'subject_end': best_right_endpoints[1],
                'query_seq': query[best_left_endpoints[0]:best_right_endpoints[0]],
                'subject_seq': subject[best_left_endpoints[1]:best_right_endpoints[1]]
            }
        return None
    
    def calculate_alignment_stats(self, query_seq, subject_seq):
        """Calcule les statistiques détaillées pour un alignement"""
        if not query_seq or not subject_seq:  # Protection contre les séquences vides
            return {
                'length': 0,
                'matches': 0,
                'identity_percent': 0,
                'mismatches': 0
            }
        
        matches = sum(q == s for q, s in zip(query_seq, subject_seq))
        length = len(query_seq)
        try:
            identity = (matches / length) * 100
        except ZeroDivisionError:
            identity = 0
        
        return {
            'length': length,
            'matches': matches,
            'identity_percent': round(identity, 2),
            'mismatches': length - matches
        }


    def format_alignment_visualization(self, query_seq, subject_seq):
        """Crée une visualisation de l'alignement avec des marqueurs de similarité"""
        match_line = ''.join('|' if q == s else ' ' for q, s in zip(query_seq, subject_seq))
        return f"Query:   {query_seq}\n         {match_line}\nSubject: {subject_seq}"

    def filter_alignments(self, alignments, min_identity=60, min_length=10):
        """Filtre les alignments selon des critères de qualité"""
        filtered_alignments = []
        for align in alignments:
            # Calculer les stats si elles n'existent pas déjà
            if 'identity_percent' not in align:
                stats = self.calculate_alignment_stats(align['query_seq'], 
                                                    align['subject_seq'])
                align.update(stats)
            
            if (align['identity_percent'] >= min_identity and 
                align['length'] >= min_length):
                filtered_alignments.append(align)
        
        return filtered_alignments

    def align(self, query, subject, min_identity=60, min_length=10):
        """Méthode principale avec sélection adaptative du mode"""
        # Sélection automatique du mode basé sur la taille des séquences
        if len(query) < 5000 and len(subject) < 10000:
            return self.sequential_align(query, subject)
        # Étape 1: Création parallèle de l'index
        start_time = time.time()
        subject_index = self.parallel_create_index(subject)
        index_time = time.time() - start_time
        
        # Étape 2: Recherche parallèle des seeds
        start_time = time.time()
        seeds = self.parallel_find_seeds(query, subject_index)
        seed_time = time.time() - start_time
        
        # Étape 3: Extension parallèle des alignements
        start_time = time.time()
        alignments = self.parallel_extend_alignments(query, subject, seeds)
        extend_time = time.time() - start_time
        
        # Éliminer les doublons basés sur les positions et séquences
        unique_alignments = []
        seen = set()
        for align in alignments:
            key = (align['query_start'], align['query_end'], 
                align['subject_start'], align['subject_end'],
                align['query_seq'], align['subject_seq'])
            if key not in seen:
                seen.add(key)
                unique_alignments.append(align)

        sorted_alignments = sorted(unique_alignments, 
                             key=lambda x: x['score'], 
                             reverse=True)

        
        # Tri final des résultats uniques
        sorted_alignments = self.filter_alignments(sorted_alignments, 
                                             min_identity=min_identity,
                                             min_length=min_length)
    
        
        for align in sorted_alignments:
            stats = self.calculate_alignment_stats(align['query_seq'], align['subject_seq'])
            align.update(stats)
        
        return {
            'alignments': sorted_alignments,
            'timings': {
                'indexing': index_time,
                'seed_finding': seed_time,
                'extension': extend_time,
                'total': index_time + seed_time + extend_time
            },
            'summary': {
                'total_alignments': len(sorted_alignments),
                'avg_identity': round(sum(a['identity_percent'] for a in sorted_alignments) / len(sorted_alignments), 2) if sorted_alignments else 0,
                'avg_length': round(sum(len(a['query_seq']) for a in sorted_alignments) / len(sorted_alignments), 2) if sorted_alignments else 0,
                'max_score': max(a['score'] for a in sorted_alignments) if sorted_alignments else 0
            }
        }
        
    def sequential_align(self, query, subject):
        """Version séquentielle de l'alignement pour comparaison"""
        start_time = time.time()
        
        # Index creation
        subject_index = defaultdict(list)
        for i in range(len(subject) - self.word_size + 1):
            word = subject[i:i + self.word_size]
            subject_index[word].append(i)
        index_time = time.time() - start_time
        
        # Seed finding
        start_time = time.time()
        seeds = []
        for i in range(len(query) - self.word_size + 1):
            word = query[i:i + self.word_size]
            if word in subject_index:
                for j in subject_index[word]:
                    seeds.append((i, j))
        seed_time = time.time() - start_time
        
        # Extension
        start_time = time.time()
        alignments = []
        for seed in seeds:
            alignment = self.extend_alignment(query, subject, seed)
            if alignment:
                alignments.append(alignment)
        extend_time = time.time() - start_time
        
        total_time = index_time + seed_time + extend_time
        
        return {
            'alignments': alignments,
            'timings': {
                'indexing': index_time,
                'seed_finding': seed_time,
                'extension': extend_time,
                'total': total_time
            }
        }

def main():
    # Initialisation avec des paramètres plus stricts
    blast = ParallelBLAST(
        word_size=4,           # Augmenté de 3 à 4 pour plus de spécificité
        match_score=2,
        mismatch_score=-2,     # Pénalité plus forte pour les mismatch
        gap_score=-3,          # Pénalité plus forte pour les gaps
        threshold=10           # Seuil plus élevé pour filtrer les alignements faibles
    )
    
    # Génération des séquences de test
    query_length = 10000     # 5x plus grand
    subject_length = 25000   # 5x plus grand

    query_seq = ''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(query_length))
    subject_seq = ''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(subject_length))
    
    print(f"\nAnalyse BLAST parallèle")
    print(f"Tailles - Requête: {query_length}, Sujet: {subject_length}")
    print(f"Nombre de processus - Indexation: {blast.n_indexers}, "
          f"Seeds: {blast.n_seed_finders}, Extension: {blast.n_extenders}")
    
    # Exécution séquentielle
    print("\nExécution séquentielle...")
    seq_results = blast.sequential_align(query_seq, subject_seq)
    seq_time = seq_results['timings']['total']
    
    # Exécution parallèle
    print("\nExécution parallèle...")
    results = blast.align(query_seq, subject_seq)
    parallel_time = results['timings']['total']
    
    # Calcul des métriques de performance
    total_processes = blast.n_indexers + blast.n_seed_finders + blast.n_extenders
    speedup = seq_time / parallel_time
    efficiency = speedup / total_processes
    
    # Affichage des performances
    print("\nMétriques de performance:")
    print(f"Temps séquentiel: {seq_time:.3f}s")
    print(f"Temps parallèle: {parallel_time:.3f}s")
    print(f"Nombre total de processus: {total_processes}")
    print(f"Accélération (speedup): {speedup:.2f}x")
    print(f"Efficacité: {(efficiency * 100):.2f}%")
    
    # Affichage des résultats détaillés pour chaque étape
    print("\nDétails par étape:")
    for stage in ['indexing', 'seed_finding', 'extension']:
        stage_speedup = seq_results['timings'][stage] / results['timings'][stage]
        print(f"{stage.capitalize()}:")
        print(f"  Temps séquentiel: {seq_results['timings'][stage]:.3f}s")
        print(f"  Temps parallèle: {results['timings'][stage]:.3f}s")
        print(f"  Accélération: {stage_speedup:.2f}x")
    
    # Affichage des performances originales
    print("\nPerformances parallèles détaillées:")
    print(f"Temps d'indexation: {results['timings']['indexing']:.3f}s")
    print(f"Temps de recherche des seeds: {results['timings']['seed_finding']:.3f}s")
    print(f"Temps d'extension: {results['timings']['extension']:.3f}s")
    print(f"Temps total: {results['timings']['total']:.3f}s")
    print(f"Nombre d'alignements trouvés: {len(results['alignments'])}")
    
    # Afficher quelques détails sur les meilleurs alignements
    if results['alignments']:
        print("\nRésumé global:")
        print(f"Nombre total d'alignements: {results['summary']['total_alignments']}")
        print(f"Identité moyenne: {results['summary']['avg_identity']}%")
        print(f"Longueur moyenne: {results['summary']['avg_length']} bp")
        print(f"Score maximum: {results['summary']['max_score']}")

    print("\nMeilleurs alignements:")
    for i, align in enumerate(results['alignments'][:3]):
        print(f"\nAlignment {i+1}:")
        print(f"Score: {align['score']}")
        print(f"Longueur: {align['length']} bp")
        print(f"Identité: {align['identity_percent']}%")
        print(f"Matches: {align['matches']}/{align['length']}")
        print("\nVisualisation:")
        print(blast.format_alignment_visualization(align['query_seq'], 
                                                align['subject_seq']))
        print(f"\nPositions:")
        print(f"Query:   {align['query_start']}-{align['query_end']}")
        print(f"Subject: {align['subject_start']}-{align['subject_end']}")

if __name__ == "__main__":
    main()

