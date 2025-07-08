import random
from collections import defaultdict
import time
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import psutil
import os

class BLAST:
    def __init__(self, word_size=3, match_score=2, mismatch_score=-1, gap_score=-2, threshold=5):
        self.word_size = word_size
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        self.threshold = threshold
        self.nucleotides = ['A', 'T', 'C', 'G']

    def generate_random_sequence(self, length):
        """Génère une séquence ADN aléatoire de la longueur spécifiée"""
        return ''.join(random.choice(self.nucleotides) for _ in range(length))

    def create_word_index(self, sequence):
        """Crée un index des mots de taille word_size dans la séquence"""
        index = defaultdict(list)
        for i in range(len(sequence) - self.word_size + 1):
            word = sequence[i:i + self.word_size]
            index[word].append(i)
        return index

    def find_seed_matches(self, query, subject_index):
        """Trouve les correspondances exactes (seeds) entre la requête et le sujet"""
        matches = []
        for i in range(len(query) - self.word_size + 1):
            word = query[i:i + self.word_size]
            if word in subject_index:
                for j in subject_index[word]:
                    matches.append((i, j))
        return matches

    def extend_alignment(self, query, subject, seed_position):
        """Étend l'alignement à partir d'une position seed"""
        query_start, subject_start = seed_position
        
        # Extension vers la gauche
        left_score = 0
        best_left_score = 0
        best_left_extension = 0
        i, j = query_start - 1, subject_start - 1
        
        while i >= 0 and j >= 0:
            if query[i] == subject[j]:
                left_score += self.match_score
            else:
                left_score += self.mismatch_score
            
            if left_score > best_left_score:
                best_left_score = left_score
                best_left_extension = query_start - i
            elif left_score < 0:
                break
                
            i -= 1
            j -= 1

        # Extension vers la droite
        right_score = 0
        best_right_score = 0
        best_right_extension = 0
        i = query_start + self.word_size
        j = subject_start + self.word_size
        
        while i < len(query) and j < len(subject):
            if query[i] == subject[j]:
                right_score += self.match_score
            else:
                right_score += self.mismatch_score
            
            if right_score > best_right_score:
                best_right_score = right_score
                best_right_extension = i - (query_start + self.word_size) + 1
            elif right_score < 0:
                break
                
            i += 1
            j += 1

        total_score = best_left_score + best_right_score + (self.word_size * self.match_score)
        
        if total_score >= self.threshold:
            q_start = query_start - best_left_extension
            q_end = query_start + self.word_size + best_right_extension
            s_start = subject_start - best_left_extension
            s_end = subject_start + self.word_size + best_right_extension
            
            return {
                'score': total_score,
                'query_start': q_start,
                'query_end': q_end,
                'subject_start': s_start,
                'subject_end': s_end,
                'query_align': query[q_start:q_end],
                'subject_align': subject[s_start:s_end]
            }
        return None

    def align(self, query, subject):
        """Effectue l'alignement BLAST entre deux séquences"""
        # Création de l'index des mots pour la séquence sujet
        subject_index = self.create_word_index(subject)
        
        # Recherche des seeds
        seed_matches = self.find_seed_matches(query, subject_index)
        
        # Extension des alignements
        alignments = []
        for seed in seed_matches:
            alignment = self.extend_alignment(query, subject, seed)
            if alignment:
                alignments.append(alignment)
        
        # Tri des alignements par score
        return sorted(alignments, key=lambda x: x['score'], reverse=True)

    def format_alignment(self, alignment):
        """Formate l'alignement pour l'affichage"""
        match_line = ''
        for i in range(len(alignment['query_align'])):
            if alignment['query_align'][i] == alignment['subject_align'][i]:
                match_line += '|'
            else:
                match_line += ' '
        
        return f"""
Score: {alignment['score']}
Query:   {alignment['query_start']}-{alignment['query_end']} {alignment['query_align']}
         {match_line}
Subject: {alignment['subject_start']}-{alignment['subject_end']} {alignment['subject_align']}
"""

def get_memory_usage():
    """Retourne l'utilisation actuelle de la mémoire en Mo"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # Conversion en Mo

def analyze_complexity():
    """Analyse les complexités temporelle et spatiale de l'algorithme BLAST"""
    blast = BLAST(word_size=3, match_score=2, mismatch_score=-1, gap_score=-2, threshold=5)
    
    sequence_sizes = list(range(100, 2100, 200))
    times = []
    memory_usages = []
    baseline_memory = get_memory_usage()  # Mémoire de base avant les tests
    
    print("Analyse des complexités temporelle et spatiale en cours...")
    total_tests = len(sequence_sizes)
    
    for i, size in enumerate(sequence_sizes):
        print(f"Progression : {i+1}/{total_tests} - Taille actuelle : {size}")
        
        # Mesures moyennes sur plusieurs essais
        trial_times = []
        trial_memories = []
        
        for trial in range(3):
            # Génération des séquences
            query = blast.generate_random_sequence(size)
            subject = blast.generate_random_sequence(size)
            
            # Mesure du temps
            start_time = time.time()
            blast.align(query, subject)
            end_time = time.time()
            
            # Mesure de la mémoire
            current_memory = get_memory_usage() - baseline_memory
            
            trial_times.append(end_time - start_time)
            trial_memories.append(current_memory)
        
        avg_time = np.mean(trial_times)
        avg_memory = np.mean(trial_memories)
        
        times.append(avg_time)
        memory_usages.append(avg_memory)
        
        print(f"Temps moyen: {avg_time:.4f} secondes")
        print(f"Utilisation mémoire: {avg_memory:.2f} Mo")
        print("-" * 50)

    # Création des graphiques
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Graphique temporel
    ax1.plot(sequence_sizes, times, 'bo-', label='Temps mesuré')
    ax1.set_xlabel('Taille des séquences')
    ax1.set_ylabel('Temps d\'exécution (secondes)')
    ax1.set_title('Complexité Temporelle')
    ax1.grid(True)
    
    # Graphique spatial (mémoire)
    ax2.plot(sequence_sizes, memory_usages, 'ro-', label='Mémoire utilisée')
    ax2.set_xlabel('Taille des séquences')
    ax2.set_ylabel('Utilisation mémoire (Mo)')
    ax2.set_title('Complexité Spatiale')
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig('blast_complexity_analysis.png')
    plt.close()
    
    # Calcul des corrélations pour la mémoire
    n = np.array(sequence_sizes) / sequence_sizes[0]
    m = np.array(memory_usages) / memory_usages[0]
    
    memory_correlations = {
        'O(n)': np.corrcoef(m, n)[0,1],
        'O(n log n)': np.corrcoef(m, n * np.log(n))[0,1],
        'O(n²)': np.corrcoef(m, n**2)[0,1]
    }
    
    return sequence_sizes, times, memory_usages, memory_correlations

if __name__ == "__main__":
    sizes, times, memories, memory_correlations = analyze_complexity()
    
    print("\nRésultats de l'analyse de complexité spatiale:")
    print("\nTaille des séquences vs Utilisation mémoire:")
    for size, memory in zip(sizes, memories):
        print(f"Taille: {size:4d} | Mémoire: {memory:.2f} Mo")
    
    print("\nCorrélations avec les complexités théoriques (mémoire):")
    for complexity, corr in memory_correlations.items():
        print(f"{complexity}: {corr:.4f}")