import random
from collections import defaultdict
import time
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

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

def analyze_time_complexity():
    """
    Analyse la complexité temporelle de l'algorithme BLAST en mesurant
    le temps d'exécution pour différentes tailles de séquences
    """
    # Paramètres pour l'analyse
    blast = BLAST(word_size=3, match_score=2, mismatch_score=-1, gap_score=-2, threshold=5)
    
    # Différentes tailles de séquences à tester
    sequence_sizes = list(range(100, 2100, 200))  # De 100 à 2000 par pas de 200
    times = []  # Pour stocker les temps d'exécution
    
    # Test pour chaque taille
    print("Analyse de la complexité temporelle en cours...")
    for size in tqdm(sequence_sizes):
        # Moyenne sur plusieurs essais pour chaque taille
        trial_times = []
        for _ in range(3):  # 3 essais par taille
            # Génération des séquences
            query = blast.generate_random_sequence(size)
            subject = blast.generate_random_sequence(size)
            
            # Mesure du temps d'exécution
            start_time = time.time()
            blast.align(query, subject)
            end_time = time.time()
            
            trial_times.append(end_time - start_time)
        
        # Calcul de la moyenne des essais
        avg_time = np.mean(trial_times)
        times.append(avg_time)
        
        # Affichage des résultats intermédiaires
        print(f"\nTaille: {size}")
        print(f"Temps moyen: {avg_time:.4f} secondes")

    # Analyse de la complexité
    # Calcul des courbes théoriques pour comparaison
    sizes_np = np.array(sequence_sizes)
    times_np = np.array(times)
    
    # Ajustement pour O(n), O(n log n) et O(n²)
    n = sizes_np / sizes_np[0]  # normalisation
    t = times_np / times_np[0]  # normalisation
    
    # Création du graphique
    plt.figure(figsize=(12, 8))
    
    # Données réelles
    plt.plot(sizes_np, times_np, 'bo-', label='Temps mesuré', linewidth=2)
    
    # Courbes théoriques
    plt.plot(sizes_np, times_np[0] * n, 'g--', label='O(n)', alpha=0.5)
    plt.plot(sizes_np, times_np[0] * n * np.log(n), 'r--', label='O(n log n)', alpha=0.5)
    plt.plot(sizes_np, times_np[0] * n**2, 'y--', label='O(n²)', alpha=0.5)
    
    plt.xlabel('Taille des séquences')
    plt.ylabel('Temps d\'exécution (secondes)')
    plt.title('Analyse de la complexité temporelle de BLAST')
    plt.legend()
    plt.grid(True)
    
    # Calcul de la corrélation avec différentes complexités
    correlations = {
        'O(n)': np.corrcoef(t, n)[0,1],
        'O(n log n)': np.corrcoef(t, n * np.log(n))[0,1],
        'O(n²)': np.corrcoef(t, n**2)[0,1]
    }
    
    print("\nAnalyse des corrélations:")
    for complexity, corr in correlations.items():
        print(f"{complexity}: {corr:.4f}")
    
    # Sauvegarde du graphique
    plt.savefig('blast_complexity.png')
    plt.close()
    
    return sequence_sizes, times, correlations

if __name__ == "__main__":
    sizes, times, correlations = analyze_time_complexity()
    
    # Affichage des résultats finaux
    print("\nRésultats de l'analyse de complexité:")
    print("\nTaille des séquences vs Temps d'exécution:")
    for size, time_taken in zip(sizes, times):
        print(f"Taille: {size:4d} | Temps: {time_taken:.4f} secondes")
    
    print("\nCorrélations avec les complexités théoriques:")
    for complexity, corr in correlations.items():
        print(f"{complexity}: {corr:.4f}")