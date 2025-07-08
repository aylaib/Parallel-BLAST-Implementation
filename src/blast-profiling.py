import random
from collections import defaultdict
import time
import cProfile
import pstats
from pstats import SortKey

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

def main_with_profiling():
    # Paramètres de l'algorithme BLAST
    blast = BLAST(word_size=3, match_score=2, mismatch_score=-1, gap_score=-2, threshold=5)
    
    # Générer des séquences de test
    query_seq = blast.generate_random_sequence(1000)  # Séquence plus longue pour mieux voir les différences
    subject_seq = blast.generate_random_sequence(2000)
    
    # Démarrer le profiling
    profiler = cProfile.Profile()
    profiler.enable()
    
    # Exécuter l'alignement
    alignments = blast.align(query_seq, subject_seq)
    
    # Arrêter le profiling
    profiler.disable()
    
    # Analyser et afficher les résultats du profiling
    stats = pstats.Stats(profiler).sort_stats(SortKey.TIME)
    
    print("\n=== Analyse détaillée des performances ===")
    print("\n1. Top 10 des fonctions les plus chronophages:")
    stats.print_stats(10)
    
    print("\n2. Analyse des appels par fonction:")
    stats.print_callers(10)

if __name__ == "__main__":
    main_with_profiling()
