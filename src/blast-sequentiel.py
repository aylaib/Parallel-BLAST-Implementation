import random
from collections import defaultdict
import time

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

def main():
    # Paramètres de l'algorithme BLAST
    blast = BLAST(word_size=3, match_score=2, mismatch_score=-1, gap_score=-2, threshold=5)
    
    # Demander la taille des séquences à l'utilisateur pour chaque séquence
    try:
        query_length = int(input("Entrez la taille de la séquence requête : "))
        subject_length = int(input("Entrez la taille de la séquence sujet : "))
        
        if query_length < 10 or subject_length < 10:
            print("Les tailles doivent être d'au moins 10 nucléotides.")
            return
    except ValueError:
        print("Veuillez entrer des nombres entiers valides.")
        return

    # Génération des séquences aléatoires avec les tailles spécifiées
    query_seq = blast.generate_random_sequence(query_length)
    subject_seq = blast.generate_random_sequence(subject_length)
    
    print("\nSéquence requête générée :", query_seq)
    print("Séquence sujet générée  :", subject_seq)
    print(f"Tailles - Requête: {len(query_seq)}, Sujet: {len(subject_seq)}")
    
    # Mesure du temps d'exécution
    start_time = time.time()
    
    # Recherche des alignements
    alignments = blast.align(query_seq, subject_seq)
    
    execution_time = time.time() - start_time
    
    # Affichage des résultats
    print(f"\nTemps d'exécution : {execution_time:.4f} secondes")
    print(f"\nNombre d'alignements trouvés : {len(alignments)}")
    
    if alignments:
        print("\nMeilleurs alignements trouvés :")
        for i, alignment in enumerate(alignments[:3], 1):  # Affiche les 3 meilleurs alignements
            print(f"\nAlignement #{i}")
            print(blast.format_alignment(alignment))
    else:
        print("\nAucun alignement significatif trouvé.")

if __name__ == "__main__":
    main()
