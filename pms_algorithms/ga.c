#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_LINES 1000
#define MAX_BUFFER 1200
#define MAX_LINE_LENGTH 100

char random_base(char* alphabet) {
    return alphabet[rand() % 4];
}

char** generate_combinations(char* alphabet, int alphabet_size) {
    /* Creating all combination for nucleotides of size 3. */
    int num_combinations = alphabet_size * alphabet_size * alphabet_size;
    char **combinations = (char **)malloc((num_combinations) * sizeof(char *));
    
    if (combinations == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    int index = 0;
    for (int i = 0; i < alphabet_size; i++) {
        for (int j = 0; j < alphabet_size; j++) {
            for (int k = 0; k < alphabet_size; k++) {
                combinations[index] = (char *)malloc((3 + 1) * sizeof(char));
                if (combinations[index] == NULL) {
                    fprintf(stderr, "Memory allocation failed\n");
                    exit(1);
                }
                combinations[index][0] = alphabet[i];
                combinations[index][1] = alphabet[j];
                combinations[index][2] = alphabet[k];
                combinations[index][3] = '\0';
                index++;
            }
        }
    }

    return combinations;
}


char* mutation(const char* motif, double mutation_prob, char* alphabet) {
    /* With given probability change motif on one with random nucleotid */
    int len = strlen(motif);
    char* mutated_motif = (char*)malloc((len + 1) * sizeof(char));
    strcpy(mutated_motif, motif);

    for (int i = 0; i < len; i++) {
        if ((double)rand() / RAND_MAX < mutation_prob) {
            mutated_motif[i] = random_base(alphabet);
        }
    }
    return mutated_motif;
}


double compute_fitness_score(char* motif, char** seqs, int seq_count) {
    /* Fitness score for motif is done by going through all kmers in input sequences and 
    finding for each sequence the smaller number of mismatch, score is the sum of this for all sequences. */
    int motif_length = strlen(motif);
    double total_score = 0.0;

    for (int i = 0; i < seq_count; i++) {
        int sequence_length = strlen(seqs[i]);
        int max_alignment_position = sequence_length - motif_length + 1;

        double best_fitness_score = 0.0;
        for (int alignment_position = 0; alignment_position < max_alignment_position; alignment_position++) {
            double current_score = 0.0;
            for (int j = 0; j < motif_length; j++) {
                if (seqs[i][alignment_position + j] == motif[j]) {
                    current_score += 1.0;
                }
            }
            double current_fitness_score = current_score / motif_length;
            if (current_fitness_score > best_fitness_score) {
                best_fitness_score = current_fitness_score;
            }
        }
        total_score += best_fitness_score;
    }
    return total_score / seq_count;
}


char* addition(const char* motif, char** seqs, int seq_count, char* alphabet) {
    /* Adding random nucleotid in the beginning and in the end of given motif and the one with higher score stays. */
    int len = strlen(motif);
    char* first = (char*)malloc((len + 2) * sizeof(char));
    char* second = (char*)malloc((len + 2) * sizeof(char));
    first[0] = random_base(alphabet);
    strcpy(first + 1, motif);
    strcpy(second, motif);
    second[len] = random_base(alphabet);
    second[len + 1] = '\0';

    double fs = compute_fitness_score(first, seqs, seq_count);
    double ss = compute_fitness_score(second, seqs, seq_count);

    char* result = fs > ss ? first : second;
    if (result == first) {
        free(second);
    } else {
        free(first);
    }
    return result;
}

char* deletion(const char* motif) {
    /* Delete last nucleotid. */
    int len = strlen(motif);
    if (len <= 1) return NULL;
    char* deleted_motif = (char*)malloc(len * sizeof(char));
    strncpy(deleted_motif, motif, len - 1);
    deleted_motif[len - 1] = '\0';
    return deleted_motif;
}

typedef struct {
    char* motif;
    double score;
} MotifScore;

int compare_motif_scores(const void* a, const void* b) {
    MotifScore* ms_a = (MotifScore*)a;
    MotifScore* ms_b = (MotifScore*)b;
    return (ms_b->score - ms_a->score) > 0 ? 1 : -1;
}

MotifScore iter_algorithm(char** seqs, int seq_count, int maxL, int maxloop, double mutation_prob, char* alphabet, int alphabet_size) {
    /*  Initialization of population with all kmers of length 3.
        For each fixed length of the motif (while loop) go through population and perform operations on it maxloop times.
        Every time check if the score of motif is better than current best, and if it is save it.
     */
    int initial_motif_count = alphabet_size*alphabet_size*alphabet_size;
    MotifScore* motifs = (MotifScore*)malloc(initial_motif_count*sizeof(MotifScore));
    
    for (int i = 0; i < initial_motif_count; i++) {
        char ** motif = generate_combinations(alphabet, alphabet_size);
        motifs[i].motif = strdup(motif[i]);
        motifs[i].score = compute_fitness_score(motif[i], seqs, seq_count);
    }

    int L = 3;
    while (L < maxL) {
        for (int i = 0; i < initial_motif_count; i++) {
            char* new_motif = addition(motifs[i].motif, seqs, seq_count, alphabet);
            double new_motif_score = compute_fitness_score(new_motif, seqs, seq_count);

            for (int loop = 0; loop < maxloop; loop++) {
                char* updated_motif = deletion(new_motif);
                if (updated_motif == NULL) {
                    free(new_motif);
                    break;
                }

                if (L > 3) {
                    char* mutated_motif = mutation(updated_motif, mutation_prob, alphabet);
                    free(updated_motif);
                    updated_motif = addition(mutated_motif, seqs, seq_count, alphabet);
                    free(mutated_motif);
                } else {
                    char* temp_motif = addition(updated_motif, seqs, seq_count, alphabet);
                    free(updated_motif);
                    updated_motif = temp_motif;
                }

                double new_upd_motif_score = compute_fitness_score(updated_motif, seqs, seq_count);
                if (new_motif_score < new_upd_motif_score) {
                    free(new_motif);
                    new_motif = updated_motif;
                    new_motif_score = new_upd_motif_score;
                } else {
                    free(updated_motif);
                }
            }
            free(motifs[i].motif);
            motifs[i].motif = new_motif;
            motifs[i].score = new_motif_score;
        }
        L++;
    }

    qsort(motifs, initial_motif_count, sizeof(MotifScore), compare_motif_scores);
    return motifs[0];
}
char **read_lines_from_file(const char *file_path, int *num_sequences) {
    FILE *file = fopen(file_path, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char **lines = (char **)malloc(MAX_LINES * sizeof(char *));
    char buffer[MAX_BUFFER];
    *num_sequences = 0;

    while (fgets(buffer, sizeof(buffer), file)) {
        buffer[strcspn(buffer, "\n")] = '\0';
        lines[*num_sequences] = strdup(buffer);
        (*num_sequences)++;
    }

    fclose(file);
    return lines;
}


char* read_alphabet(const char *file_name, int *num_chars) {
    FILE *file = fopen(file_name, "r");
    if (!file) {
        perror("Failed to open file for azbuka\n");
        return NULL;
    }
    char *line = (char *)malloc(MAX_LINE_LENGTH * sizeof(char));

    if (fgets(line, MAX_LINE_LENGTH, file) != NULL) {
        *num_chars = 0;
        while (line[*num_chars] != '\0' && line[*num_chars] != '\n') {
            (*num_chars)++;
        }
        return line;
    } else {
        free(line); 
        return NULL;
    }
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    if (argc < 5) {
        fprintf(stderr, "Argumenti: <duzina_motiva> <datoteka_sa_sekvencama> <verovatnoca_mutacija [0,1]> <broj_iteracija> [<datoteka_sa_azbukom>]\n");
        return 1;
    }

    int k = atoi(argv[1]);
    int maxloop = atoi(argv[4]);
    double mutation_prob = strtod(argv[3], NULL);

    int num_sequences;
    char **sequences = read_lines_from_file(argv[2], &num_sequences);

    char *alphabet = (char *)malloc(25*sizeof(char));
    int alphabet_size;

    if (argc == 6) {
        alphabet = read_alphabet(argv[5], &alphabet_size);
    } else{
        alphabet = strdup("acgt");
        alphabet_size = 4;
    }

    clock_t start = clock();
    MotifScore best_motif = iter_algorithm(sequences, num_sequences, k, maxloop, mutation_prob, alphabet, alphabet_size);
    clock_t end = clock();
    printf("Vreme: %fs\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("Najbolji motiv: %s, Ocena: %.2f\n", best_motif.motif, best_motif.score);


    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i]);
    }
    free(sequences);
    free(best_motif.motif);
    free(alphabet);
    return 0;
}

