#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define MAX_LINES 1000
#define MAX_BUFFER 1200
#define MAX_KMERS 100000000
#define NUM_SEQ 20
#define MAX_LINE_LENGTH 100

typedef struct {
    char *motif;
    int count;
} MotifResult;

int hamming_distance(const char *s1, const char *s2, int length) {
    int distance = 0;
    for (int i = 0; i < length; i++) {
        if (s1[i] != s2[i]) {
            distance++;
        }
    }
    return distance;
}
void generate_all_kmers(char **kmers, int k, int index, char *current_kmer, int *kmer_count, char* alphabet, int alphabet_size) {
    if (index == k) {
        current_kmer[index] = '\0';
        strcpy(kmers[*kmer_count], current_kmer);
        (*kmer_count)++;
        return;
    }

    for (int i = 0; i < alphabet_size; i++) {
        current_kmer[index] = alphabet[i];
        generate_all_kmers(kmers, k, index + 1, current_kmer, kmer_count, alphabet, alphabet_size);
    }
}

char **init_kmers(int k, int *num_kmers, int alphabet_size) {
    *num_kmers = pow(alphabet_size, k);
    char **kmers = malloc(*num_kmers * sizeof(char *));
    for (int i = 0; i < *num_kmers; i++) {
        kmers[i] = malloc((k + 1) * sizeof(char));
    }
    return kmers;
}

void free_kmers(char **kmers, int num_kmers) {
    for (int i = 0; i < num_kmers; i++) {
        free(kmers[i]);
    }
    free(kmers);
}

MotifResult motif_finding_with_mismatches(char **sequences, int num_sequences, int k, int d, char* alphabet, int alphabet_size) {
    /*
        Generate all possible combinations of the nucleotides and for each check if it has neighbours
        in all of the input sequences.
     */
    int num_kmers;
    char **all_kmers = init_kmers(k, &num_kmers, alphabet_size);
    char current_kmer[k + 1];
    int kmer_count = 0;
    generate_all_kmers(all_kmers, k, 0, current_kmer, &kmer_count, alphabet, alphabet_size);

    MotifResult best_motif_result = {NULL, 0};

    for (int i = 0; i < num_kmers; i++) {
        char *candidate = all_kmers[i];
        int *count = (int *)malloc(num_sequences*sizeof(int));
        for (int s = 0; s<num_sequences; s++){
            count[s] = 0;
        }

        for (int j = 0; j < num_sequences; j++) {
            char *seq = sequences[j];
            int seq_len = strlen(seq);

            for (int m = 0; m <= seq_len - k; m++) {
                char kmer[k + 1];
                strncpy(kmer, seq + m, k);
                kmer[k] = '\0';

                if (hamming_distance(candidate, kmer, k) <= d) {
                    count[j]=1;
                }
            }
        }
        bool all_seq = true;
        for (int w = 0; w < num_sequences; w++) {
            if (count[w] == 0) {
                all_seq = false;
            }
        }
        if (all_seq) {
            best_motif_result.motif = strdup(candidate);
        }
        free(count);
    }

    free_kmers(all_kmers, num_kmers);
    return best_motif_result;
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

    if (argc < 4) {
        fprintf(stderr, "Argumenti: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> [<datoteka_sa_azbukom>]\n");
        return 1;
    }

    int l = atoi(argv[1]);
    int d = atoi(argv[2]);

    if (l<0 || d<0){
        fprintf(stderr, "Argumenti duzina motiva i dozvoljene mutacije moraju biti veci od 0.\n");
        return 1;
    }

    int num_sequences;
    char **sequences = read_lines_from_file(argv[3], &num_sequences);

    char *alphabet = (char *)malloc(25*sizeof(char));
    int alphabet_size;

    if (argc == 5) {
        alphabet = read_alphabet(argv[4], &alphabet_size);
    } else{
        alphabet = strdup("acgt");
        alphabet_size = 4;
    }

    clock_t start = clock();
    MotifResult result = motif_finding_with_mismatches(sequences, num_sequences, l, d, alphabet, alphabet_size);
    clock_t end = clock();

    printf("Motiv: %s\n", result.motif);
    printf("Vreme: %lfs\n", (double)(end - start) / CLOCKS_PER_SEC);

    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i]);
    }
    free(sequences);
    free(alphabet);
    return 0;
}
