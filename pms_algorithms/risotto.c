#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#define MAX_MOTIF_LENGTH 25
#define MAX_VALID_MOTIFS 10
#define HASH_SIZE 100000
#define MAX_LINES 1000
#define MAX_BUFFER 1200
#define MAX_LINE_LENGTH 100


typedef struct Entry {
    char *key;
    int value;
    struct Entry *next;
} Entry;

typedef struct {
    Entry **buckets;
    int size;
} HashTable;

int hamming_distance(const char *seq1, const char *seq2, int length) {
    int distance = 0;
    for (int i = 0; i < length; i++) {
        if (seq1[i] != seq2[i]) {
            distance++;
        }
    }
    return distance;
}

bool is_valid(const char *motif, char **sequences, int num_sequences, int quorum, int max_mismatches) {
    /* Checking if quorum is satisfied. */
    bool *is_there = (bool *)malloc(num_sequences * sizeof(bool));
    memset(is_there, false, num_sequences * sizeof(bool));
    int mot_len = strlen(motif);
    for (int j = 0; j < num_sequences; j++) {
        for (int i = 0; i <= strlen(sequences[j]) - mot_len; i++) {
            if (hamming_distance(motif, &sequences[j][i], mot_len) <= max_mismatches) {
                is_there[j] = true;
                break;
            }
        }
    }

    int count = 0;
    for (int j = 0; j < num_sequences; j++) {
        if (is_there[j]) {
            count++;
        }
    }
    free(is_there);
    return count >= quorum;
}

unsigned int hash_fun(const char *key, int table_size) {
    unsigned long hash = 5381;
    int c;
    while ((c = *key++)) {
        hash = ((hash << 5) + hash) + c;
    }
    return hash % table_size;
}

void insert(HashTable *table, const char *key, int value) {
    unsigned int bucket = hash_fun(key, table->size);
    Entry *current = table->buckets[bucket];
    while (current != NULL) {
        if (strcmp(current->key, key) == 0) {
            current->value = value;
            return;
        }
        current = current->next;
    }

    Entry *new_entry = malloc(sizeof(Entry));
    new_entry->key = strdup(key);
    new_entry->value = value;
    new_entry->next = table->buckets[bucket];
    table->buckets[bucket] = new_entry;
}

void free_table(HashTable *table) {
    for (int i = 0; i < table->size; i++) {
        Entry *entry = table->buckets[i];
        while (entry != NULL) {
            Entry *temp = entry;
            entry = entry->next;
            free(temp->key);
            free(temp);
        }
    }
    free(table->buckets);
    free(table);
}

HashTable* create_table(int size) {
    HashTable *table = malloc(sizeof(HashTable));
    table->size = size;
    table->buckets = malloc(size * sizeof(Entry*));

    for (int i = 0; i < size; i++) {
        table->buckets[i] = NULL;
    }

    return table;
}


int lookup(HashTable *table, const char *key) {
    unsigned int bucket = hash_fun(key, table->size);

    Entry *entry = table->buckets[bucket];
    while (entry != NULL && strcmp(entry->key, key) != 0) {
        entry = entry->next;
    }

    if (entry == NULL) {
        return -1;
    } else {
        return entry->value;
    }
}

void extract_single_motif(const char *motif, char **sequences, int num_sequences, int quorum, int max_mismatches, int k_min, int k_max, char **valid_motifs, int *valid_motif_count, HashTable *max_ext, char* alphabet, int alphabet_size) {
    
    /*
        Recurive function that in lexicographical order adds nucleotides and checks if the new motifs are valid.
        In the max_ect table it stores the value of maximal extensibility for already visited suffixes and uses that
        info for earlier stopping if the egde doesn't lead to motif.
    */
    
    char motif_alpha[MAX_MOTIF_LENGTH];
    for (int a = 0; a < alphabet_size; a++) {
        snprintf(motif_alpha, sizeof(motif_alpha), "%s%c", motif, alphabet[a]);
        char *x = motif_alpha;

        while (strlen(x) > 0 && lookup(max_ext, x) == -1) {
            x++;
        }

        if (strlen(x) > 0 && lookup(max_ext, x) + strlen(motif_alpha) < k_min) {
            int value = lookup(max_ext, motif_alpha);
            insert(max_ext, motif_alpha, value);
            continue;
        }
        if (is_valid(motif_alpha, sequences, num_sequences, quorum, max_mismatches)) {
            if (strlen(motif_alpha) >= k_min) {
                valid_motifs[(*valid_motif_count)++] = strdup(motif_alpha);
            }
            if (strlen(motif_alpha) < k_max) {
                extract_single_motif(motif_alpha, sequences, num_sequences, quorum, max_mismatches, k_min, k_max, valid_motifs, valid_motif_count, max_ext, alphabet, alphabet_size);
            } else {
                insert(max_ext, motif_alpha, INT_MAX);
            }
        } else {
            if (strlen(motif_alpha) < k_min) {
                insert(max_ext, motif_alpha, 0);
            } else{ 
                char motif_prefix[MAX_MOTIF_LENGTH];
                int j = 0;
                for (j = 0; j < (k_min - 1); j++){
                    motif_prefix[j] = motif_alpha[j];
                }
                motif_prefix[j] = '\0';
                if (lookup(max_ext, motif_prefix) < strlen(motif_alpha) - (k_min - 1)) {
                    insert(max_ext, motif_prefix, strlen(motif_alpha) - (k_min - 1));
                }
            }
        }
    }

    /*Each node gets max_ext value that is taken from its child with maximum value and increased by 1*/
    if (strlen(motif) < k_min - 1) {
        int max_child = -1;
        for (int a = 0; a < alphabet_size; a++) {
            char child[25];
            char c = alphabet[a];
            strcpy(child, motif);
            char temp[2];
            temp[0] = c;
            temp[1] = '\0';
            strcat(child, temp);
            if (max_child < lookup(max_ext, child)){
                max_child = lookup(max_ext, child);
            }
        }
        insert(max_ext, motif, max_child + 1);
    }
}

char **read_lines_from_file(const char *file_path, int *num_lines) {
    FILE *file = fopen(file_path, "r");
    if (!file) {
        perror("Failed to open file");
        return NULL;
    }

    char **lines = (char **)malloc(MAX_LINES * sizeof(char *));
    char buffer[MAX_BUFFER];
    *num_lines = 0;

    while (fgets(buffer, sizeof(buffer), file)) {
        lines[*num_lines] = strdup(buffer);
        (*num_lines)++;
    }
    fclose(file);
    return lines;
}


char* read_alphabet(const char *file_name, int *num_chars) {
    FILE *file = fopen(file_name, "r");
    if (!file) {
        perror("Failed to open file");
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

    char *valid_motifs[MAX_VALID_MOTIFS];
    int valid_motif_count = 0;
    HashTable *max_ext = create_table(HASH_SIZE);

    clock_t start = clock();
    extract_single_motif("", sequences, num_sequences, num_sequences, d, l, l, valid_motifs, &valid_motif_count, max_ext, alphabet, alphabet_size);
    clock_t end = clock();

    printf("Pronadjeni motivi: \n");
    for (int i = 0; i < valid_motif_count; i++) {
        printf("%s\n", valid_motifs[i]);
        free(valid_motifs[i]);
    }
    printf("Vreme: %lfs\n", (double)(end - start) / CLOCKS_PER_SEC);

    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i]);
    }
    free(sequences);
    free_table(max_ext);
    free(alphabet);
    return 0;
}