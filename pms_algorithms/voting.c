#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define TABLE_SIZE 40000003
#define MAX_VARIANTS 100000
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

unsigned int hash(const char *key, int table_size) {
    unsigned long hash = 5381;
    int c;
    while ((c = *key++)) {
        hash = ((hash << 5) + hash) + c;
    }
    return hash % table_size;
}

HashTable* create_table(int size) {
    HashTable *table = malloc(sizeof(HashTable));
    table->size = size;
    table->buckets = calloc(size, sizeof(Entry*));
    return table;
}

void insert(HashTable *table, const char *key, int value) {
    unsigned int bucket = hash(key, table->size);
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

int lookup(HashTable *table, const char *key) {
    unsigned int bucket = hash(key, table->size);
    Entry *current = table->buckets[bucket];
    while (current != NULL) {
        if (strcmp(current->key, key) == 0) {
            return current->value;
        }
        current = current->next;
    }
    return -1;
}

void increment(HashTable *table, HashTable *R, const char *key, int i) {
    unsigned int bucket = hash(key, table->size);
    Entry *current = table->buckets[bucket];
    if (i==0 && lookup(R, key) != 0){
        if (current == NULL){
            Entry *new_entry = malloc(sizeof(Entry));
            new_entry->key = strdup(key);
            new_entry->value = 1;
            new_entry->next = table->buckets[bucket];
            table->buckets[bucket] = new_entry;
            insert(R, key, i);
        }
    }else if (lookup(R, key) != i){
        while (current != NULL) {
            if (strcmp(current->key, key) == 0) {
                current->value++;
                insert(R, key, i);
                return;
            }
            current = current->next;
        }
    }
}

void free_table(HashTable *table) {
    for (int i = 0; i < table->size; i++) {
        Entry *current = table->buckets[i];
        while (current != NULL) {
            Entry *tmp = current;
            current = current->next;
            free(tmp->key);
            free(tmp);
        }
    }
    free(table->buckets);
    free(table);
}

char* slice(char *sequence, int start, int length) {
    char *substring = malloc((length + 1) * sizeof(char));
    if (substring == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    strncpy(substring, &sequence[start], length);
    substring[length] = '\0';
    return substring;
}

void generate_neighbors(char *str, char *neighbor, int pos, int d, int length, char ***results, int *result_count, int *result_capacity, char* alphabet) {
    /*
        Recursive function that finds all neighbours until differences are less or equal to d.
        There are two cases, when recursive function is called but we keep the character on position pos, and d stays the same.
        Second case is when the letter is changed with all 3 different nucleotids and recursive call is made with d-1.
    */
    
    if (d < 0) {
        return;
    }
    if (pos == length) {
        if (*result_count >= *result_capacity) {
            *result_capacity *= 4;
            *results = realloc(*results, *result_capacity * sizeof(char *));
            if (*results == NULL) {
                perror("Failed to realloc");
                exit(EXIT_FAILURE);
            }
        }
        (*results)[*result_count] = strdup(neighbor);
        if ((*results)[*result_count] == NULL) {
            perror("Failed to strdup");
            exit(EXIT_FAILURE);
        }
        (*result_count)++;
        return;
    }

    neighbor[pos] = str[pos];
    generate_neighbors(str, neighbor, pos + 1, d, length, results, result_count, result_capacity, alphabet);

    for (int i = 0; alphabet[i] != '\0'; i++) {
        if (alphabet[i] != str[pos]) {
            neighbor[pos] = alphabet[i];
            generate_neighbors(str, neighbor, pos + 1, d - 1, length, results, result_count, result_capacity, alphabet);
        }
    }
}

void get_variants(char *s, int d, char ***variants, int *count, char* alphabet) {
    /*
        Allocate memory for neighbours of lmer s and call recursive function for finding them all.
     */
    int length = strlen(s);
    int capacity = MAX_VARIANTS;
    *variants = malloc(capacity * sizeof(char *));
    if (*variants == NULL) {
        perror("Failed to malloc");
        exit(EXIT_FAILURE);
    }
    char *neighbor = malloc(length + 1);
    if (neighbor == NULL) {
        perror("Failed to malloc");
        exit(EXIT_FAILURE);
    }
    neighbor[length] = '\0';

    *count = 0;
    generate_neighbors(s, neighbor, 0, d, length, variants, count, &capacity, alphabet);

    free(neighbor);
}

void voting_algorthm(char **sequences, int num_sequences, int l, int d, char* alphabet) {
    /*
        Initialization of hash tables, V to count votes and R to keep track of which sequence voted.
        For each lmer from the input sequences, find all its neighbors and increment the result for all its neighbors in table V.
        It is not allowed for more than one l-mer from the same string to cast more than one vote for one neighbor.
        After voting step is done, next step is checking votes and printing found motifs.
    */
    HashTable *V = create_table(TABLE_SIZE);
    HashTable *R = create_table(TABLE_SIZE);

    for (int i = 0; i < num_sequences; i++) {
        for (int j = 0; j <= strlen(sequences[i]) - l; j++) {
            char *lmer = slice(sequences[i], j, l);
            char **variants;
            int variant_count;
            get_variants(lmer, d, &variants, &variant_count, alphabet);

            for (int k = 0; k < variant_count; k++) {
                char *variant = variants[k];
                increment(V, R, variant, i);
                // if (lookup(R, variant) != i) {
                //     increment(V, variant, i);
                //     insert(R, variant, i);
                // }
            }
            for (int k = 0; k < variant_count; k++) {
                free(variants[k]);
            }
            free(lmer);
            free(variants);
        }
    }

    for (int i = 0; i < TABLE_SIZE; i++) {
        Entry *entry = V->buckets[i];
        if (entry == NULL) {
            continue;
        }
        while (entry != NULL) {
            if (entry->value == num_sequences){
                printf("Pronadjen motiv: %s\n", entry->key);
            }
            entry = entry->next;
        }
    }
    free_table(V);
    free_table(R);
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
        fprintf(stderr, "Argumenti duzina motiva i broj dozvoljenih mutacije moraju biti veci od 0.\n");
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
    voting_algorthm(sequences, num_sequences, l, d, alphabet);
    clock_t end = clock();

    printf("Vreme: %lfs\n", (double)(end - start) / CLOCKS_PER_SEC);
    
    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i]);
    }
    free(sequences);
    free(alphabet);
    return 0;
}
