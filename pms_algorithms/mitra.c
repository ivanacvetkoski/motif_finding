#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define MAX_LINES 1000
#define MAX_BUFFER 1200
#define MAX_SEQ_LEN 650
#define MAX_LINE_LENGTH 100

typedef struct Dictionary{
    int *table;
    int size;
    struct Dictionary* next;
} Dictionary;

typedef struct Node {
    int id;
    int distance;
    char *motif;
} Node;


typedef struct MotifNode {
    char *motif;
    struct MotifNode *next;
} MotifNode;


MotifNode* create_motif_node(const char *motif) {
    MotifNode *new_node = (MotifNode *)malloc(sizeof(MotifNode));
    if (new_node == NULL) {
        printf("Memory allocation failed create node1\n");
        exit(1);
    }
    new_node->motif = strdup(motif);
    if (new_node->motif == NULL) {
        printf("Memory allocation failed create node2\n");
        free(new_node);
        exit(1);
    }
    new_node->next = NULL;
    return new_node;
}

void add_motif(MotifNode **head, const char *motif) {
    MotifNode *new_node = create_motif_node(motif);
    new_node->next = *head;
    *head = new_node;
}

void free_motifs(MotifNode *head) {
    MotifNode *current = head;
    while (current != NULL) {
        MotifNode *next = current->next;
        free(current->motif);
        free(current);
        current = next;
    }
}

Dictionary *make_index_list(char **seq, int num_sequences, int l) {
    int total_length = MAX_SEQ_LEN*num_sequences;
    Dictionary *ind_dict = (Dictionary*)malloc(sizeof(Dictionary));
    ind_dict->table = (int *)calloc(total_length, sizeof(int));

    if (ind_dict->table == NULL) {
        printf("Memory allocation failed make index list\n");
        exit(1);
    }
    ind_dict->size = total_length;
    ind_dict->next = NULL;

    return ind_dict;
}

Node* make_node(int distance, char *motif) {
    Node *node = (Node *)malloc(sizeof(Node));
    if (node == NULL) {
        printf("Memory allocation failed make node1\n");
        exit(1);
    }
    node->distance = distance;
    node->motif = strdup(motif);
    if (node->motif == NULL) {
        printf("Memory allocation failed make node2\n");
        free(node);
        exit(1);
    }
    return node;
}


Node *get_next_node(Node *node, int position, int l, char **sequences, int seq_len, char *motif, char e) {
    /*
        Creating child node from parent by concatening char e to motif.
        Each lmer in input sequences has unique position (sequence index and position in t) so the lmers are saved with those
        two number instead of whole strings.
    */

    size_t new_len = strlen(motif) + 2;
    char *new_motif = (char *)malloc(new_len);
    if (new_motif == NULL) {
        printf("Memory allocation failed get next 2\n");
        exit(1);
    }
    strncpy(new_motif, motif, strlen(motif));
    new_motif[strlen(motif)] = e;
    new_motif[strlen(motif) + 1] = '\0';

    Node *next_node = make_node(node->distance + 1, new_motif);
    return next_node;
}

int hamming_distance(char *str1, char *str2) {
    int count = 0;
    for (int i = 0; i < strlen(str1); i++) {
        if (str1[i] != str2[i]) {
            count++;
        }
    }
    return count;
}


void virtual_dfs(Node *start, Dictionary *tables, int l, char **sequences, int num_sequences, int k, int d, MotifNode **motifs, char* alphabet, int alphabet_size) {
    if (start->distance == l) {
        add_motif(motifs, start->motif);
        return;
    }

    for (int i = 0; i < alphabet_size; i++) {

        Node *next_node = get_next_node(start, start->distance, l, sequences, num_sequences, start->motif, alphabet[i]);

        Dictionary *updated_table = (Dictionary*)malloc(sizeof(Dictionary));
        updated_table->table = (int *)malloc(MAX_SEQ_LEN*num_sequences*sizeof(int));
        memcpy(updated_table->table, tables->table, tables->size * sizeof(int));
        updated_table->size = tables->size;
        updated_table->next = NULL;
        int count = 0;

        for (int j = 0; j < updated_table->size; j++) {
            int whole = j / MAX_SEQ_LEN;
            int rem = j - MAX_SEQ_LEN * whole;
            char *lmer = (char *) malloc((l+1) * sizeof(char));
            int len_sq = strlen(sequences[whole]);
            if (rem > len_sq - l + 1){
                j += MAX_SEQ_LEN - rem - 1;
                free(lmer);
                continue;
            }
            strncpy(lmer, &sequences[whole][rem], l);
            lmer[l] = '\0';
            if (lmer[start->distance] != alphabet[i]) {
                int tree = updated_table->table[j];
                updated_table->table[j]++;
            }
            if (updated_table->table[j] <= d) {
                count++;
            }
            free(lmer);
        }

        if (count >= k) {
            virtual_dfs(next_node, updated_table, l, sequences, num_sequences, k, d, motifs, alphabet, alphabet_size);
        }

        free(updated_table->table);
        free(updated_table);
        free(next_node->motif);
        free(next_node);
    }
}

void filter_motifs(MotifNode **motifs, char **sequences, int num_sequences, int l, int d) {
    /* Go through all possible motifs and filter ones that don't have neighbours in each input sequences. */
    MotifNode* motif_list = (*motifs);
    while (motif_list) {
        int matches = 0;
        for (int j = 0; j < num_sequences; j++) {
            int s_len = strlen(sequences[j]);
            for (int k = 0; k < s_len - l + 1; k++) {
                char* current_substring = (char *) malloc((l+1)*sizeof(char));
                strncpy(current_substring, &sequences[j][k], l);
                current_substring[l] = '\0';

                int distance = hamming_distance(motif_list->motif, current_substring);
                if (distance <= d) {
                    matches++;
                    break;
                }
            }
        }

        if (matches == num_sequences) {
            printf("Motiv: %s\n", motif_list->motif);
        }
        motif_list = motif_list->next;
    }

}

void find_motifs(int l, int d, char **sequences, int num_sequences, char*alphabet, int alphabet_size) {
    /* 
        Create root node from which will traveling begun.
        Perform DFS and filter motifs.
     */
    Node *root_node = make_node(0, "");
    Dictionary *list_tables = make_index_list(sequences, num_sequences, l);

    list_tables->next = NULL;
    MotifNode *motifs = NULL;
    
    virtual_dfs(root_node, list_tables, l, sequences, num_sequences, num_sequences, d, &motifs, alphabet, alphabet_size);
    filter_motifs(&motifs, sequences, num_sequences, l, d);
    
    free_motifs(motifs);
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

    clock_t start = clock();
    find_motifs(l, d, sequences, num_sequences, alphabet, alphabet_size);
    clock_t end = clock();

    printf("Vreme: %lfs\n", (double)(end - start) / CLOCKS_PER_SEC);

    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i]);
    }
    free(sequences);
    return 0;
}
