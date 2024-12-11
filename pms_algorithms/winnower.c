#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define MAX_SEQUENCES 20
#define MAX_SEQ_LENGTH 650
#define MAX_NEIGHBOURS 50000
#define HASH_TABLE_SIZE 15000
#define HASH_SET_SIZE 100
#define CAPACITY 200000
#define MAX_LINES 1000
#define MAX_BUFFER 1200
#define MAX_LINE_LENGTH 100


typedef struct Node {
    char *motif;
    int sequence_index;
    int position;
    unsigned int hash_value;
    struct Node **neighbours;
    int neigh_count;
    int neigh_capacity;
} Node;

typedef struct HashMap {
    Node **table;
    int size;
    int r_size;
} HashMap;

typedef struct {
    Node **nodes;
    int size;
} Clique;

typedef struct {
    char **data;
    int size;
    int capacity;
} ConsensusSet;

ConsensusSet* create_consensus_set(int capacity) {
    ConsensusSet *set = (ConsensusSet *)malloc(sizeof(ConsensusSet));
    set->data = (char **)malloc(capacity * sizeof(char *));
    set->size = 0;
    set->capacity = capacity;
    return set;
}

void add_consensus(ConsensusSet *set, char *consensus) {
    if (set->size >= set->capacity) {
        set->capacity *= 2;
        set->data = (char **)realloc(set->data, set->capacity * sizeof(char *));
    }
    set->data[set->size] = strdup(consensus);
    set->size++;
}

int contains_consensus(ConsensusSet *set, char *consensus) {
    for (int i = 0; i < set->size; i++) {
        if (strcmp(set->data[i], consensus) == 0) {
            return 1;
        }
    }
    return 0;
}

void free_consensus_set(ConsensusSet *set) {
    for (int i = 0; i < set->size; i++) {
        free(set->data[i]);
    }
    free(set->data);
    free(set);
}

void compute_consensus_motif(Node **clique, int clique_size, char *consensus, int l, char* alphabet, int alphabet_size) {
    /* Find k-mer from list of k-mers that have on each position nucleotide that appears on most in that position in all k-mers */
    int **counts = (int **)malloc(l * sizeof(int *));
    for (int i = 0; i < l; i++) {
        counts[i] = (int *)calloc(alphabet_size, sizeof(int));
    }
    for (int i = 0; i < clique_size; i++) {
        char *sequence = clique[i]->motif;
        for (int j = 0; j < l; j++) {
            for (int r = 0; r < alphabet_size; r++){
                if (sequence[j] == alphabet[r]){
                    counts[j][r]++;
                }
            }
        }
    }

    for (int i = 0; i < l; i++) {
        int max_count = 0;
        char consensus_base = 'a';
        for (int j = 0; j < alphabet_size; j++) {
            if (counts[i][j] > max_count) {
                max_count = counts[i][j];
                consensus_base = alphabet[j];
            }
        }
        consensus[i] = consensus_base;
    }
    consensus[l] = '\0';
}

int hash_func(int sequence_index, int position) {
    return sequence_index*MAX_SEQ_LENGTH + position;
}

Node* create_node(const char *motif, int sequence_index, int position) {
    Node *new_node = (Node *)malloc(sizeof(Node));
    new_node->motif = strdup(motif);
    new_node->sequence_index = sequence_index;
    new_node->position = position;
    new_node->hash_value = hash_func(sequence_index, position);
    new_node->neigh_capacity = 20;
    new_node->neighbours = (Node**)malloc(new_node->neigh_capacity * sizeof(Node*));
    new_node->neigh_count = 0;
    return new_node;
}

void add_neighbor(Node *node,  Node *node2) {
    if (node->neigh_count == node->neigh_capacity) {
        node->neigh_capacity *= 2;
        node->neighbours = (Node**)realloc(node->neighbours, node->neigh_capacity * sizeof(Node*));
    }
    node->neighbours[node->neigh_count++] = node2;
}


HashMap* create_hash_map() {
    HashMap *hashmap = (HashMap *)malloc(sizeof(HashMap));
    hashmap->table = (Node **)malloc(HASH_TABLE_SIZE * sizeof(Node *));
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        hashmap->table[i] = NULL;
    }
    hashmap->size = HASH_TABLE_SIZE;
    hashmap->r_size = 0;
    return hashmap;
}

void add_to_hash_map(HashMap *map, Node *node) {
    unsigned int h = node->hash_value;
    map->table[h] = node;
    map->r_size++;
}

Node* get_from_hash_map(HashMap *map, unsigned int hash_value) {
    unsigned int h = hash_value;
    Node *temp = map->table[h];
    if (temp != NULL) {
        if (temp->hash_value == hash_value) {
            return temp;
        }
    }
    return NULL;
}


void delete_first_element(Node **array, int *size) {
    if (*size <= 0) {
        return;
    }
    for (int i = 1; i < *size; i++) {
        array[i - 1] = array[i];
    }

    (*size)--;
}


int hamming_distance(const char *str1, const char *str2) {
    int distance = 0;
    while (*str1 && *str2) {
        if (*str1 != *str2) {
            distance++;
        }
        str1++;
        str2++;
    }
    return distance;
}

Node* find_or_add_node(HashMap *hashmap, const char *motif, int sequence_index, int position) {
    /* All k-mers have unique positions in input, and it can be represented through index of input sequence and position in it. */
    int hash_id = hash_func(sequence_index, position);
    Node *node = get_from_hash_map(hashmap, hash_id);
    if (node){
        return node;
    }
    Node *new_node = create_node(motif, sequence_index, position);
    add_to_hash_map(hashmap, new_node);

   return new_node;
}


void add_edge(HashMap *hashmap, const char *motif1, int seq_idx1, int pos1, const char *motif2, int seq_idx2, int pos2) {
    Node *node1 = find_or_add_node(hashmap, motif1, seq_idx1, pos1);
    Node *node2 = find_or_add_node(hashmap, motif2, seq_idx2, pos2);

    add_neighbor(node1, node2);
    add_neighbor(node2, node1);
}


void construct_hashmap(HashMap *hashmap, char *sequences[], int num_sequences, int l, int d) {
    /*
        Creating hashmap that will store all k-mers as keys and its neigbours as value.
        Condition for two kmers to be neighbours is that they differ on 2d or less positions.
    */
    for (int i = 0; i < num_sequences; i++) {
        for (int j = i + 1; j < num_sequences; j++) {
            for (int k = 0; k <= strlen(sequences[i]) - l; k++) {
                char motif1[l + 1];
                strncpy(motif1, sequences[i] + k, l);
                motif1[l] = '\0';
                for (int m = 0; m <= strlen(sequences[j]) - l; m++) {
                    char motif2[l + 1];
                    strncpy(motif2, sequences[j] + m, l);
                    motif2[l] = '\0';
                    if (hamming_distance(motif1, motif2) <= 2 * d) {
                        add_edge(hashmap, motif1, i, k, motif2, j, m);
                    }
                }
            }
        }
    }
}

void remove_edge(Node *node1, Node *node2) {
    /* Go through all neighbours of node1 and delete node2 in that list, and the other way around. */
    int found = 0;
    for (int i = 0; i < node1->neigh_count; i++) {
        if (node1->neighbours[i]->hash_value == node2->hash_value && strcmp(node1->neighbours[i]->motif, node2->motif)==0) {
            found = 1;
        }
        if (found && i < node1->neigh_count - 1) {
            node1->neighbours[i] = node1->neighbours[i + 1];
        }
    }
    if (found) {
        node1->neigh_count--;
    }

    found = 0;
    for (int i = 0; i < node2->neigh_count; i++) {
        if (node2->neighbours[i]->hash_value == node1->hash_value && strcmp(node2->neighbours[i]->motif, node1->motif)==0) {
            found = 1;
        }
        if (found && i < node2->neigh_count - 1) {
            node2->neighbours[i] = node2->neighbours[i + 1];
        }
    }
    if (found) {
        node2->neigh_count--;
    }
}

int compare_ints(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}

int* intersect_sets(int *set1, int size1, int *set2, int size2, int *intersection_size) {
    int *intersection = (int *)malloc((size1 < size2 ? size1 : size2) * sizeof(int));

    qsort(set1, size1, sizeof(int), compare_ints);
    qsort(set2, size2, sizeof(int), compare_ints);

    int i = 0, j = 0;
    int index = 0;
    while (i < size1 && j < size2) {
        if (set1[i] < set2[j]) {
            i++;
        } else if (set1[i] > set2[j]) {
            j++;
        } else {
            intersection[index++] = set1[i];
            i++;
            j++;
        }
    }

    *intersection_size = index;
    intersection = (int *)realloc(intersection, index * sizeof(int));

    return intersection;
}

void winnower_k2(HashMap *map, int q) {
    /*
        Checking for all egdes if number of common neighbours for nodes that edge connect is smaller than q-2.
        If it is we remove that edge because that edge cannot be in the clique.
     */
    int remove_count = 0;
    int capacity = CAPACITY;
    Node **edges_to_remove = (Node **)malloc(capacity * sizeof(Node *));
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        Node *node_first = map->table[i];
        if (node_first){
            int fr_neigh_hashes[node_first->neigh_count];
            for (int p=0; p<node_first->neigh_count; p++){
                fr_neigh_hashes[p] = node_first->neighbours[p]->hash_value;
            }
            for (int j = 0; j < node_first->neigh_count; j++) {
                Node* node_second = node_first->neighbours[j];

                if (!node_second) {
                    continue;
                }
                int sc_neigh_hashes[node_second->neigh_count];
                for (int p=0; p<node_second->neigh_count; p++){
                    sc_neigh_hashes[p] = node_second->neighbours[p]->hash_value;
                }
                
                int triangle_count;
                int *intersection = intersect_sets(fr_neigh_hashes, node_first->neigh_count, sc_neigh_hashes, node_second->neigh_count, &triangle_count);

                if (triangle_count < q - 2) {
                    if (remove_count >= capacity) {
                        capacity*=2;
                        Node **temp = (Node **)realloc(edges_to_remove, capacity * sizeof(Node *));
                        edges_to_remove = temp;
                    }
                    edges_to_remove[remove_count] = (Node *)malloc(sizeof(Node));
                    edges_to_remove[remove_count + 1] = (Node *)malloc(sizeof(Node));
                    edges_to_remove[remove_count] = node_first;
                    edges_to_remove[remove_count + 1] = node_second;
                    remove_count += 2;
                }
            }
        }
    }

    for (int i = 0; i < remove_count; i += 2) {
        remove_edge(edges_to_remove[i], edges_to_remove[i + 1]);
    }
    
    free(edges_to_remove);
}


void winnower_k3(HashMap *map, int q) {
    /*
        Checking for all egdes if number of cliques of size 4 that they build with all of their common neigbours is smaller than q-3, 
        we remove that edge because that edge cannot be in the clique.
     */
    int remove_count = 0;
    int capacity = CAPACITY;
    Node **edges_to_remove = (Node **)malloc(capacity * sizeof(Node *));
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        Node *node_first = map->table[i];
        if (node_first){
            int fr_neigh_hashes[node_first->neigh_count];
            for (int p=0; p<node_first->neigh_count; p++){
                fr_neigh_hashes[p] = node_first->neighbours[p]->hash_value;
            }
            for (int j = 0; j < node_first->neigh_count; j++) {
                Node* node_second = node_first->neighbours[j];

                if (!node_second) {
                    continue;
                }
                int sc_neigh_hashes[node_second->neigh_count];
                for (int p=0; p<node_second->neigh_count; p++){
                    sc_neigh_hashes[p] = node_second->neighbours[p]->hash_value;
                }

                int triangle_count;
                int *intersection = intersect_sets(fr_neigh_hashes, node_first->neigh_count, sc_neigh_hashes, node_second->neigh_count, &triangle_count);
                int counter_4 = 0;
                if (triangle_count > q - 3){
                    for (int r = 0; r < triangle_count; r++){
                        Node* node_thr = map->table[intersection[r]];
                        int th_neigh_hashes[node_thr->neigh_count];
                        for (int p=0; p<node_thr->neigh_count; p++){
                            th_neigh_hashes[p] = node_thr->neighbours[p]->hash_value;
                        }

                        int four;
                        int * cliq = intersect_sets(th_neigh_hashes, node_thr->neigh_count, intersection, triangle_count, &four);

                        if (four >= q - 3) {
                            counter_4++;
                        }
                    }
                }
                if (counter_4 < q - 2) {
                    if (remove_count >= capacity) {
                        capacity*=2;
                        Node **temp = (Node **)realloc(edges_to_remove, capacity * sizeof(Node *));
                        edges_to_remove = temp;
                    }
                    edges_to_remove[remove_count] = (Node *)malloc(sizeof(Node));
                    edges_to_remove[remove_count + 1] = (Node *)malloc(sizeof(Node));
                    edges_to_remove[remove_count] = node_first;
                    edges_to_remove[remove_count + 1] = node_second;
                    remove_count += 2;
                }
            }
        }
    }

    printf("Uklanjanje %d grana\n", remove_count / 2);
    for (int i = 0; i < remove_count; i += 2) {
        remove_edge(edges_to_remove[i], edges_to_remove[i + 1]);
    }
    free(edges_to_remove);
}

void winnower_k4(HashMap *map, int q) {
    /*
        Checking for all egdes if number of cliques of size 5 that they build with all of their common neigbours is smaller than q-3, 
        we remove that edge because that edge cannot be in the clique.
     */
    int remove_count = 0;
    int capacity = CAPACITY;
    Node **edges_to_remove = (Node **)malloc(capacity * sizeof(Node *));

    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        Node *node_first = map->table[i];
        if (node_first) {
            int fr_neigh_hashes[node_first->neigh_count];
            for (int p = 0; p < node_first->neigh_count; p++) {
                fr_neigh_hashes[p] = node_first->neighbours[p]->hash_value;
            }

            for (int j = 0; j < node_first->neigh_count; j++) {
                Node *node_second = node_first->neighbours[j];
                if (!node_second) continue;

                int sc_neigh_hashes[node_second->neigh_count];
                for (int p = 0; p < node_second->neigh_count; p++) {
                    sc_neigh_hashes[p] = node_second->neighbours[p]->hash_value;
                }

                int triangle_count;
                int *triangle_intersection = intersect_sets(fr_neigh_hashes, node_first->neigh_count, sc_neigh_hashes, node_second->neigh_count, &triangle_count);
                int counter_4 = 0;

                if (triangle_count >= q - 3) {
                    for (int r = 0; r < triangle_count; r++) {
                        Node *node_third = map->table[triangle_intersection[r]];
                        int th_neigh_hashes[node_third->neigh_count];
                        for (int p = 0; p < node_third->neigh_count; p++) {
                            th_neigh_hashes[p] = node_third->neighbours[p]->hash_value;
                        }

                        int quad_count;
                        int *quad_intersection = intersect_sets(th_neigh_hashes, node_third->neigh_count, triangle_intersection, triangle_count, &quad_count);

                        if (quad_count >= q - 3) {
                            counter_4++;
                        }

                        free(quad_intersection);
                    }
                }

                if (counter_4 < q - 2) {
                    if (remove_count >= capacity) {
                        capacity *= 2;
                        Node **temp = (Node **)realloc(edges_to_remove, capacity * sizeof(Node *));
                        edges_to_remove = temp;
                    }
                    edges_to_remove[remove_count] = node_first;
                    edges_to_remove[remove_count + 1] = node_second;
                    remove_count += 2;
                }

                free(triangle_intersection);
            }
        }
    }

    printf("Uklanjanje %d grana\n", remove_count / 2);
    for (int i = 0; i < remove_count; i += 2) {
        remove_edge(edges_to_remove[i], edges_to_remove[i + 1]);
    }

    free(edges_to_remove);
}



void bron_kerbosch(Clique *R, Node **P, int p_size, Node **X, int x_size, int k, ConsensusSet *unique_consensus_set, int x_cap, int l, char* alphabet, int alphabet_size) {
    /* Recursive algorithm for checking if set of nodes form clique by checking all combinations. */
    
    if (p_size == 0 && x_size == 0) {
        if (R->size == k) {
            char consensus[l + 1];
            compute_consensus_motif(R->nodes, R->size, consensus, l, alphabet, alphabet_size);

            if (!contains_consensus(unique_consensus_set, consensus)) {
                add_consensus(unique_consensus_set, consensus);
            }
        }
    }

    while (p_size > 0){
        Node *v = P[0];
        delete_first_element(P, &p_size);

        R->nodes[R->size++] = v;

        Node **new_P = (Node **)malloc(p_size * sizeof(Node *));
        Node **new_X = (Node **)malloc(x_cap * sizeof(Node *));
        int new_p_size = 0;
        int new_x_size = 0;
        for (int j = 0; j < p_size; j++) {
            if (P[j] != NULL && v->hash_value != P[j]->hash_value) {
                for (int l = 0; l < v->neigh_count; l++) {
                    if (v->neighbours[l]->hash_value == P[j]->hash_value) {
                        new_P[new_p_size++] = P[j];
                        break;
                    }
                }
            }
        }
        for (int j = 0; j < x_size; j++) {
            if (X[j] != NULL && v->hash_value != X[j]->hash_value) {
                for (int t = 0; t < v->neigh_count; t++) {
                    if (v->neighbours[t]->hash_value == X[j]->hash_value) {
                        new_X[new_x_size++] = X[j];
                        break;
                    }
                }
            }
        }

        bron_kerbosch(R, new_P, new_p_size, new_X, new_x_size, k, unique_consensus_set, x_cap, l, alphabet, alphabet_size);
        R->size--;
        X[x_size++] = v;
        free(new_P);
        free(new_X);
    }
}


void find_cliques(HashMap *hashmap, int k, int l, char* alphabet, int alphabet_size) {
    /*
        We initialize values for Bron Kerbosch algorithm and create consensus set so that one clique will be identified
        as same as the clique with the same nodes but different order.
    */
    Clique R;
    R.nodes = (Node **)malloc(MAX_SEQUENCES * sizeof(Node *));
    for (int i = 0; i < MAX_SEQUENCES; i++) {
        R.nodes[i] = (Node *)malloc(sizeof(Node));
    }
    R.size = 0;

    Node **P = (Node **)malloc(hashmap->size * sizeof(Node *));
    Node **X = (Node **)malloc(hashmap->size * sizeof(Node *));


    int index = 0;
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        Node *node = hashmap->table[i];
        if (node!=NULL){
            P[index++] = node;
        }
    }

    int p_size = index;
    int x_size = 0;
    ConsensusSet *unique_consensus_set = create_consensus_set(HASH_SET_SIZE);
    bron_kerbosch(&R, P, p_size, X, x_size, k, unique_consensus_set, p_size, l, alphabet, alphabet_size);
    for (int i=0; i<unique_consensus_set->size; i++){
        printf("Found motif: %s\n", unique_consensus_set->data[i]);
    }

    free(P);
    free(X);
    free_consensus_set(unique_consensus_set);
}

char **read_lines_from_file(const char *file_path, int *num_lines) {
    FILE *file = fopen(file_path, "r");
    if (!file) {
        perror("Failed to open file!\n");
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

    if (argc < 5) {
        fprintf(stderr, "Argumenti: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> <k-uslov odsecanja> [<datoteka_sa_azbukom>]\n");
        return 1;
    }

    int l = atoi(argv[1]);
    int d = atoi(argv[2]);
    int k = atoi(argv[4]);

    if (l<0 || d<0){
        fprintf(stderr, "Argumenti duzina motiva i dozvoljene mutacije moraju biti veci od 0.\n");
        return 1;
    }

    int num_sequences;
    char **sequences = read_lines_from_file(argv[3], &num_sequences);

    char *alphabet = (char *)malloc(25*sizeof(char));
    int alphabet_size;

    if (argc == 6) {
        alphabet = read_alphabet(argv[5], &alphabet_size);
    } else{
        alphabet = strdup("acgt");
        alphabet_size = 4;
    }

    clock_t start_creating = clock();
    HashMap *hashmap = create_hash_map();
    construct_hashmap(hashmap, sequences, num_sequences, l, d);
    clock_t end_creating = clock();

    printf("Made all conections for %d nodes.\n", hashmap->r_size);

    if (k==2){
        printf("Winnower k=2\n");
        clock_t start_w2 = clock();
        for (int i=0; i<4; i++){
            winnower_k2(hashmap, num_sequences);
        }

        find_cliques(hashmap, num_sequences, l, alphabet, alphabet_size);
        clock_t end_w2 = clock();

        printf("Vreme Winnower k=2: %lf\n", (double)(end_creating - start_creating) / CLOCKS_PER_SEC + (double)(end_w2 - start_w2) / CLOCKS_PER_SEC);

    }else if(k==3){
        printf("Winnower k=3\n");
        clock_t start_w3 = clock();
        for (int i=0; i<2; i++){
            winnower_k3(hashmap, num_sequences);
        }
        clock_t end_w3 = clock();

        clock_t start_find_clique2 = clock();
        find_cliques(hashmap, num_sequences, l, alphabet, alphabet_size);
        clock_t end_find_clique2 = clock();

        printf("Vreme Winnower k=3: %lf\n", (double)(end_creating - start_creating) / CLOCKS_PER_SEC + 
        (double)(end_w3 - start_w3) / CLOCKS_PER_SEC + (double)(end_find_clique2 - start_find_clique2) / CLOCKS_PER_SEC);
    }else if(k==4){
        printf("Winnower k=4\n");
        clock_t start_w4 = clock();
        for (int i=0; i<1; i++){
            winnower_k4(hashmap, num_sequences);
        }
        clock_t end_w4 = clock();

        clock_t start_find_clique4 = clock();
        find_cliques(hashmap, num_sequences, l, alphabet, alphabet_size);
        clock_t end_find_clique4 = clock();

        printf("Vreme Winnower k=4: %lf\n", (double)(end_creating - start_creating) / CLOCKS_PER_SEC + 
        (double)(end_w4 - start_w4) / CLOCKS_PER_SEC + (double)(end_find_clique4 - start_find_clique4) / CLOCKS_PER_SEC);
    }else{
        fprintf(stderr, "Argument k mora biti 2, 3 ili 4.\n");
        return 1;
    }
    
    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i]);
    }
    free(sequences);
    free(hashmap);
    free(alphabet);
    return 0;
}