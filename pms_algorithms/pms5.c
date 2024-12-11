#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define MAX_LEN 20
#define MAX_VISITED 100000
#define MAX_Q_SIZE 1000000
#define Q_TRESHOLD 8000
#define MAX_LINES 1000
#define MAX_BUFFER 1200
#define SET_SIZE 100
#define N 5
#define MAX_LINE_LENGTH 100


typedef struct Node {
    int position;
    char value[MAX_LEN];
    struct Node** children;
    int num_children;
    int depth;
} Node;

typedef struct Queue {
    Node** items;
    int front;
    int rear;
    int size;
    int capacity;
} Queue;

typedef struct Set {
    char** items;
    int size;
    int capacity;
} Set;

Set* create_set(int capacity) {
    Set* set = (Set*)malloc(sizeof(Set));
    set->items = (char**)malloc(capacity * sizeof(char*));
    set->size = 0;
    set->capacity = capacity;
    return set;
}

void add_to_set(Set* set, const char* item) {
    if (set->size == set->capacity) {
        set->capacity *= 2;
        set->items = (char**)realloc(set->items, set->capacity * sizeof(char*));
    }
    set->items[set->size] = strdup(item);
    set->size++;
}

bool set_contains(Set* set, const char* item) {
    for (int i = 0; i < set->size; i++) {
        if (strcmp(set->items[i], item) == 0) {
            return true;
        }
    }
    return false;
}

Set* set_union(Set* set1, Set* set2) {
    Set* result = create_set(set1->size + set2->size);
    for (int i = 0; i < set1->size; i++) {
        add_to_set(result, set1->items[i]);
    }
    for (int i = 0; i < set2->size; i++) {
        if (!set_contains(result, set2->items[i])) {
            add_to_set(result, set2->items[i]);
        }
    }
    return result;
}

Set* set_intersection(Set* set1, Set* set2) {
    Set* result = create_set(set1->size);
    for (int i = 0; i < set1->size; i++) {
        if (set_contains(set2, set1->items[i])) {
            add_to_set(result, set1->items[i]);
        }
    }
    return result;
}

void free_set(Set* set) {
    for (int i = 0; i < set->size; i++) {
        free(set->items[i]);
    }
    free(set->items);
    free(set);
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

Queue* createQueue(int capacity) {
    Queue* queue = (Queue*)malloc(sizeof(Queue));
    queue->capacity = capacity;
    queue->front = queue->size = 0;
    queue->rear = capacity - 1;
    queue->items = (Node**)malloc(queue->capacity * sizeof(Node*));
    return queue;
}

bool isQueueFull(Queue* queue) {
    return (queue->size == queue->capacity);
}

bool isQueueEmpty(Queue* queue) {
    return (queue->size == 0);
}

void enqueue(Queue* queue, Node* item) {
    if (isQueueFull(queue)) {
        return;
    }
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->items[queue->rear] = item;
    queue->size = queue->size + 1;
}

Node* dequeue(Queue* queue) {
    if (isQueueEmpty(queue)) {
        return NULL;
    }
    Node* item = queue->items[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size = queue->size - 1;
    return item;
}
Node* createNode(char* value, int position, int depth) {
    Node* node = (Node*)malloc(sizeof(Node));
    strcpy(node->value, value);
    node->position = position;
    node->depth = depth;
    node->num_children = 0;
    node->children = NULL;
    return node;
}

void addChild(Node* parent, Node* child) {
    parent->num_children++;
    parent->children = (Node**)realloc(parent->children, parent->num_children * sizeof(Node*));
    parent->children[parent->num_children - 1] = child;
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

void calculate_n_values(char *x, char *y, char *z, int l, int n_values[l][N]) {
    /* Counting types of differences between three k-mers */
    for (int p = 0; p < l; p++) {
        int n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0;
        for (int i = p; i < l; i++) {
            if (x[i] == y[i] && y[i] == z[i]) {
                n1++;
            } else if (x[i] == y[i] && y[i] != z[i]) {
                n2++;
            } else if (x[i] == z[i] && z[i] != y[i]) {
                n3++;
            } else if (y[i] == z[i] && z[i] != x[i]) {
                n4++;
            } else {
                n5++;
            }
        }
        n_values[p][0] = n1;
        n_values[p][1] = n2;
        n_values[p][2] = n3;
        n_values[p][3] = n4;
        n_values[p][4] = n5;
    }
}


void freeTree(Node* node) {
    for (int i = 0; i < node->num_children; i++) {
        freeTree(node->children[i]);
    }
    free(node->children);
    free(node);
}


bool is_visited(const char* value, char visited[MAX_VISITED][MAX_LEN]) {
    for (int i = 0; i < MAX_VISITED; i++) {
        if (strcmp(visited[i], value) == 0) {
            return true;
        }
    }
    return false;
}

void mark_visited(const char* value, char visited[MAX_VISITED][MAX_LEN]) {
    for (int i = 0; i < MAX_VISITED; i++) {
        if (visited[i][0] == '\0') {
            strncpy(visited[i], value, MAX_LEN - 1);
            visited[i][MAX_LEN - 1] = '\0';
            break;
        }
    }
}


void generate_tree_bfs(Node* root, int max_depth, char* alphabet, int alphabet_size) {
    /* 
        Storing nodes in queue, iterate through nucleotides and create new child node,
        that will differ from parent in position that is the same as depth that child is at.
     */
    Queue* queue = createQueue(MAX_Q_SIZE);
    enqueue(queue, root);

    while (!isQueueEmpty(queue)) {
        Node* current_node = dequeue(queue);
        if (current_node->depth == max_depth) {
            continue;
        }
        int n = strlen(current_node->value);
        char *x = (char *) malloc(MAX_LEN * sizeof(char));
        strcpy(x, current_node->value);

        for (int p1 = current_node->position + 1; p1 <= n; p1++) {
            for (int s = 0; s < alphabet_size; s++) {
                if (x[p1 - 1] != alphabet[s]) {
                    char *new_x = (char *) malloc(MAX_LEN * sizeof(char));
                    strcpy(new_x, x);
                    new_x[p1 - 1] = alphabet[s];

                    Node* new_node = createNode(new_x, p1, current_node->depth + 1);
                    addChild(current_node, new_node);
                    enqueue(queue, new_node);
                    free(new_x);
                }
            }
        }
        free(x);
    }

    while (!isQueueEmpty(queue)) { //
        Node* node_to_free = dequeue(queue);//
        freeTree(node_to_free);//
    }//

    free(queue->items);
    free(queue);
}


Node* make_T_neigh(char* x, int d, char* alphabet, int alphabet_size) {
    /* Generate tree with all neighbours for k-mer x. */
    Node* root = createNode(x, 0, 0);
    generate_tree_bfs(root, d, alphabet, alphabet_size);
    return root;
}

int distance(char* kmer, char** sequences, int num_sequences, int j) {
    int n = strlen(kmer);
    int maxim = 0;

    for (int s_idx = j; s_idx < num_sequences; s_idx++) {
        char* s = sequences[s_idx];
        int minim = n;

        for (int i = 0; i <= strlen(s) - n; i++) {
            char *substr = slice(s, i, n);
            int dis = hamming_distance(kmer, substr);
            if (dis < minim) {
                minim = dis;
            }
        }

        if (minim > maxim) {
            maxim = minim;
        }
    }

    return maxim;
}


Set* fullprune(char* x, char* y, char* z, int d, Node* root, bool******** ilp_table) {
    /*
        Start DFS for the given tree, check if each node differes from nodes x, y and z on less than d positions each.
        If so the node is common neighbour and it is saved.
        Then we check for value from table where all linear combinations are saved if node has potential to have a 
        descendant that can be neighbour of all three nodes, if the answer is no, we perform backtracking.
    */
    int l = strlen(x);
    Set* q = create_set(100);
    int d_xy = hamming_distance(x, y);
    int d_xz = hamming_distance(x, z);

    int n_values[MAX_LEN][N];
    calculate_n_values(x, y, z, l, n_values);

    Queue* stack = createQueue(1000);
    enqueue(stack, root);

    while (!isQueueEmpty(stack)) {
        Node* node = dequeue(stack);
        int p = node->position;
        char* t = node->value;

        int d_xt = hamming_distance(x, t);
        int d_yt = hamming_distance(y, t);
        int d_zt = hamming_distance(z, t);

        if (d_xt <= d && d_yt <= d && d_zt <= d) {
            add_to_set(q, t);
        }

        char* x1 = (char *) malloc(MAX_LEN * sizeof(char));
        char* y1 = (char *) malloc(MAX_LEN * sizeof(char));
        char* z1 = (char *) malloc(MAX_LEN * sizeof(char));
        char* t1 = (char *) malloc(MAX_LEN * sizeof(char));

        strncpy(x1, x, p);
        x1[p] = '\0';
        strncpy(y1, y, p);
        y1[p] = '\0';
        strncpy(z1, z, p);
        z1[p] = '\0';
        strncpy(t1, t, p);
        t1[p] = '\0';

        int d_x1t1 = hamming_distance(x1, t1);
        int d_y1t1 = hamming_distance(y1, t1);
        int d_z1t1 = hamming_distance(z1, t1);

        if (d_y1t1 > d || d_z1t1 > d) {
            free(x1); //
            free(y1);//
            free(z1);//
            free(t1);//
            continue;
        }

        int p1 = d - d_x1t1;
        int p2 = d - d_y1t1;
        int p3 = d - d_z1t1;

        if (p < strlen(t)) {
            int n1 = n_values[p][0];
            int n2 = n_values[p][1];
            int n3 = n_values[p][2];
            int n4 = n_values[p][3];
            int n5 = n_values[p][4];
            bool ilp_t = ilp_table[n1][n2][n3][n4][n5][p1][p2][p3];
            if (!ilp_t) {
                free(x1);//
                free(y1);//
                free(z1);//
                free(t1);//
                continue;//
            } else {
                for (int i = 0; i < node->num_children; i++) {
                    enqueue(stack, node->children[i]);
                }
            }
        }
        free(x1);
        free(y1);
        free(z1);
        free(t1);
    }

    while (!isQueueEmpty(stack)) {//
        Node* node_to_free = dequeue(stack); //
        freeTree(node_to_free); //
    }//

    free(stack->items);
    free(stack);
    return q;
}



bool******** allocate_8d_array(int l, int d) {
    bool******** array = (bool********)malloc((l + 1) * sizeof(bool*******));
    for (int i = 0; i <= l; i++) {
        array[i] = (bool*******)malloc((l + 1) * sizeof(bool******));
        for (int j = 0; j <= l; j++) {
            array[i][j] = (bool******)malloc((l + 1) * sizeof(bool*****));
            for (int k = 0; k <= l; k++) {
                array[i][j][k] = (bool*****)malloc((l + 1) * sizeof(bool****));
                for (int m = 0; m <= l; m++) {
                    array[i][j][k][m] = (bool****)malloc((l + 1) * sizeof(bool***));
                    for (int n = 0; n <= l; n++) {
                        array[i][j][k][m][n] = (bool***)malloc((d + 1) * sizeof(bool**));
                        for (int o = 0; o <= d; o++) {
                            array[i][j][k][m][n][o] = (bool**)malloc((d + 1) * sizeof(bool*));
                            for (int p = 0; p <= d; p++) {
                                array[i][j][k][m][n][o][p] = (bool*)malloc((d + 1) * sizeof(bool));
                            }
                        }
                    }
                }
            }
        }
    }
    return array;
}

void free_8d_array(bool******** array, int l, int d) {
    for (int i = 0; i <= l; i++) {
        for (int j = 0; j <= l; j++) {
            for (int k = 0; k <= l; k++) {
                for (int m = 0; m <= l; m++) {
                    for (int n = 0; n <= l; n++) {
                        for (int o = 0; o <= d; o++) {
                            for (int p = 0; p <= d; p++) {
                                free(array[i][j][k][m][n][o][p]);
                            }
                            free(array[i][j][k][m][n][o]);
                        }
                        free(array[i][j][k][m][n]);
                    }
                    free(array[i][j][k][m]);
                }
                free(array[i][j][k]);
            }
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}


bool******** read_table_from_file(int l, int d, char* filename) {
    /*
        Read previously created table with values for all combinations of ns and ps.
        Table has 8 int numbers and one bool -> 3 2 0 0 1 2 2 1 true
    */
    bool******** table = allocate_8d_array(l, d);
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file");
        return NULL;
    }

    int n1, n2, n3, n4, n5, p1, p2, p3;
    char value[6];
    while (fscanf(file, "%d %d %d %d %d %d %d %d %s\n", &n1, &n2, &n3, &n4, &n5, &p1, &p2, &p3, value) == 9) {
        table[n1][n2][n3][n4][n5][p1][p2][p3] = (strcmp(value, "true") == 0);
    }

    fclose(file);
    return table;
}

Set* pms5(char** sequences, int num_sequences, int l, int d, char* filename, char* alphabet, int alphabet_size) {
    /*
        Iterating through k-mers in first sequences and k-mers from all pairs of rest of the sequences.
        For each pair of three k-mers use fullprune algorithm to find common neighbours.
        Store union of neighbours from all k-mers from paired sequences and fixed k-mer from first one in set Q.
        Check if the set is smaller than trashold, if so check if there is real motif, if not do the intersection of set Q
        with set form previous pair of sequences. This is done in iterations for all k-mers in first sequences.
     */
    Set* res_mot = create_set(SET_SIZE);
    if (num_sequences % 2 == 0) {
        sequences[num_sequences] = sequences[num_sequences - 1];
        num_sequences++;
    }
    int p = (num_sequences - 1) / 2;
    char* s1 = sequences[0];
    Set* q1 = create_set(SET_SIZE);
    bool******** ilp_table = read_table_from_file(l, d, filename);

    for (int i = 0; i < strlen(s1) - l + 1; i++) {
        char *x = (char *)malloc(MAX_LEN*sizeof(char));
        strncpy(x, s1 + i, l);
        x[l] = '\0';
        Node* root = make_T_neigh(x, d, alphabet, alphabet_size);
        int k = 0;
        for (k = 0; k < p; k++) {
            Set* q = create_set(100);
            char* s2 = sequences[2 * k + 1];
            char* s3 = sequences[2 * k + 2];
            for (int j = 0; j < strlen(s2) - l + 1; j++) {
                char *y = (char *)malloc(MAX_LEN*sizeof(char));
                strncpy(y, s2 + j, l);
                y[l] = '\0';
                for (int r = 0; r < strlen(s3) - l + 1; r++) {
                    char *z = (char *)malloc(MAX_LEN*sizeof(char));
                    strncpy(z, s3 + r, l);
                    z[l] = '\0';
                    Set* neighbors_of_all = fullprune(x, y, z, d, root, ilp_table);
                    Set* temp_q = set_union(q, neighbors_of_all);
                    free_set(q);
                    q = temp_q;
                    free_set(neighbors_of_all);
                    free(z);
                }
                free(y);
            }
            if (k == 0) {
                free_set(q1);
                q1 = q;
            } else {
                Set* temp_q1 = set_intersection(q1, q);
                free_set(q1);
                q1 = temp_q1;
                free_set(q);
            }
            if (q1->size < Q_TRESHOLD) {
                break;
            }
        }
        for (int s = 0; s < q1->size; s++) {
            if (distance(q1->items[s], sequences, num_sequences, k*2+3) <= d) {
                add_to_set(res_mot, q1->items[s]);
                printf("Pronadjen motiv: %s\n", q1->items[s]);
            }
        }
        free(x);
        freeTree(root);
    }
    free_set(q1);
    free_8d_array(ilp_table, l, d);
    return res_mot;
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

    if (argc < 5) {
        fprintf(stderr, "Argumenti: <duzina_motiva> <broj_mutacija> <datoteka_sa_sekvencama> <datoteka_ilr_tabela> [<datoteka_sa_azbukom>]\n");
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

    if (argc == 6) {
        alphabet = read_alphabet(argv[5], &alphabet_size);
    } else{
        alphabet = strdup("acgt");
        alphabet_size = 4;
    }

    clock_t start = clock();
    pms5(sequences, num_sequences, l, d, argv[4], alphabet, alphabet_size);
    clock_t end = clock();

    printf("Vreme: %lfs\n", (double)(end - start) / CLOCKS_PER_SEC);

    for (int i = 0; i < num_sequences; i++) {
        free(sequences[i]);
    }
    free(sequences);
    free(alphabet);
    return 0;
}