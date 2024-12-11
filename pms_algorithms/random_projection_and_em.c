#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>


#define MAX_LINES 1000
#define MAX_BUFFER 1200
#define EM_ITER 5
#define MAX_SEQ_LEN 650
#define MAX_MOT_LEN 20
#define MAX_LINE_LENGTH 100

int num_buckets = 0;
int num_lmers = 0;

typedef struct {
    char *value;
    int count;
} Sequence;

typedef struct {
    int *positions;
    int length;
} Projection;

typedef struct {
    int size;
    char **lmers;
} Bucket;

void create_lmers(char **sequences, int num_sequences, int l, char ***lmers) {
    *lmers = (char **)malloc((MAX_SEQ_LEN) * MAX_MOT_LEN * sizeof(char *));
    int index = 0;
    for (int i = 0; i < num_sequences; i++) {
        int s_l = strlen(sequences[i]) - 1;
        for (int j = 0; j < s_l - l + 1; j++) {
            (*lmers)[index] = (char *)malloc((l + 1) * sizeof(char));
            strncpy((*lmers)[index], sequences[i] + j, l);
            (*lmers)[index][l] = '\0';
            index++;
            num_lmers++;
        }
    }
}

Projection random_projection(int l, int k) {
    Projection proj;
    proj.length = k;
    proj.positions = (int *)malloc(k * sizeof(int));
    srand(time(NULL));
    int count = 0;
    while (count < k) {
        int unique = 1;
        int pos = rand() % (l - 1);
        for (int i = 0; i < count; i++){
            if (proj.positions[i] == pos){
                unique = 0;
                break;
            }
        }
        if (unique) {
            proj.positions[count] = pos;
            count++;
        }
    }
    
    return proj;
}


Bucket *hash_lmers(char **lmers, int num_lmers, Projection proj) {
    int max_buckets = num_lmers;
    Bucket *buckets = (Bucket *)malloc(max_buckets * sizeof(Bucket));
    for (int i = 0; i < num_lmers; i++) {
        if (!lmers[i]){
            continue;
        }
        char *key = (char *)malloc((proj.length) * sizeof(char));
        for (int j = 0; j < proj.length; j++) {
            key[j] = lmers[i][proj.positions[j]];
        }
        int found = 0;
        for (int j = 0; j < num_buckets; j++) {
            bool mat = true;
            for (int k = 0; k < proj.length; k++){
                if (key[k] != buckets[j].lmers[0][proj.positions[k]]){
                    mat = false;
                    break;
                }
            }
            if (mat) {
                buckets[j].lmers = (char **)realloc(buckets[j].lmers, (buckets[j].size + 1) * sizeof(char *));
                buckets[j].lmers[buckets[j].size] = lmers[i];
                buckets[j].size++;
                found = 1;
                break;
            }
        }
        if (!found) {
            buckets[num_buckets].lmers = (char **)malloc(sizeof(char *));
            buckets[num_buckets].lmers[0] = strdup(lmers[i]);
            buckets[num_buckets].size = 1;
            num_buckets++;
        }
        free(key);
    }
    return buckets;
}

double **allocate_2d_array(int rows, int cols) {
    double **array = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        array[i] = (double *)malloc(cols * sizeof(double));
    }
    return array;
}


void initialize_pwm_from_bucket(Bucket bucket, int l, double** pwm, char* alphabet, int alphabet_size, double* background) {

    for (int i = 0; i < alphabet_size; i++) {
        for (int j = 0; j < l; j++) {
            pwm[i][j] = background[i];
        }
    }

    int counts[alphabet_size];
    for (int i = 0; i < alphabet_size; i++){
        counts[alphabet[i]] = i;
    }

    for (int i = 0; i < bucket.size; i++) {
        char *lmer = strdup(bucket.lmers[i]);
        for (int j = 0; j < l; j++) {
            pwm[counts[(unsigned char)lmer[j]]][j] += 1.0;
        }
    }

    for (int j = 0; j < l; j++) {
        double sum = 0.0;
        for (int i = 0; i < alphabet_size; i++) {
            sum += pwm[i][j];
        }
        if (sum > 0) {
            for (int i = 0; i < alphabet_size; i++) {
                pwm[i][j] /= sum;
            }
        }
    }

}

void initialize_pwm_from_motifs(char** motifs, int l, double** pwm, int num_seq, char* alphabet, int alphabet_size, double* background) {
    for (int i = 0; i < alphabet_size; i++) {
        for (int j = 0; j < l; j++) {
            pwm[i][j] = background[i];
        }
    }

    int counts[alphabet_size];
    for (int i = 0; i < alphabet_size; i++){
        counts[alphabet[i]] = i;
    }

    for (int i = 0; i < num_seq; i++) {
        char *lmer = strdup(motifs[i]);
        for (int j = 0; j < l; j++) {
            pwm[counts[(unsigned char)lmer[j]]][j] += 1.0;
        }
    }

    for (int j = 0; j < l; j++) {
        double sum = 0.0;
        for (int i = 0; i < alphabet_size; i++) {
            sum += pwm[i][j];
        }
        if (sum > 0) {
            for (int i = 0; i < alphabet_size; i++) {
                pwm[i][j] /= sum;
            }
        }
    }
}


void free_2d_array(double **array, int rows) {
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}
void normalize(double *array, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += array[i];
    }
    for (int i = 0; i < size; i++) {
        array[i] /= sum;
    }
}

void e_step(char **sequences, int num_seqs, int num_lmers, double **pwm, int motif_length, double* lr, char* alphabet, int alphabet_size, double* background) {

    for (int i = 0; i < alphabet_size; i++){
        for (int m=0; m < motif_length; m++){
            double w_alpha = 0.0;
            double w_all = 0.0;
            for (int j = 0; j < num_seqs; j++) {
                for (int k = 0; k < strlen(sequences[j]) - motif_length; k++) {
                    if ((char)sequences[j][k + m] == alphabet[i]) {
                        w_alpha += lr[MAX_SEQ_LEN*j+k];
                    }
                    w_all += lr[MAX_SEQ_LEN*j+k];
                }
            }
            pwm[i][m] = w_alpha/w_all;
        }
    }
}

void m_step(char **sequences, int num_seqs, int num_lmers, double **pwm, int motif_length, double *lr, char* alphabet, int alphabet_size, double* background) {
    for (int i = 0; i < num_seqs; i++) {
        for (int j = 0; j < strlen(sequences[i]) - motif_length; j++) {
            double likelihood_ratio_w = 1.0;
            double likelihood_ratio_b = 1.0;
            int p;
            for (int m = 0; m < motif_length; m++){
                for (int r = 0; r < alphabet_size; r++){
                    if (sequences[i][j+m] == alphabet[r]){
                        p=r;
                        break;
                    }
                }
                likelihood_ratio_w *= pwm[p][m];
                likelihood_ratio_b *= background[p];
            }

            lr[MAX_SEQ_LEN*i + j] = likelihood_ratio_w/likelihood_ratio_b;
        }
    }
}

void em_algorithm(char **sequences, int num_seqs, int num_lmers, double **pwm, int motif_length, int max_iter, double tol, char* alphabet, int alphabet_size, double* background) {
    double *lr = (double *)malloc(MAX_MOT_LEN * MAX_SEQ_LEN *sizeof(double));

    for (int iteration = 0; iteration < max_iter; iteration++) {
        double **old_pwm = allocate_2d_array(alphabet_size, motif_length);
        for (int i = 0; i < alphabet_size; i++) {
            for (int j = 0; j < motif_length; j++) {
                old_pwm[i][j] = pwm[i][j];
            }
        }

        if (iteration != 0){
            e_step(sequences, num_seqs, num_lmers, pwm, motif_length, lr, alphabet, alphabet_size, background);
            m_step(sequences, num_seqs, num_lmers, pwm, motif_length, lr, alphabet, alphabet_size, background);
        }else{
            m_step(sequences, num_seqs, num_lmers, pwm, motif_length, lr, alphabet, alphabet_size, background);
        }

        if (iteration != 0){
            double norm_diff = 0.0;
            for (int i = 0; i < alphabet_size; i++) {
                for (int j = 0; j < motif_length; j++) {
                    norm_diff += (pwm[i][j] - old_pwm[i][j]) * (pwm[i][j] - old_pwm[i][j]);
                }
            }
            norm_diff = sqrt(norm_diff);

            free_2d_array(old_pwm, alphabet_size);
            if (norm_diff < tol) {
                break;
            }
        }
    }
    free(lr);
}

void compute_consensus_motif(char **motifs, int num_seqs, int l, char* consensus, char* alphabet, int alphabet_size, double* background) {
    int (*counts)[alphabet_size] = malloc(l * sizeof(*counts));
    for(int i = 0; i < l; i++) {
        for(int j = 0; j < alphabet_size; j++) {
            counts[i][j] = 0;
        }
    }
    for (int i = 0; i < num_seqs; i++) {
        for (int j = 0; j < l; j++) {
            for (int r = 0; r < alphabet_size; r++){
                if (motifs[i][j] == alphabet[r]){
                    counts[j][r]++;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < l; i++) {
        int max_count = 0;
        char consensus_base = alphabet[0];
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


void find_best_motifs(char **sequences, int num_sequences, double **pwm, char **best_motifs, int motif_length, char* alphabet, int alphabet_size){
    char *current_best_motif = (char *)malloc((motif_length+1) * sizeof(char));
    double best_prob = 0.0;
    for (int i = 0; i < num_sequences; i++){
        for (int j = 0; j < strlen(sequences[i]) - motif_length + 1; j++) {
            double current_prob = 1.0;
            for (int k = 0; k < motif_length; k++) {
                int p;
                for (int r = 0; r < alphabet_size; r++){
                    if (sequences[i][j+k] == alphabet[r]){
                        p=r;
                        break;
                    }
                }
                current_prob *= pwm[p][k];
            }
            if (best_prob < current_prob){
                best_prob = current_prob;
                strncpy(current_best_motif, sequences[i] + j, motif_length);
                current_best_motif[motif_length] = '\0';
            }
        }
        strncpy(best_motifs[i], current_best_motif, motif_length);
        best_motifs[i][motif_length] = '\0';
    }
    free(current_best_motif);

}


double calculate_likelihood_ratio(char **best_motifs, double **S_pwm, int num_motifs, int motif_length, char* alphabet, int alphabet_size, double* background){
    double likelihood_ratio_Sb = 1.0;
    for (int i = 0; i < num_motifs; i++) {
        double likelihood_ratio_w = 1.0;
        double likelihood_ratio_b = 1.0;
        int p;
        for (int m=0; m < motif_length; m++){
            for (int r = 0; r < alphabet_size; r++){
                if (best_motifs[i][m] == alphabet[r]){
                    p=r;
                    break;
                }
            }
            likelihood_ratio_w *= S_pwm[p][m];
            likelihood_ratio_b *= background[p];
        }
        likelihood_ratio_Sb *= likelihood_ratio_w/likelihood_ratio_b;
    }
    return likelihood_ratio_Sb;
}

void projection_algorithm(char **sequences, int num_sequences, int l, int k, int s, int max_trials, char* alphabet, int alphabet_size) {
    char **lmers;
    // Generate all lmers from all input sequences
    create_lmers(sequences, num_sequences, l, &lmers);

    double best_likelihood_ratio = 0.0;
    char* consensus = (char *)malloc((l+1)*sizeof(char));
    double background[MAX_LINE_LENGTH];
    for (int i = 0; i < alphabet_size; i++){
        background[i] = 1.0/alphabet_size;
    }
    for (int trial = 0; trial < max_trials; trial++) {
        // Pick random k indexes, that will be our projection
        Projection proj = random_projection(l, k);
        // Make buckets by the projection and group all lmers in them
        Bucket *buckets = hash_lmers(lmers, num_lmers, proj);

        // For each bucket that has more than s lmers in it do EM algorithm
        for (int i = 0; i < num_buckets; i++) {
            if (buckets[i].size >= s) {
                double **pwm = allocate_2d_array(alphabet_size, l);
                // Make PWM matrix using lmers from the bucket
                initialize_pwm_from_bucket(buckets[i], l, pwm, alphabet, alphabet_size, background);
                // Do E and M steps EM_ITER times 
                em_algorithm(sequences, num_sequences, num_lmers, pwm, l, EM_ITER, 1e-4, alphabet, alphabet_size, background);

                // From last PWM matrix EM algorithm made find best kmers from each sequence
                char **best_motifs = (char **)malloc(num_sequences * sizeof(char *));
                for (int j = 0; j < num_sequences; j++){
                    best_motifs[j] = (char *)malloc(MAX_MOT_LEN * sizeof(char));
                }
                find_best_motifs(sequences, num_sequences, pwm, best_motifs, l, alphabet, alphabet_size);

                // Make new PWM matrix from best kmers
                double **S_pwm = allocate_2d_array(alphabet_size, l);
                initialize_pwm_from_motifs(best_motifs, l, S_pwm, num_sequences, alphabet, alphabet_size, background);

                // Calculate likelihood for this best motifs and if its best from all buckets save consensus kmer as result
                double likelihood_ratio = calculate_likelihood_ratio(best_motifs, S_pwm, num_sequences, l, alphabet, alphabet_size, background);
                if (likelihood_ratio > best_likelihood_ratio){
                    best_likelihood_ratio = likelihood_ratio;
                    compute_consensus_motif(best_motifs, num_sequences, l, consensus, alphabet, alphabet_size, background);
                }
                free_2d_array(pwm, alphabet_size);
                for (int j = 0; j < num_sequences; j++){
                    free(best_motifs[j]);
                }
                free(best_motifs);
            }
        }
        num_buckets = 0;
    }
    
    printf("Pronadjen motiv: %s\n", consensus);

    for (int i = 0; i < num_lmers; i++) {
        free(lmers[i]);
    }
}

char **read_lines_from_file(const char *file_path, int *num_lines) {
    FILE *file = fopen(file_path, "r");
    if (!file) {
        perror("Failed to open file");
        exit(1);
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

    if (argc < 6) {
        fprintf(stderr, "Argumenti: <duzina_motiva> <broj_proj> <s-filtriranje_korpi> <iteracije> <datoteka_sa_sekvencama> [<datoteka_sa_azbukom>]\n");
        return 1;
    }

    int l = atoi(argv[1]);
    int k = atoi(argv[2]);
    int s = atoi(argv[3]);
    int iter = atoi(argv[4]);


    if (l<0){
        fprintf(stderr, "Argumenti duzina motiva i dozvoljene mutacije moraju biti veci od 0.\n");
        return 1;
    }

    int num_sequences;
    char **sequences = read_lines_from_file(argv[5], &num_sequences);

        char *alphabet = (char *)malloc(25*sizeof(char));
    int alphabet_size;

    if (argc == 7) {
        alphabet = read_alphabet(argv[6], &alphabet_size);
    } else{
        alphabet = strdup("acgt");
        alphabet_size = 4;
    }

    clock_t start = clock();
    projection_algorithm(sequences, num_sequences, l, k, s, iter, alphabet, alphabet_size);
    clock_t end = clock();

    printf("TIME: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
    free(alphabet);
    free(sequences);
    return 0;
}
