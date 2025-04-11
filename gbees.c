// gbees.c, https://github.com/bhanson10/gbees
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <lapacke.h>

#include "gbees.h"

#define TOL 1E-8
#define FNV_OFFSET 14695981039346656037ULL
#define FNV_PRIME 1099511628211ULL

void exit_nomem(const char* error_string) {
    fprintf(stderr, "Error: memory allocation during %s.\n", error_string);
    exit(EXIT_FAILURE);
}

Meas Meas_create(int dim, const char* M_DIR, const char* M_FILE) {
    char M_PATH[256];
    snprintf(M_PATH, sizeof(M_PATH), "%s/%s", M_DIR, M_FILE);

    FILE* m_file = fopen(M_PATH, "r");
    if (m_file == NULL) {
        fprintf(stderr, "Error: could not open file %s", M_PATH);
        exit(EXIT_FAILURE);
    }

    Meas M;
    M.dim = dim; 
    M.mean = malloc(dim * sizeof(double));
    M.cov = malloc(dim * sizeof(double*));
    if (M.mean == NULL || M.cov == NULL) {
        const char* error_string = "measurement creation"; 
        exit_nomem(error_string); 
    }
    for (int i = 0; i < dim; i++) {
        M.cov[i] = malloc(dim * sizeof(double));
        if (M.cov[i] == NULL) {
            const char* error_string = "measurement creation"; 
            exit_nomem(error_string); 
        }
    }

    char line[256];
    char* token;
    int count = 0;

    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // mean vector
    token = strtok(line, " ");
    while (token != NULL && count < dim) { // read mean vector
        M.mean[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0;

    // Read covariance matrix
    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    for (int i = 0; i < dim; i++) { // read covariance matrix
        fgets(line, sizeof(line), m_file);
        token = strtok(line, " ");
        while (token != NULL && count < dim) {
            M.cov[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
        count = 0;
    }

    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // read T value
    M.T = strtod(line, NULL);

    fclose(m_file);
    return M;
}

void Meas_free(Meas* M) {
    if (M->mean != NULL) {
        free(M->mean);
    }
    if (M->cov != NULL) {
        for (int i = 0; i < M->dim; i++) {
            if (M->cov[i] != NULL) {
                free(M->cov[i]);
            }
        }
        free(M->cov);
    }
    return; 
}

Grid Grid_create(int dim, double t0, double thresh, Meas M, double* factor, bool ROTATE){
    Grid G; 
    G.dim = dim; 
    G.thresh = thresh; 
    G.t = t0; 
    G.dt = DBL_MAX; 
    G.center = malloc(dim * sizeof(double));
    G.dx = malloc(dim * sizeof(double));
    G.factor = malloc(dim * sizeof(double));
    if (G.center == NULL || G.dx == NULL || G.factor == NULL){
        const char* error_string = "grid creation"; 
        exit_nomem(error_string); 
    }
    for(int i = 0; i < dim; i++){
        G.center[i] = M.mean[i]; 
        G.factor[i] = factor[i]; 
    }
    G.hi_bound = DBL_MAX; 
    G.lo_bound = -DBL_MAX; 

    double A[dim * dim]; int count = 0; 
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            A[count] = M.cov[i][j]; 
            count++;
        }
    }


    G.R = (double**)malloc(dim * sizeof(double*));
    G.Rt = (double**)malloc(dim * sizeof(double*));
    if((G.R == NULL)||(G.Rt == NULL)){
        const char* error_string = "grid creation"; 
        exit_nomem(error_string);
    }
    for (int i = 0; i < dim; i++) {
        G.R[i] = (double*)malloc(dim * sizeof(double));
        G.Rt[i] = (double*)malloc(dim * sizeof(double));
        if((G.R[i] == NULL)||(G.Rt[i] == NULL)){
            const char* error_string = "grid creation"; 
            exit_nomem(error_string);
        }
    }

    if(ROTATE){
        double eigvals[dim]; // Array to hold eigenvalues
        int lda = dim;       // Leading dimension of the matrix
        int info;            // LAPACK info code
        info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', dim, A, lda, eigvals);
        if (info != 0) {
            perror("LAPACK failed to compute eigenvalues.\n");
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                G.R[i][j]  = A[j * dim + i];
                G.Rt[j][i] = A[j * dim + i];
            }
        }

        // Calculating dx based on new rotated grid
        double RtCov[G.dim][G.dim];
        for (int i = 0; i < G.dim; i++) {
            for (int j = 0; j < G.dim; j++) {
                RtCov[i][j] = 0;
                for (int k = 0; k < G.dim; k++) {
                    RtCov[i][j] += G.Rt[i][k] * M.cov[k][j];
                }
            }
        }
        double RtCovR[G.dim][G.dim];
        for (int i = 0; i < G.dim; i++) {
            for (int j = 0; j < G.dim; j++) {
                RtCovR[i][j] = 0;
                for (int k = 0; k < G.dim; k++) {
                    RtCovR[i][j] += RtCov[i][k] * G.R[k][j];
                }
            }
        }
        for(int i = 0; i < G.dim; i++){
            G.dx[i] = sqrt(RtCovR[i][i]) / (2.0 * G.factor[i]); 
        }
    }else{
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++){
                if(i == j){
                    G.R[i][j]  = 1;
                    G.Rt[i][j] = 1;
                }else{
                    G.R[i][j]  = 0;
                    G.Rt[i][j] = 0;
                }
            }
        }

        // Calculating dx based on original measurement covariance grid
        for(int i = 0; i < G.dim; i++){
            G.dx[i] = sqrt(M.cov[i][i]) / (2.0 * G.factor[i]); 
        }
    }

    return G;
}

void Grid_free(Grid* G) {
    // Free the allocated memory
    free(G->center);
    free(G->dx);
    free(G->factor);
    for (int i = 0; i < G->dim; i++){
        free(G->R[i]); free(G->Rt[i]);
    }
    free(G->R); free(G->Rt); 
            
    // Set pointers to NULL to avoid dangling pointers
    G->center = NULL;
    G->dx = NULL;
    G->factor = NULL;
    G->R = NULL;
    G->Rt = NULL;
    return; 
}

Traj Traj_create(int n, double* coef){
    Traj T; 
    T.coef = malloc(n * sizeof(double));
    if (T.coef == NULL) {
        const char* error_string = "trajectory creation"; 
        exit_nomem(error_string); 
    }
    for(int i = 0; i < n; i++){
        T.coef[i] = coef[i]; 
    }
    return T;
}

void Traj_free(Traj* T) {
    if (T->coef != NULL) {
        free(T->coef);
        T->coef = NULL;
    }
    return;
}

HashTableEntry* HashTableEntry_create(int dim, int* state, double prob, double J){
    // Allocate space for hash table struct.
    HashTableEntry* entry = (HashTableEntry*)malloc(sizeof(HashTableEntry));
    if (entry == NULL) {
        const char* error_string = "hash table entry creation"; 
        exit_nomem(error_string);
    }

    entry->state = malloc(dim * sizeof(int));
    entry->prob = prob;
    entry->v = (double*)malloc(dim * sizeof(double));
    entry->ctu = (double*)malloc(dim * sizeof(double));
    entry->i_nodes = (HashTableEntry**)malloc(dim * sizeof(HashTableEntry*));
    entry->k_nodes = (HashTableEntry**)malloc(dim * sizeof(HashTableEntry*));
    entry->dcu = 0; 
    entry->cfl_dt = 0; 
    entry->new_f = 0; 
    entry->ik_f = 0; 
    entry->bound_val = J; 

    // Check for memory allocation failure
    if(entry->state == NULL || entry->v == NULL || entry->ctu == NULL || entry->state == NULL || entry->i_nodes == NULL || entry->k_nodes == NULL){
        if (entry->state) free(entry->state);
        if (entry->v) free(entry->v);
        if (entry->ctu) free(entry->ctu);
        if (entry->i_nodes) free(entry->i_nodes);
        if (entry->k_nodes) free(entry->k_nodes);
        const char* error_string = "hash table entry attributes creation"; 
        exit_nomem(error_string);
    }

    for (int i = 0; i < dim; i++) {
        entry->state[i] = state[i];
        entry->v[i] = 0.0;
        entry->ctu[i] = 0.0;
        entry->i_nodes[i] = NULL;
        entry->k_nodes[i] = NULL;
    }
    entry->next = NULL; 
    return entry;
}

HashTable* HashTable_create(int CAPACITY) {
    // Allocate space for hash table struct.
    HashTable* P = malloc(sizeof(HashTable));
    if (P == NULL) {
        const char* error_string = "hash table creation"; 
        exit_nomem(error_string); 
    }
    P->a_count = 0;
    P->tot_count = 0;
    P->capacity = CAPACITY;

    P->entries = (HashTableEntry**)malloc(P->capacity * sizeof(HashTableEntry*));
    if (P->entries == NULL) {
        free(P); // error, free table before we return!
        const char* error_string = "hash table entries creation"; 
        exit_nomem(error_string); 
    }
    for(int i = 0; i < P->capacity; i++){
        P->entries[i] = NULL; 
    }
    return P;
}

void HashTable_free(HashTable* P) {
    // First free allocated keys.
    for (size_t i = 0; i < P->capacity; i++) {
        HashTableEntry* temp; 
        HashTableEntry* head = P->entries[i]; 
        while(head != NULL){
            temp = head; 
            head = head->next; 
            free(temp->state);
            free(temp->v);
            free(temp->ctu);
            free(temp->i_nodes);
            free(temp->k_nodes);
            free(temp); 
        }
    }

    // Then free entries array and table itself.
    free(P->entries);
    free(P);
    return; 
}

/*==============================================================================
                        NON-MEMBER FUNCTION DEFINITIONS
==============================================================================*/
// Return 64-bit FNV-1a hash for state (NUL-terminated). See description:
// https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function
uint64_t FNV1a(int* state, int dim) {
    uint64_t hash = FNV_OFFSET;
    for (int p = 0; p < dim; p++) {
        hash ^= (uint64_t)(state[p]);
        hash *= FNV_PRIME;
    }
    return hash;
}

#define ROTATE(x, k) (((x) << (k)) | ((x) >> (64 - (k))))

// Mix function as defined in Bob Jenkins' lookup3 hash
#define MIX(a, b, c) \
{ \
  a -= c;  a ^= ROTATE(c, 43);  c += b; \
  b -= a;  b ^= ROTATE(a, 39);  a += c; \
  c -= b;  c ^= ROTATE(b, 29);  b += a; \
  a -= c;  a ^= ROTATE(c, 34);  c += b; \
  b -= a;  b ^= ROTATE(a, 28);  a += c; \
  c -= b;  c ^= ROTATE(b, 20);  b += a; \
}

// Final mix to ensure all bits are affected
#define FINAL(a, b, c) \
{ \
  c ^= b; c -= ROTATE(b, 14); \
  a ^= c; a -= ROTATE(c, 11); \
  b ^= a; b -= ROTATE(a, 25); \
  c ^= b; c -= ROTATE(b, 16); \
  a ^= c; a -= ROTATE(c, 4);  \
  b ^= a; b -= ROTATE(a, 14); \
  c ^= b; c -= ROTATE(b, 24); \
}

uint64_t lookup3(const int* state, size_t dim) {
    uint64_t a, b, c;
    a = b = c = 0xdeadbeefdeadbeefLL + (dim << 3);

    while (dim > 2) {
        a += (uint64_t)state[0];
        b += (uint64_t)state[1];
        c += (uint64_t)state[2];
        MIX(a, b, c);
        state += 3;  // Move to the next 3 integers
        dim -= 3;
    }

    // Handle the remaining elements
    switch (dim) {
        case 2:
            b += (uint64_t)state[1];
            a += (uint64_t)state[0];
            break;
        case 1:
            a += (uint64_t)state[0];
            break;
        case 0:  // Nothing to add
            break;
    }
    
    FINAL(a, b, c);
    return c;  // Return the final hash value
}

#define get16bits(d) ((uint32_t)((d) & 0xFFFF))

uint32_t SuperFastHash(const int* state, int dim) {
    uint32_t hash = (uint32_t)dim, tmp;
    int rem;

    if (dim <= 0 || state == NULL) return 0;

    rem = dim & 1;  // Check remainder to see if there is one 32-bit int left
    dim >>= 1;      // Each loop step processes 2 * 16-bit units (one 32-bit integer)

    /* Main loop */
    for (; dim > 0; dim--) {
        hash  += get16bits(state[0]);   // Lower 16 bits of first 32-bit integer
        tmp    = (get16bits(state[0] >> 16) << 11) ^ hash;  // Upper 16 bits shifted and XORed
        hash   = (hash << 16) ^ tmp;
        state++;  // Move to the next 32-bit integer
        hash  += hash >> 11;
    }

    /* Handle end cases */
    if (rem) {  // If there's one 32-bit integer left to process
        hash += get16bits(state[0]);
        hash ^= hash << 16;
        hash ^= (int8_t)(state[0] >> 16) << 18;  // Process the higher 16 bits
        hash += hash >> 11;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

uint64_t DJBX33A(const int* state, int dim) {
    uint64_t hash = 5381;  // Initial value for the DJBX33A hash function

    for (int i = 0; i < dim; ++i) {
        // Process the 32-bit integer as 4 separate bytes
        hash = ((hash << 5) + hash) + (state[i] & 0xFF);         // Add the lowest byte
        hash = ((hash << 5) + hash) + ((state[i] >> 8) & 0xFF);  // Add the next byte
        hash = ((hash << 5) + hash) + ((state[i] >> 16) & 0xFF); // Add the next byte
        hash = ((hash << 5) + hash) + ((state[i] >> 24) & 0xFF); // Add the highest byte
    }

    return hash;  // Return the 64-bit hash value
}

uint64_t MurmurHash(const int* state, int dim) {
    const uint64_t seed = 0;  // Initial seed value
    const uint64_t m = 0xc6a4a7935bd19ea9ULL;
    const int r = 47;

    uint64_t hash = seed ^ (dim * m);  // Seed hash with the dimension

    for (int i = 0; i < dim; ++i) {
        uint64_t k = (uint64_t)state[i];  // Convert int to uint64_t for processing

        // Mix the bits of the key
        k *= m;
        k ^= k >> r;
        k *= m;

        // Combine the hash with the mixed key
        hash ^= k;
        hash *= m;  // Multiply by the constant
    }

    // Finalization
    hash ^= hash >> r;
    hash *= m;
    hash ^= hash >> r;

    return hash;  // Return the 64-bit hash value
}

uint64_t BuzHash(const int* state, int dim) {
    uint64_t hash = 0;  // Initialize the hash value
    uint64_t prime = 0x5bd1e995;  // A prime number used in the hash function

    for (int i = 0; i < dim; ++i) {
        hash ^= (uint64_t)state[i];      // Combine the current integer with the hash
        hash *= prime;                   // Multiply by a prime number
        hash ^= hash >> 47;              // Mix the bits
    }

    // Finalization step
    hash *= prime;                       // Final multiplication
    hash ^= hash >> 47;                  // Final mixing step

    return hash;  // Return the 64-bit hash value
}

bool same_state(int* state1, int* state2, int dim){
    for(int i = 0; i < dim; i++){
        if(state1[i] != state2[i]){
            return false; 
        }
    }
    return true; 
}

double mc(double th){ // MC flux limiter
    double min1 = fmin((1 + th)/2.0, 2.0);
    return fmax(0.0, fmin(min1, 2*th)); 
}

void inv_mat(double** mat, double* inverse, int size) {
    int i, j, k;
    double ratio;
    double* augmented = (double*)malloc(size * size * 2 * sizeof(double));
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            augmented[i * 2 * size + j] = mat[i][j];
            augmented[i * 2 * size + (j + size)] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (i = 0; i < size; i++) {
        if (augmented[i * 2 * size + i] == 0) {
            perror("Error: matrix inversion error, zero pivot element");
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < size; j++) {
            if (i != j) {
                ratio = augmented[j * 2 * size + i] / augmented[i * 2 * size + i];
                for (k = 0; k < 2 * size; k++) {
                    augmented[j * 2 * size + k] -= ratio * augmented[i * 2 * size + k];
                }
            }
        }
    }
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            inverse[i * size + j] = augmented[i * 2 * size + (j + size)] / augmented[i * 2 * size + i];
        }
    }
    free(augmented);
    return; 
}

void mul_mat_vec(double* matrix, double* vector, double* result, int size) {
    int i, j;
    for (i = 0; i < size; i++) {
        result[i] = 0;
        for (j = 0; j < size; j++) {
            result[i] += matrix[i * size + j] * vector[j];
        }
    }
    return; 
}

double* matrix_vector_multiply(int n, int m, double** matrix, double* vector){
    double* result = malloc(n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        result[i] = 0.0;  // Initialize each element of the result vector
        for (int j = 0; j < m; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    
    return result; 
}

double dot_product(double* vec1, double* vec2, int size) {
    int i;
    double result = 0;
    for (i = 0; i < size; i++) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

double gauss_probability(int dim, double* x, Meas M){ // MC Calculate gaussian probability at state x given mean and covariance
    double* M_inv = (double*)malloc(dim * dim * sizeof(double));
    double* M_inv_x = (double*)malloc(dim * sizeof(double));
    double exp_result;
    
    double diff[dim];
    for(int i = 0; i < dim; i++){
        diff[i] = x[i] - M.mean[i]; 
    }
    inv_mat(M.cov, M_inv, dim);
    mul_mat_vec(M_inv, diff, M_inv_x, dim);
    double dot_prod = dot_product(diff, M_inv_x, dim);
    exp_result = exp(-0.5 * dot_prod);

    free(M_inv);
    free(M_inv_x);

    return exp_result;
}

/*==============================================================================
                    BINARY SEARCH TREE FUNCTION DEFINITIONS
==============================================================================*/
void HashTableEntry_insert(HashTable* P, int* state, double prob, int dim, double J) {
    // uint64_t hash = FNV1a(state, dim);
    // uint64_t hash = lookup3(state, dim);
    // uint64_t hash = SuperFastHash(state, dim);
    // uint64_t hash = DJBX33A(state, dim);
    // uint64_t hash = MurmurHash(state, dim);
    uint64_t hash = BuzHash(state, dim);
    size_t index = (size_t)(hash & (uint64_t)(P->capacity - 1));
    HashTableEntry* current = P->entries[index]; 

    while(current != NULL){
        if(same_state(state, current->state, dim)){
            return; 
        }
        current = current->next; 
    }

    HashTableEntry* entry = HashTableEntry_create(dim, state, prob, J); 
    entry->next = P->entries[index]; 
    P->entries[index] = entry; 
    return;
}

void HashTable_insert(HashTable* P, Grid* G, Traj T, int* state, double prob, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    if(BOUNDS){
        double x[G->dim];
        for(int k = 0; k < G->dim; k++){
            x[k] = G->dx[k]*state[k] + G->center[k];
        }

        double J = BOUND_f(x, T.coef); 

        if(J >= G->lo_bound && J <= G->hi_bound){
            HashTableEntry_insert(P, state, prob, G->dim, J);
        }
        return;
    }
    HashTableEntry_insert(P, state, prob, G->dim, 0.0);
    return; 
}

HashTableEntry* HashTable_search(HashTable* P, int* state, int dim) {
    // uint64_t hash = FNV1a(state, dim);
    // uint64_t hash = lookup3(state, dim);
    // uint64_t hash = SuperFastHash(state, dim);
    // uint64_t hash = DJBX33A(state, dim);
    // uint64_t hash = MurmurHash(state, dim);
    uint64_t hash = BuzHash(state, dim);
    size_t index = (size_t)(hash & (uint64_t)(P->capacity - 1));

    HashTableEntry* entry = P->entries[index]; 
    while(entry != NULL){
        if(same_state(state, entry->state, dim)){
            return entry; 
        }
        entry = entry->next; 
    }
    return NULL;
}

void HashTable_delete(HashTable* P, int* state, int dim){
    // uint64_t hash = FNV1a(state, dim);
    // uint64_t hash = lookup3(state, dim);
    // uint64_t hash = SuperFastHash(state, dim);
    // uint64_t hash = DJBX33A(state, dim);
    // uint64_t hash = MurmurHash(state, dim);
    uint64_t hash = BuzHash(state, dim);
    size_t index = (size_t)(hash & (uint64_t)(P->capacity - 1));

    HashTableEntry* current = P->entries[index]; 
    HashTableEntry* previous = NULL; 
    while(current != NULL){
        if(same_state(current->state, state, dim)){
            if(previous == NULL){
                P->entries[index] = current->next; 
            }else{
                previous->next = current->next; 
            }
            free(current->state);
            free(current->v);
            free(current->ctu);
            free(current->i_nodes);
            free(current->k_nodes);
            free(current); 
            return; 
        }
        previous = current; 
        current = current->next; 
    }
    return; 
}

int get_size(HashTable* P){
    int size = 0; 
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                size += 1; 
                last = last->next; 
            }
        }
    }
    return size;
}

/*==============================================================================
                        GBEES FUNCTION DEFINITIONS
==============================================================================*/
void initialize_adv(void (*f)(double*, double*, double, double*), HashTable* P, Grid* G, Traj T, bool TV){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                if((last->new_f == 0)||(TV == true)){
                    double x[G->dim];
                    for(int i = 0; i < G->dim; i++){
                        x[i] = G->dx[i]*last->state[i]; 
                    }
                    
                    double xk[G->dim]; 
                    double* xR;
                    double advR[G->dim]; 
                    double* adv; 
                    double sum = 0;
                    for(int i = 0; i < G->dim; i++){
                        memcpy(xk, x, G->dim * sizeof(double));
                        xk[i] += G->dx[i] / 2; 
                        xR = matrix_vector_multiply(G->dim, G->dim, G->R, xk); 
                        for(int j = 0; j < G->dim; j++) xR[j] += G->center[j];
                        (*f)(advR, xR, G->t, T.coef);
                        adv = matrix_vector_multiply(G->dim, G->dim, G->Rt, advR); 
                        last->v[i] = adv[i];
                        sum += fabs(last->v[i]) / G->dx[i];

                        free(xR);  
                        free(adv);
                    }
                    last->new_f = 1;
                    last->cfl_dt = 1.0/sum;
                }
                last = last->next; 
            }
        }
    }
    return;
}

void initialize_ik_nodes(HashTable* P, Grid* G){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                if(last->ik_f == 0){
                    for(int i = 0; i < G->dim; i++){
                        // Initializing i, k nodes
                        int i_state[G->dim]; memcpy(i_state, last->state, G->dim * sizeof(int)); i_state[i] = i_state[i] - 1;  
                        int k_state[G->dim]; memcpy(k_state, last->state, G->dim * sizeof(int)); k_state[i] = k_state[i] + 1;
                        HashTableEntry* i_node = HashTable_search(P, i_state, G->dim); HashTableEntry* k_node = HashTable_search(P, k_state, G->dim);
                        last->i_nodes[i] = i_node; last->k_nodes[i] = k_node;

                        if(i_node != NULL){
                            i_node->k_nodes[i] = last; 
                        }
                        if(k_node != NULL){
                            k_node->i_nodes[i] = last; 
                        }
                    }
                    last->ik_f = 1; 
                }
                last = last->next; 
            }
        }
    }
    return; 
}

void recursive_loop(HashTable* P, Grid* G, Meas M, Traj T, int level, int* current_state, double* current_state_vec, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    if (level == G->dim) {
        double* current_state_vecR = matrix_vector_multiply(G->dim, G->dim, G->R, current_state_vec); // B.L.H.
        for(int j = 0; j < G->dim; j++) current_state_vecR[j] += G->center[j];
        double prob = gauss_probability(M.dim, current_state_vecR, M);
        HashTable_insert(P, G, T, current_state, prob, BOUNDS, BOUND_f);

        free(current_state_vecR);
        return;
    }

    int start = (int) round(-6 * G->factor[level]); 
    int end   = (int) round(6 * G->factor[level]); 
    for (int i = start; i <= end; i++) {
        current_state[level] = i;
        current_state_vec[level] = i * G->dx[level];
        recursive_loop(P, G, M, T, level + 1, current_state, current_state_vec, BOUNDS, BOUND_f);
    }
    return; 
}

void initialize_grid(void (*f)(double*, double*, double, double*), HashTable* P, Grid* G, Meas M, Traj T, bool TV, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    int current_state[G->dim]; double current_state_vec[G->dim];
    recursive_loop(P, G, M, T, 0, current_state, current_state_vec, BOUNDS, BOUND_f);
    initialize_adv(f, P, G, T, TV);
    initialize_ik_nodes(P, G);
    return; 
}

void set_bounds(HashTable* P, Grid* G){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                G->lo_bound = fmin(G->lo_bound, last->bound_val); 
                G->hi_bound = fmax(G->hi_bound, last->bound_val); 
                last = last->next;
            }
        }
    }
    return;
}

void get_sum(HashTable* P, double* prob_sum){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                *prob_sum += last->prob;
                last = last->next;
            }
        }
    }
    return;
}

void divide_sum(HashTable* P, double prob_sum, Grid* G){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                last->prob /= prob_sum; 
                if(last->prob >= G->thresh){
                    P->a_count += 1; 
                }
                P->tot_count += 1; 
                last = last->next; 
            }
        }
    }
    return;
}

void normalize_tree(HashTable* P, Grid* G){
    P->a_count = 0; P->tot_count = 0; 
    double prob_sum = 0;
    get_sum(P, &prob_sum);
    divide_sum(P, prob_sum, G);
    return; 
}

char* concat_file(const char* str1, const char* str2, int num1, const char* str3, int num2) {
    int num1_len = snprintf(NULL, 0, "%d", num1);
    int num2_len = snprintf(NULL, 0, "%d", num2);
    char* str4 = ".txt"; 
    size_t total_len = strlen(str1) + strlen(str2) + num1_len + strlen(str3) + num2_len + strlen(str4) + 1; 
    char* result = (char*)malloc(total_len * sizeof(char));
    if (result == NULL) {
        perror("Error: memory allocation failure during concat_file");
        exit(EXIT_FAILURE);
    }
    strcpy(result, str1);
    strcat(result, str2);
    char num1_str[num1_len + 1];
    sprintf(num1_str, "%d", num1);
    strcat(result, num1_str);
    strcat(result, str3);
    char num2_str[num2_len + 1];
    sprintf(num2_str, "%d", num2);
    strcat(result, num2_str);
    strcat(result, str4);

    return result;
}

void write_cells(FILE* myfile, HashTable* P, Grid G){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                if(last->prob >= G.thresh){
                    fprintf(myfile, "%.10e", last->prob);
                    
                    double x[G.dim];
                    for (int i = 0; i < G.dim; i++){
                        x[i] = G.dx[i] * last->state[i];
                    }
                    double* xR = matrix_vector_multiply(G.dim, G.dim, G.R, x); 
                    for(int i = 0; i < G.dim; i++){
                        xR[i] += G.center[i]; 
                        fprintf(myfile, " %.10e", xR[i]);
                    }
                    fprintf(myfile, "\n");
                    free(xR); 
                }
                last = last->next; 
            }
        }
    }
    return;
}

void record_pdf(HashTable* P, const char* FILE_NAME, Grid G, const double t){
    FILE* file = fopen(FILE_NAME, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: could not open file %s", FILE_NAME);
        exit(EXIT_FAILURE);
    }
    fprintf(file, "%lf\n", t);
    write_cells(file, P, G);
    fclose(file);
    return; 
}

void create_neighbors(HashTable* P, Grid G, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                if(last->prob >= G.thresh){
                    double current_v[G.dim];
                    memcpy(current_v, last->v, G.dim * sizeof(double));
                    int current_state[G.dim]; 
                    memcpy(current_state, last->state, G.dim * sizeof(int));
                    int new_state[G.dim];
                    for (int i = 0; i < G.dim; i++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        // Checking Forward Faces
                        if(current_v[i] > 0){
                            if(last->k_nodes[i] == NULL){
                                new_state[i] += 1;
                                HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);

                                // Checking Edges
                                for (int j = 0; j < G.dim; j++){
                                    memcpy(new_state, current_state, G.dim * sizeof(int));
                                    new_state[i] += 1;
                                    if(j != i){
                                        if(current_v[j] > 0){
                                            new_state[j] += 1;
                                            HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);
                                        }else if (current_v[j] < 0){
                                            new_state[j] -= 1;
                                            HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);
                                        }
                                    }
                                }
                            }else{
                                // Checking Edges
                                for (int j = 0; j < G.dim; j++){
                                    memcpy(new_state, current_state, G.dim * sizeof(int));
                                    new_state[i] += 1;
                                    if(j != i){
                                        if(current_v[j] > 0){
                                            if(last->k_nodes[i]->k_nodes[j] == NULL){
                                                new_state[j] += 1;
                                                HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);
                                            }
                                        }else if (current_v[j] < 0){
                                            if(last->k_nodes[i]->i_nodes[j] == NULL){
                                                new_state[j] -= 1;
                                                HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);
                                            }
                                        }
                                    }
                                }
                            }
                        // Checking Backward Faces
                        }else if(current_v[i] < 0){
                            if(last->i_nodes[i] == NULL){
                                new_state[i] -= 1;
                                HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);

                                // Checking Edges
                                for (int j = 0; j < G.dim; j++){
                                    memcpy(new_state, current_state, G.dim * sizeof(int));
                                    new_state[i] -= 1;
                                    if(j != i){
                                        if(current_v[j] > 0){
                                            new_state[j] += 1;
                                            HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);
                                        }else if (current_v[j] < 0){
                                            new_state[j] -= 1;
                                            HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);
                                        }
                                    }
                                }
                            }else{
                                // Checking Edges
                                for (int j = 0; j < G.dim; j++){
                                    memcpy(new_state, current_state, G.dim * sizeof(int));
                                    new_state[i] -= 1;
                                    if(j != i){
                                        if(current_v[j] > 0){
                                            if(last->i_nodes[i]->k_nodes[j] == NULL){
                                                new_state[j] += 1;
                                                HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f); 
                                            }
                                        }else if (current_v[j] < 0){
                                            if(last->i_nodes[i]->i_nodes[j] == NULL){
                                                new_state[j] -= 1;
                                                HashTable_insert(P, &G, T, new_state, 0.0, BOUNDS, BOUND_f);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                last = last->next; 
            }
        }
    }
    return; 
}

void grow_tree(void (*f)(double*, double*, double, double*), HashTable* P, Grid G, Traj T, bool TV, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    create_neighbors(P, G, T, BOUNDS, BOUND_f);
    initialize_adv(f, P, &G, T, TV);
    initialize_ik_nodes(P, &G);
    return; 
}

void check_cfl_condition(HashTable* P, Grid* G){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                G->dt = fmin(G->dt, last->cfl_dt);
                last = last->next; 
            }
        }
    }
    return;
}

void update_dcu(HashTable* P, Grid G){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                last->dcu = 0; 
                HashTableEntry* i_node; HashTableEntry* k_node; 
                for(int i = 0; i < G.dim; i++){
                    last->ctu[i] = 0.0; 
                    i_node = last->i_nodes[i]; k_node = last->k_nodes[i];

                    double dcu_p = 0; double dcu_m = 0; 

                    if(k_node != NULL){
                        dcu_p = fmax(last->v[i], 0.0) * last->prob + fmin(last->v[i], 0.0) * k_node->prob;
                    }else{
                        dcu_p = fmax(last->v[i], 0.0) * last->prob;
                    }
                    if(i_node != NULL){
                        dcu_m = fmax(i_node->v[i], 0.0)*i_node->prob + fmin(i_node->v[i], 0.0)*last->prob;
                    }

                    last->dcu -= (G.dt/G.dx[i])*(dcu_p-dcu_m);
                }
                last = last->next;
            }
        }
    }
    return; 
}

void update_ctu(HashTable* P, Grid G){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                for(int i = 0; i < G.dim; i++){
                    HashTableEntry* i_node = last->i_nodes[i];
                    HashTableEntry* j_node; HashTableEntry* p_node;
                    if(i_node!=NULL){
                        double F = G.dt*(last->prob-i_node->prob)/(2*G.dx[i]);
                        for(int j = 0; j < G.dim; j++){
                            if (j!=i){
                                j_node = last->i_nodes[j];
                                p_node = i_node->i_nodes[j];

                                last->ctu[j] -= fmax(i_node->v[i], 0.0) * fmax(last->v[j], 0.0) * F;
                                i_node->ctu[j] -= fmin(i_node->v[i], 0.0) * fmax(i_node->v[j], 0.0) * F;

                                if(j_node!=NULL){
                                    j_node->ctu[j] -= fmax(i_node->v[i], 0.0) * fmin(j_node->v[j], 0.0) * F;
                                }
                                if(p_node!=NULL){
                                    p_node->ctu[j] -= fmin(i_node->v[i], 0.0) * fmin(p_node->v[j], 0.0) * F;
                                }
                            }
                        }

                        // High-Resolution Correction Terms
                        double th;
                        if (i_node->v[i]>0){
                            HashTableEntry* i_i_node = i_node->i_nodes[i];
                            if(i_i_node != NULL){
                                th = (i_node->prob - i_i_node->prob)/(last->prob - i_node->prob);
                            }else{
                                th = (i_node->prob)/(last->prob - i_node->prob); 
                            }
                        }else{
                            HashTableEntry* k_node = last->k_nodes[i]; 
                            if(k_node != NULL){
                                th = (k_node->prob - last->prob)/(last->prob - i_node->prob);
                            }else{
                                th = (-last->prob)/(last->prob - i_node->prob);
                            }
                        }
                        i_node->ctu[i] += fabs(i_node->v[i])*(G.dx[i]/G.dt - fabs(i_node->v[i]))*F*mc(th);
                    }
                }
                last = last->next; 
            }
        }
    }
    return; 
}

void godunov_method(HashTable* P, Grid G){
    update_dcu(P, G); 
    update_ctu(P, G);
    return; 
}

void update_prob(HashTable* P, Grid G){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                last->prob += last->dcu;
                for(int i = 0; i < G.dim; i++){
                    HashTableEntry* i_node = last->i_nodes[i]; 
                    if(i_node != NULL){
                        last->prob -= (G.dt/G.dx[i])*(last->ctu[i] - i_node->ctu[i]);
                    }else{
                        last->prob -= (G.dt/G.dx[i])*(last->ctu[i]);
                    }
                }
                last->prob = fmax(last->prob, 0.0); 
                last = last->next;
            }
        }
    }
    return; 
}

void mark_cells(HashTable* P, Grid G, double* del_probs, int** del_states, int* idx_s){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                last->ik_f = 0; bool DELETE = true;
                if(last->prob < G.thresh){
                    for(int i = 0; i < G.dim; i++){
                        // Looking at Backwards Node
                        HashTableEntry* i_node = last->i_nodes[i]; 

                        if(i_node != NULL){
                            if ((i_node->v[i]>0)&&(i_node->prob >= G.thresh)){
                                DELETE = false;
                                break;
                            }else{
                                for (int j = 0; j < G.dim; j++){
                                    if(j!=i){
                                        HashTableEntry* i_i_node = i_node->i_nodes[j]; 
                                        if(i_i_node != NULL){
                                            if ((i_i_node->v[i]>0)&&(i_i_node->v[j]>0)&&(i_i_node->prob >= G.thresh)){
                                                DELETE = false;
                                                break;
                                            }
                                        }
                                        HashTableEntry* i_k_node = i_node->k_nodes[j]; 
                                        if(i_k_node != NULL){
                                            if ((i_k_node->v[i]>0)&&(i_k_node->v[j]<0)&&(i_k_node->prob >= G.thresh)){
                                                DELETE = false;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        // Looking at Forwards Node
                        HashTableEntry* k_node = last->k_nodes[i]; 

                        if(k_node != NULL){
                            if ((k_node->v[i]<0)&&(k_node->prob >= G.thresh)){
                                DELETE = false;
                                break;
                            }else{
                                for (int j = 0; j < G.dim; j++){
                                    if(j!=i){
                                        HashTableEntry* k_i_node = k_node->i_nodes[j];
                                        if(k_i_node != NULL){
                                            if ((k_i_node->v[i]<0)&&(k_i_node->v[j]>0)&&(k_i_node->prob >= G.thresh)){
                                                DELETE = false;
                                                break;
                                            }
                                        }
                                        
                                        HashTableEntry* k_k_node = k_node->k_nodes[j];
                                        if(k_k_node != NULL){
                                            if ((k_k_node->v[i]<0)&&(k_k_node->v[j]<0)&&(k_k_node->prob >= G.thresh)){
                                                DELETE = false;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if(DELETE){
                        del_probs[*idx_s] = last->prob; 
                        del_states[*idx_s] = last->state; 
                        *idx_s += 1; 
                    }
                }
                last = last->next; 
            }
        }
    }
    return; 
}

#ifdef __linux__ 
  int compare_indices(const void* a, const void* b, void* del_probs) {
#else  
  int compare_indices(void* del_probs, const void* a, const void* b) {
#endif  
    const double* double_list = (double*)del_probs;
    int idx1 = *(const int*)a;
    int idx2 = *(const int*)b;

    if (double_list[idx1] < double_list[idx2])
        return -1;
    else if (double_list[idx1] > double_list[idx2])
        return 1;
    else
        return 0;
}

void sort_by_double(double* del_probs, int** del_states, size_t n) {
    int* indices = malloc(n * sizeof(int));
    if (indices == NULL) {
        const char* error_string = "sorted indices creation"; 
        exit_nomem(error_string); 
    }

    for (size_t i = 0; i < n; i++) {
        indices[i] = i;
    }

    // Sort the indices based on the corresponding doubles    
    #ifdef __linux__ 
        qsort_r(indices, n, sizeof(int), compare_indices, del_probs);
    #else  
        qsort_r(indices, n, sizeof(int), del_probs, compare_indices);  
    #endif 

    // Create a temporary array to hold the sorted uint64_t's
    double* sorted_del_probs = (double*)malloc(n * sizeof(double));
    int** sorted_del_states = (int**)malloc(n * sizeof(int*));
    if ((sorted_del_probs == NULL)||(sorted_del_states == NULL)) {
        free(indices);
        const char* error_string = "sorted keys creation"; 
        exit_nomem(error_string); 
    }

    for (size_t i = 0; i < n; i++) {
        sorted_del_probs[i] = del_probs[indices[i]];
        sorted_del_states[i] = del_states[indices[i]];
    }

    for (size_t i = 0; i < n; i++) {
        del_probs[i] = sorted_del_probs[i];
        del_states[i] = sorted_del_states[i];
    }

    free(indices);
    free(sorted_del_probs);
    free(sorted_del_states);
    return; 
}

void delete_cells(HashTable* P, Grid G, double* del_probs, int** del_states, int idx){
    double del_prob_sum = 0; 
    for(int i = 0; i < idx; i++){
        del_prob_sum += del_probs[i]; 
        if(del_probs[i]/(1-del_prob_sum) < G.thresh){
            HashTable_delete(P, del_states[i], G.dim);
        }else{
            break;
        }
    }
    return; 
}

void prune_tree(HashTable* P, Grid G){
    double* del_probs = malloc(P->tot_count * sizeof(double));
    int** del_states = (int**)malloc(P->tot_count * sizeof(int*));    
    if ((del_probs == NULL)||(del_states == NULL)){
        const char* error_string = "del_probs/del_keys creation"; 
        exit_nomem(error_string); 
    }
    int idx = 0; 
    mark_cells(P, G, del_probs, del_states, &idx);
    sort_by_double(del_probs, del_states, idx);
    delete_cells(P, G, del_probs, del_states, idx);
    initialize_ik_nodes(P, &G);
    free(del_probs); 
    free(del_states); 
    return; 
}

void meas_up_recursive(void (*h)(double*, double*, double, double*), HashTable* P, Grid G, Meas M, Traj T){
    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                double x[G.dim];
                for(int i = 0; i < G.dim; i++){
                    x[i] = G.dx[i]*last->state[i];
                }

                double* xR = matrix_vector_multiply(G.dim, G.dim, G.R, x); 
                for(int i = 0; i < G.dim; i++) xR[i] += G.center[i];
                double y[M.dim];
                (*h)(y, xR, G.t, T.coef); 

                double prob = gauss_probability(M.dim, y, M);   
                last->prob *= prob;
                last = last->next; 

                free(xR); 
            }
        }
    }
    return;
}


void record_collisions(HashTable* P, const char* FILE_NAME){
    FILE* myfile = fopen(FILE_NAME, "w");
    if (myfile == NULL) {
        fprintf(stderr, "Error: could not open file %s", FILE_NAME);
        exit(EXIT_FAILURE);
    }

    for(int idx = 0; idx < P->capacity; idx++){
        HashTableEntry* head = P->entries[idx];
        int entry_count = 0; 
        if(head != NULL){
            HashTableEntry* last = head;
            while(last != NULL){
                entry_count += 1; 
                last = last->next; 
            }
        }
        fprintf(myfile, "%d\n", entry_count);
    }

    fclose(myfile);
    return;
}

void run_gbees(void (*f)(double*, double*, double, double*), void (*h)(double*, double*, double, double*), double (*BOUND_f)(double*, double*), Grid G, Meas M, Traj T, char* P_DIR, char* M_DIR, int NUM_DIST, int NUM_MEAS, int DEL_STEP, int OUTPUT_FREQ, int CAPACITY, int DIM_h, bool OUTPUT, bool RECORD, bool MEASURE, bool BOUNDS, bool COLLISIONS, bool TV){
    char* P_PATH; char* C_PATH; 
    double RECORD_TIME;      

    HashTable* P = HashTable_create(CAPACITY); 

    printf("Initializing distribution...\n\n");

    initialize_grid(f, P, &G, M, T, TV, BOUNDS, BOUND_f); 
    if(BOUNDS){G.lo_bound = DBL_MAX; G.hi_bound = -DBL_MAX; set_bounds(P, &G);} 
    normalize_tree(P, &G); 

    printf("Entering time marching...\n\n");

    clock_t start = clock(); 
    for(int nm = 0; nm < NUM_MEAS; nm++){
        printf("Timestep: %d-0, Program time: %f s, Sim. time: %f", nm, ((double)(clock()-start))/CLOCKS_PER_SEC, G.t); 
        printf(" TU, Active/Total Cells: %zu/%zu\n", P->a_count, P->tot_count); 
        if(RECORD){P_PATH = concat_file(P_DIR, "/P", nm, "/pdf_", 0); record_pdf(P, P_PATH, G, G.t); free(P_PATH);};
        if(COLLISIONS){C_PATH = concat_file(P_DIR, "/C", nm, "/c_", 0); record_collisions(P, C_PATH); free(C_PATH);}

        double mt = 0; int record_count = 1; int step_count = 1; double rt; 
        while(fabs(mt - M.T) > TOL) { // time between discrete measurements

            rt = 0; RECORD_TIME = M.T/(NUM_DIST-1);      
            while (rt < RECORD_TIME) { // time between PDF recordings

                grow_tree(f, P, G, T, TV, BOUNDS, BOUND_f);
                check_cfl_condition(P, &G); 
                G.dt = fmin(G.dt, RECORD_TIME - rt);
                rt += G.dt; G.t += G.dt; 
                godunov_method(P, G);
                update_prob(P, G);
                normalize_tree(P, &G); 

                if (step_count % DEL_STEP == 0) { // deletion procedure
                    prune_tree(P, G);
                    normalize_tree(P, &G);
                }

                if ((OUTPUT) && (step_count % OUTPUT_FREQ == 0)) { // print size to terminal
                    printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f", nm, step_count, ((double)(clock()-start))/CLOCKS_PER_SEC, G.t); 
                    printf(" TU, Active/Total Cells: %zu/%zu\n", P->a_count, P->tot_count); 
                }

                if(COLLISIONS){C_PATH = concat_file(P_DIR, "/C", nm, "/c_", step_count); record_collisions(P, C_PATH); free(C_PATH);}

                step_count += 1; G.dt = DBL_MAX;
            }
            
            if (((step_count-1) % OUTPUT_FREQ != 0)||(!OUTPUT)){ // print size to terminal  
                printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f", nm, step_count - 1, ((double)(clock()-start))/CLOCKS_PER_SEC, G.t); 
                printf(" TU, Active/Total Cells: %zu/%zu\n", P->a_count, P->tot_count); 
            }
            
            if (RECORD) { // record PDF
                P_PATH = concat_file(P_DIR, "/P", nm, "/pdf_", record_count); 
                record_pdf(P, P_PATH, G, G.t);
                record_count += 1;
                free(P_PATH);
            }

            if (COLLISIONS) { // record collisions
                C_PATH = concat_file(P_DIR, "/C", nm, "/c_", record_count); 
                record_collisions(P, C_PATH);
                record_count += 1;
                free(C_PATH);
            }

            mt += rt;
        }

        if ((MEASURE) && (nm < NUM_MEAS - 1)) { // peform discrete measurement update

            printf("\nPERFORMING BAYESIAN UPDATE AT: %f TU...\n\n", G.t);

            // freeing previous measurement 
            Meas_free(&M);                                                    

            // reading new measurement
            char* M_FILE_NAME = "measurement";                                             
            char* M_FILE_EXT = ".txt"; 
            int length = snprintf(NULL, 0, "%s%d%s", M_FILE_NAME, nm+1, M_FILE_EXT) + 1;
            char* M_FILE = (char*)malloc(length);
            if (M_FILE == NULL) {
                const char* error_string = "discrete measurement update"; 
                exit_nomem(error_string); 
            }
            snprintf(M_FILE, length, "%s%d%s", M_FILE_NAME, nm+1, M_FILE_EXT);
            M = Meas_create(DIM_h, M_DIR, M_FILE);                                       
            free(M_FILE); 

            // performing discrete update
            meas_up_recursive(h, P, G, M, T);                                     
            normalize_tree(P, &G); 
            prune_tree(P, G);
            normalize_tree(P, &G); 
        }
    }
    
    Meas_free(&M); 
    Grid_free(&G); 
    Traj_free(&T);
    HashTable_free(P);

    printf("Time marching complete.\n");

    return;
}
