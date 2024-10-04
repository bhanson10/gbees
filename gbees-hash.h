// gbees-hash.h, https://github.com/bhanson10/gbees-hash
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef GBEES_HASH_H
#define GBEES_HASH_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

/*==============================================================================
                            STRUCTURE DEFINITIONS
==============================================================================*/
typedef struct Meas {
    int dim; 
    double *mean;
    double **cov;
    double T;
} Meas;

Meas Meas_create(int dim, const char* M_DIR, const char* M_FILE);

void Meas_free(Meas *M);

typedef struct Grid {
    int dim; 
    double thresh;
    double t; 
    double dt;
    double *center;
    double *dx;
    double hi_bound;
    double lo_bound;
} Grid;

Grid Grid_create(int dim, double thresh, double* center, double* dx);

void Grid_free(Grid* G);

typedef struct Traj {
    double *coef;
} Traj;

Traj Traj_create(int n, double* coef);

void Traj_free(Traj* T);

typedef struct HashTableEntry HashTableEntry;

struct HashTableEntry { 
    int *state;
    double prob;
    double *v;
    double *ctu;
    HashTableEntry **i_nodes;
    HashTableEntry **k_nodes;
    double dcu;
    double cfl_dt;
    int new_f;
    int ik_f;
    double bound_val; 
    HashTableEntry* next; 
};

typedef struct HashTable HashTable;

struct HashTable { 
    HashTableEntry **entries; 
    size_t capacity; 
    size_t a_count; 
    size_t tot_count; 
};

HashTable* HashTable_create(int CAPACITY);

HashTableEntry* HashTableEntry_create(int dim, int* state, double prob, double J);

void HashTable_free(HashTable* P);

/*==============================================================================
                        NON-MEMBER FUNCTION DEFINITIONS
==============================================================================*/
void exit_nomem(const char* error_string);

uint64_t FNV1a(int* state, int dim);

bool same_state(int* state1, int* state2, int dim);

double mc(double th);

void inv_mat(double** mat, double* inverse, int size);

void mul_mat_vec(double* matrix, double* vector, double* result, int size);

double dot_product(double* vec1, double* vec2, int size);

double gauss_probability(int dim, double* x, Meas M);

/*==============================================================================
                    HASH TABLE FUNCTION DEFINITIONS
==============================================================================*/
void HashTableEntry_insert(HashTable* P, int* state, double prob, int dim, double J);

void HashTable_insert(HashTable* P, Grid* G, Traj T, int* state, double prob, bool BOUNDS, double (*BOUND_f)(double*, double*));

HashTableEntry* HashTable_search(HashTable* P, int* state, int dim);

void HashTable_delete(HashTable* P, int* state, int dim); 

int get_size(HashTable* P);

/*==============================================================================
                        GBEES FUNCTION DEFINITIONS
==============================================================================*/
void initialize_adv(void (*f)(double*, double*, double, double*, double*), HashTable* P, Grid* G, Traj T);

void initialize_ik_nodes(HashTable* P, Grid* G);

void recursive_loop(HashTable* P, Grid* G, Meas M, Traj T, int level, int* current_state, double* current_state_vec, bool BOUNDS, double (*BOUND_f)(double*, double*));

void initialize_grid(void (*f)(double*, double*, double, double*, double*), HashTable* P, Grid* G, Meas M, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*));

void set_bounds(HashTable* P, Grid* G);

void get_sum(HashTable* P, double* prob_sum);

void divide_sum(HashTable* P, double prob_sum, Grid* G);

void normalize_tree(HashTable* P, Grid* G);

char* concat_file(const char* str1, const char* str2, int num1, const char* str3, int num2);

void write_cells(FILE* myfile, HashTable* P, Grid G);

void record_pdf(HashTable* P, const char* FILE_NAME, Grid G, const double t);

void create_neighbors(HashTable* P, Grid G, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*));

void grow_tree(void (*f)(double*, double*, double, double*, double*), HashTable* P, Grid G, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*));

void check_cfl_condition(HashTable* P, Grid* G);

void update_dcu(HashTable* P, Grid G);

void update_ctu(HashTable* P, Grid G);

void godunov_method(HashTable* P, Grid G);

void update_prob(HashTable* P, Grid G);

void mark_cells(HashTable* P, Grid G, double* del_probs, int** del_states, int* idx_s);

#ifdef __linux__ 
  int compare_indices(const void *a, const void *b, void *del_probs);
#else  
  int compare_indices(void *del_probs, const void *a, const void *b);
#endif  

void sort_by_double(double* del_probs, int** del_states, size_t n);

void delete_cells(HashTable* P, Grid G, double* del_probs, int** del_states, int idx);

void prune_tree(HashTable* P, Grid G);

void meas_up_recursive(void (*h)(double*, double*, double, double*, double*), HashTable* P, Grid G, Meas M, Traj T);

void record_collisions(HashTable* P, const char* FILE_NAME); 

void run_gbees(void (*f)(double*, double*, double, double*, double*), void (*h)(double*, double*, double, double*, double*), double (*BOUND_f)(double*, double*), Grid G, Meas M, Traj T, char* P_DIR, char* M_DIR, int NUM_DIST, int NUM_MEAS, int DEL_STEP, int OUTPUT_FREQ, int CAPACITY, int DIM_h, bool OUTPUT, bool RECORD, bool MEASURE, bool BOUNDS, bool COLLISIONS);

#endif // GBEES_HASH_H
