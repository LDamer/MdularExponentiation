#ifndef ___TO_IMPLEMENT___
#define ___TO_IMPLEMENT___

#include <gmp.h>

int get_random_exponent(mpz_t e, int bits, gmp_randstate_t randomState);
int exponentiate_binary(mpz_t result, mpz_t g, mpz_t e, mpz_t modulus, long *count_S, long *count_M);
int exponentiate_k_ary(mpz_t result, mpz_t g, mpz_t e, mpz_t modulus, int k, long *count_S, long *count_M);
int exponentiate_sliding_window(mpz_t result, mpz_t g, mpz_t e, mpz_t modulus, int k, long *count_S, long *count_M);

size_t get_fb_LUT_size(size_t maximal_exponent_length);
size_t get_imp_fb_LUT_size(size_t maximal_exponent_length, unsigned window_size);
mpz_t *allocate_fixed_base_LUT(size_t lut_size);
void free_fixed_base_LUT(mpz_t *lut, size_t lut_size);
void precompute_fixed_base_LUT(mpz_t *lut, mpz_t g, mpz_t modulus, size_t lut_size);
void precompute_improved_fixed_base_LUT(mpz_t *lut, mpz_t g, mpz_t modulus, size_t lut_size, unsigned window_size);
int exponentiate_fixed_based(mpz_t result, mpz_t e, mpz_t modulus, mpz_t *lut, size_t maximal_exponent_length, long *count_S, long *count_M);
int exponentiate_improved_fixed_based(mpz_t result, mpz_t e, mpz_t modulus, mpz_t *lut, size_t maximal_exponent_length, unsigned window_size, long *count_S, long *count_M);


int ecc_double_add(mpz_t resultX, mpz_t resultY, mpz_t a, mpz_t b, mpz_t p, mpz_t q, mpz_t inputX, mpz_t inputY, mpz_t factor, long *count_D, long *count_A);
int ecc_naf_double_add(mpz_t resultX, mpz_t resultY, mpz_t a, mpz_t b, mpz_t p, mpz_t q, mpz_t inputX, mpz_t inputY, mpz_t factor, long *count_D, long *count_A);
#endif
