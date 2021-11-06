#include <gmp.h> /* GNU Multi Precision library */
#include "to_implement.h"
#include "ecc_functions.h"

int get_random_exponent(mpz_t e, int bits, gmp_randstate_t randomState)
{
    mpz_urandomb(e, randomState, bits);
    return 0;
}


int exponentiate_binary(mpz_t result, mpz_t g, mpz_t e, mpz_t modulus, long *count_S, long *count_M)
{
    unsigned l = mpz_sizeinbase(e, 2);
    mpz_set(result, g);
    for(int i = l-2; i >= 0; i--){
        mpz_mul(result, result, result);
        mpz_mod(result, result, modulus);
        (*count_S)++;
        if(mpz_tstbit(e, i)){
            mpz_mul(result, result, g);
            mpz_mod(result, result, modulus);
            (*count_M)++;
        }
    }
	return 0; 
}


int exponentiate_k_ary(mpz_t result, mpz_t g, mpz_t e, mpz_t modulus, int k, long *count_S, long *count_M)
{
    //l stores the number of bits
    unsigned l = mpz_sizeinbase(e,2);
    //t stores the number of blocks of size k
    unsigned t = (l+1)/k;
    //SETUP LUT
    int array_size;
    //determine needed LUT size: 2^k entries....but without the first entry 2^0=1 because this would be unnecessary.=> array_size = 2^k-1
    //we catch this case while computation, to avoid performing a TLU for multiplying with 1.
    switch(k){
        case 3: array_size=8;break;
        case 4: array_size=16;break;
        case 5: array_size=32;break;
        case 6: array_size=64;break;
        case 7: array_size=128;break;
        case 8: array_size=256;break;
        case 9: array_size=512;break;
        case 10: array_size=1024;break;
        default: printf("\n\nK is too big!\n\n");break;
    }
    mpz_t lut[array_size];
    //initialize trivial entries in LUT
    mpz_init(lut[0]);
    mpz_set_ui(lut[0], 1);
    mpz_init(lut[1]);
    mpz_set(lut[1], g);
    //precomputation of the whole LUT
    for(int i = 2; i < array_size; i++){
        mpz_init(lut[i]);
        mpz_mul(lut[i], lut[i-1], g);
        mpz_mod(lut[i], lut[i], modulus);
        //gmp_printf("%d: %Zd\n", i, lut[i]);
        //first round is g*g, therefore a squaring. Rest is multiplication
        (*count_M)++;
    }
    //first round of for-loop is g*g but we counted a multiplication.....correcting that now
    (*count_M)--;
    (*count_S)++;
    //index will store the current block E_i of the exponent
    mpz_t index;
    mpz_init(index);
    //read first Block E_t of exponent and store it in index(step 2 in the algorithm in the script)
    for(int i = 0; i < k; i++){
        if(mpz_tstbit(e,(t+1)*k-1-i)){
            mpz_setbit(index, k-i-1);
        }
    }
    //gmp_printf("E_t = %Zd\n", index);
    //use the first Block from above and initialize 'result' with g**E_t from the LUT instead of squaring/multiplying with 1.
    mpz_set(result, lut[mpz_get_ui(index)]);
    mpz_set_ui(index, 0);
    //compute rest of the exponent using the LUT(step 3 in the script)
    for(int i = t; i > 0; i--){
        //Iterate over thing, we have to do k times... 1.) determine the next k bit an store in 'index' 2.) square k times.
        for(int j = 0; j < k; j++){
            //Sqaure to create space
            mpz_mul(result, result, result);
            mpz_mod(result, result, modulus);
            (*count_S)++;
            //test the next corresponding bit
            if(mpz_tstbit(e, t*k-1-j+(i-t)*k)){
                mpz_setbit(index, k-j-1);
            }
        }
        //gmp_printf("E_i = %Zd\n", index);
        //multiply with the corresponding entry of the LUT using the counted bits stored in index
        //avoid multiplying with 1.
        //printf("E_i = %d\n", mpz_get_ui(index));
        if(mpz_get_ui(index) > 0){
            mpz_mul(result, result, lut[mpz_get_ui(index)]);
            mpz_mod(result, result, modulus);
            (*count_M)++;
        }else{
            //printf("\n exponent should be 0...but actually ist 0");
            //gmp_printf("%Zd\n", index);
        }
        mpz_set_ui(index, 0);
    }
    //free the array
    for(int i = 0; i < array_size; i++){
        mpz_clear(lut[i]);
    }
    //free the variable 'index'
    mpz_clear(index);
	//gmp_printf("\nreturn=%Zd\n", result);
	return 0; 
}


int exponentiate_sliding_window(mpz_t result, mpz_t g, mpz_t e, mpz_t modulus, int k, long *count_S, long *count_M)
{
    //l stores the number of bits
    unsigned l = mpz_sizeinbase(e,2);
    //SETUP LUT
    int array_size;
    //determine needed LUT size: 2^k entries....but without the first entry 2^0=1 because this would be unnecessary.=> array_size = 2^k-1
    //we catch this case while computation, to avoid performing a TLU for multiplying with 1.
    switch(k){
        case 3: array_size=4;break;
        case 4: array_size=8;break;
        case 5: array_size=16;break;
        case 6: array_size=32;break;
        case 7: array_size=64;break;
        case 8: array_size=128;break;
        case 9: array_size=256;break;
        case 10: array_size=512;break;
        default: printf("\n\nK is too big!\n\n");break;
    }
    mpz_t lut[array_size];
    //initialize trivial entries in LUT
    mpz_init(lut[0]);
    mpz_set(lut[0], g);

    mpz_t g_squared;
    mpz_init(g_squared);
    mpz_mul(g_squared, g, g);
    mpz_mod(g_squared, g_squared, modulus);
    (*count_S)++;
    //precomputation of the whole LUT
    for(int i = 1; i < array_size; i++){
        mpz_init(lut[i]);
        mpz_mul(lut[i], lut[i-1], g_squared);
        mpz_mod(lut[i], lut[i], modulus);
        //gmp_printf("%d: %Zd\n", i, lut[i]);
        //first round is g*g, therefore a squaring. Rest is multiplication
        (*count_M)++;
    }
    mpz_clear(g_squared);
    //index will store the current window of the exponent
    mpz_t index;
    mpz_init(index);
    mpz_set_ui(index, 0);
    int last_bit_index;
    //compute rest of the exponent using the sliding window algorithm.
    for(int i = l-1; i >= 0; i--){
        if(mpz_tstbit(e, i) == 0){
            //printf("we found a 0 at %d and squared\n", i);
            //we need to square and skip to the next bit
            mpz_mul(result, result, result);
            mpz_mod(result, result, modulus);
            //gmp_printf("result is now = %Zd\n", result);
            (*count_S)++;
            continue;
        }else{
            //determine the position of the last 1-bit in the maximum possible range
            //printf("we found a 1 at %d\n", i);
            last_bit_index = i;
            for(int j = i-k+1; j < i; j++){
                //printf("watching at %d\n", j);
                if(mpz_tstbit(e, j) == 1){
                    //we found the very last 1 bit in the window
                    last_bit_index = j;
                    break;
                }
            }
            //printf("the window goes to %d\n", last_bit_index);
            //read the window bits and store it in index
            for(int j = last_bit_index; j < i+1; j++){
                if(mpz_tstbit(e, j)){
                   // printf("set index to 1 at pos %d\n", j-last_bit_index);
                    mpz_setbit(index, j-last_bit_index);
                }
            }
            if(i < l-1){
            //always true but not in the first round...wen can use set() and dont need to sqaure and mul -> see else branch
            //now we need to create space for the next window
                for(int j = 0; j < i-last_bit_index+1; j++){
                    mpz_mul(result, result, result);
                    mpz_mod(result, result, modulus);
                    (*count_S)++;
                }
                //printf("we squared %d times\n", i-last_bit_index+1);
                //now we can multiply with the corresponding LUT entry
                mpz_mul(result, result, lut[mpz_get_ui(index) / 2]);
                mpz_mod(result, result, modulus);
                (*count_M)++;
                //gmp_printf("we multiplied with %Zd\n", lut[mpz_get_ui(index) / 2]);
                //gmp_printf("result is now = %Zd\n", result);
            }else{
                mpz_set(result, lut[mpz_get_ui(index) / 2]);

            }
            //reset index variable for futher computations
            mpz_set_ui(index, 0);
            //decrease i because we computed multiple bits
            i = last_bit_index;
        }
    }
    //free the array
    for(int i = 0; i < array_size; i++){
        mpz_clear(lut[i]);
    }
    //free the variable 'index'
    mpz_clear(index);
    //gmp_printf("result = %Zd\n", result);
	//gmp_printf("\nreturn=%Zd\n", result);
	return 0; 
}

size_t get_fb_LUT_size(size_t maximal_exponent_length) {
	//highly efficient!
	return maximal_exponent_length;
}

size_t get_imp_fb_LUT_size(size_t maximal_exponent_length, unsigned window_size) {   
    if((maximal_exponent_length)%window_size == 0){
       return maximal_exponent_length/window_size;
    }else{
        return (maximal_exponent_length)/window_size + 1;
    }

    return 0;
}

mpz_t *allocate_fixed_base_LUT(size_t lut_size) {
    mpz_t * lut = malloc(lut_size * sizeof(mpz_t));
    for(int i = 0; i < lut_size; i++){
        mpz_init(lut[i]);
        //mpz_set_ui(lut[i], i);
    }
    return lut;
}

void free_fixed_base_LUT(mpz_t *lut, size_t lut_size) {
    for(int i = 0; i < lut_size; i++){
        mpz_clear(lut[i]);
    }
}

void precompute_fixed_base_LUT(mpz_t *lut, mpz_t g, mpz_t n, size_t lut_size) {
     mpz_set(*lut, g);
     for(int i = 1; i < lut_size; i++){
        mpz_mul(*(lut+i), *(lut+i-1), *(lut+i-1));
        mpz_mod(*(lut+i),*(lut+i), n);
     }
}

void precompute_improved_fixed_base_LUT(mpz_t *lut, mpz_t g, mpz_t modulus, size_t lut_size, unsigned window_size) {
    
    mpz_set(*lut, g);
    for(int i = 1; i < lut_size; i++){
        //square it
        mpz_mul(*(lut+i), *(lut+i-1), *(lut+i-1));
        mpz_mod(*(lut+i),*(lut+i), modulus);
        for(int j = 0; j < window_size-1; j++){
            //sqaure it the neccessary times according to b or window_size.
            mpz_mul(*(lut+i), *(lut+i), *(lut+i));
            mpz_mod(*(lut+i),*(lut+i), modulus);
        }
     }
}

int exponentiate_fixed_based(mpz_t result, mpz_t e, mpz_t modulus, mpz_t *lut, size_t maximal_exponent_length, long *count_S, long *count_M) {
	//first one is "geschenkt"
	mpz_set(result, lut[maximal_exponent_length-1]);
	//compute the rest
	for(int i = maximal_exponent_length-2; i >= 0; i--){
        if(mpz_tstbit(e, i)){
            // found a one, so we need to multiply...
            mpz_mul(result, result, lut[i]);
            mpz_mod(result, result, modulus);
            (*count_M)++;
        }
	}
	return 0; 
}

int exponentiate_improved_fixed_based(mpz_t result, mpz_t e, mpz_t modulus, mpz_t *lut, size_t maximal_exponent_length, unsigned window_size, long *count_S, long *count_M) {

    mpz_set_ui(result, 1); //A
    mpz_t B;
    mpz_init(B);
    mpz_set_ui(B, 1);//B
    int t;
    int t = maximal_exponent_length/window_size;
    if((maximal_exponent_length)%window_size == 0){
       t -= 1;
    }
	int b = 1 << (window_size);
    //compute e in base b.
    mpz_t * eInBaseB = malloc((t+1) * sizeof(mpz_t));
    mpz_t mask;
    mpz_init(mask);
    mpz_set_ui(mask, 1);
    for(int i = 0; i < window_size-1; i++){
        mpz_mul_2exp (mask, mask, 1);
        mpz_setbit(mask, 0);
    }
    mpz_mul_2exp (mask, mask, t*window_size);
    //go through e with the mask
    for(int i = 0; i < t+1; i++){
        mpz_init(*(eInBaseB+i));
        mpz_and(*(eInBaseB+i), e, mask);
        mpz_tdiv_q_2exp(*(eInBaseB+i), *(eInBaseB+i), (t-i)*window_size);
        mpz_tdiv_q_2exp(mask, mask, window_size);
    }
    //actual algorithm from script
    for(int j = b-1; j > 0; j--){
        for(int i = 0; i <= t; i++){
       	    if(mpz_cmp_ui(*(eInBaseB+i), j) == 0){
                mpz_mul(B, B, lut[t-i]);
                mpz_mod(B, B, modulus);
                (*count_M)++;
            }
        }
        mpz_mul(result, result, B);
        mpz_mod(result, result, modulus);
        (*count_M)++;
    }
    mpz_clear(B);
    mpz_clear(mask);
    free(eInBaseB);
return 0;
}

int ecc_double_add(mpz_t resultX, mpz_t resultY, mpz_t a, mpz_t b, mpz_t p, mpz_t q, mpz_t inputX, mpz_t inputY, mpz_t factor, long *count_D, long *count_A)
{
	unsigned l = mpz_sizeinbase(factor, 2);
	//reduce factor mod group oder q
	mpz_mod(factor, factor, q);
	//set result. Do not need to compute first bit of the factor.
    	mpz_set(resultX, inputX);
   	mpz_set(resultY, inputY);
    	for(int i = l-2; i >= 0; i--){
        //mpz_mul(result, result, result);
        ecc_op_double(resultX, resultY, a, b, p, resultX, resultY);
        //mpz_mod(result, result, modulus);
        (*count_D)++;
        if(mpz_tstbit(factor, i)){
            ecc_op_add(resultX, resultY, a, b, p, resultX, resultY, inputX, inputY);
            //mpz_mod(result, result, modulus);
            (*count_A)++;
        }
    }
	return 0;
}

int ecc_naf_double_add(mpz_t resultX, mpz_t resultY, mpz_t a, mpz_t b, mpz_t p, mpz_t q, mpz_t inputX, mpz_t inputY, mpz_t factor, long *count_D, long *count_A)
{
	mpz_mod(factor, factor, q);
	unsigned l = mpz_sizeinbase(factor, 2);
	//store naf representation in array. potentially one bit longer
	int naf[l+1];
	//array initializes NOT with all zeros, but sometimes with random values....(??)
	//if the first bit is 1 by accident, but the naf would compute a zero, the naf repr. would be wrong
	//therefore, initialize naf[l] = 0; additionally
	naf[l] = 0;
	int i = 0;
    while(mpz_cmp_ui(factor, 0) > 0){
        if(mpz_tstbit(factor, 0)){
            naf[i] = 2-(2*mpz_tstbit(factor, 1) + mpz_tstbit(factor, 0));
            if(naf[i] == -1){
                //printf("read -1\n");
                mpz_add_ui(factor, factor, 1);
            }else if(naf[i] == 1){
                //printf("read 1\n");
                mpz_sub_ui(factor, factor, 1);
                //gmp_printf("factor = %Zd\n", factor);
            }
        }else{
            //printf("read 0\n");
            naf[i] = 0;
        }
        mpz_fdiv_q_ui(factor, factor, 2);
        i++;
	}
	//double and add with naf representation
	//fiqure out where the first one is
	//printf("NAF REP:\n");
    for(int i = 0; i < l+1; i++){
        //printf("%d: %d\n", i, naf[i]);
    }
	int starting_point;
	if(naf[l] == 1){
        //printf("starting at l-1 \n");
        starting_point = l-1;
	}else{
	    //printf("starting at l-2 \n");
        starting_point = l-2;
	}
	//store inverse point
	mpz_t invInputY;
	mpz_init(invInputY);
    mpz_ui_sub(invInputY, 0, inputY);
	mpz_mod(invInputY, invInputY, p);

	mpz_set(resultX, inputX);
    mpz_set(resultY, inputY);
    for(int i = starting_point; i >= 0; i--){
        //printf("lookin at bit %d : %d\n", i, naf[i]);
        ecc_op_double(resultX, resultY, a, b, p, resultX, resultY);
        (*count_D)++;
        if(naf[i] == 1){
            ecc_op_add(resultX, resultY, a, b, p, resultX, resultY, inputX, inputY);
            (*count_A)++;
        }else if(naf[i] == -1){
            ecc_op_add(resultX, resultY, a, b, p, resultX, resultY, inputX, invInputY);
            (*count_A)++;
        }
    }
	return 0; 
}
