#include <iostream>
#include <immintrin.h>
#include <assert.h>
//#include <avxintrin.h>
//#include <emmintrin.h>
#include <math.h>

void do_block(int MAT_SIZE, int M, int N, int K, double* A, double* B, double* C){

    for(int i = 0; i < M; i = i+6){
        for(int j = 0; j < N; j = j+8){
            __m256d mA0,mA1,mB0,mB1;

            __m256d result0_0  = _mm256_set1_pd(0);
            __m256d result1_0  = _mm256_set1_pd(0);
            __m256d result2_0  = _mm256_set1_pd(0);
            __m256d result3_0  = _mm256_set1_pd(0);
            __m256d result4_0  = _mm256_set1_pd(0);
            __m256d result5_0  = _mm256_set1_pd(0);

            __m256d result0_1  = _mm256_set1_pd(0);
            __m256d result1_1  = _mm256_set1_pd(0);
            __m256d result2_1  = _mm256_set1_pd(0);
            __m256d result3_1  = _mm256_set1_pd(0);
            __m256d result4_1  = _mm256_set1_pd(0);
            __m256d result5_1  = _mm256_set1_pd(0);

            for(int k = 0; k < K; k++){

                mB0 = _mm256_loadu_pd(&B[MAT_SIZE * k + j + 4 * 0]);
                mB1 = _mm256_loadu_pd(&B[MAT_SIZE * k + j + 4 * 1]);


                mA0 = _mm256_set1_pd(A[k + (i + 0) * MAT_SIZE]);    // Load float @ A's col k, row m+0 into reg
                mA1 = _mm256_set1_pd(A[k + (i + 1) * MAT_SIZE]);    // Load float @ A's col k, row m+1

                result0_0 = _mm256_add_pd(result0_0,_mm256_mul_pd(mB0,mA0)); // Getting one column of C
                result0_1 = _mm256_add_pd(result0_1,_mm256_mul_pd(mB1,mA0)); // Getting one column of C appended to the previous
                result1_0 = _mm256_add_pd(result1_0,_mm256_mul_pd(mB0,mA1)); // Getting one column of C
                result1_1 = _mm256_add_pd(result1_1,_mm256_mul_pd(mB1,mA1)); // Getting one column of C

                mA0 = _mm256_set1_pd(A[k + (i + 2) * MAT_SIZE]);    // Load float @ A's col k, row m+0 into reg
                mA1 = _mm256_set1_pd(A[k + (i + 3) * MAT_SIZE]);    // Load float @ A's col k, row m+1


                result2_0 = _mm256_add_pd(result2_0,_mm256_mul_pd(mB0,mA0)); // Getting one column of C
                result2_1 = _mm256_add_pd(result2_1,_mm256_mul_pd(mB1,mA0)); // Getting one column of C appended to the previous
                result3_0 = _mm256_add_pd(result3_0,_mm256_mul_pd(mB0,mA1)); // Getting one column of C
                result3_1 = _mm256_add_pd(result3_1,_mm256_mul_pd(mB1,mA1)); // Getting one column of C

                mA0 = _mm256_set1_pd(A[k + (i + 4) * MAT_SIZE]);    // Load float @ A's col k, row m+0 into reg
                mA1 = _mm256_set1_pd(A[k + (i + 5) * MAT_SIZE]);    // Load float @ A's col k, row m+1


                result4_0 = _mm256_add_pd(result4_0,_mm256_mul_pd(mB0,mA0)); // Getting one column of C
                result4_1 = _mm256_add_pd(result4_1,_mm256_mul_pd(mB1,mA0)); // Getting one column of C appended to the previous
                result5_0 = _mm256_add_pd(result5_0,_mm256_mul_pd(mB0,mA1)); // Getting one column of C
                result5_1 = _mm256_add_pd(result5_1,_mm256_mul_pd(mB1,mA1)); // Getting one column of C
            }


            mA0 =  _mm256_loadu_pd(&C[(i + 0)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result0_0);_mm256_storeu_pd(&C[(i + 0)*MAT_SIZE+j + 0*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 1)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result1_0);_mm256_storeu_pd(&C[(i + 1)*MAT_SIZE+j + 0*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 2)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result2_0);_mm256_storeu_pd(&C[(i + 2)*MAT_SIZE+j + 0*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 3)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result3_0);_mm256_storeu_pd(&C[(i + 3)*MAT_SIZE+j + 0*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 4)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result4_0);_mm256_storeu_pd(&C[(i + 4)*MAT_SIZE+j + 0*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 5)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result5_0);_mm256_storeu_pd(&C[(i + 5)*MAT_SIZE+j + 0*4],mA0);

            mA0 =  _mm256_loadu_pd(&C[(i + 0)*MAT_SIZE+j + 1*4]);mA0 = _mm256_add_pd(mA0,result0_1);_mm256_storeu_pd(&C[(i + 0)*MAT_SIZE+j + 1*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 1)*MAT_SIZE+j + 1*4]);mA0 = _mm256_add_pd(mA0,result1_1);_mm256_storeu_pd(&C[(i + 1)*MAT_SIZE+j + 1*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 2)*MAT_SIZE+j + 1*4]);mA0 = _mm256_add_pd(mA0,result2_1);_mm256_storeu_pd(&C[(i + 2)*MAT_SIZE+j + 1*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 3)*MAT_SIZE+j + 1*4]);mA0 = _mm256_add_pd(mA0,result3_1);_mm256_storeu_pd(&C[(i + 3)*MAT_SIZE+j + 1*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 4)*MAT_SIZE+j + 1*4]);mA0 = _mm256_add_pd(mA0,result4_1);_mm256_storeu_pd(&C[(i + 4)*MAT_SIZE+j + 1*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 5)*MAT_SIZE+j + 1*4]);mA0 = _mm256_add_pd(mA0,result5_1);_mm256_storeu_pd(&C[(i + 5)*MAT_SIZE+j + 1*4],mA0);






        }
    }
}
int main() {

    srand((unsigned int)0x100);

    const int MAT_SIZE = 24;
    double A[MAT_SIZE*MAT_SIZE] __attribute__((aligned(64)));
    double B[MAT_SIZE*MAT_SIZE] __attribute__((aligned(256)));
    double C[MAT_SIZE*MAT_SIZE]__attribute__((aligned(64)));

    double * vec4 = new double[MAT_SIZE*MAT_SIZE];


    for(int i = 0; i < MAT_SIZE*MAT_SIZE; i++){
        A[i] = float(rand()%100) / 100.0;//drand48();double(i%MAT_SIZE);//double(rand()%100) / 100.0;//drand48();;
        B[i] = float(rand()%100) / 100.0;//drand48();
	    C[i] = 0;
	    vec4[i] = 0;
    }

//    for(int i = 0; i < MAT_SIZE*MAT_SIZE; i++){
//        std::cout << "STarting = " << i << "\n";
//         *((__m256d*) (&C[i]));
//        std::cout << "Ending = " << i << "\n";
//    }
    int BLOCK_SIZE = 28;

    int stride = MAT_SIZE/BLOCK_SIZE;




    /* For each block-row of A */
    for (int i = 0; i < stride; i ++) {
        /* For each block-column of B */
        for (int j = 0; j < stride; j ++) {
            /* Accumulate block dgemms into block of C */
            for (int k = 0; k < stride; k++) {
                /* Correct block dimensions if block "goes off edge of" the matrix */
                /* Perform individual block dgemm */
                do_block(MAT_SIZE, BLOCK_SIZE, BLOCK_SIZE , BLOCK_SIZE, A + i * MAT_SIZE*BLOCK_SIZE + k*BLOCK_SIZE,
                         B + k * MAT_SIZE*BLOCK_SIZE + j*BLOCK_SIZE, C + i * MAT_SIZE*BLOCK_SIZE + j*BLOCK_SIZE);
            }
        }
    }
//    assert(false);




    for(int i = 0; i <MAT_SIZE; i++){
        for(int j = 0; j <MAT_SIZE; j++){
            for(int k = 0; k <MAT_SIZE; k++){
                vec4[i*MAT_SIZE + k] += A[i*MAT_SIZE+j]*B[j*MAT_SIZE + k];
            }

        }
    }
    bool check = true;
    for(int i = 0; i < MAT_SIZE*MAT_SIZE; i++){
        if(fabs(vec4[i] - C[i]) > 1E-6){
            std::cout << vec4[i] << " " << C[i] << "\n";
            check = false;
        }
    }
    std::cout << check << "\n";
//    delete vec2;
//    delete vec3;
    delete vec4;
}
