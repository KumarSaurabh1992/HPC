//
// Created by maksbh on 2/16/19.
//
#include "../include/DGEMM.h"
void multiply_naive(int N, double * A, double *B, double *C){
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < N; k++) {
                C[i + N * j] += A[i + N * k] * B[k + N * j];
            }
        }
    }
}

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
void multiply4x4(int MAT_SIZE, int M, int N, int K, double* A, double* B, double* C ){

    for(int i = 0; i < M;i+=4){
        for(int j = 0; j < N;j=j+4){
            __m256d mA0, mA1,mB0;
            __m256d result0_0  = _mm256_set1_pd(0);
            __m256d result1_0  = _mm256_set1_pd(0);
            __m256d result2_0  = _mm256_set1_pd(0);
            __m256d result3_0  = _mm256_set1_pd(0);
            for(int k = 0; k < K; k++){
                mB0 = _mm256_loadu_pd(&B[MAT_SIZE * k + j + 4*0]);

                mA0 = _mm256_set1_pd(A[k + (i + 0) * MAT_SIZE]);    // Load float @ A's col k, row m+0 into reg
                mA1 = _mm256_set1_pd(A[k + (i + 1) * MAT_SIZE]);    // Load float @ A's col k, row m+1

                result0_0 = _mm256_add_pd(result0_0,_mm256_mul_pd(mB0,mA0)); // Getting one column of C
                result1_0 = _mm256_add_pd(result1_0,_mm256_mul_pd(mB0,mA1)); // Getting one column of C

                mA0 = _mm256_set1_pd(A[k + (i + 2) * MAT_SIZE]);    // Load float @ A's col k, row m+0 into reg
                mA1 = _mm256_set1_pd(A[k + (i + 3) * MAT_SIZE]);    // Load float @ A's col k, row m+1


                result2_0 = _mm256_add_pd(result2_0,_mm256_mul_pd(mB0,mA0)); // Getting one column of C
                result3_0 = _mm256_add_pd(result3_0,_mm256_mul_pd(mB0,mA1)); // Getting one column of C
            }

            mA0 =  _mm256_loadu_pd(&C[(i + 0)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result0_0);_mm256_storeu_pd(&C[(i + 0)*MAT_SIZE+j + 0*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 1)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result1_0);_mm256_storeu_pd(&C[(i + 1)*MAT_SIZE+j + 0*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 2)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result2_0);_mm256_storeu_pd(&C[(i + 2)*MAT_SIZE+j + 0*4],mA0);
            mA0 =  _mm256_loadu_pd(&C[(i + 3)*MAT_SIZE+j + 0*4]);mA0 = _mm256_add_pd(mA0,result3_0);_mm256_storeu_pd(&C[(i + 3)*MAT_SIZE+j + 0*4],mA0);



        }
    }
}

/** A is rowStride stride A(i,k)**/
void packA(double * PackA, const double * A, const int MatSize, const int BlockSizeI, const int BlockSizeK,
        const int rowStride, const int  starti, const int startk ){

    const int realBlockSize = std::min(BlockSizeI,MatSize - starti);
    const int num_blocks_complete = realBlockSize/rowStride;
    const int num_blocks_incomplete = realBlockSize%rowStride;

    for(int count = 0; count < num_blocks_complete; count++){
        for(int i = 0; i < rowStride; i++){
            for(int k = 0; k < std::min(BlockSizeK,MatSize - startk);k++){
                int packed_id = count * BlockSizeK*rowStride + k*rowStride + i;
                int mat_id = (startk*MatSize +  starti) + (count * rowStride + i) + k * MatSize ;
                PackA[packed_id] =
                        /*A[(starti + count * rowStride + i)*MatSize + startk + k];*/
                    A[mat_id ];
                    //  std::cout << "( "  << packed_id << " to , " << mat_id << ") = " << " " << A[mat_id]  << "\n";

            }
        }
    }

    if(num_blocks_incomplete > 0){
    //     for(int k = 0; k < std::min(BlockSizeK,MatSize - startk);k++){
    //         for(int i = 0; i < num_blocks_incomplete;i++){
    //              int packed_id = num_blocks_complete * BlockSizeI*rowStride + i*rowStride + k;
    //             int mat_id = (startk*MatSize +  starti) + (num_blocks_complete * rowStride + k) * MatSize + i ;
    //             PackA[packed_id] =
    //                     /*A[(starti + count * rowStride + i)*MatSize + startk + k];*/
    //                 A[mat_id];
    //                 //  std::cout << "In( "  << packed_id << " to , " << mat_id << ") = " << " " << A[mat_id]  << "\n";
    //         }
    //         for(int k = num_blocks_incomplete; k < rowStride; k++){
    //             PackA[num_blocks_complete * BlockSizeI*rowStride + i*rowStride + k] = 0;
    //         }
        // }
       
            for(int k = 0; k < std::min(BlockSizeK,MatSize - startk); k++){
                 for(int i = 0; i < num_blocks_incomplete;i++){
                int packed_id = num_blocks_complete * BlockSizeK*rowStride + k*rowStride + i;
                int mat_id = (startk*MatSize +  starti) + (num_blocks_complete * rowStride + i) + k * MatSize ;
                PackA[packed_id] =
                        A[mat_id];
            }
            for(int i = num_blocks_incomplete; i < std::min(BlockSizeK,MatSize - startk); i++){
                PackA[num_blocks_complete * BlockSizeK*rowStride + k*rowStride + i] = 0;
            }
        }
    }
        }





/** B is colStride stride B(k,j)**/
void packB(double * PackB, const double * B, const int MatSize, const int BlockSizeJ, const int BlockSizeK,
           const int colStride, const int  startj, const int startk ){
       
    const int realBlockSize = std::min(BlockSizeJ,MatSize - startj);
    const int num_blocks_complete = realBlockSize/colStride;
    const int num_blocks_incomplete = realBlockSize%colStride;

    for(int count = 0; count < num_blocks_complete; count++){

            for(int k = 0; k < std::min(BlockSizeK,MatSize - startk);k++){
                for(int j = 0; j < colStride; j++){
                int packed_id = count * BlockSizeK*colStride + k*colStride + j;
                int mat_id = (startj*MatSize +  startk) + ( j + colStride*count) * MatSize + k ;
                PackB[packed_id] = B[mat_id];

            //    std::cout << "( "  << packed_id << " to , " << mat_id << ")\n";

            }
        }
    }
    
    if(num_blocks_incomplete > 0)
    for(int k = 0; k < std::min(BlockSizeK,MatSize - startk);k++){
        int packed_id;
        for(int j = 0; j < num_blocks_incomplete; j++){
            packed_id = num_blocks_complete * BlockSizeK*colStride + k*colStride + j;
            int mat_id = (startj*MatSize +  startk) + ( j + colStride*num_blocks_complete) * MatSize + k ;
            // std::cout << "Non zeros at " << packed_id << " " << PackB[packed_id] << " "<< B[mat_id] << " "<< i << "\n";
            PackB[packed_id] = B[mat_id];

        }
        for(int j = num_blocks_incomplete; j < colStride; j++){
        
	        // std::cout << "Zeros at id =  " << packed_id << "\n";
            PackB[num_blocks_complete * BlockSizeK*colStride + k*colStride + j] = 0;
        }
    }
}

