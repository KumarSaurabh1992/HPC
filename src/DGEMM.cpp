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

void dgemm_row(const int N, const double * A, const double * B, double *C ){

    const int BS_J = 24;
    const int BS_K = 24;
    const int BS_I = 24;
    const int BC_J = 8;
    const int BC_I = 6;

     double *Apack = (double *)_mm_malloc(BS_I * BS_K * sizeof(double), 64);
    double *Bpack = (double *)_mm_malloc(BS_J * BS_K * sizeof(double), 64);

     for (int j_loop = 0; j_loop < N; j_loop += BS_J) {
    for (int k_loop = 0; k_loop < N; k_loop += BS_K) {
      packB(Bpack,B,N,BS_J, BS_K, BC_J, j_loop, k_loop);
      for (int i_loop = 0; i_loop < N; i_loop += BS_I) {
        packA(Apack,A,N,BS_I, BS_K, BC_I, i_loop, k_loop);
        int jc_loc = 0;
        for (int jc = j_loop; jc < j_loop +std::min(BS_J, N - j_loop); jc += BC_J, jc_loc += 1) {
          double *Blocal = &Bpack[jc_loc * BS_K * BC_J];
          int ic_loc = 0;
          for (int ic = i_loop; ic < i_loop +std::min(BS_I, N - i_loop); ic += BC_I, ic_loc += 1) {
            double *Alocal = &Apack[ic_loc * BS_K * BC_I];
            __m256d result0_0 = _mm256_setzero_pd();
            __m256d result1_0 = _mm256_setzero_pd();
            __m256d result2_0 = _mm256_setzero_pd();
            __m256d result3_0 = _mm256_setzero_pd();
            __m256d result4_0 = _mm256_setzero_pd();
            __m256d result5_0 = _mm256_setzero_pd();
            __m256d result0_1 = _mm256_setzero_pd();
            __m256d result1_1 = _mm256_setzero_pd();
            __m256d result2_1 = _mm256_setzero_pd();
            __m256d result3_1 = _mm256_setzero_pd();
            __m256d result4_1 = _mm256_setzero_pd();
            __m256d result5_1 = _mm256_setzero_pd();
            __m256d mA0, mA1, mB0, mB1;
            for (int k_loc = 0; k_loc <std::min(BS_K, N - k_loop); k_loc += 1) {

               
              mB0 = _mm256_loadu_pd(&Blocal[k_loc * BC_J + 0]);
              mB1 = _mm256_loadu_pd(&Blocal[k_loc * BC_J + 4]);
               
              mA0 = _mm256_set1_pd(Alocal[k_loc * BC_I + 0]);
              mA1 = _mm256_set1_pd(Alocal[k_loc * BC_I + 1]);
              result0_0 = _mm256_add_pd(result0_0, _mm256_mul_pd(mA0, mB0 ));
              result0_1 = _mm256_add_pd(result0_1, _mm256_mul_pd(mA0, mB1 ));
              result1_0 = _mm256_add_pd(result1_0, _mm256_mul_pd(mA1, mB0 ));
              result1_1 = _mm256_add_pd(result1_1, _mm256_mul_pd(mA1, mB1 ));

              mA0 = _mm256_set1_pd(Alocal[k_loc * BC_I + 2]);
              mA1 = _mm256_set1_pd(Alocal[k_loc * BC_I + 3]);
              result2_0 = _mm256_add_pd(result2_0, _mm256_mul_pd(mA0, mB0 ));
              result2_1 = _mm256_add_pd(result2_1, _mm256_mul_pd(mA0, mB1 ));
              result3_0 = _mm256_add_pd(result3_0, _mm256_mul_pd(mA1, mB0 ));
              result3_1 = _mm256_add_pd(result3_1, _mm256_mul_pd(mA1, mB1 ));

              mA0 = _mm256_set1_pd(Alocal[k_loc * BC_I + 4]);
              mA1 = _mm256_set1_pd(Alocal[k_loc * BC_I + 5]);
              result4_0 = _mm256_add_pd(result4_0, _mm256_mul_pd(mA0, mB0 ));
              result4_1 = _mm256_add_pd(result4_1, _mm256_mul_pd(mA0, mB1 ));
              result5_0 = _mm256_add_pd(result5_0, _mm256_mul_pd(mA1, mB0 ));
              result5_1 = _mm256_add_pd(result5_1, _mm256_mul_pd(mA1, mB1 ));
              
              
            }
            
            double buffer[BC_I * BC_J];// __attribute__((aligned(32)));
            _mm256_storeu_pd(&buffer[0 * BC_J + 0], result0_0);
            _mm256_storeu_pd(&buffer[0 * BC_J + 4], result0_1);
            _mm256_storeu_pd(&buffer[1 * BC_J + 0], result1_0);
            _mm256_storeu_pd(&buffer[1 * BC_J + 4], result1_1);
            _mm256_storeu_pd(&buffer[2 * BC_J + 0], result2_0);
            _mm256_storeu_pd(&buffer[2 * BC_J + 4], result2_1);
            _mm256_storeu_pd(&buffer[3 * BC_J + 0], result3_0);
            _mm256_storeu_pd(&buffer[3 * BC_J + 4], result3_1);
            _mm256_storeu_pd(&buffer[4 * BC_J + 0], result4_0);
            _mm256_storeu_pd(&buffer[4 * BC_J + 4], result4_1);
            _mm256_storeu_pd(&buffer[5 * BC_J + 0], result5_0);
            _mm256_storeu_pd(&buffer[5 * BC_J + 4], result5_1);
            
            for (int ii = 0; ii <std::min(BC_I, N - ic); ii++) {
              for (int jj = 0; jj <std::min(BC_J, N - jc); jj++) {
                C[ic + ii + N * (jc + jj)] += buffer[ii * BC_J + jj];
              }
            }
          }
        }
      }
    }
  }
  
  _mm_free(Apack);
  _mm_free(Bpack);

}

void dgemm_col(const int N, const double * A, const double * B, double *C ){

    const int BS_J = 96;
    const int BS_K = 48;
    const int BS_I = 24;
    const int BC_J = 6;
    const int BC_I = 8;

     double *Apack = (double *)_mm_malloc(BS_I * BS_K * sizeof(double), 64);
    double *Bpack = (double *)_mm_malloc(BS_J * BS_K * sizeof(double), 64);

     for (int j_loop = 0; j_loop < N; j_loop += BS_J) {
    for (int k_loop = 0; k_loop < N; k_loop += BS_K) {
      // pack B
      packB(Bpack,B,N,BS_J, BS_K, BC_J, j_loop, k_loop);
      for (int i_loop = 0; i_loop < N; i_loop += BS_I) {
        // pack A
        packA(Apack,A,N,BS_I, BS_K, BC_I, i_loop, k_loop);
        int jc_loc = 0;
        for (int jc = j_loop; jc < j_loop +std::min(BS_J, N - j_loop); jc += BC_J, jc_loc += 1) {
          double *Blocal = &Bpack[jc_loc * BS_K * BC_J];
          int ic_loc = 0;
          for (int ic = i_loop; ic < i_loop +std::min(BS_I, N - i_loop); ic += BC_I, ic_loc += 1) {
            double *Alocal = &Apack[ic_loc * BS_K * BC_I];
            __m256d result0_0 = _mm256_setzero_pd();
            __m256d result1_0 = _mm256_setzero_pd();
            __m256d result2_0 = _mm256_setzero_pd();
            __m256d result3_0 = _mm256_setzero_pd();
            __m256d result4_0 = _mm256_setzero_pd();
            __m256d result5_0 = _mm256_setzero_pd();
            __m256d result0_1 = _mm256_setzero_pd();
            __m256d result1_1 = _mm256_setzero_pd();
            __m256d result2_1 = _mm256_setzero_pd();
            __m256d result3_1 = _mm256_setzero_pd();
            __m256d result4_1 = _mm256_setzero_pd();
            __m256d result5_1 = _mm256_setzero_pd();
            __m256d mA0, mA1, mB0, mB1;
            for (int k_loc = 0; k_loc <std::min(BS_K, N - k_loop); k_loc += 1) {

              mB0 = _mm256_loadu_pd(&Alocal[k_loc * BC_I + 0]);
              mB1 = _mm256_loadu_pd(&Alocal[k_loc * BC_I + 4]);
            
               
              mA0 = _mm256_set1_pd(Blocal[k_loc * BC_J + 0]);
              mA1 = _mm256_set1_pd(Blocal[k_loc * BC_J + 1]);
              result0_0 = _mm256_add_pd(result0_0, _mm256_mul_pd(mA0, mB0 ));
              result0_1 = _mm256_add_pd(result0_1, _mm256_mul_pd(mA0, mB1 ));
              result1_0 = _mm256_add_pd(result1_0, _mm256_mul_pd(mA1, mB0 ));
              result1_1 = _mm256_add_pd(result1_1, _mm256_mul_pd(mA1, mB1 ));

              mA0 = _mm256_set1_pd(Blocal[k_loc * BC_J + 2]);
              mA1 = _mm256_set1_pd(Blocal[k_loc * BC_J + 3]);
              result2_0 = _mm256_add_pd(result2_0, _mm256_mul_pd(mA0, mB0 ));
              result2_1 = _mm256_add_pd(result2_1, _mm256_mul_pd(mA0, mB1 ));
              result3_0 = _mm256_add_pd(result3_0, _mm256_mul_pd(mA1, mB0 ));
              result3_1 = _mm256_add_pd(result3_1, _mm256_mul_pd(mA1, mB1 ));


              mA0 = _mm256_set1_pd(Blocal[k_loc * BC_J + 4]);
              mA1 = _mm256_set1_pd(Blocal[k_loc * BC_J + 5]);
              result4_0 = _mm256_add_pd(result4_0, _mm256_mul_pd(mA0, mB0 ));
              result4_1 = _mm256_add_pd(result4_1, _mm256_mul_pd(mA0, mB1 ));
              result5_0 = _mm256_add_pd(result5_0, _mm256_mul_pd(mA1, mB0 ));
              result5_1 = _mm256_add_pd(result5_1, _mm256_mul_pd(mA1, mB1 ));
            }
            
            double buffer[BC_I * BC_J] ;//__attribute__((aligned(32)));
            _mm256_storeu_pd(&buffer[0 * BC_I + 0], result0_0);
            _mm256_storeu_pd(&buffer[0 * BC_I + 4], result0_1);
            _mm256_storeu_pd(&buffer[1 * BC_I + 0], result1_0);
            _mm256_storeu_pd(&buffer[1 * BC_I + 4], result1_1);
            _mm256_storeu_pd(&buffer[2 * BC_I + 0], result2_0);
            _mm256_storeu_pd(&buffer[2 * BC_I + 4], result2_1);
            _mm256_storeu_pd(&buffer[3 * BC_I + 0], result3_0);
            _mm256_storeu_pd(&buffer[3 * BC_I + 4], result3_1);
            _mm256_storeu_pd(&buffer[4 * BC_I + 0], result4_0);
            _mm256_storeu_pd(&buffer[4 * BC_I + 4], result4_1);
            _mm256_storeu_pd(&buffer[5 * BC_I + 0], result5_0);
            _mm256_storeu_pd(&buffer[5 * BC_I + 4], result5_1);
        
            for (int ii = 0; ii <std::min(BC_I, N - ic); ii++) {
              for (int jj = 0; jj <std::min(BC_J, N - jc); jj++) {
                C[ic + ii + N * (jc + jj)] += buffer[ii  + jj*BC_I];
              }
            }
          }
        }
      }
    }
  }
  
  _mm_free(Apack);
  _mm_free(Bpack);

}

