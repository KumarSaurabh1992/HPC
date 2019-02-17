#include "include/DGEMM.h"
using namespace std;



int main() {

    srand((unsigned int)0x100);

    const int MAT_SIZE = 5;
    const int N = MAT_SIZE;
    double A[MAT_SIZE*MAT_SIZE];// __attribute__((aligned(64)));
    double B[MAT_SIZE*MAT_SIZE];// __attribute__((aligned(256)));
    double C[MAT_SIZE*MAT_SIZE];//__attribute__((aligned(64)));

    double * vec4 = new double[MAT_SIZE*MAT_SIZE];

    const int BS_J = 24;
    const int BS_K = 24;
    const int BS_I = 24;
    const int BC_J = 8;
    const int BC_I = 6;
    double *Apack = (double *)_mm_malloc(BS_I * BS_K * sizeof(double), 64);
    double *Bpack = (double *)_mm_malloc(BS_J * BS_K * sizeof(double), 64);
    for(int i = 0; i < MAT_SIZE*MAT_SIZE; i++){
        A[i] = float(rand()%100) / 100.0;//drand48();double(i%MAT_SIZE);//double(rand()%100) / 100.0;//drand48();;
        B[i] = float(rand()%100) / 100.0;//drand48();
	    C[i] = 0;
	    vec4[i] = 0;
    }
    std::cout << "-------------------B-------------------\n";
    for(int i = 0; i < MAT_SIZE; i++){
        for(int j = 0; j < MAT_SIZE; j++){
            std::cout << B[i + j*MAT_SIZE] << " " ;
        }
        std::cout << "\n";
    }
    std::cout << "-------------------A-------------------\n";
    for(int i = 0; i < MAT_SIZE; i++){
        for(int j = 0; j < MAT_SIZE; j++){
            std::cout << A[i + j*MAT_SIZE] << " " ;
        }
        std::cout << "\n";
    }
        std::cout << "-------------------Over-------------------\n";



     for (int jb = 0; jb < N; jb += BS_J) {
    for (int kb = 0; kb < N; kb += BS_K) {
      // pack B
      packB(Bpack,B,N,BS_J, BS_K, BC_J, jb, kb);
      for (int ib = 0; ib < N; ib += BS_I) {
        // pack A
        packA(Apack,A,N,BS_I, BS_K, BC_I, ib, kb);
        int jc_loc = 0;
        for (int jc = jb; jc < jb + min(BS_J, N - jb); jc += BC_J, jc_loc += 1) {
          double *Blocal = &Bpack[jc_loc * BS_K * BC_J];
          int ic_loc = 0;
          for (int ic = ib; ic < ib + min(BS_I, N - ib); ic += BC_I, ic_loc += 1) {
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
            for (int k_loc = 0; k_loc < min(BS_K, N - kb); k_loc += 1) {

               
              // Cleanup
              mB0 = _mm256_loadu_pd(&Blocal[k_loc * BC_J + 0]);
              mB1 = _mm256_loadu_pd(&Blocal[k_loc * BC_J + 4]);
               std::cout << "B0 = "  << mB0[0] << " "<< mB0[1] << " "<< mB0[2] << " "<< mB0[3] << "\n";
               std::cout << "B1 = " << mB1[0] << " "<< mB1[1] << " "<< mB1[2] << " "<< mB1[3] << "\n";

               std::cout << "Alocal = "  << Alocal[0] << " "<< Alocal[1] << " "<< Alocal[2] << " "<< Alocal[3] 
               <<  " " << Alocal[4] << " "<< Alocal[5] << "\n";
               
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
            
            double buffer[BC_I * BC_J] __attribute__((aligned(32)));
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
            for (int ii = 0; ii < min(BC_I, N - ic); ii++) {
              for (int jj = 0; jj < min(BC_J, N - jc); jj++) {
                //   std::cout << "Value at " << 
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


    for(int i = 0; i <MAT_SIZE; i++){
        for(int j = 0; j <MAT_SIZE; j++){
            for(int k = 0; k <MAT_SIZE; k++){
                vec4[i + k*MAT_SIZE] += A[i+j*MAT_SIZE]*B[j + MAT_SIZE* k];
            }

        }
    }
//     std::cout << "------------------Correct -----------------------------------\n";
//     for(int i = 0; i <MAT_SIZE; i++){
//         for(int j = 0; j <MAT_SIZE; j++){
//             std::cout << vec4[i+MAT_SIZE*j] << " ";
//         }
//         std::cout << "\n";
// }
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
