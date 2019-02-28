#include "include/DGEMM.h"
using namespace std;


int main() {

    srand((unsigned int)0x100);

    const int MAT_SIZE = 200;
    const int N = MAT_SIZE;
    double A[MAT_SIZE*MAT_SIZE];// __attribute__((aligned(64)));
    double B[MAT_SIZE*MAT_SIZE];// __attribute__((aligned(256)));
    double C[MAT_SIZE*MAT_SIZE];//__attribute__((aligned(64)));

    double * vec4 = new double[MAT_SIZE*MAT_SIZE];

    
   
    for(int i = 0; i < MAT_SIZE*MAT_SIZE; i++){
        A[i] = double(rand()%100) / 100.0;//drand48();double(i%MAT_SIZE);//double(rand()%100) / 100.0;//drand48();;
        B[i] = double(rand()%100) / 100.0;//drand48();
	    C[i] = 0;
	    vec4[i] = 0;
    }

    // std::cout << "-------------------B-------------------\n";
    // for(int i = 0; i < MAT_SIZE; i++){
    //     for(int j = 0; j < MAT_SIZE; j++){
    //         std::cout << B[i + j*MAT_SIZE] << " " ;
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << "-------------------A-------------------\n";
    // for(int i = 0; i < MAT_SIZE; i++){
    //     for(int j = 0; j < MAT_SIZE; j++){
    //         std::cout << A[i + j*MAT_SIZE] << " " ;
    //     }
    //     std::cout << "\n";
    // }
    //     std::cout << "-------------------Over-------------------\n";


    dgemm_col(MAT_SIZE,A,B,C);

    for(int i = 0; i <MAT_SIZE; i++){
        for(int j = 0; j <MAT_SIZE; j++){
            for(int k = 0; k <MAT_SIZE; k++){
                vec4[i + k*MAT_SIZE] += A[i+j*MAT_SIZE]*B[j + MAT_SIZE* k];
            }

        }
    }

//         for(int i = 0; i <MAT_SIZE; i++){
//         for(int j = 0; j <MAT_SIZE; j++){
//             std::cout << vec4[i+MAT_SIZE*j] << " ";
//         }
//         std::cout << "\n";
// }
//     std::cout << "------------------Correct -----------------------------------\n";
//     for(int i = 0; i <MAT_SIZE; i++){
//         for(int j = 0; j <MAT_SIZE; j++){
//             std::cout << C[i+MAT_SIZE*j] << " ";
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
