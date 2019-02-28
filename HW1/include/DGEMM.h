//
// Created by maksbh on 2/16/19.
//

#include <iostream>
#include <immintrin.h>
#include <assert.h>
//#include <avxintrin.h>
//#include <emmintrin.h>
#include <math.h>
#ifndef TUTORIALS_DGEMM_H
#define TUTORIALS_DGEMM_H

void multiply_naive(int N, double * A, double *B, double *C);

void do_block(int MAT_SIZE, int M, int N, int K, double* A, double* B, double* C);
void multiply4x4(int MAT_SIZE, int M, int N, int K, double* A, double* B, double* C );

void packA(double * PackA, const double * A, const int MatSize, const int BlockSizeI, const int BlockSizeK,
           const int rowStride, const int  starti, const int startk);

void packB(double * PackB, const double * B, const int MatSize, const int BlockSizeJ, const int BlockSizeK,
           const int columnStride, const int  startj, const int startk );

void dgemm_row(const int N, const double * A, const double * B, double *C );

void dgemm_col(const int N, const double * A, const double * B, double *C );

#endif //TUTORIALS_DGEMM_H
