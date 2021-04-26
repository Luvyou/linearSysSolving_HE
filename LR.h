/*
 * @Descripttion: 
 * @version: 
 * @Author: LvU
 * @Date: 2020-10-13 10:53:25
 * @LastEditors: LvU
 * @LastEditTime: 2020-10-15 09:14:52
 */
#ifndef LR_H_
#define LR_H_

#include "../src/HEAAN.h"
#include "VectorUtils.h"
#include <fstream>

class LR{
public:
    static void leastSqaure(double** A, double* b, double * x, int size);
    static void leastSqaureHE(long logq, long logp, long logn, long dim_m, double** A, double* b, double * x, ofstream& outfile);
    static void leastSqaureHE(long logq, long logp, long logn, long dim_m, ofstream& outfile);
    static void leastSqaureHE(long logq, long logp, long logn, long dim_m, long models, ofstream& outfile);
    
};

#endif