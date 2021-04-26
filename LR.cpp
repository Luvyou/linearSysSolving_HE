/*
 * @Descripttion: 
 * @version: 
 * @Author: LvU
 * @Date: 2020-10-13 10:53:26
 * @LastEditors: LvU
 * @LastEditTime: 2020-10-15 19:08:16
 */

#include "LR.h"

void LR::leastSqaureHE(long logq, long logp, long logn, long dim_m, double** AA,double* bb, double* xx, ofstream& outfile){
    srand(time(NULL));
	SetNumThreads(16);
	TimeUtils timeutils1;
	TimeUtils timeutils2;
	TimeUtils timeccr;
	TimeUtils timeccb;

	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addLeftRotKeys(secretKey);
	SchemeAlgo algo(scheme);

	long n = (1 << logn);
	double** A = new double*[dim_m];
	double* b = new double[n];
    double** x_m = new double*[dim_m];
	double* x = new double[dim_m];

	cout << endl << "********************************************************" << endl;
	cout << "print Params:.............." << endl;
	cout << "logq:=" << logq << " logp:=" << logp << " dimension:=" << dim_m << endl;
	cout << "loading Matrix A...." << endl;
	cout << "A:=" << endl;
	StringUtils::showVec(A[0], 121);
	for(long i = 0; i < dim_m; ++ i){
		A[i] = new double[n];
		for(long k = 0; k < dim_m; ++ k) A[i][k] = AA[i][k];
		for(long j = dim_m; j < n; ++ j) A[i][j] = 0;
		StringUtils::showVec(A[i], n);
	}
	for(long k = 0; k < dim_m; ++ k) b[k] = bb[k];
	for(long j = dim_m; j < n; ++ j) b[j] = 0;

	outfile << endl << "********************************************************" << endl;
	outfile << "print Params:.............." << endl;
	outfile << "logq:=" << logq << " logp:=" << logp << " dimension:=" << dim_m << endl;
	outfile << "loading Matrix A...." << endl;
	outfile << "A:=" << endl;
	for(long i = 0; i < dim_m; ++ i){
		for(long j = dim_m; j < n; ++ j) A[i][j] = 0;
		StringUtils::writeVec(A[i], n, outfile);
	}
	outfile << "loading vector b...." << endl;
	for(long i = dim_m; i < n; ++ i) b[i] = 0;
	outfile << "b：= " << endl; 
	//StringUtils::showVec(b, n);
	StringUtils::writeVec(b, n, outfile);
    
	timeutils2.start("leastSqaure..........");
    LR::leastSqaure(A, b, x, dim_m);
	timeutils2.stop("finish leastSqaure..........");
	outfile << "Time consuming:= " << timeutils2.timeElapsed << "ms" << endl;

    
	
	timeutils1.start("leastSqaureHE..........");

	Ciphertext *c_A = new Ciphertext[dim_m];
	Ciphertext **c_C = new Ciphertext *[dim_m-1];
	Ciphertext c_b;
	
	for(long i = 0; i < dim_m; i ++){
	    x_m[i] = new double[dim_m];
	}
	for(long i = 1; i < dim_m; ++i ){
		c_C[i-1] = new Ciphertext[i];
	}
	scheme.encrypt(c_b, b, n, logp, logq);
	for(long i = 0; i < dim_m; ++ i){
		scheme.encrypt(c_A[i], A[i], n, logp, logq);
		//delete [] A[i];
	}
	//delete [] A;
	
	
	Ciphertext *len = new Ciphertext[dim_m];
	Ciphertext *len_inv = new Ciphertext[dim_m];
	double c = 0;
	complex<double> **lenm = new complex<double> *[dim_m];
    complex<double> **lenm_inv = new complex<double> *[dim_m];
	for(long i = 0; i < dim_m; ++ i){
		lenm_inv[i] = new complex<double>[n];
	}
	

	scheme.square(len[0], c_A[0]);
	scheme.reScaleByAndEqual(len[0], logp);
	algo.partialSlotsSumAndEqual(len[0], n);
	c = EvaluatorUtils::randomReal(10);
	scheme.multByConstAndEqual(len[0], c, logp);
	scheme.reScaleByAndEqual(len[0], logp);
	lenm[0] = scheme.decrypt(secretKey, len[0]);
	//StringUtils::showVec(lenm[0], n);

	
	//NTL_EXEC_RANGE(n, first, last)
	for(long i = 0; i < n; ++ i){
		lenm_inv[0][i].real(1./lenm[0][0].real());
	}
    //NTL_EXEC_RANGE_END

	//StringUtils::showVec(lenm_inv[0], n);
	scheme.encrypt(len_inv[0], lenm_inv[0], n, logp, logq);
	scheme.multByConstAndEqual(len_inv[0], c, logp);
	lenm_inv[0] = scheme.decrypt(secretKey ,len_inv[0]);
	//StringUtils::showVec(lenm_inv[0], n);
	scheme.reScaleByAndEqual(len_inv[0], logp);
	

    Ciphertext v;
	Ciphertext temp;
	for(long i = 1; i <= dim_m-1; ++ i){
		outfile << "start calculating c_R and c_B[" << i << "]" << endl;
		//timeccb.start(" calculating c_B ");
		
		v.copy(c_A[i]);
		for(long j = 0; j <= i-1; ++ j){
			if(len_inv[j].logq < c_A[i].logq)
				scheme.modDownToAndEqual(c_A[i], len_inv[j].logq);  
			else 
				scheme.modDownToAndEqual(len_inv[j], c_A[i].logq);
			scheme.mult(temp, len_inv[j], c_A[i]);
			scheme.reScaleByAndEqual(temp, logp);

            if(temp.logq < c_A[j].logq)
				scheme.modDownToAndEqual(c_A[j], temp.logq); 
			else 
				scheme.modDownToAndEqual(temp, c_A[j].logq);
			scheme.mult(c_C[i-1][j], temp, c_A[j]);
			scheme.reScaleByAndEqual(c_C[i-1][j], logp);
			algo.partialSlotsSumAndEqual(c_C[i-1][j], n);

			if(c_C[i-1][j].logq < c_A[j].logq)
				scheme.modDownToAndEqual(c_A[j], c_C[i-1][j].logq);
			  
			else 
				scheme.modDownToAndEqual(c_C[i-1][j], c_A[j].logq);
			scheme.mult(temp, c_C[i-1][j], c_A[j]);   
			scheme.reScaleByAndEqual(temp, logp);
            
			if(v.logq < temp.logq)
				scheme.modDownToAndEqual(temp, v.logq);
			else 
				scheme.modDownToAndEqual(v, temp.logq);
			scheme.subAndEqual(v, temp);
			
		}
		//timeccb.stop("finish calculating c_B ");
		//outfile << "calculating c_B[" << i << "] time: " << timeccb.timeElapsed << "ms" << endl;
		
		c_A[i].copy(v);
		complex<double> *BB = scheme.decrypt(secretKey, c_A[i]);
		outfile << "c_B[" << i << "]:=" << endl;
		StringUtils::writeVec(BB, n, outfile);
		scheme.square(len[i], c_A[i]); 
		scheme.reScaleByAndEqual(len[i], logp);
		algo.partialSlotsSumAndEqual(len[i], n);
		c = EvaluatorUtils::randomReal(10);
	    scheme.multByConstAndEqual(len[i], c, logp);
		scheme.reScaleByAndEqual(len[i], logp);
		lenm[i] = scheme.decrypt(secretKey, len[i]);

		//NTL_EXEC_RANGE(n, first, last)
		for(long t = 0; t < n; ++ t){
			lenm_inv[i][t].real(1./lenm[i][t].real());
		}
        //NTL_EXEC_RANGE_END

		scheme.encrypt(len_inv[i], lenm_inv[i], n, logp, logq);
		scheme.multByConstAndEqual(len_inv[i], c, logp);
		scheme.reScaleByAndEqual(len_inv[i], logp);

        //outfile << "start calcucating c_R[ ][" << i << "]" << endl;
		for(long j = 0; j <= i-1; ++j){
			
			//timeccr.start("calcucating c_R");
			scheme.negateAndEqual(c_C[i-1][j]);
			for(long k = j+1; k <= i-1; ++k){

				if(c_C[k-1][j].logq < c_C[i-1][k].logq)
					scheme.modDownToAndEqual(c_C[i-1][k], c_C[k-1][j].logq);
			    else 
				scheme.modDownToAndEqual(c_C[k-1][j], c_C[i-1][k].logq);
				scheme.mult(temp, c_C[k-1][j], c_C[i-1][k]);   
				scheme.reScaleByAndEqual(temp, logp);
				if(c_C[i-1][j].logq < temp.logq)
			       scheme.modDownToAndEqual(temp, c_C[i-1][j].logq);
			    else if(c_C[i-1][j].logq > temp.logq)
			       scheme.modDownToAndEqual(c_C[i-1][j], temp.logq);
				scheme.subAndEqual(c_C[i-1][j], temp);
			}
			//timeccr.stop("finish calcucating c_R");
		}
			
	}

	outfile << "computing AX = b;" << endl;
	timeutils2.start("part2 calculating x");

	//NTL_EXEC_RANGE(dim_m, first, last)
	for(long i = 0; i < dim_m; ++i){
		c_A[i].logq > c_b.logq ? scheme.modDownToAndEqual(c_A[i], c_b.logq) : scheme.modDownToAndEqual(c_b, c_A[i].logq);
        scheme.multAndEqual(c_A[i], c_b);
		scheme.reScaleByAndEqual(c_A[i], logp);
		algo.partialSlotsSumAndEqual(c_A[i], n);
		c_A[i].logq > len_inv[i].logq ? scheme.modDownToAndEqual(c_A[i], len_inv[i].logq) : scheme.modDownToAndEqual(len_inv[i], c_A[i].logq);
		scheme.multAndEqual(c_A[i], len_inv[i]);
		scheme.reScaleByAndEqual(c_A[i], logp);
	}
	//NTL_EXEC_RANGE_END

	//NTL_EXEC_RANGE(dim_m, first, last) 此处不可并行
	for(int i = 0; i < dim_m; ++i){
		Ciphertext temp;
		for(int j = i + 1; j < dim_m; ++ j){
			c_C[j-1][i].logq < c_A[j].logq ? scheme.modDownToAndEqual(c_A[j], c_C[j-1][i].logq) : scheme.modDownToAndEqual(c_C[j-1][i], c_A[j].logq);
			scheme.mult(temp, c_C[j-1][i], c_A[j]);
			scheme.reScaleByAndEqual(temp, logp);
			c_A[i].logq < temp.logq ? scheme.modDownToAndEqual(temp, c_A[i].logq) : scheme.modDownToAndEqual(c_A[i], temp.logq);
			scheme.addAndEqual(c_A[i], temp);
		}
	}
	//NTL_EXEC_RANGE_END

	timeutils1.stop("finish leastSqaureHE..........");
	timeutils2.stop("finish part2, calculating x");
	outfile << "get result part2 x time: " << timeutils2.timeElapsed << "ms" << endl;
	
	complex<double>* mm5;
	for(int i = 0; i < dim_m; ++ i){
		mm5 = scheme.decrypt(secretKey, c_A[i]);
		outfile << "---------------------" << endl;
		for(int z = 0; z < dim_m; ++ z){
			x_m[i][z] = mm5[z].real();
			outfile << "m LR : " << i << " :" << x[i];
			outfile << "  d LR : " << i << " :" << x_m[i][z];
			outfile << "  e LR : " << i << " :" << x[i]-x_m[i][z] << endl;
		}
		xx[i] = x_m[i][0];
		outfile << "---------------------" << endl;
	}
	delete [] mm5;
	
	outfile << "X.logq=:";
	for(int i = 0 ; i < dim_m; ++ i){
		outfile << c_A[i].logq << " ";
	}
	outfile << "    consumeMUL:" << (logq - c_A[0].logq)/logp << endl;
	outfile << "HE Time consuming:= " << timeutils1.timeElapsed << "ms" << endl;
	outfile << endl << "********************************************************" << endl;
	outfile << "finish HE" << endl;
	
	delete [] len;
	delete [] len_inv;
	//delete [] x;
	//delete [] b;
	for(long i = 0; i < dim_m; ++i){
		delete [] A[i];
		delete [] x_m[i];
		delete [] lenm[i];
		delete [] lenm_inv[i];
	}
	delete [] A;
	delete [] b;
	delete [] x_m;
	delete [] lenm;
	delete [] lenm_inv;
	for(long i =0 ; i < dim_m-1; ++i) delete []c_C[i];
	delete [] c_C;
}
void LR::leastSqaureHE(long logq, long logp, long logn, long dim_m, ofstream& outfile){
    srand(time(NULL));
	SetNumThreads(16);
	TimeUtils timeutils1;
	TimeUtils timeutils2;
	TimeUtils timeccr;
	TimeUtils timeccb;

	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addLeftRotKeys(secretKey);
	SchemeAlgo algo(scheme);

	long n = (1 << logn);

	
	double **A = new double *[dim_m];
	double *b = new double[n];
    double* x = new double[dim_m];
    double** x_m = new double*[dim_m];


	outfile << endl << "********************************************************" << endl;
	outfile << "print Params:.............." << endl;
	outfile << "logq:=" << logq << " logp:=" << logp << " dimension:=" << dim_m << endl;
	outfile << "generating Matrix A...." << endl;
	//outfile << "A:=" << endl;
	for(long i = 0; i < dim_m; ++ i){
		A[i] = EvaluatorUtils::randomRealArray(n, 100.);
		for(long j = dim_m; j < n; ++ j) A[i][j] = 0;
		//StringUtils::showVec(A[i], n);
		//StringUtils::writeVec(A[i], n, outfile);
	}
	outfile << "generating vector b...." << endl;
	b = EvaluatorUtils::randomRealArray(n, 100);
	for(long i = dim_m; i < n; ++ i) b[i] = 0;
	//outfile << "b：= " << endl; 
	//StringUtils::showVec(b, n);
	//StringUtils::writeVec(b, n, outfile);
    
	timeutils2.start("leastSqaure..........");
    LR::leastSqaure(A, b, x, dim_m);
	timeutils2.stop("finish leastSqaure..........");
	outfile << "Time consuming:= " << timeutils2.timeElapsed << "ms" << endl;

    
	
	timeutils1.start("leastSqaureHE..........");

	Ciphertext *c_A = new Ciphertext[dim_m];
	Ciphertext **c_C = new Ciphertext *[dim_m-1];
	Ciphertext c_b;
	
	for(long i = 0; i < dim_m; i ++){
	    x_m[i] = new double[dim_m];
	}
	for(long i = 1; i < dim_m; ++i ){
		c_C[i-1] = new Ciphertext[i];
	}
	scheme.encrypt(c_b, b, n, logp, logq);
	for(long i = 0; i < dim_m; ++ i){
		scheme.encrypt(c_A[i], A[i], n, logp, logq);
		delete [] A[i];
	}
	delete [] A;
	
	
	Ciphertext *len = new Ciphertext[dim_m];
	Ciphertext *len_inv = new Ciphertext[dim_m];
	double c = 0;
	complex<double> **lenm = new complex<double> *[dim_m];
    complex<double> **lenm_inv = new complex<double> *[dim_m];
	for(long i = 0; i < dim_m; ++ i){
		lenm_inv[i] = new complex<double>[n];
	}
	

	scheme.square(len[0], c_A[0]);
	scheme.reScaleByAndEqual(len[0], logp);
	algo.partialSlotsSumAndEqual(len[0], n);
	c = EvaluatorUtils::randomReal(10);
	scheme.multByConstAndEqual(len[0], c, logp);
	scheme.reScaleByAndEqual(len[0], logp);
	lenm[0] = scheme.decrypt(secretKey, len[0]);
	//StringUtils::showVec(lenm[0], n);

	
	//NTL_EXEC_RANGE(n, first, last)
	for(long i = 0; i < n; ++ i){
		lenm_inv[0][i].real(1./lenm[0][0].real());
	}
    //NTL_EXEC_RANGE_END

	//StringUtils::showVec(lenm_inv[0], n);
	scheme.encrypt(len_inv[0], lenm_inv[0], n, logp, logq);
	scheme.multByConstAndEqual(len_inv[0], c, logp);
	lenm_inv[0] = scheme.decrypt(secretKey ,len_inv[0]);
	//StringUtils::showVec(lenm_inv[0], n);
	scheme.reScaleByAndEqual(len_inv[0], logp);
	

    Ciphertext v;
	Ciphertext temp;
	for(long i = 1; i <= dim_m-1; ++ i){
		outfile << "start calculating c_R and c_B[" << i << "]" << endl;
		//timeccb.start(" calculating c_B ");
		
		v.copy(c_A[i]);
		for(long j = 0; j <= i-1; ++ j){
			if(len_inv[j].logq < c_A[i].logq)
				scheme.modDownToAndEqual(c_A[i], len_inv[j].logq);  
			else 
				scheme.modDownToAndEqual(len_inv[j], c_A[i].logq);
			scheme.mult(temp, len_inv[j], c_A[i]);
			scheme.reScaleByAndEqual(temp, logp);

            if(temp.logq < c_A[j].logq)
				scheme.modDownToAndEqual(c_A[j], temp.logq); 
			else 
				scheme.modDownToAndEqual(temp, c_A[j].logq);
			scheme.mult(c_C[i-1][j], temp, c_A[j]);
			scheme.reScaleByAndEqual(c_C[i-1][j], logp);
			algo.partialSlotsSumAndEqual(c_C[i-1][j], n);

			if(c_C[i-1][j].logq < c_A[j].logq)
				scheme.modDownToAndEqual(c_A[j], c_C[i-1][j].logq);
			  
			else 
				scheme.modDownToAndEqual(c_C[i-1][j], c_A[j].logq);
			scheme.mult(temp, c_C[i-1][j], c_A[j]);   
			scheme.reScaleByAndEqual(temp, logp);
            
			if(v.logq < temp.logq)
				scheme.modDownToAndEqual(temp, v.logq);
			else 
				scheme.modDownToAndEqual(v, temp.logq);
			scheme.subAndEqual(v, temp);
			
		}
		//timeccb.stop("finish calculating c_B ");
		//outfile << "calculating c_B[" << i << "] time: " << timeccb.timeElapsed << "ms" << endl;
		
		c_A[i].copy(v);
		scheme.square(len[i], c_A[i]); 
		scheme.reScaleByAndEqual(len[i], logp);
		algo.partialSlotsSumAndEqual(len[i], n);
		c = EvaluatorUtils::randomReal(10);
	    scheme.multByConstAndEqual(len[i], c, logp);
		scheme.reScaleByAndEqual(len[i], logp);
		lenm[i] = scheme.decrypt(secretKey, len[i]);

		//NTL_EXEC_RANGE(n, first, last)
		for(long t = 0; t < n; ++ t){
			lenm_inv[i][t].real(1./lenm[i][t].real());
		}
        //NTL_EXEC_RANGE_END

		scheme.encrypt(len_inv[i], lenm_inv[i], n, logp, logq);
		scheme.multByConstAndEqual(len_inv[i], c, logp);
		scheme.reScaleByAndEqual(len_inv[i], logp);

        //outfile << "start calcucating c_R[ ][" << i << "]" << endl;
		for(long j = 0; j <= i-1; ++j){
			
			//timeccr.start("calcucating c_R");
			scheme.negateAndEqual(c_C[i-1][j]);
			for(long k = j+1; k <= i-1; ++k){

				if(c_C[k-1][j].logq < c_C[i-1][k].logq)
					scheme.modDownToAndEqual(c_C[i-1][k], c_C[k-1][j].logq);
			    else 
				scheme.modDownToAndEqual(c_C[k-1][j], c_C[i-1][k].logq);
				scheme.mult(temp, c_C[k-1][j], c_C[i-1][k]);   
				scheme.reScaleByAndEqual(temp, logp);
				if(c_C[i-1][j].logq < temp.logq)
			       scheme.modDownToAndEqual(temp, c_C[i-1][j].logq);
			    else if(c_C[i-1][j].logq > temp.logq)
			       scheme.modDownToAndEqual(c_C[i-1][j], temp.logq);
				scheme.subAndEqual(c_C[i-1][j], temp);
			}
			//timeccr.stop("finish calcucating c_R");
		}
			
	}

	outfile << "computing AX = b;" << endl;
	timeutils2.start("part2 calculating x");

	//NTL_EXEC_RANGE(dim_m, first, last)
	for(long i = 0; i < dim_m; ++i){
		c_A[i].logq > c_b.logq ? scheme.modDownToAndEqual(c_A[i], c_b.logq) : scheme.modDownToAndEqual(c_b, c_A[i].logq);
        scheme.multAndEqual(c_A[i], c_b);
		scheme.reScaleByAndEqual(c_A[i], logp);
		algo.partialSlotsSumAndEqual(c_A[i], n);
		c_A[i].logq > len_inv[i].logq ? scheme.modDownToAndEqual(c_A[i], len_inv[i].logq) : scheme.modDownToAndEqual(len_inv[i], c_A[i].logq);
		scheme.multAndEqual(c_A[i], len_inv[i]);
		scheme.reScaleByAndEqual(c_A[i], logp);
	}
	//NTL_EXEC_RANGE_END

	//NTL_EXEC_RANGE(dim_m, first, last) 此处不可并行
	for(int i = 0; i < dim_m; ++i){
		Ciphertext temp;
		for(int j = i + 1; j < dim_m; ++ j){
			c_C[j-1][i].logq < c_A[j].logq ? scheme.modDownToAndEqual(c_A[j], c_C[j-1][i].logq) : scheme.modDownToAndEqual(c_C[j-1][i], c_A[j].logq);
			scheme.mult(temp, c_C[j-1][i], c_A[j]);
			scheme.reScaleByAndEqual(temp, logp);
			c_A[i].logq < temp.logq ? scheme.modDownToAndEqual(temp, c_A[i].logq) : scheme.modDownToAndEqual(c_A[i], temp.logq);
			scheme.addAndEqual(c_A[i], temp);
		}
	}
	//NTL_EXEC_RANGE_END

	timeutils1.stop("finish leastSqaureHE..........");
	timeutils2.stop("finish part2, calculating x");
	outfile << "get result part2 x time: " << timeutils2.timeElapsed << "ms" << endl;
	
	complex<double>* mm5;
	for(int i = 0; i < dim_m; ++ i){
		mm5 = scheme.decrypt(secretKey, c_A[i]);
		outfile << "---------------------" << endl;
		for(int z = 0; z < dim_m; ++ z){
			x_m[i][z] = mm5[z].real();
			outfile << "m LR : " << i << " :" << x[i];
			outfile << "  d LR : " << i << " :" << x_m[i][z];
			outfile << "  e LR : " << i << " :" << x[i]-x_m[i][z] << endl;
		}
		outfile << "---------------------" << endl;
	}
	delete [] mm5;
	
	outfile << "X.logq=:";
	for(int i = 0 ; i < dim_m; ++ i){
		outfile << c_A[i].logq << " ";
	}
	outfile << "    consumeMUL:" << (logq - c_A[0].logq)/logp << endl;
	outfile << "HE Time consuming:= " << timeutils1.timeElapsed << "ms" << endl;
	outfile << endl << "********************************************************" << endl;
	outfile << "finish HE" << endl;
	
	delete [] len;
	delete [] len_inv;
	delete [] x;
	delete [] b;
	for(long i = 0; i < dim_m; ++i){
		delete [] x_m[i];
		delete [] lenm[i];
		delete [] lenm_inv[i];
	}
	delete [] x_m;
	delete [] lenm;
	delete [] lenm_inv;
	for(long i =0 ; i < dim_m-1; ++i) delete []c_C[i];
	delete [] c_C;
}

void LR::leastSqaureHE(long logq, long logp, long logn, long dim_m, long models, ofstream& outfile){
    srand(time(NULL));
	SetNumThreads(16);
	TimeUtils timeutils1;
	TimeUtils timeutils2;
	TimeUtils timeccr;
	TimeUtils timeccb;

	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addLeftRotKeys(secretKey);
	SchemeAlgo algo(scheme);

	long n = (1 << logn);

	
	double **A = new double *[dim_m];
	double *b = new double[n];
	double **sA = new double *[dim_m];
	double *sb = new double[n];
	//double *b = new double[n]{1, 24, 5, 0, 1, 67, 7, 0};
    double* x = new double[models*dim_m];
    double** x_m = new double*[dim_m];

	//A[0] = new double[n]{2, 2, 2, 1, 1, 0.5, 5.5, 0};
	//A[1] = new double[n]{0.5, 0.5, 7.6, 0, 2, 3, 3.7, 1};
	

	outfile << endl << "********************************************************" << endl;
	outfile << "print Params:.............." << endl;
	outfile << "logq:=" << logq << " logp:=" 
			<< logp << " dimension:=" << dim_m << " models:= " << models << endl;
	outfile << "generating Matrix A...." << endl;
	//outfile << "A:=" << endl;
	for(long i = 0; i < dim_m; ++ i){
		A[i] = EvaluatorUtils::randomRealArray(n, 100.);
		for(long j = models*dim_m; j < n; ++ j) A[i][j] = 0;
		//StringUtils::showVec(A[i], n);
		//StringUtils::writeVec(A[i], n, outfile);
	}
	outfile << "generating vector b...." << endl;
	b = EvaluatorUtils::randomRealArray(n, 100);
	for(long i = models*dim_m; i < n; ++ i) b[i] = 0;
	//outfile << "b：= " << endl; 
	//StringUtils::showVec(b, n);
	//StringUtils::writeVec(b, n, outfile);

	

    
	
	timeutils1.start("leastSqaureHE..........");

	Ciphertext *c_A = new Ciphertext[dim_m];
	Ciphertext **c_C = new Ciphertext *[dim_m-1];
	Ciphertext c_b;
	
	for(long i = 0; i < dim_m; i ++){
	    x_m[i] = new double[models*dim_m];
	}
	for(long i = 1; i < dim_m; ++i ){
		c_C[i-1] = new Ciphertext[i];
	}
	scheme.encrypt(c_b, b, n, logp, logq);
	for(long i = 0; i < dim_m; ++ i){
		scheme.encrypt(c_A[i], A[i], n, logp, logq);
		//delete [] A[i];
	}
	//delete [] A;
	cout << "error 1!" << endl;
	
	Ciphertext *len = new Ciphertext[dim_m];
	Ciphertext *len_inv = new Ciphertext[dim_m];
	double c = 0;
	complex<double> **lenm = new complex<double> *[dim_m];
    complex<double> **lenm_inv = new complex<double> *[dim_m];
	for(long i = 0; i < dim_m; ++ i){
		lenm_inv[i] = new complex<double>[n];
	}
	cout << "error 2!" << endl;

	scheme.square(len[0], c_A[0]);
	scheme.reScaleByAndEqual(len[0], logp);
	cout << "error 30!" << endl;
	algo.partialSlotsSumAndEqual(len[0], models, n);
	cout << "error 31!" << endl;
	c = EvaluatorUtils::randomReal(10);
	scheme.multByConstAndEqual(len[0], c, logp);
	scheme.reScaleByAndEqual(len[0], logp);
	lenm[0] = scheme.decrypt(secretKey, len[0]);
	//StringUtils::showVec(lenm[0], n);
    cout << "error 3!" << endl;
	
	//NTL_EXEC_RANGE(n, first, last)
	for(long i = 0; i < n; ++ i){
		lenm_inv[0][i].real(1./lenm[0][i].real());
	}
    //NTL_EXEC_RANGE_END

	//StringUtils::showVec(lenm_inv[0], n);
	scheme.encrypt(len_inv[0], lenm_inv[0], n, logp, logq);
	scheme.multByConstAndEqual(len_inv[0], c, logp);
	lenm_inv[0] = scheme.decrypt(secretKey ,len_inv[0]);
	//StringUtils::showVec(lenm_inv[0], n);
	scheme.reScaleByAndEqual(len_inv[0], logp);
	

    Ciphertext v;
	Ciphertext temp;
	for(long i = 1; i <= dim_m-1; ++ i){
		outfile << "start calculating c_R and c_B[" << i << "]" << endl;
		//timeccb.start(" calculating c_B ");
		
		v.copy(c_A[i]);
		for(long j = 0; j <= i-1; ++ j){
			if(len_inv[j].logq < c_A[i].logq)
				scheme.modDownToAndEqual(c_A[i], len_inv[j].logq);  
			else 
				scheme.modDownToAndEqual(len_inv[j], c_A[i].logq);
			scheme.mult(temp, len_inv[j], c_A[i]);
			scheme.reScaleByAndEqual(temp, logp);

            if(temp.logq < c_A[j].logq)
				scheme.modDownToAndEqual(c_A[j], temp.logq); 
			else 
				scheme.modDownToAndEqual(temp, c_A[j].logq);
			scheme.mult(c_C[i-1][j], temp, c_A[j]);
			scheme.reScaleByAndEqual(c_C[i-1][j], logp);
			algo.partialSlotsSumAndEqual(c_C[i-1][j], models, n);

			if(c_C[i-1][j].logq < c_A[j].logq)
				scheme.modDownToAndEqual(c_A[j], c_C[i-1][j].logq);
			  
			else 
				scheme.modDownToAndEqual(c_C[i-1][j], c_A[j].logq);
			scheme.mult(temp, c_C[i-1][j], c_A[j]);   
			scheme.reScaleByAndEqual(temp, logp);
            
			if(v.logq < temp.logq)
				scheme.modDownToAndEqual(temp, v.logq);
			else 
				scheme.modDownToAndEqual(v, temp.logq);
			scheme.subAndEqual(v, temp);
			
		}
		//timeccb.stop("finish calculating c_B ");
		//outfile << "calculating c_B[" << i << "] time: " << timeccb.timeElapsed << "ms" << endl;
		
		c_A[i].copy(v);
		scheme.square(len[i], c_A[i]); 
		scheme.reScaleByAndEqual(len[i], logp);
		algo.partialSlotsSumAndEqual(len[i], models, n);
		c = EvaluatorUtils::randomReal(10);
	    scheme.multByConstAndEqual(len[i], c, logp);
		scheme.reScaleByAndEqual(len[i], logp);
		lenm[i] = scheme.decrypt(secretKey, len[i]);

		//NTL_EXEC_RANGE(n, first, last)
		for(long t = 0; t < n; ++ t){
			lenm_inv[i][t].real(1./lenm[i][t].real());
		}
        //NTL_EXEC_RANGE_END

		scheme.encrypt(len_inv[i], lenm_inv[i], n, logp, logq);
		scheme.multByConstAndEqual(len_inv[i], c, logp);
		scheme.reScaleByAndEqual(len_inv[i], logp);

        //outfile << "start calcucating c_R[ ][" << i << "]" << endl;
		for(long j = 0; j <= i-1; ++j){
			
			//timeccr.start("calcucating c_R");
			scheme.negateAndEqual(c_C[i-1][j]);
			for(long k = j+1; k <= i-1; ++k){

				if(c_C[k-1][j].logq < c_C[i-1][k].logq)
					scheme.modDownToAndEqual(c_C[i-1][k], c_C[k-1][j].logq);
			    else 
				scheme.modDownToAndEqual(c_C[k-1][j], c_C[i-1][k].logq);
				scheme.mult(temp, c_C[k-1][j], c_C[i-1][k]);   
				scheme.reScaleByAndEqual(temp, logp);
				if(c_C[i-1][j].logq < temp.logq)
			       scheme.modDownToAndEqual(temp, c_C[i-1][j].logq);
			    else if(c_C[i-1][j].logq > temp.logq)
			       scheme.modDownToAndEqual(c_C[i-1][j], temp.logq);
				scheme.subAndEqual(c_C[i-1][j], temp);
			}
			//timeccr.stop("finish calcucating c_R");
		}
			
	}

	outfile << "computing AX = b;" << endl;
	timeutils2.start("part2 calculating x");

	//NTL_EXEC_RANGE(dim_m, first, last)
	for(long i = 0; i < dim_m; ++i){
		c_A[i].logq > c_b.logq ? scheme.modDownToAndEqual(c_A[i], c_b.logq) : scheme.modDownToAndEqual(c_b, c_A[i].logq);
        scheme.multAndEqual(c_A[i], c_b);
		scheme.reScaleByAndEqual(c_A[i], logp);
		algo.partialSlotsSumAndEqual(c_A[i], models, n);
		c_A[i].logq > len_inv[i].logq ? scheme.modDownToAndEqual(c_A[i], len_inv[i].logq) : scheme.modDownToAndEqual(len_inv[i], c_A[i].logq);
		scheme.multAndEqual(c_A[i], len_inv[i]);
		scheme.reScaleByAndEqual(c_A[i], logp);
	}
	//NTL_EXEC_RANGE_END

	//NTL_EXEC_RANGE(dim_m, first, last) 此处不可并行
	for(int i = 0; i < dim_m; ++i){
		Ciphertext temp;
		for(int j = i + 1; j < dim_m; ++ j){
			c_C[j-1][i].logq < c_A[j].logq ? scheme.modDownToAndEqual(c_A[j], c_C[j-1][i].logq) : scheme.modDownToAndEqual(c_C[j-1][i], c_A[j].logq);
			scheme.mult(temp, c_C[j-1][i], c_A[j]);
			scheme.reScaleByAndEqual(temp, logp);
			c_A[i].logq < temp.logq ? scheme.modDownToAndEqual(temp, c_A[i].logq) : scheme.modDownToAndEqual(c_A[i], temp.logq);
			scheme.addAndEqual(c_A[i], temp);
		}
	}
	//NTL_EXEC_RANGE_END

	timeutils1.stop("finish leastSqaureHE..........");
	timeutils2.stop("finish part2, calculating x");
	outfile << "get result part2 x time: " << timeutils2.timeElapsed << "ms" << endl;
	complex<double>* mm5;
	for(long mod = 0; mod < models; ++ mod){
		for(long i = 0; i < dim_m; ++ i){
			sA[i] = new double[dim_m];
			for(long s = 0, t = mod; s < dim_m; ++ s, t += models) 
				sA[i][s] = A[i][t];
		}
		for(long s = 0, t = mod; s < dim_m; ++ s, t += models) sb[s] = b[t];
		timeutils2.start("leastSqaure..........");
    	LR::leastSqaure(sA, sb, x, dim_m);
		timeutils2.stop("finish leastSqaure..........");
		outfile << "Time consuming:= " << timeutils2.timeElapsed << "ms" << endl;
		outfile << "models " << mod << endl;
		for(int i = 0; i < dim_m; ++ i){
			mm5 = scheme.decrypt(secretKey, c_A[i]);
			outfile << "---------------------" << endl;
			for(int z = mod; z < n; z += models){
				x_m[i][z] = mm5[z].real();
				outfile << "m LR : " << i << " :" << x[i];
				outfile << "  d LR : " << i << " :" << x_m[i][z];
				outfile << "  e LR : " << i << " :" << x[i]-x_m[i][z] << endl;
			}
			outfile << "---------------------" << endl;
		}
	}
	
	delete [] mm5;
	
	outfile << "X.logq=:";
	for(int i = 0 ; i < dim_m; ++ i){
		outfile << c_A[i].logq << " ";
	}
	outfile << "    consumeMUL:" << (logq - c_A[0].logq)/logp << endl;
	outfile << "HE Time consuming:= " << timeutils1.timeElapsed << "ms" << endl;
	outfile << endl << "********************************************************" << endl;
	outfile << "finish HE" << endl;
	
	for(long i = 0; i < dim_m; ++ i){
		delete [] A[i];
	}
	delete [] A;
	delete [] len;
	delete [] len_inv;
	delete [] x;
	delete [] b;
	for(long i = 0; i < dim_m; ++i){
		delete [] x_m[i];
		delete [] lenm[i];
		delete [] lenm_inv[i];
	}
	delete [] x_m;
	delete [] lenm;
	delete [] lenm_inv;
	for(long i =0 ; i < dim_m-1; ++i) delete []c_C[i];
	delete [] c_C;
	outfile << "exit..." << endl;
}

void LR::leastSqaure(double** A, double* b, double * x, int size){
    double** B = new double*[size];
    double** C = new double*[size];
    double** R = new double*[size];
    double* len = new double[size];
    double* v = new double[size];
    double* len_inv = new double[size];

    for(int i = 0; i < size; ++ i){
        B[i] = new double[size];
        C[i] = new double[size];
        R[i] = new double[size];
        for(int j = 0; j < size; ++ j){
            B[i][j] = 0;
            C[i][j] = 0;
            R[i][j] = 0;
        }
    }
    for(int j = 0; j < size; ++ j){
            v[j] = 0;
        }

    VectorUtils::vecCopy(B[0], A[0], size);
    cout << "B[0]:=";
    VectorUtils::showVec(B[0], size);
    VectorUtils::vecDot(len[0], B[0], B[0], size);
    cout << "len[0]:=" << len[0] << endl;
    len_inv[0] = 1 / len[0];
    cout << "len_inv[0]:=" << len_inv[0] << endl;

    double* temp = new double[size]; 
    for(int i = 1; i <= size - 1; ++ i){
        VectorUtils::vecCopy(v, A[i], size);
        for(int j = 0; j <= i-1; ++ j){
            VectorUtils::vecDot(C[j][i], A[i], B[j], size);
            C[j][i] *= len_inv[j];
            cout << "c[" << j << "," << i << "]:=" <<  C[j][i] << endl;
            VectorUtils::vecMulConst(temp, B[j], C[j][i], size);
            //VectorUtils::showVec(B[j], size);
            //VectorUtils::showVec(temp, size);
            VectorUtils::vecNegate(temp, temp, size);
            //VectorUtils::showVec(temp, size);
            VectorUtils::vecAddVec(v, temp, v, size);
            //VectorUtils::showVec(v, size);
        }
        for(int j = 0; j <= i-1; ++ j){
            R[j][i] = -C[j][i];
            for(int k = j+1; k <= i-1; ++ k){
                R[j][i] = R[j][i] - R[j][k] * C[k][i];
            }
            cout << "R[" << j << "," << i << "]:=" <<  R[j][i] << endl;
        }
        VectorUtils::vecCopy(B[i], v, size);
        cout << "B[" << i << "]:=";
        VectorUtils::showVec(B[i], size);
        VectorUtils::vecDot(len[i], B[i], B[i], size);
        len_inv[i] = 1 / len[i];
        cout << "len[" << i << "]:=" << len[i] << endl;
        cout << "len_inv[" << i << "]:=" << len_inv[i] << endl;
    }
    delete [] temp;
    
    for(int i = 0; i < size; ++ i){
        VectorUtils::vecDot(x[i], B[i], b, size);
        x[i] *= len_inv[i];
    }
    for(int i = 0; i < size; ++ i){
        double temp = x[i];
        for(int j = i+1; j < size; ++ j){
            temp += x[j] * R[i][j];
        }
        x[i] = temp;
    }
	VectorUtils::showVec(x, size);
	for(int i = 0; i < size; ++ i){
        delete [] B[i];
        delete [] C[i];
        delete [] R[i];
        
    }
    delete [] B;
    delete [] C;
    delete [] R;
	delete [] len;
    delete [] v;
    delete [] len_inv;
	
};
