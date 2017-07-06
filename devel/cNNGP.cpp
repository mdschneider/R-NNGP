#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//Description: update B and F.
void updateBF(double *B, double *F, double *c, double *C, double *D, double *d, int *nnIndxLU, int *CIndx, int n, double *thetaAlpha, int covModel, double *bk, double nuMax){
    
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';

  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuMax));
  int threadID = 0;

  double phi, nu, alpha;
  phi = thetaAlpha[0];
  nu = 0;
  alpha = thetaAlpha[1];
  if(covModel == 2){//matern
    nu = thetaAlpha[2];
  }
  
#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID)
#endif
    for(i = 0; i < n; i++){
      threadID = omp_get_thread_num();
      if(i > 0){
	for(k = 0; k < nnIndxLU[n+i]; k++){
	  c[nnIndxLU[i]+k] = spCor(d[nnIndxLU[i]+k], phi, nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= k; l++){
	    C[CIndx[i]+l*nnIndxLU[n+i]+k] = spCor(D[CIndx[i]+l*nnIndxLU[n+i]+k], phi, nu, covModel, &bk[threadID*nb]); 
	    if(l == k){
	      C[CIndx[i]+l*nnIndxLU[n+i]+k] += alpha;
	    }
	  }
	}
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[CIndx[i]], &nnIndxLU[n+i], &c[nnIndxLU[i]], &inc, &zero, &B[nnIndxLU[i]], &inc);
	F[i] = 1.0 + alpha - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[nnIndxLU[i]], &inc);
      }else{
	B[i] = 0;
	F[i] = 1.0 + alpha;
      }
    }

}

extern "C" {
  
  SEXP cNNGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coords_r,
	     SEXP X0_r, SEXP n0_r, SEXP coords0_r, SEXP nnIndx0_r, SEXP g_r, SEXP thetaAlpha_r, 
	     SEXP m_r, SEXP sigmaSqIG_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r){
    
    int h, i, j, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *upper = "U";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *rside = "R";
    char const *lside = "L";
    
    //get args
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    double *coords = REAL(coords_r);

    double *X0 = REAL(X0_r);
    int n0 = INTEGER(n0_r)[0];
    double *coords0 = REAL(coords0_r);
    int *nnIndx0 = INTEGER(nnIndx0_r);
    int g = INTEGER(g_r)[0];
    double *thetaAlpha = REAL(thetaAlpha_r); //nTheta x g

    int m = INTEGER(m_r)[0];
    double sigmaSqIGa = REAL(sigmaSqIG_r)[0]; double sigmaSqIGb = REAL(sigmaSqIG_r)[1];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("NNGP Conjugate model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("------------\n");
      Rprintf("Priors and hyperpriors:\n");
      Rprintf("\tbeta flat.\n");
      Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);
      Rprintf("------------\n");
      if(corName == "matern"){
	Rprintf("Considering %i combinations of phi, nu, and alpha.\n\n", g);
      }else{
	Rprintf("Considering %i combinations of phi and alpha.\n\n", g);
      }
      if(n0 > 0){
	Rprintf("Predicting at %i locations.\n\n", n0);
      }
      Rprintf("------------\n");
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
    
    //parameters
    int nTheta = 2; //phi, alpha
    
    if(corName == "matern"){
      nTheta++;//nu
    }
    

    //allocated for the nearest neighbor index vector (note, first location has no neighbors).
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    int *nnIndx = (int *) R_alloc(nIndx, sizeof(int));
    double *d = (double *) R_alloc(nIndx, sizeof(double));
    int *nnIndxLU = (int *) R_alloc(2*n, sizeof(int)); //first column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but will simplifying some parallelization).

    //make the neighbor index
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\tBuilding neighbor index\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }
    
    mkNNIndx(n, m, coords, nnIndx, d, nnIndxLU);
       
    //other stuff
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(n, sizeof(double));
    double *c =(double *) R_alloc(nIndx, sizeof(double));
    
    int *CIndx = (int *) R_alloc(2*n, sizeof(int)); //index for D and C.
    for(i = 0, j = 0; i < n; i++){//zero should never be accessed
      j += nnIndxLU[n+i]*nnIndxLU[n+i];
      if(i == 0){
	CIndx[n+i] = 0;
	CIndx[i] = 0;
      }else{
	CIndx[n+i] = nnIndxLU[n+i]*nnIndxLU[n+i]; 
	CIndx[i] = CIndx[n+i-1] + CIndx[i-1];
      }
    }
    
    double *C = (double *) R_alloc(j, sizeof(double)); zeros(C, j);
    double *D = (double *) R_alloc(j, sizeof(double)); zeros(D, j);
    
    for(i = 0; i < n; i++){
      for(k = 0; k < nnIndxLU[n+i]; k++){   
	for(l = 0; l <= k; l++){
	  D[CIndx[i]+l*nnIndxLU[n+i]+k] = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
	}
      }
    }

    //prediction
    int mm = m*m;
    double *tmp_m = new double[m];
    double *D0 = new double[n0*mm];
    double *d0 = new double[n0*m];
    double *C0 = new double[n0*mm];
    double *c0 = new double[n0*m];
    double *tmp_mn0 = new double[n0*m];
    
    for(i = 0; i < n0; i++){
      for(k = 0; k < m; k++){
	d0[i*m+k] = dist2(coords0[i], coords0[n0+i], coords[nnIndx0[k*n0+i]], coords[n+nnIndx0[k*n0+i]]);
	for(l = 0; l < m; l++){
	  D0[i*mm+k*m+l] = dist2(coords[nnIndx0[k*n0+i]], coords[n+nnIndx0[k*n0+i]], coords[nnIndx0[l*n0+i]], coords[n+nnIndx0[l*n0+i]]);
	}
      }
    }
    

    // //return stuff  
    // SEXP betaSamples_r, thetaSamples_r;
    // PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++;
    // PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nSamples)); nProtect++; 

    //other stuff
    double status = 0;
    int pp = p*p;
    int mp = m*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *beta = (double *) R_alloc(g*p, sizeof(double));
    double *ab = (double *) R_alloc(g*nTheta, sizeof(double));
    double *y0Hat = (double *) R_alloc(g*n0, sizeof(double));
    double *y0HatVar = (double *) R_alloc(g*n0, sizeof(double));
    double *tmp_mp = (double *) R_alloc(mp, sizeof(double));

    //get max nu
    double *bk = NULL;
    double nuMax = 0;
    
    if(corName == "matern"){
      for(i = 0; i < g; i++){
	if(thetaAlpha[i*nTheta+2] > nuMax){
	  nuMax = thetaAlpha[i*nTheta+2];
	}
      }
      
      bk = (double *) R_alloc(nThreads*(static_cast<int>(1.0+nuMax)), sizeof(double));
    }
    
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    GetRNGstate();
    
    for(k = 0; k < g; k++){

      //update B and F
      updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, &thetaAlpha[k*nTheta], covModel, bk, nuMax);

      //estimation
      for(i = 0; i < p; i++){
	tmp_p[i] = Q(B, F, &X[n*i], y, n, nnIndx, nnIndxLU); //X'C^{-1}y
	for(j = 0; j <= i; j++){
	  tmp_pp[j*p+i] = Q(B, F, &X[n*j], &X[n*i], n, nnIndx, nnIndxLU);//X'C^{-1}X
	}
      }

      //V = inv(B), V = inv(tmp_pp), V = tmp_pp after dpotri
      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}

      //g = solve(B, v), g = solve(tmp_pp, tmp_p), g = beta
      F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, &beta[k*p], &inc);

      //a 
      ab[k*2] = sigmaSqIGa + 1.0*n/2.0;
      
      //b
      ab[k*2+1] = sigmaSqIGb + 0.5*(Q(B, F, y, y, n, nnIndx, nnIndxLU) - ddot_(&p, &beta[k*p], &inc, tmp_p, &inc));
    
      //prediction
      for(i = 0; i < n0; i++){
	//stopped here
	//make z
	for(j = 0; j < m; j++){
	  c0[i*m+j] = exp(-phiAlpha[k]*d0[i*m+j]);
	}
	
	//Make M
	for(j = 0; j < mm; j++){
	  C0[i*mm+j] = exp(-phiAlpha[k]*D0[i*mm+j]);
	}
	
	for(j = 0; j < m; j++){
	  C0[i*mm+j*m+j] += phiAlpha[g+k]; 
	}
	
	//make w
	dpotrf_(&lower, &m, &C0[i*mm], &m, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
	dpotri_(&lower, &m, &C0[i*mm], &m, &info); if(info != 0){cout << "c++ error: dpotri failed" << endl;}
	dsymv_(&lower, &m, &one, &C0[i*mm], &m, &c0[i*m], &inc, &zero, &tmp_mn0[i*m], &inc);
	
	//make hat(y)
	for(j = 0; j < m; j++){
	  tmp_m[j] = y[nnIndx0[j*n0+i]] - ddot_(&p, &X[nnIndx0[j*n0+i]], &n, &beta[k*p], &inc);
	}
	
	yHat[k*n0+i] = ddot_(&p, &X0[i], &n0, &beta[k*p], &inc) + ddot_(&m, &tmp_mn0[i*m], &inc, tmp_m, &inc);
	
	//make u
	for(j = 0; j < m; j++){
	  dcopy_(&p, &X[nnIndx0[j*n0+i]], &n, &tmp_mp[j], &m);
	}
	
	dgemv_(&ytran, &m, &p, &one, tmp_mp, &m, &tmp_mn0[i*m], &inc, &zero, tmp_p, &inc);
	
	for(j = 0; j < p; j++){
	  tmp_p[j] = X0[j*n0+i] - tmp_p[j];
	}
	
	//make v_y and var(y)
	dsymv_(&lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc);
	
	yHatVar[k*n0+i] = 2.0*ab[k*2+1] * (ddot_(&p, tmp_p, &inc, tmp_p2, &inc) + 1.0 + phiAlpha[g+k] - ddot_(&m, &tmp_mn0[i*m], &inc, &c0[i*m], &inc))/(2.0*ab[k*2]-1.0);
	
      }
      // if(thetaUpdate){
	
      // 	thetaUpdate = false;
	
      // 	///////////////
      // 	//update beta 
      // 	///////////////
      // 	for(i = 0; i < p; i++){
      // 	  tmp_p[i] = Q(B, F, &X[n*i], y, n, nnIndx, nnIndxLU);
      // 	  for(j = 0; j <= i; j++){
      // 	    tmp_pp[j*p+i] = Q(B, F, &X[n*j], &X[n*i], n, nnIndx, nnIndxLU);
      // 	  }
      // 	}
	
      // 	F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      // 	F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
      // 	F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc);
      // 	F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      // }
      
      // mvrnorm(beta, tmp_p2, tmp_pp, p);
      
      // ///////////////
      // //update theta
      // ///////////////
      // F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc);
      // F77_NAME(daxpy)(&n, &negOne, y, &inc, tmp_n, &inc);
      
      // //current    
      // logDetCurrent = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, theta, tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nuUnifb);
      
      // QCurrent = Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU);

      // logPostCurrent = -0.5*logDetCurrent - 0.5*QCurrent;
      // logPostCurrent += log(theta[phiIndx] - phiUnifa) + log(phiUnifb - theta[phiIndx]); 
      // logPostCurrent += -1.0*(1.0+sigmaSqIGa)*log(theta[sigmaSqIndx])-sigmaSqIGb/theta[sigmaSqIndx]+log(theta[sigmaSqIndx]);
      // logPostCurrent += -1.0*(1.0+tauSqIGa)*log(theta[tauSqIndx])-tauSqIGb/theta[tauSqIndx]+log(theta[tauSqIndx]);

      //  if(corName == "matern"){
      // 	 logPostCurrent += log(theta[nuIndx] - nuUnifa) + log(nuUnifb - theta[nuIndx]); 
      //  }
      
      // //candidate
      // thetaCand[phiIndx] = logitInv(rnorm(logit(theta[phiIndx], phiUnifa, phiUnifb), tuning[phiIndx]), phiUnifa, phiUnifb);
      // thetaCand[sigmaSqIndx] = exp(rnorm(log(theta[sigmaSqIndx]), tuning[sigmaSqIndx]));
      // thetaCand[tauSqIndx] = exp(rnorm(log(theta[tauSqIndx]), tuning[tauSqIndx]));

      // if(corName == "matern"){
      // 	thetaCand[nuIndx] = logitInv(rnorm(logit(theta[nuIndx], nuUnifa, nuUnifb), tuning[nuIndx]), nuUnifa, nuUnifb);
      // }
      
      // //update B and F
      // logDetCand = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, thetaCand, tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nuUnifb);
      
      // QCand = Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU);
      
      // logPostCand = -0.5*logDetCand - 0.5*QCand;
      // logPostCand += log(thetaCand[phiIndx] - phiUnifa) + log(phiUnifb - thetaCand[phiIndx]); 
      // logPostCand += -1.0*(1.0+sigmaSqIGa)*log(thetaCand[sigmaSqIndx])-sigmaSqIGb/thetaCand[sigmaSqIndx]+log(thetaCand[sigmaSqIndx]);
      // logPostCand += -1.0*(1.0+tauSqIGa)*log(thetaCand[tauSqIndx])-tauSqIGb/thetaCand[tauSqIndx]+log(thetaCand[tauSqIndx]);

      //  if(corName == "matern"){
      // 	 logPostCand += log(thetaCand[nuIndx] - nuUnifa) + log(nuUnifb - thetaCand[nuIndx]); 
      //  }
      
      // if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){
      // 	thetaUpdate = true;
      // 	dcopy_(&nTheta, thetaCand, &inc, theta, &inc);
      // 	accept++;
      // 	batchAccept++;
      // }

      // //save samples
      // F77_NAME(dcopy)(&p, beta, &inc, &REAL(betaSamples_r)[s*p], &inc);
      // F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[s*nTheta], &inc);

      // //report
     //report
   
	if(verbose){
	  if(corName == "matern"){
	    Rprintf("Update for phi=%.5f, nu=%.5f, and alpha=%.5f\n", thetaAlpha[k*nTheta], thetaAlpha[k*nTheta+2], thetaAlpha[k*nTheta+1]);
	  }else{
	    Rprintf("Update for phi=%.5f and alpha=%.5f\n", thetaAlpha[k*nTheta], thetaAlpha[k*nTheta+1]);
	  }
    
          #ifdef Win32
      	  R_FlushConsole();
          #endif
      	}

      
      R_CheckUserInterrupt();
    }
    
    PutRNGstate();

    // //make return object
    // SEXP result_r, resultName_r;
    // int nResultListObjs = 2;

    // PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    // PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    // SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    // SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.samples")); 
    
    // SET_VECTOR_ELT(result_r, 1, thetaSamples_r);
    // SET_VECTOR_ELT(resultName_r, 1, mkChar("p.theta.samples")); 

    // namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(m_r);
  }
}
