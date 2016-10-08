///
/// C++ suport code for Prabhakaran, Sandhya, Sudhir Raman, Julia E. Vogt, and Volker Roth. "Automatic model selection in archetype analysis." In Joint DAGM (German Association for Pattern Recognition) and OAGM Symposium, pp. 458-467. Springer Berlin Heidelberg, 2012.
///
/// Code authors - Prof. Dr. Volker Roth and Sandhya Prabhakaran,
///                Department of Mathematics and Computer Science
///                University of Basel, Switzerland
///
/// February 2012.
///


#include "my_mat.h"
#include "my_mat_bool_scale.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "MultiTask_Lasso.h"



extern "C"{
  int MultiTaskLasso(int* initIndex, int*N, int*K, int* dim, double* Amat, int* task_labels, double* y,  double* beta_path, int* LOGREG, int* infty, double* kappa, double* errtol, double* stepsize, int* iterN, int* modulo, int* update_groups,int *max_iter, int*  max_groups){

    
  bool bool_infty = false;
  if(*infty)
    bool_infty = true;
 
  bool bool_LOGREG = false;
  if(*LOGREG)
    bool_LOGREG = true;
      
  vector<int> s_0(2,0);
  vector<vector<int> > start(*K,s_0);

  My_Vector targets(*N);

  My_Matrix A(*N,*dim);


  int idx =0;
  for(int k = 0; k<*K; k++){
    
    start[k][0] = idx;
    start[k][1] = idx-1;
    for(int i =0; i< *N; i++){
      if(task_labels[i] == k){
	
	for(int j =0; j<*dim; j++)
	  A[idx][j] = Amat[i* (*dim)+j]; // encode in R as "as.double(t(A))"
	targets[idx] = y[i];
	start[k][1]++; // increment stop indicator
	idx++;
      }
    }
  }
     

 
  vector<int> init_idx(2);
  init_idx[0] = initIndex[0];
  init_idx[1] = initIndex[1];
 
  

  L_ONE_INFTY LASSO(A, targets, start, bool_infty,  *max_groups, init_idx);
  
  LASSO.SetParams(  *kappa, *errtol,  *stepsize,  *iterN,  *modulo,  *update_groups, *max_iter);
  
  LASSO.SetLOGREG(bool_LOGREG);  
  
  LASSO.optimize();

  vector<vector<My_Vector> > beta = LASSO.get_beta_path();
  
  vector<vector<int> > work_index = LASSO.get_work_index_path();

  int i_idx = *iterN -1;
  vector<My_Vector> beta_last = beta[i_idx];
  
 
  for(int j=0; j< (*dim); j++){
    for(int k =0; k< (*K); k++){
      beta_path[j * (*K) + k] = 0;
    }
  }
  
  for(int j=0; j<beta_last.size();j++){
    int j_idx = work_index[i_idx][j];
    for(int k =0; k< (*K); k++){
      beta_path[ j_idx * (*K) +k] = beta_last[j].el(k); // load into R as "matrix(..., nrow = dim, ncol= K, byrow=TRUE)"
    }
  }
  
  return 0;


}
}
