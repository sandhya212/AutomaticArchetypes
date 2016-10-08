#include "MultiTask_Lasso.h"



struct sortPerm {                                 /* struct can have may fields */
        double dPart;                             /* at least one will be used  */
        int iPart;                             /* for the sort key           */
};

sortPerm* sort_vect;

int CStruct(const void* e1, const void* e2)     /* struct comparison func */
{
  if(((sortPerm *)e1)->dPart - ((sortPerm *)e2)->dPart <0)
    return 1;
  else
    return (-1);
}


int CStructIncreasing(const void* e1, const void* e2)     /* struct comparison func */
{
  if(((sortPerm *)e1)->dPart - ((sortPerm *)e2)->dPart <0)
    return -1;
  else
    return 1;
}


int L_ONE_INFTY::ensure_constraint_active(double& lambda,  // current lagrange para
					  vector<int>& work_index, // original indices of groups in working set
					  const double& current_kappa,
					  vector<My_Vector>& beta, // regression weights
					  My_Vector& residualV,
					  double step
					  ){


  double maxg;
  int maxgi;
  double gn;
  int up_index;
  vector<int>::iterator v_iter;

  int J = J_global;
  int workJ = beta.size();

  while(lambda < 0){

  
   
    maxg = -1;
    maxgi = -1;
    My_Matrix_bool_scale G;
    
   

    for(int j=0; j< J;j++){
      
      v_iter = find(work_index.begin(), work_index.end(),j);
      if(v_iter !=  work_index.end())
	continue;

   
      gn = Norm(GRAD_J_RESIDUAL(residualV, j ));
      
    
      
    
      if(gn > maxg){
	maxg = gn;
	maxgi = j;
      }
    }
  
    
    up_index = maxgi; // the group with the largest gradient magnitude
    cout << "L_1-constraint not active.  Updating group: " <<  up_index << " with "<< maxg << endl;
    
    if(maxg < 1e-9){
      cout << "NN LS solution is feasible!!!" << endl; 
   
      return -1;
    }
  

    // update variable sets and beta

    My_Vector empty(K);
    work_index.push_back(up_index);
    beta.push_back(empty);
    

     // new call to optimizer

    if(infty)
      lambda = BSparse(work_index,  current_kappa, beta,  residualV, step);
    else{
      if(ell2){
	lambda = BSparse_L2(work_index,  current_kappa, beta,  residualV, step);
      }else{
	lambda = BSparse_LP(work_index,  current_kappa, beta,  residualV, step, Norm_p, Norm_gamma);
      }
    }
  }
  return 1;
}


My_Vector  L_ONE_INFTY::A_multiTask_times_x(const My_Matrix& X, const vector<vector<int> >& StartStop, const vector<int>& work_index, const vector<My_Vector>& beta){

  int work_J = work_index.size();
  
  My_Vector yhat(X.Nrows());
  for(int k =0; k< K; k++){ // loop over tasks
   
    for(int i = StartStop[k][0]; i<=StartStop[k][1];i++){
      double sum = 0.0;
      for(int j =0; j<work_J; j++)
	sum += X.el(i,work_index[j])*beta[j][k];
      yhat.el(i) = sum;
    }
  }
  return yhat;
  
}

My_Vector L_ONE_INFTY::Y_HAT(const vector<int>& work_index, const vector<My_Vector>& beta){

  double res;
 
  int work_J = work_index.size();
  
  My_Vector yhat(N_tot);

  for(int k =0; k< K; k++){ // loop over tasks
   
    for(int i = StartStop_global[k][0]; i<=StartStop_global[k][1];i++){

   
      double sum = 0.0;
      for(int j =0; j<work_J; j++)
	sum += X_global.el(i,work_index[j])*beta[j][k];
      yhat.el(i) = sum;
    }
  }
  if(LOGREG){
    for(int k=0; k<N_tot;k++){
      res =  1.0/(1.0+exp(-yhat.el(k)));
      if(res < 1e-7)
	res =1e-7;
      if(res > 1-1e-7)
	res =1-1e-7;
      yhat.el(k) = res;
    }
  }
  return yhat;
}

My_Vector L_ONE_INFTY::GRAD_J_RESIDUAL(const My_Vector& RESIDUAL,  const int& j ){
   
  My_Vector grad_j(K);


  int start =0;

  for(int k =0; k< K; k++){ // loop over tasks
    double sum = 0.0;
    for(int i = StartStop_global[k][0]; i<=StartStop_global[k][1];i++){
	sum += X_global.el(i,j)*RESIDUAL[i];
    }
    grad_j[k] = -sum;
  }
  return grad_j;
}



double L_ONE_INFTY::Likelihood(const vector<int>& work_index, const vector<My_Vector>& beta,My_Vector&  y_minus_pi){

  
  My_Vector xb = Y_HAT(work_index, beta);
 
  y_minus_pi = y-xb;
  double likeli = 0;
  
  if(LOGREG){
    for(int k=0; k<N_tot;k++){
      likeli +=  y.el(k) * log(xb.el(k)) + (1-y.el(k))*log(1-xb.el(k));
    }
  }
  else{
    likeli = -pow(Norm(y_minus_pi),2);
  }
  return likeli/N_tot;
}


////// Projection l_1/2 ///////////////////////////
double L_ONE_INFTY::BSparse_L2(vector<int>& work_index, 
				  const double& current_kappa,
				  vector<My_Vector>& beta, 
				  My_Vector&  y_minus_pi,
				  double s){
  
 
  int work_J = work_index.size();

  My_Vector M(work_J);
  My_Vector NormB(work_J);
  int nzn;

  vector<bool> activeSet(work_J,true); 
 
 
 
  int d =0;
  for(int j =0; j<work_J; j++){
    d += beta[j].Length();
  }
 
  vector<My_Vector> b(work_J);



  double s_buff = s;
  double sumNorm;
  double gn;
  int iter;
  double err = 1e100;
 
  double  lambda_avg;

  iter = 0;

 
  vector<My_Vector> grad(work_J);
  double likeli_new;
 

  cout << "work_index: " <<' ';
  for(int j =0; j< work_J; j++)
    cout << work_index[j] << ' ';
  cout << endl;

  double  likeli_old = Likelihood(work_index, beta, y_minus_pi);
  My_Vector  last_y_minus_pi =  y_minus_pi;
  vector<My_Vector> last_beta = beta;
  
  vector<My_Vector> u(work_J);
 
  int count =0;

  cout << "starting iterations.... "<< likeli_old <<endl;
  while(err > errtol &&  iter < max_iter){
    
  
    

    iter++;

    y_minus_pi = y-Y_HAT(work_index, beta);
  
    for(int j =0; j<work_J; j++){
      // grad[j] = GRAD_J(work_index, beta,j);
      grad[j] = GRAD_J_RESIDUAL(y_minus_pi, work_index[j] );
    }
     
    
    for(int j =0; j< work_J; j++){
      activeSet[j] = true;
    }   

   
    double maxSum =0;
  
   
  
    for(int j =0; j<work_J; j++)
      b[j] = beta[j] - grad[j]*s; 
  
   
   
    vector<My_Vector> beta_old = beta;
  

    // check constraint
    maxSum =0;
    for(int j =0; j< work_J; j++){
      maxSum += Norm(b[j]);
    }

    /*
    if(likeli_old >  -1e-1){

      cout <<"Likelihood > -1e-1" << endl;
	  
      return -1;
    }
    */
    if(maxSum <  current_kappa){
   
    
      beta = b;

      if(count%1000 == 0){
	  
	gn =0.;
	for(int j =0; j< work_J; j++){
	  gn += Norm(grad[j]);
	}

	//cout <<"   CONSTRAINT IS NOT  ACTIVE !!! "<< maxSum <<' '<<  current_kappa<< endl;


	likeli_new = Likelihood(work_index, beta,  y_minus_pi);
   

	//cout << 	likeli_new << ' ' << gn << endl;
	
	if( likeli_new - likeli_old <  -1e-4 && lambda_is_active){

	  cout  <<"============================ "<< "reset: " << iter <<' ' <<   likeli_old <<' ' <<  likeli_new <<' '<<  likeli_new-likeli_old;
	  // reset 
	  s /=2;
	  
	  cout <<  " new stepsize: " << s << endl;
	  beta = last_beta;
	  y_minus_pi =  last_y_minus_pi;
	  iter = 0;
	  
	  
	  continue;
	} 
	
	
	if(gn < 1e-7){
	  
	  cout <<"Error: constraint not active and gradient norm  < 1e-7" << endl;
	  
	  return -1;
	}
      }
      
      iter = 0;
      count++;
   
      continue;
    }
    
    //////////// project

   
    nzn = work_J;
    bool isNonNegative = false;
    while(isNonNegative == false){

      sumNorm =0;
      
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == false){
	  M.el(j) =0;
	  continue;
	}
	NormB.el(j) =  Norm(b[j]);
	sumNorm += 	NormB.el(j);
      }
      
      isNonNegative = true;
      double sumM =0;
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == false)
	  continue;
	
	M.el(j) = 	NormB.el(j)  + (current_kappa-sumNorm)/nzn; 
      
	sumM += M.el(j);

	if(M.el(j) < -1e-11)
	  isNonNegative = false;
      }

   
    
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == false)
	  continue;
	if(M.el(j) <1e-9){
	  activeSet[j] = false;
	
	  nzn--;
	
	
	  M.el(j) = 0;
	    
	}
      }
    }

    
  
    for(int j =0; j< work_J; j++){
    
      if(activeSet[j] == false){
	beta[j] = 0;
	continue;
      }
      beta[j] = b[j] * (M.el(j)/NormB.el(j));
    }
    
  
  
    

    ////////////
    
    
    My_Vector lambda(work_J);
    My_Vector v(work_J);
    lambda_avg =0;
    int nznC =0;

    // error 
    
    err = 100;

  
    
    if(iter%modulo == 0){
    


      y_minus_pi = y-Y_HAT(work_index, beta);
      
      for(int j =0; j<work_J; j++){
	grad[j] = GRAD_J_RESIDUAL(y_minus_pi, work_index[j] );
      }
     
      /*
      for(int j =0; j<work_J; j++){
	grad[j] = GRAD_J(work_index, beta,j);
      }
      */

    
      likeli_new = Likelihood(work_index, beta,y_minus_pi);

      
      cout<<"                               LOG-LIKELIHOOD: "
	  <<   likeli_old <<' ' 
	  <<  likeli_new <<  ' ' << lambda_is_active <<" ### "<< endl;



      double maxB =0;
      for(int j =0; j< work_J; j++){
	maxB +=  Norm(beta[j]);

      }
      cout << "maxB: " << maxB << ' '<<current_kappa<<endl;
    
      
      if( likeli_new - likeli_old <  -1e-4 && lambda_is_active){

	cout  <<"============================ "<< "reset: " << iter <<' ' <<   likeli_old <<' ' <<  likeli_new <<' '<<  likeli_new-likeli_old;
	// reset 
	s /=2;
	
	cout <<  " new stepsize: " << s << endl;
	beta = last_beta;
	y_minus_pi = last_y_minus_pi;
	iter = 0;
	
	
	continue;
      }

      last_beta = beta;
      last_y_minus_pi = y_minus_pi;
      likeli_old =  likeli_new;
  
      nznC =0;
      
      double max_lambda = 0.;
      for(int j =0; j< work_J;j++){

	 lambda[j] = Norm(grad[j]);
	 
	
	if(activeSet[j] == true){

	  lambda_avg +=  lambda[j];
	  nznC++;
	}

	if(lambda[j] > max_lambda)
	  max_lambda = lambda[j];
	

      }
     
      lambda_avg  /= nznC;
      double err_act =0;
    
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == false )
	  continue;
        err_act += fabs( lambda[j] -lambda_avg);

	
      }
      
       err_act /= (lambda_avg*nznC);
    
    


       cout << "err_act: "<< err_act << ' ' << lambda_avg << ' '<< nznC << ' ' << iter << endl;
       for(int j =0; j< work_J; j++){
	 if(activeSet[j] == true)
	   cout <<  Norm(grad[j]) <<'\t';
       }      
       cout << endl;
      
      double err_inact=0;
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == true)
	  continue;


	double v = Norm(grad[j])/lambda_avg;
	
	if(v-1>0)
	  err_inact += (v-1);
      }
    
      
      
      cout << "err_inact: "<< err_inact << endl;
      
      err =  err_act +  err_inact;

     
    }
  
  }
  

  double maxB =0;
  for(int j =0; j< work_J; j++){
    maxB += Norm(beta[j]);
  }
  cout << "maxB: " << maxB << ' '<<current_kappa<<endl;


  // adjust work_index and beta
   
  for (int j = 0; j < work_J; j++) {
   
     if( activeSet[j] == false){
    
      
       work_index.erase(work_index.begin() + j);
       beta.erase(beta.begin() + j);
       activeSet.erase(activeSet.begin() + j);
       work_J--;
       j--;
     }
  }

  
 
  return  lambda_avg;
 }



////// Projection l_1/infty ///////////////////////////

double L_ONE_INFTY::BSparse(vector<int>& work_index, 
				  const double& current_kappa,
				  vector<My_Vector>& beta, 
				  My_Vector&  y_minus_pi,
				  double s){
  
 
  int work_J = work_index.size();
 
  vector<bool> activeSet(work_J,true); 
 
 
 
  int d =0;
  for(int j =0; j<work_J; j++){
    d += beta[j].Length();
  }
 
  vector<My_Vector> b(work_J);



  double s_buff = s;
  double sumNorm;
  double gn;
  int iter;
  double err = 1e100;
 
  double  lambda_avg;

  iter = 0;

 
  vector<My_Vector> grad(work_J);
  double likeli_new;
 

  cout << "work_index: " <<' ';
  for(int j =0; j< work_J; j++)
    cout << work_index[j] << ' ';
  cout << endl;

  double  likeli_old = Likelihood(work_index, beta, y_minus_pi);
  My_Vector  last_y_minus_pi =  y_minus_pi;
  vector<My_Vector> last_beta = beta;
  
  vector<My_Vector> u(work_J);
 
  int count =0;

  cout << "starting iterations.... "<< likeli_old <<endl;
  while(err > errtol &&  iter < max_iter){
    
  
    

    iter++;

    y_minus_pi = y-Y_HAT(work_index, beta);
      
    for(int j =0; j<work_J; j++){
      grad[j] = GRAD_J_RESIDUAL(y_minus_pi, work_index[j] );
    }
    /*
    for(int j =0; j<work_J; j++){
      grad[j] = GRAD_J(work_index, beta,j);
    }
    */
    
    for(int j =0; j< work_J; j++){
      activeSet[j] = true;
    }   

   
    double maxSum =0;
  
   
  
    for(int j =0; j<work_J; j++)
      b[j] = beta[j] - grad[j]*s; 
  
   
   
    vector<My_Vector> beta_old = beta;
  

    // check constraint
    maxSum =0;
    for(int j =0; j< work_J; j++){
      maxSum += Norm_infty(b[j]);
    }

 
   
    if(maxSum <  current_kappa){
   
    
      beta = b;

      if(count%1000 == 0){
	  
	gn =0.;
	for(int j =0; j< work_J; j++){
	  gn += Norm(grad[j]);
	}

	cout <<"   CONSTRAINT IS NOT  ACTIVE !!! "<< maxSum <<' '<<  current_kappa<< endl;


	likeli_new = Likelihood(work_index, beta,  y_minus_pi);
   

	cout << 	likeli_new << ' ' << gn << endl;
	
	if( likeli_new - likeli_old <  -1e-4 && lambda_is_active){

	  cout  <<"============================ "<< "reset: " << iter <<' ' <<   likeli_old <<' ' <<  likeli_new <<' '<<  likeli_new-likeli_old;
	  // reset 
	  s /=2;
	  
	  cout <<  " new stepsize: " << s << endl;
	  beta = last_beta;
	  y_minus_pi =  last_y_minus_pi;
	  iter = 0;
	  
	  
	  continue;
	} 
	
	
	if(gn < 1e-6){
	  
	  cout <<"Error: constraint not active and gradient norm  < 1e-6" << endl;
	  
	  return -1;
	}
      }
      
      iter = 0;
      count++;
   
      continue;
    }
    
 
    
    // now project on l1_infty-ball

  
   
    sortPerm* R_t = new sortPerm[d+work_J];
    int rc =0;
    sumNorm =0;
    for(int j =0; j< work_J; j++){
      u[j] = ABS_SORT(b[j]);  // sort the absolute values of bV in each group
     
    
      R_t[rc].dPart = 0;
      R_t[rc].iPart = j;
      rc++;

      double sum_u =0.0;

      for(int i = 0; i < beta[j].Length()-1; i++){
	sum_u +=  u[j][i];


	R_t[rc].dPart = sum_u - (i+1) *  u[j][i+1];
	R_t[rc].iPart = j;
	rc++;
      }
      sum_u +=  u[j][beta[j].Length()-1];
     

      R_t[rc].dPart = sum_u;
      R_t[rc].iPart = j;
      rc++;
     
      sumNorm +=   Norm_infty(b[j]);
      
    }
  
 
  

    qsort(R_t,d+work_J,sizeof(struct sortPerm),CStructIncreasing);
  
    
    vector<int> j_occurences(work_J,0);
  
    int t_break = 0;
    double n_prev =  sumNorm;
    double sum_G,sum_G_old;
    double norm, next_norm;
    int idx = 0;

    while(R_t[idx].dPart < 1e-12){
      j_occurences[R_t[idx].iPart]++;
      idx++;
    }
  
    idx--;

    sum_G = 0.0;
    for(int j =0; j<work_J;j++){
    
      sum_G += (-1.0/j_occurences[j]);
    }
    norm = sumNorm;
    next_norm = n_prev +     sum_G * (R_t[idx+1].dPart);
  

    vector<bool> j_flag(work_J,true);

    while(next_norm > current_kappa){
    
      idx++;
      norm = next_norm;
      j_occurences[R_t[idx].iPart]++;

    

      sum_G = 0.0;
      for(int j =0; j<work_J;j++){
	int nj = beta[j].Length();
	if( j_flag[j] == true && j_occurences[j] <= nj)
	  sum_G += (-1.0/j_occurences[j]);
	else{
	  j_flag[j] = false;
	  j_occurences[j] = nj;
	}

      }

     
      next_norm = norm +  sum_G * (R_t[idx+1].dPart- R_t[idx].dPart);
 
    }
  
    double theta = R_t[idx].dPart + (current_kappa-norm)/sum_G;
   

    vector<double> mu(work_J);
    for(int j =0; j<work_J;j++){
      double sum_u = 0.;
     
            
      for(int i =0; i<j_occurences[j]; i++)
	sum_u += u[j][i];
      

      mu[j] = (sum_u-theta)/(j_occurences[j]);
      if( mu[j] <0)
	mu[j] =0;
    }

    for(int j =0; j< work_J; j++){
     
      if( mu[j] < 1e-12){
	activeSet[j] = false;
	beta[j] = 0;
      }
      else{

	for(int i = 0; i<beta[j].Length(); i++){
	  if(fabs(b[j][i]) >= mu[j]){
	    beta[j][i] =  mu[j];
	    if(b[j][i] < 0)
	      beta[j][i] = -mu[j];
	  }
	  else{
	    beta[j][i] = b[j][i]; 
	  }
	 
	}

      }
     
    }
 
   
  
    delete [] R_t;
   
   
    My_Vector lambda(work_J);
    My_Vector v(work_J);
    lambda_avg =0;
    int nznC =0;

    // error 
    
    err = 100;

  
    
    if(iter%modulo == 0){
    

      y_minus_pi = y-Y_HAT(work_index, beta);
      
      for(int j =0; j<work_J; j++){
	grad[j] = GRAD_J_RESIDUAL(y_minus_pi, work_index[j] );
      }
      /*
      for(int j =0; j<work_J; j++){
	grad[j] = GRAD_J(work_index, beta,j);
      }
      */

    
      likeli_new = Likelihood(work_index, beta,y_minus_pi);

      
      cout<<"                               LOG-LIKELIHOOD: "
	  <<   likeli_old <<' ' 
	  <<  likeli_new <<  ' ' << lambda_is_active <<" ### "<< endl;



      double maxB =0;
      for(int j =0; j< work_J; j++){
	maxB +=  Norm_infty(beta[j]);

      }
      cout << "maxB: " << maxB << ' '<<current_kappa<<endl;
    
      
      if( likeli_new - likeli_old <  -1e-4 && lambda_is_active){

	cout  <<"============================ "<< "reset: " << iter <<' ' <<   likeli_old <<' ' <<  likeli_new <<' '<<  likeli_new-likeli_old;
	// reset 
	s /=2;
	
	cout <<  " new stepsize: " << s << endl;
	beta = last_beta;
	y_minus_pi = last_y_minus_pi;
	iter = 0;
	
	
	continue;
      }

      last_beta = beta;
      last_y_minus_pi = y_minus_pi;
      likeli_old =  likeli_new;
  
      nznC =0;
      
      double max_lambda = 0.;
      for(int j =0; j< work_J;j++){

	 lambda[j] = Norm_L1(grad[j]);
	 
	
	if(activeSet[j] == true){

	  lambda_avg +=  lambda[j];
	  nznC++;
	}

	if(lambda[j] > max_lambda)
	  max_lambda = lambda[j];
	

      }
     
      lambda_avg  /= nznC;
      double err_act =0;
    
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == false )
	  continue;
        err_act += fabs( lambda[j] -lambda_avg);

	
      }
      
       err_act /= (lambda_avg*nznC);
    
    


      cout << "err_act: "<< err_act << ' ' << lambda_avg << endl;
      
      
      double err_inact=0;
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == true)
	  continue;

	
	double v = Norm_L1(grad[j])/lambda_avg;
	
	if(v-1>0)
	  err_inact += (v-1);
      }
      
      
      cout << "err_inact: "<< err_inact << endl;
      
      err =  err_act +  err_inact;

     
    }
  
  }
  

  double maxB =0;
  for(int j =0; j< work_J; j++){
    maxB += Norm_infty(beta[j]);
  }
  cout << "maxB: " << maxB << ' '<<current_kappa<<endl;


  // adjust work_index and beta
   
  for (int j = 0; j < work_J; j++) {
   
     if( activeSet[j] == false){
    
      
       work_index.erase(work_index.begin() + j);
       beta.erase(beta.begin() + j);
       activeSet.erase(activeSet.begin() + j);
       work_J--;
       j--;
     }
  }

  
 
  return  lambda_avg;
 }




////// Projection l_1/p ///////////////////////////
double L_ONE_INFTY::BSparse_LP(vector<int>& work_index, 
				  const double& current_kappa,
				  vector<My_Vector>& beta, 
				  My_Vector&  y_minus_pi,
				  double s,double p,double gamma){


  int work_J = work_index.size();

  My_Vector M(work_J);
  My_Vector NormB(work_J);
  My_Vector pNormB(work_J);
  int nzn;

  vector<bool> activeSet(work_J,true); 

 
 
 
  int d =0;
  for(int j =0; j<work_J; j++){
    d += beta[j].Length();
  }
 
  vector<My_Vector> b(work_J);
  vector<My_Vector> b2(work_J);


  double s_buff = s;
  double sumNorm;
  double sumNorm2;
  double psumNorm2;
  double gn;
  int iter;
  double err = 1e100;
 
  double  lambda_avg;

  iter = 0;

 
  vector<My_Vector> grad(work_J);
  vector<My_Vector> grad2(work_J);
  double likeli_new;
 

  cout << "work_index: " <<' ';
  for(int j =0; j< work_J; j++)
    cout << work_index[j] << ' ';
  cout << endl;

  double  likeli_old = Likelihood(work_index, beta, y_minus_pi);
  My_Vector  last_y_minus_pi =  y_minus_pi;
  vector<My_Vector> last_beta = beta;
  
  vector<My_Vector> u(work_J);
 
  int count =0;

  cout << "starting iterations.... "<< likeli_old <<endl;
  while(err > errtol &&  iter < max_iter){
    

    

    iter++;

    y_minus_pi = y-Y_HAT(work_index, beta);
  
    for(int j =0; j<work_J; j++){
      // grad[j] = GRAD_J(work_index, beta,j);
      grad[j] = GRAD_J_RESIDUAL(y_minus_pi, work_index[j] );
    }
     
    
    for(int j =0; j< work_J; j++){
      activeSet[j] = true;
    }   

   
    double maxSum =0;
  
   
  
    for(int j =0; j<work_J; j++){
      b[j] = beta[j] - grad[j]*s; 
    }
 


    vector<My_Vector> beta_old = beta;
  

    // check constraint
    maxSum =0;
    for(int j =0; j< work_J; j++){
      maxSum += pNorm(b[j],p);
    }
      
  

    if(maxSum <  current_kappa){
    
  
      beta = b;

    
      if(count%1000 == 0){
	gn =0.;
	for(int j =0; j< work_J; j++){
	  gn += Norm(grad[j]);
	}

	cout <<"   CONSTRAINT IS NOT  ACTIVE !!! "<< maxSum <<' '<<  current_kappa<< endl;


	likeli_new = Likelihood(work_index, beta,  y_minus_pi);
   

	cout << 	likeli_new << ' ' << gn << endl;
	
	if( likeli_new - likeli_old <  -1e-4 && lambda_is_active){

	  cout  <<"============================ "<< "reset: " << iter <<' ' <<   likeli_old <<' ' <<  likeli_new <<' '<<  likeli_new-likeli_old;
	  // reset 
	  s /=2;
	  
	  cout <<  " new stepsize: " << s << endl;
	  beta = last_beta;
	  y_minus_pi =  last_y_minus_pi;
	  iter = 0;
	  
	  
	  continue;
	} 
	
	
	if(gn < 1e-7){
	  
	  cout <<"Error: constraint not active and gradient norm  < 1e-7" << endl;
	  
	  return -1;
	}
      }
      
      iter = 0;
      count++;
   
      continue;
    }
   
    //////////// project




    nzn = work_J;
    int nzn_buffer = nzn;

  
   
    vector<My_Vector> old_beta = beta;
   
    double sumnorm=0;


    
    double Lagrange_mu;
   
    nzn = work_J;
      
    
    double mu_min=0;
    double mu_max;
    double r = 1./((2./p)-1);
    double pstar= 1+1./r;
  
    old_beta=beta;

    vector<My_Vector> tau = beta;
    vector<double> lambdagroup(work_J);
    vector<double> lambdasum(work_J);
    double sum;
      
      
      
    mu_max = 1.5;
      
   
    Lagrange_mu = ((mu_min+mu_max)/2e2);
      
      
      
    sumnorm=0;
      
    double diffbeta = 100;
      
    //adapt mu
    while(fabs(sumnorm-current_kappa)> 1e-3){

      // we now hav a new mu  => initialize beta, activeSet and nzn
    

      for(int j =0; j< work_J; j++){
	
	activeSet[j] = true;
	
	for(int i =0; i< beta[j].Length(); i++){
	  
	  if( fabs(beta[j][i]) < 1e-4 ) 
	    beta[j][i] = 1e-4;
	  
	}
	nzn = nzn_buffer;
      }
      

      diffbeta = 100;

      while(diffbeta>1e-6 ){

	old_beta = beta;

	sumnorm=0;
	for(int j =0; j< work_J; j++){	 
	  if(activeSet[j]==false){
	    continue;
	  }
	  sumnorm = sumnorm+pNorm(beta[j],p);
	 
	}

	


	diffbeta=0;
	
	for(int j =0; j< work_J; j++){
	  
	  if(activeSet[j]==false){
	    beta[j] = 0;
	    continue;
	  }
	  
	  sum=0;
	  for(int i =0; i< beta[j].Length(); i++){
	    sum=sum + pow(fabs(beta[j][i]),2./pstar);
	  }
	  
	  lambdasum[j] = sum;
	  lambdagroup[j] = pNorm(beta[j],p)/sumnorm;
	  
	  
	  for(int i =0; i< beta[j].Length(); i++){
	    tau[j][i]=pow(fabs(beta[j][i]),2.0-(2.0/pstar))/pow(lambdasum[j],pstar-1);
	    beta[j][i]=(1.0/(1.0+(Lagrange_mu/(tau[j][i]*lambdagroup[j]))))*b[j][i];
	  }
	  
	  
	  
	  
	  if( Norm_infty(beta[j]) < 1e-5){
	    activeSet[j] = false;
	    beta[j] = 0;
	    nzn--;
	  }

	  diffbeta=diffbeta+Norm_infty(beta[j]-old_beta[j]);
	  
	  
	 
	}
	
	diffbeta=diffbeta/work_J;
	//	cout << diffbeta <<endl;

	//end of inner while loop
      }
      // we now have a new beta


      	    
	
      sumnorm=0;
      for(int j =0; j< work_J; j++){
	if(activeSet[j]==false){
	  continue;
	}
	sumnorm= sumnorm+pNorm(beta[j],p);
      }
	
    
    
      if(sumnorm<current_kappa){
	//take smaller mu
	mu_max=Lagrange_mu;
	Lagrange_mu=((mu_max+mu_min)/2);

      }else{
	//take larger mu
	 mu_min=Lagrange_mu;
	 Lagrange_mu=((mu_min+mu_max)/2);
	
      }
      
      
       
     
      
    }
    // cout << Lagrange_mu <<' '<< sumnorm<<' '<< current_kappa<<endl;
    // cout << " ======================================================"<< endl;

      





    ////////////
    
    
    My_Vector lambda(work_J);
    My_Vector v(work_J);
    lambda_avg =0;
    int nznC =0;

    // error 
    
    err = 100;
    double help=iter%modulo;


    if(iter % modulo == 0){


      y_minus_pi = y-Y_HAT(work_index, beta);
      



      for(int j =0; j<work_J; j++){
	grad[j] = GRAD_J_RESIDUAL(y_minus_pi, work_index[j] );

      }


    
      likeli_new = Likelihood(work_index, beta,y_minus_pi);

      
      cout<<"                               LOG-LIKELIHOOD: "
	  <<   likeli_old <<' ' 
	  <<  likeli_new <<  ' ' << lambda_is_active <<" ### "<< endl;




      double maxB =0;
      for(int j =0; j< work_J; j++){
	maxB +=  pNorm(beta[j],p);

      }
      cout << "maxB: " << maxB << ' '<<current_kappa<< " we have " << nzn <<" active groups" << endl;
  

      if( likeli_new - likeli_old <  -1e-4 && lambda_is_active){

	cout  <<"============================ "<< "reset: " << iter <<' ' <<   likeli_old <<' ' <<  likeli_new <<' '<<  likeli_new-likeli_old;
	// reset 
	s /=2;
	
	cout <<  " new stepsize: " << s << endl;
	beta = last_beta;
	y_minus_pi = last_y_minus_pi;
	iter = 0;
	

	
	continue;
      }

      last_beta = beta;
      last_y_minus_pi = y_minus_pi;
      likeli_old =  likeli_new;
  
      nznC =0;
      
      double max_lambda = 0.;
      for(int j =0; j< work_J;j++){

	 lambda[j] = gammaNorm(grad[j],gamma);
 
	
	if(activeSet[j] == true){

	  lambda_avg +=  lambda[j];
	  nznC++;
	}

	if(lambda[j] > max_lambda)
	  max_lambda = lambda[j];
	

      }
     
      lambda_avg  /= nznC;


      double err_act =0;
    
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == false )
	  continue;
        err_act += fabs( lambda[j] -lambda_avg);

	
      }
      
       err_act /= (lambda_avg*nznC);
    
    


       cout << "err_act: "<< err_act << ' ' << lambda_avg << ' ' <<   Lagrange_mu <<endl;
      
      
      double err_inact=0;
      for(int j =0; j< work_J; j++){
	if(activeSet[j] == true)
	  continue;

	
	double v = gammaNorm(grad[j],gamma)/lambda_avg;
   
       
	

	if(v-1>0)
	  err_inact += (v-1);
      }
      
      
      cout << "err_inact: "<< err_inact << endl;

      
      err =  err_act +  err_inact;
 cout << "err: "<< err << endl;
     
    }


  }
  

  double maxB =0;
  for(int j =0; j< work_J; j++){
    maxB += pNorm(beta[j],p);
  }
  cout << "maxB: " << maxB << ' '<<current_kappa<<endl;

  // adjust work_index and beta
   
  for (int j = 0; j < work_J; j++) {
   
     if( activeSet[j] == false){
    
      
       work_index.erase(work_index.begin() + j);
       beta.erase(beta.begin() + j);
       activeSet.erase(activeSet.begin() + j);
       work_J--;
       j--;
     }
  }

  
 
  return  lambda_avg;
 }



int L_ONE_INFTY::optimize_current( vector<int>& work_index,  const double& current_kappa, vector<My_Vector>& beta,  double& lambda, My_Vector& residualV){


  
   // first optimization

   lambda_is_active = false;

   if(lambda>0)
     lambda_is_active = true;

    
    if(infty)
      lambda = BSparse(work_index,  current_kappa, beta,  residualV, stepsize/2);
    else{
      if(ell2){
	lambda = BSparse_L2(work_index,  current_kappa, beta,  residualV, stepsize/2);
      }else{
	lambda = BSparse_LP(work_index,  current_kappa, beta,  residualV, stepsize/2, Norm_p, Norm_gamma);
      }
    } 

   
   cout <<"first lambda: " << lambda << endl;
   
   if(lambda>0)
     lambda_is_active = true;
   
   int constraint_is_active;
   if(lambda < 0){
     
     // ensure active constraint
     double step =     stepsize/4.0;
     constraint_is_active = 
       ensure_constraint_active(lambda,
				work_index, 
				current_kappa, 
				beta, 
				residualV,
				step
				);
     
     
     if( constraint_is_active == -1)
       return 0;
     
     cout <<"now we should have a positive lambda: " << lambda << endl;
   }
   
   if(lambda>0)
     lambda_is_active = true;
 
   
   // now iterate until zero lagrangian violation
   
   int max_err_index;
   double err_inact,max_err_inact;
   vector<int>::iterator v_iter;
   vector<bool> guessed_active_old;
   
 
     
   max_err_inact = 1e100;
   int j_counter =0;
  
   
  while( max_err_inact > errtol){

    vector<bool> guessed_active(J_global, false);
    
     max_err_index = -1;
     err_inact = -1;
     max_err_inact=-1;

     My_Matrix_bool_scale G;
     double v;
     j_counter =0;
     
     for(int j =0; j< J_global; j++){
       v_iter = find(work_index.begin(), work_index.end(),j);
       if(v_iter !=  work_index.end())
	 continue;


       if(work_index.size()>1 && guessed_active_old.size() > 0)
	 if(guessed_active_old[j] == false)
	   continue;

      
       
     
       if(infty){
	 v = Norm_L1(GRAD_J_RESIDUAL(residualV, j ))/ lambda;
       }else{
	 if(ell2){ // l_2
	   v = Norm(GRAD_J_RESIDUAL(residualV, j ))/ lambda;	
	 }else{// l_p
	   v = gammaNorm(GRAD_J_RESIDUAL(residualV, j ),Norm_gamma)/ lambda;
	 }
       }
     
      
      
     
       if(v-1>0){
	 err_inact = v-1;

	 if(err_inact > max_err_inact){
	   max_err_inact = err_inact;
	   max_err_index =j;
	 }

	 sort_vect[ j_counter].dPart = err_inact;
	 sort_vect[ j_counter].iPart = j;
	 j_counter++;
	 guessed_active[j] = true;
       }
     }
     
    
     guessed_active_old =  guessed_active;

     cout <<"max violation of inactive groups: "<<  max_err_inact<< ' '<< lambda << endl;
    if(max_err_inact<0)
      break;
    
    qsort(sort_vect, j_counter,sizeof(struct sortPerm),CStruct);
  
     // update
  
    for(int k =0; k< update_groups; k++){
      if(k == j_counter)
	break;
      cout <<" Updating group: "<<   sort_vect[k].iPart  << endl;

      My_Vector empty(K);
      work_index.push_back(sort_vect[k].iPart);
      beta.push_back(empty);
    
    }
       
    
       // optimize
       
  
     
    if(infty)
      lambda = BSparse(work_index,  current_kappa, beta,  residualV, stepsize);
    else{
      if(ell2){
	lambda = BSparse_L2(work_index,  current_kappa, beta,  residualV, stepsize);
      }else{
	lambda = BSparse_LP(work_index,  current_kappa, beta,  residualV, stepsize, Norm_p, Norm_gamma);
      }
    }

     // cout << "lambda: " << lambda << endl;
   
     
    
     
     if(lambda <0){
       constraint_is_active = 
	 ensure_constraint_active(lambda,
				  work_index,
				  current_kappa, 
				  beta,   
				  residualV,
				  stepsize);
       
       if( constraint_is_active == -1)
	 return  1;
     }
     
  }
   
 

  cout << "we have " <<work_index.size() << " active groups."<<endl;

 
  
   return 1;
}

int L_ONE_INFTY::optimize(){

 
  sort_vect = new sortPerm[J_global];
 

  double act_kappa =0;
  double kappa_step = kappa/iterN;


  int first_group = init_idx[0]; 
  vector<int> work_index(1,first_group);
  
  My_Vector first_beta(K);
  first_beta =  0.25*(kappa_step/sqrt(K));
  for(int k =0; k<K;k++)
    first_beta.el(k) = rand()/(double)RAND_MAX;

  
  int second_group = init_idx[1]; // the 1st group
  work_index.push_back(second_group);
  
 

  int n = y.Length();
 

  vector<My_Vector> beta(work_index.size(),first_beta);

  
  // check constraint
  double maxSum =0;
  for(int j =0; j< work_index.size(); j++){
    // cout << beta[j] << endl;
    maxSum += Norm(beta[j]);
  }
  cout << "initial norm: "<< maxSum << endl;
  for(int j =0; j< work_index.size(); j++){
    beta[j] = beta[j] * (0.4/ maxSum);
    beta[j] = beta[j] * kappa_step;
  }
  maxSum =0;
  for(int j =0; j< work_index.size(); j++){
    // cout << beta[j] << endl;
    maxSum += Norm(beta[j]);
  }
  cout << "initial norm: "<< maxSum << endl;
  double lambda = -1;
  My_Vector residualV(n);

 

  for(int iter =0; iter <iterN; iter++){

    
     act_kappa += kappa_step;
  
     optimize_current( work_index,
		       act_kappa, 
		       beta,
		       lambda,
		       residualV);

     work_index_path.push_back(work_index);
     beta_path.push_back(beta);
     kappa_path.push_back(act_kappa);
     lambda_path.push_back(lambda);
     /*
     for(int i =0; i<beta.size();i++)
       cout << beta[i] <<"|";
     cout << endl;
     */
    
     // uniqueness
     vector<int>::iterator v_iter;
     vector<int> non_unique(1,-1);
     bool do_unique = false;
     if(do_unique){
       for(int j =0; j< J_global; j++){
	 v_iter = find(work_index.begin(), work_index.end(),j);
	 if(v_iter !=  work_index.end())
	   continue;
	 double v;
	

	 if(infty){
	   v = Norm_L1(GRAD_J_RESIDUAL(residualV, j ))/ lambda;
	 }else{
	   if(ell2){ // l_2
	     v = Norm(GRAD_J_RESIDUAL(residualV, j ))/ lambda;	
	   }else{// l_p
	     v = gammaNorm(GRAD_J_RESIDUAL(residualV, j ),Norm_gamma)/ lambda;
	   }
	 }


	 if(v>1-errtol){
	   non_unique.push_back(j);
	 }
       }
       cout << "Incomplete: ";
       for(int j=1;j< non_unique.size();j++){
	 cout << non_unique[j] <<' ';
       }
       cout << endl;
     }
     non_unique_path.push_back(non_unique);
     
     if( yV.Length() >0){
       double balErr, thresh;
       
       if(LOGREG){
	 vector<double> ErrVec = balancedERROR(work_index, beta);
	 balErr = ErrVec[0];
	 thresh = ErrVec[1];
	 cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>   Balanced Error: "<< balErr << endl;
       }else{
	 balErr = squaredERROR(work_index, beta);
	 thresh = -1;
	 cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>   Squared Error: "<< balErr << endl;
       }
       balErr_path.push_back(balErr);
       bestThresh_path.push_back(thresh);
     }

     if(max_groups > 0 && work_index.size() > max_groups){
       for(int j = iter+1; j<iterN; j++){
	 balErr_path.push_back(1e100);
	 work_index_path.push_back(work_index);
	 beta_path.push_back(beta);
	 kappa_path.push_back(act_kappa);
	 lambda_path.push_back(lambda);
	 non_unique_path.push_back(non_unique);
       }
       break;
     }

  }

  delete [] sort_vect;
}

vector<double> L_ONE_INFTY::balancedERROR(const vector<int>& work_index, const vector<My_Vector>& beta){

  vector<double> ErrVec(2,1e100);
  int nV = yV.Length();
  if(nV <1)
    return ErrVec;

  double thresh,best_thresh;
 
 
  
  My_Vector yVhat(nV);
  My_Vector res(nV);
  int work_J = work_index.size();
 
 
  yVhat = A_multiTask_times_x(XV_global, StartStopV_global, work_index, beta);

  for(int k=0; k<nV;k++){
    res.el(k) =  1.0/(1.0+exp(-yVhat.el(k)));
  }
  double balERR_min = 2.0;
  double balERR = 1.0;
  int t_iter = 25;
  thresh = 0.05;
  best_thresh = thresh;
  double t_inc = (0.95-thresh)/t_iter;
  
  for(int t =0; t< t_iter; t++){
    int t1 =0, f1 =0, t2 =0, f2 =0;
    for(int k=0; k<nV;k++){
      //res =  1.0/(1.0+exp(-yVhat.el(k)));
      //cout << res << " " << yV.el(k) << endl;
      if(yV.el(k) > 0.5){
	if( res.el(k) >  thresh)
	  t1++;
	else
	  f1++;
      }
      else{
	if(res.el(k) <=  thresh)
	  t2++;
	else
	  f2++;
      }
    }
    balERR = 0.5*((double)f1/(t1+f1) + (double)f2/(t2+f2));
    if( balERR_min > balERR){
      balERR_min =  balERR;
      best_thresh = thresh;
    }
    thresh += t_inc;
  }

   
  ErrVec[0] =  balERR_min;
  ErrVec[1] =  best_thresh;
  return  ErrVec;
  
  
 
  
}

/*
 vector<double> L_ONE_INFTY::balancedERROR_perTask(const vector<int>& work_index, const vector<My_Vector>& beta){
   vector<double> E_0(K,1e100)
   
   vector< vector<double> >ErrVec(2,E_0);
   int nT = yV.Length();
   if(nT <1)
     return ErrVec;
   
  double thresh,best_thresh;
 
 
  
  My_Vector yhat(nT);
  My_Vector res(nT);
  int work_J = work_index.size();
  for(int j =0; j<work_J; j++)
    yhat = yhat + A_shifted_times_x(XV_global[work_index[j]],beta[j]);
 
  for(int k=0; k<nT;k++){
    res.el(k) =  1.0/(1.0+exp(-yhat.el(k)));
  }
  double balERR_min = 2.0;
  double balERR = 1.0;
  int t_iter = 25;
  
  double t_inc = (0.95-thresh)/t_iter;
  double res_d;
  int XTr =  XT[0].Nrows();
   
  for(int k=0; k<K;k++){
  
    thresh = 0.05;
    best_thresh = thresh;

    for(int t =0; t< t_iter; t++){
      int t1 =0, f1 =0, t2 =0, f2 =0;
      for(int j=0; j<XTr;j++){
	res_d =   res.el(k*XTr +j) // 1.0/(1.0+exp(-yThat.el(k*XTr +j)));
	  
	  if(yV.el(k*XTr +j) > 0.5){
	    if( res_d >  thresh)
	      t1++;
	    else
	      f1++;
	  }
	  else{
	    if(res_d <=  thresh)
	      t2++;
	    else
	      f2++;
	  }
      }
      
      balERR = 0.5*((double)f1/(t1+f1) + (double)f2/(t2+f2));
      if( balERR_min > balERR){
	balERR_min =  balERR;
	best_thresh = thresh;
      }
      thresh += t_inc;
    }
    ERR[k] = balERR_min;
    THR[k] = best_thresh; 
  }
  ErrVec[0] = ERR;
  ErrVec[1] = THR;
  return  ErrVec;
 }
*/

double L_ONE_INFTY::squaredERROR(const vector<int>& work_index, const vector<My_Vector>& beta){

  int nV = yV.Length();
  if(nV <1)
    return 1e100;

  double thresh,best_thresh;
  double res;
 
  
  My_Vector yVhat(nV);
 
 
  yVhat = A_multiTask_times_x(XV_global, StartStopV_global, work_index, beta);

 
    double sum = 0;
    for(int k=0; k<nV;k++){
     
      sum += pow(yV.el(k)- yVhat.el(k),2);
     
    }
    
    return sum/double(nV);
 
}

/*
double L_ONE_INFTY::correlationCoeff(const vector<int>& work_index, const vector<My_Vector>& beta){

  int nT = yV.Length();
  if(nT <1)
    return 1e100;


  My_Vector yhat(nT);
  int work_J = work_index.size();
  for(int j =0; j<work_J; j++)
    yhat = yhat + A_shifted_times_x(XV_global[work_index[j]],beta[j]);
 
  
  
  double sum_sq_x = 0;
  double sum_sq_y = 0;
  double sum_coproduct = 0;
  double mean_x = yhat.el(0);
  double mean_y = yV.el(0);
  double sweep,delta_x,delta_y;
  for(int i=1;i<nV;i++){
    sweep = i / (i+1.0);
    delta_x =  yhat.el(i) - mean_x;
    delta_y = yV.el(i) - mean_y;
    sum_sq_x += delta_x * delta_x * sweep;
    sum_sq_y += delta_y * delta_y * sweep;
    sum_coproduct += delta_x * delta_y * sweep;
    mean_x += delta_x / (i+1);
    mean_y += delta_y / (i+1) ;
  }
  double pop_sd_x = sqrt( sum_sq_x / nV);
  double pop_sd_y = sqrt( sum_sq_y / nV );
  double cov_x_y = sum_coproduct / nV;
  double correlation = cov_x_y / (pop_sd_x * pop_sd_y);
  return correlation;
}
*/
double  L_ONE_INFTY::predict_with_best_beta(const My_Matrix & XT,  const My_Vector & targetsT, const vector<vector <int> >& StartStopT){

   int nT = targetsT.Length();
   if(nT <1 || balErr_path.size()<1)
     return -1;

   double thresh;
   double res;
   
   double err_best = balErr_path[0];
   thresh =  bestThresh_path[0];
   vector<My_Vector> beta_best = beta_path[0];
   vector<int> work_index_best = work_index_path[0];

   for(int k =1; k<balErr_path.size();k++)
     if(balErr_path[k] < err_best){
        err_best = balErr_path[k];
	beta_best = beta_path[k];
	work_index_best = work_index_path[k];
	thresh =  bestThresh_path[k];
     }
   My_Vector yThat(nT);
 
   yThat = A_multiTask_times_x(XT, StartStopT, work_index_best, beta_best);

   if(LOGREG){
     int t1 =0, f1 =0, t2 =0, f2 =0;
     for(int k=0; k<nT;k++){
       res =  1.0/(1.0+exp(-yThat.el(k)));
       
       if(targetsT.el(k) > 0.5){
	 if( res >  thresh)
	   t1++;
	 else
	   f1++;
       }
       else{
	 if(res <=  thresh)
	   t2++;
	 else
	   f2++;
       }
     }
     double balERR = 0.5*((double)f1/(t1+f1+1e-6) + (double)f2/(t2+f2+1e-6));
     
     return balERR;
   }
   else{
     
     double sum = 0;
     for(int k=0; k<nT;k++){
     
       sum += pow(targetsT.el(k)- yThat.el(k),2);
      
     }
    
     return sum/double(nT);
   }  
}

vector<double>  L_ONE_INFTY::predict_with_best_beta_perTask(const My_Matrix & XT,  const My_Vector & targetsT, const vector<vector <int> >& StartStopT){

  
   vector<double> balERRvec(K,-1);
   int nT = targetsT.Length();
   if(nT <1 || balErr_path.size()<1)
     return  balERRvec;

   double thresh =  bestThresh_path[0];;
   double res;

  
   
   double err_best = balErr_path[0];
   vector<My_Vector> beta_best = beta_path[0];
   vector<int> work_index_best = work_index_path[0];

   for(int k =1; k<balErr_path.size();k++)
     if(balErr_path[k] < err_best){
        err_best = balErr_path[k];
	beta_best = beta_path[k];
	work_index_best = work_index_path[k];
	thresh =  bestThresh_path[k];
     }
   My_Vector yThat(nT);
 
   yThat = A_multiTask_times_x(XT, StartStopT, work_index_best, beta_best);

   if(LOGREG){
   
   
     for(int k=0; k<K;k++){
       int t1 =0, f1 =0, t2 =0, f2 =0;
       for(int j=StartStopT[k][0]; j<=StartStopT[k][1];j++){
	 res =  1.0/(1.0+exp(-yThat.el(j)));
       
	 if(targetsT.el(j) > 0.5){
	   if( res >  thresh)
	     t1++;
	   else
	     f1++;
	 }
	 else{
	   if(res <=  thresh)
	     t2++;
	   else
	     f2++;
	 }
       }
       double balERR = 0.5*((double)f1/(t1+f1+1e-6) + (double)f2/(t2+f2+1e-6));
       balERRvec[k] = balERR;
     }
     return balERRvec;
   }
   else{
    
   
     for(int k=0; k<K;k++){
     
       double sum = 0;
       for(int j=StartStopT[k][0]; j<=StartStopT[k][1];j++){

	 sum += pow(targetsT.el(j)- yThat.el(j),2);
	
       }
      
       balERRvec[k] =  sum/double(StartStopT[k][1]-StartStopT[k][0]+1);
     }
     return balERRvec;
   }  
}
void  L_ONE_INFTY::Rplot(){
  
    ofstream outWorkIndex("Idx.asc",ios::out);
      for(int i=0; i< work_index_path.size(); i++){
      int s = work_index_path[i].size();
      
      for(int j =0; j<s-1; j++){
	//cout << work_index_path[i][j]<<'\t';
	outWorkIndex << work_index_path[i][j]<<'\t';
      }
      outWorkIndex << work_index_path[i][s-1]<< endl;
      //cout << work_index_path[i][s-1]<< endl;
    } 
    outWorkIndex.close();
    
    ofstream outBetaNorm("Norm.asc",ios::out);
    for(int i=0; i<beta_path.size();i++){
      for(int j =0; j<beta_path[i].size()-1;j++){
	if(infty)
	  outBetaNorm << Norm_infty(beta_path[i][j])<<'\t';
	else{
	  if(ell2){
	    outBetaNorm << Norm(beta_path[i][j])<<'\t';
	  }else{
	    outBetaNorm << pNorm(beta_path[i][j], Norm_p)<<'\t';
	  }
	}
      }
      if(infty)
	outBetaNorm << Norm_infty(beta_path[i][beta_path[i].size()-1])<< endl;
      else{
	if(ell2){
	  outBetaNorm << Norm(beta_path[i][beta_path[i].size()-1])<< endl;
	}else{
	  outBetaNorm << pNorm(beta_path[i][beta_path[i].size()-1], Norm_p)<< endl;
	}
      }
    } 
    outBetaNorm.close();
    
    ofstream outErr("BalErr.asc",ios::out);
    for(int j =0; j<balErr_path.size()-1;j++)
      outErr << balErr_path[j] <<'\t';
    outErr << balErr_path[balErr_path.size()-1]<<endl;
    outErr.close();
    
    ofstream outKappa("Kappa.asc",ios::out);
    for(int j =0; j<kappa_path.size()-1;j++)
      outKappa << kappa_path[j] <<'\t';
    outKappa << kappa_path[kappa_path.size()-1]<<endl;
    outKappa.close(); 
    
    ofstream outInfty("Infty.asc",ios::out);
    outInfty << infty << endl;
    outInfty.close();

    ofstream outEll2("Ell2.asc",ios::out);
    outEll2 << ell2 << endl;
    outEll2.close();
    
    system("R CMD BATCH ReadOut.R");
}
