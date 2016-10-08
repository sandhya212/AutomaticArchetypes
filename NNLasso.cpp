#include "NNLasso.h"

int NNLasso_small(vector<double>& y, Matrix<double>& X, double kap, double eps, Matrix<double>& XTranX,Matrix<double> &beta ){

    int d = X.cols();
    int n = y.size();
    double kappa = kap;
    double epsilon = eps;
    Matrix<double> r(n,1); //check if this is a vector
    //r = y
   // Matrix<double> beta(X.cols(),1,0); send from call
    Matrix<double> y_hat(y.size(),1);/////check

    int nonzero, l1 = 0;
    int j = -1;
    Matrix<double> correl(X.cols(),1); ///////check
    
    while(l1 < (kappa-epsilon) && (nonzero < min(d,n)+1)){      // && (*(Norm(r),Norm(r)) > 1e-10)){ 
	if (j < 0){
	    correl = t(X) * r; //check
	}	
	else{
	    correl = correl - epsilon * XTranX(_,j);
	}// end if
        if (max(correl) < 1e-10){
	  break;
	}
	j = maxind(correl);
	beta[j] = beta[j] + epsilon;
	r = r - epsilon * X(_,j);
	y_hat = y_hat + epsilon * X(_,j);
	l1 = l1 + epsilon;
	
    } //end while

    return 0;

}// end function

//vector<vector<double> > NNLasso(vector<double>& y, vector<vector<double> >& X, double& kappa, double& epsilon,vector<vector<double> > XTran){
///*My_Vector NNLasso(My_Vector& y, My_Matrix& X, double& kappa, double& epsilon,My_Matrix& XTranX){

//     int d = X.Ncols();
//     int n = y.Length();
//     My_Vector r(y);
//     My_Vector beta(d);
//     My_Vector yhat(n);
// 
//     int nonzero,l1 = 0;
//     int j = -1;
//    // vector<double> correl(d,0);
//     double correl[d];
//     while(l1 < (kappa-epsilon) && (nonzero < min_element(d,n)+1) && (*(Norm(r),Norm(r)) > 1e-10)){ 
// 	if (j < 0){
// 	    correl = Transpose(X) * r; //check
// 	}
// 	else{
// 	  //  correl = correl - epsilon * XTranX[j];//wrong
//         correl -=  *(epsilon,XTranX[j]);//wrong
// 	}//end if
// 	double* first = correl;
// 	double* last = correl + d -1;
// 	if (max_element(first,last) < 1e-10){
// 	  break;
// 	}
// 
// 	//j = correlmax
// 	beta[j] = beta[j] + epsilon;
// 	r = r - *(epsilon, X[j]); //wrong
// 	yhat = yhat + *(epsilon, X[j]); //wrong
// 	l1 = l1 + epsilon;
// 	//nonzero = 
// 
//   
//     }//end while
//     return beta;*/
//}


// My_Vector NNLassolarge(My_Vector& y, My_Matrix& X, double& kappa, double& epsilon){
// 
//     int d = X.Ncols();
//     int n = y.Length();
//     My_Vector r(y);
//     My_Vector beta(d);
//     My_Vector yhat(n);
// 
//     int nonzero,l1 = 0;
//     int j = -1;
//     //My_Vector correl(d);
//     vector<double> correl(d,0);
//   
//     while(l1 < (kappa-epsilon) && (nonzero < min_element(d,n)+1) && (*(Norm(r),Norm(r)) > 1e-10)){ 
// 
//          correl = Transpose(X) * r;
// //          if (correl.Max() < 1e-10){
// // 	    break;
// // 	 }
// 	 //j = correlmax
// 	 beta[j] = beta[j] + epsilon;
// 	 r = r - *(epsilon, X[j]); //wrong
// 	 yhat = yhat + *(epsilon, X[j]); //wrong
// 	 l1 = l1 + epsilon;	
// 
//     }//end while
//     return beta;
// }