#ifndef L_ONE_INFTY_HEADER
#define L_ONE_INFTY_HEADER


#include "my_mat.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>


class My_Vector;
class My_Matrix;

using namespace std;

static My_Matrix DUMMY_MAT(0,0);
static My_Vector DUMMY_INT_VEC(0);
static vector<vector<int> > DUMMY_INT_MAT(0); 

static int my_init_indices[] = {0,1};
static vector<int>  init_idx_default (my_init_indices, my_init_indices + sizeof(my_init_indices) / sizeof(int) );




//typedef std::vector<My_Matrix> VecOfMyMat; 

class L_ONE_INFTY {
  
  friend class My_Matrix;
  friend class My_Vector;
  
  
 public:
  
 L_ONE_INFTY( const My_Matrix & X,  const My_Vector & targets,  const vector<vector <int> >& StartStop):X_global(X),y(targets),StartStop_global(StartStop),XV_global(DUMMY_MAT), yV(DUMMY_INT_VEC),StartStopV_global( DUMMY_INT_MAT){
    J_global = X_global.Ncols();
    kappa = 30;
    errtol = 1e-2;
    stepsize = 2e-3;
    iterN = 50;
    modulo = 1000;
    update_groups = 2;
    infty = true;
    max_groups = -1;
    LOGREG =  true;
    max_iter = 200000;
    N_tot = y.Length();
    K = StartStop_global.size();
    init_idx =  init_idx_default;
  };

  L_ONE_INFTY( const My_Matrix & X,  const My_Vector & targets,  const vector<vector <int> >& StartStop,  bool _infty,   int _max_groups, vector<int> init_idx_): X_global(X),y(targets),StartStop_global(StartStop),XV_global(DUMMY_MAT), yV(DUMMY_INT_VEC),StartStopV_global( DUMMY_INT_MAT),  infty(_infty), ell2(!_infty){
    J_global = X_global.Ncols();
    kappa = 30;
    errtol = 1e-2;
    stepsize = 2e-3;
    iterN = 50;
    modulo = 1000;
    update_groups = 2;
    max_groups = _max_groups;
   
    LOGREG =  true;
    max_iter = 200000;
    N_tot = y.Length();
    K = StartStop_global.size();
    init_idx =  init_idx_;
  };

 L_ONE_INFTY( const My_Matrix & X,  const My_Vector & targets,  const vector<vector <int> >& StartStop, 
	      const My_Matrix & XV,  const My_Vector & targetsV, const vector<vector <int> >& StartStopV)
   :X_global(X),y(targets),StartStop_global(StartStop),XV_global(XV),yV(targetsV),StartStopV_global(StartStopV){
    J_global = X_global.Ncols();
    kappa = 30;
    errtol = 1e-2;
    stepsize = 2e-3;
    iterN = 50;
    modulo = 1000;
    update_groups = 2;
    infty = true;
    ell2 = false;
    max_groups = -1;
    LOGREG =  true;
    max_iter = 200000;
    N_tot = y.Length();
    K = StartStop_global.size();
    init_idx =  init_idx_default;
  };

 L_ONE_INFTY( const My_Matrix & X,  const My_Vector & targets,  const vector<vector <int> >& StartStop,  
	      const My_Matrix & XV,  const My_Vector & targetsV, const vector<vector <int> >& StartStopV,
	      bool _infty, bool _ell2)
   :X_global(X),y(targets),StartStop_global(StartStop),XV_global(XV),yV(targetsV), StartStopV_global(StartStopV), infty(_infty), ell2(_ell2){
    J_global = X_global.Ncols();
    kappa = 30;
    errtol = 1e-2;
    stepsize = 2e-3;
    iterN = 50;
    modulo = 1000;
    update_groups = 2;
    max_groups = -1;
    LOGREG =  true;
    max_iter = 200000;
    N_tot = y.Length();
    K = StartStop_global.size();
    init_idx =  init_idx_default;
  };

 L_ONE_INFTY( const My_Matrix & X,  const My_Vector & targets,  const vector<vector <int> >& StartStop,  
	      const My_Matrix & XV,  const My_Vector & targetsV, const vector<vector <int> >& StartStopV,
	      bool _infty,  bool _ell2, int _max_groups)
   :X_global(X),y(targets),StartStop_global(StartStop),XV_global(XV),yV(targetsV),  StartStopV_global(StartStopV), infty(_infty), ell2(_ell2){
    J_global = X_global.Ncols();
    kappa = 30;
    errtol = 1e-2;
    stepsize = 2e-3;
    iterN = 50;
    modulo = 1000;
    update_groups = 2;
    max_groups = _max_groups;
    LOGREG =  true;
    max_iter = 200000;
    N_tot = y.Length();
    K = StartStop_global.size();
    init_idx =  init_idx_default;
  };

 L_ONE_INFTY( const My_Matrix & X,  const My_Vector & targets,  const vector<vector <int> >& StartStop,  
	      const My_Matrix & XV,  const My_Vector & targetsV, const vector<vector <int> >& StartStopV,
	      bool _infty,  bool _ell2, int _max_groups, double _Norm_p, double _Norm_gamma)
   :X_global(X),y(targets),StartStop_global(StartStop),XV_global(XV),yV(targetsV),  StartStopV_global(StartStopV), infty(_infty), ell2(_ell2), Norm_p(_Norm_p), Norm_gamma(_Norm_gamma){
    J_global = X_global.Ncols();
    kappa = 30;
    errtol = 1e-2;
    stepsize = 2e-3;
    iterN = 50;
    modulo = 1000;
    update_groups = 2;
    max_groups = _max_groups;
    LOGREG =  true;
    max_iter = 200000;
    N_tot = y.Length();
    K = StartStop_global.size();
    init_idx =  init_idx_default;
  };
  

  
  ~L_ONE_INFTY(){};   
  void SetParams(double kappa_new, double errtol_new,  double stepsize_new,  int iterN_new,  int modulo_new,  int update_groups_new, int m_iter, double _Norm_p = 1.0, double _Norm_gamma = 1.0){
    kappa =  kappa_new;
    errtol  = errtol_new;
    stepsize  = stepsize_new;
    iterN  = iterN_new;
    modulo  = modulo_new; 
    update_groups  = update_groups_new;
    max_iter = m_iter;
    Norm_p = _Norm_p;
    Norm_gamma = _Norm_gamma; 
  };

  void SetLOGREG(bool logreg){
    LOGREG = logreg;
  }
  vector<vector<int> > get_work_index_path(){return  work_index_path;};
  vector<vector<My_Vector> > get_beta_path(){return beta_path;}
  vector<double> get_kappa_path(){return kappa_path;}
  vector<double> get_lambda_path(){return lambda_path;}
  vector<vector<int> >get_non_unique_path(){return non_unique_path;};
  vector<double> get_balErr_path(){return balErr_path;};

  int optimize();
  void Rplot();
  double predict_with_best_beta(const My_Matrix & XT,  const My_Vector & targetsT,  const vector<vector <int> >& StartStopT);
  vector<double> predict_with_best_beta_perTask(const My_Matrix & XT,  const My_Vector & targetsT, const vector<vector <int> >& StartStopT);



 protected:
  int ensure_constraint_active(double& lambda,  // current lagrange para
			       vector<int>& work_index, // indices of groups in working set
			       const double& current_kappa,
			       vector<My_Vector>& beta, // regression weights
			       My_Vector& residualV,
			       double step
			       );


  int optimize_current( vector<int>& work_index,  
			const double& current_kappa, 
			vector<My_Vector>& beta,  
			double& lambda, 
			My_Vector& residualV);

 

  //for other exponential family models  adjust the following 3 functions...
  My_Vector Y_HAT(const vector<int>& work_index, const vector<My_Vector>& beta);

  My_Vector GRAD_J_RESIDUAL(const My_Vector& RESIDUAL,  const int& j );
    
  double Likelihood(const vector<int>& work_index, const vector<My_Vector>& beta, My_Vector&  y_minus_pi);



  
  double BSparse(vector<int>& work_index,  
		 const double& current_kappa,
		 vector<My_Vector>& beta, 
		 My_Vector&  y_minus_pi,
		 double s);

  double BSparse_L2(vector<int>& work_index,  
		       const double& current_kappa,
		       vector<My_Vector>& beta, 
		       My_Vector&  y_minus_pi,
		       double s);

  double BSparse_LP(vector<int>& work_index, 
		    const double& current_kappa,
		    vector<My_Vector>& beta, 
		    My_Vector&  y_minus_pi,
		    double s,double p,double gamma);

  vector<double> balancedERROR(const vector<int>& work_index, const vector<My_Vector>& beta);
  
  double squaredERROR(const vector<int>& work_index, const vector<My_Vector>& beta);

  My_Vector A_multiTask_times_x(const My_Matrix& X, const vector<vector<int> >& StartStop, const vector<int>& work_index, const vector<My_Vector>& beta);
  

    
private:
  const My_Matrix& X_global;
  const My_Vector& y;
  const vector<vector <int> >& StartStop_global;
  const My_Matrix& XV_global;
  const My_Vector& yV;
  const vector<vector <int> >& StartStopV_global;

  int J_global;

  int K;
  int N_tot;

  double kappa;
  double errtol;
  double stepsize;
  int iterN;
  int modulo; 
  int update_groups;
  
  vector< vector<int> >work_index_path;
  vector < vector<My_Vector> > beta_path;
  vector<double> kappa_path;
  vector<double> lambda_path;
  vector< vector<int> >non_unique_path;
  vector<double> balErr_path;
  vector<double> bestThresh_path;
  bool  lambda_is_active;

  bool infty;
  int max_groups;
  bool LOGREG;
  int  max_iter;
  double  Norm_p, Norm_gamma;
  bool ell2;

  vector<int> init_idx;

};
#endif
