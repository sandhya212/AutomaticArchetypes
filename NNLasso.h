//#include "my_mat.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>


// #include "scythestat/rng.h"
// #include "scythestat/rng/mersenne.h"
// #include "scythestat/distributions.h"
// #include "scythestat/ide.h"
// #include "scythestat/la.h"
#include "scythestat/matrix.h"
#include "scythestat/smath.h"
#include "scythestat/stat.h"
#include "scythestat/optimize.h"
// #include "scythestat/lapack.h"

//class My_Vector;
//class My_Matrix;
using namespace scythe;
using namespace std;

//vector<vector<double> > NNLasso(vector<double>& y, vector<vector<double> >& X, double& kappa, double& epsilon,vector<vector<double> > XTran);
//  My_Vector A_multiTask_times_x(const My_Matrix& X, const vector<vector<int> >& StartStop, const vector<int>& work_index, const vector<My_Vector>& beta);
// y,X,kappa, epsilon, XtX =  t(X) %*% X
//My_Vector NNLassolarge(My_Vector& y, My_Matrix& X, double& kappa, double& epsilon);

int NNLasso_small(vector<double>& y, Matrix<double> & X, double kappa, double epsilon, Matrix<double> & XTranX,Matrix<double> &A );

// int local_Analysis(int min_overlap, int K, int  nSample, int nT, int max_reads_in_window, int WindowStart, 
// int WindowIncrement, vector<vector<int> >& WindowStartStop,  int max_seq_length,int min_seq_start, const vector<int>& Positions_Start, 
//  const vector<vector<int> >& Reads,  bool have_quality_scores, const vector<string> & quality_scores, const vector<string>& IDs,  string prefix, string  SQ_string,  
// bool have_true_haplotypes,   const vector<vector<int> >& trueHaplos,    const  vector<int>& trueHaplo_Positions_Start,   const   vector<string>& trueHaplo_IDs,  
// vector<int>& foundClusters,  double entropy_threshold, double entropy_fraction,  const vector<int>& entropy_select){
// 
// 
//  local_Analysis( min_overlap, K, nSample, nT, max_reads_in_window, WindowStart,  local_window_size, WindowStartStop,  
//  max_sequence_length, min_seq_start, Positions_Start, Reads, have_quality_scores, quality_scores, IDs,  
// prefix,SQ_string, have_true_haplotypes, trueHaplos, trueHaplo_Positions_Start, trueHaplo_IDs, foundClusters,  entropy_threshold, entropy_fraction, entropy_select);