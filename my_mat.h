#ifndef MY_MAT1_H
#define MY_MAT1_H

#include <iostream>
#include "my_mat_bool.h"
#include "my_mat_bool_scale.h"
#include "my_vec.h"
#include <vector>

using namespace std;

class My_Matrix_bool;
class My_Matrix_bool_scale;
class My_Vector;

class My_Matrix {

  friend class My_Matrix_bool;
  friend class My_Matrix_bool_scale;
  friend class My_Matrix_bool_scale_sparse;
  friend class My_Vector;

  friend ostream &operator<<(ostream &, const My_Matrix &); 
  friend istream &operator>>(istream &, My_Matrix &);
  
  friend  My_Matrix Transpose(const My_Matrix& a);
  friend  My_Vector A_shifted_times_x(const My_Matrix & A, const My_Vector & x); /* computes Ax on a "shifted" matrix A (multi-task learning)*/

 friend  My_Vector A_shifted_T_times_x(const My_Matrix & A,  const My_Vector & x); // computes A^t x on a "shifted" matrix A (multi-task learning)


public:
  My_Matrix();
  My_Matrix(const My_Matrix& M);
  My_Matrix(int, int);                     
  ~My_Matrix();                                     
  double el(int r, int c) const { return ptr[r][c]; };
  double & el(int r, int c) { return ptr[r][c]; };
  
  void ReDimension(int, int);
  void Enlarge(int, int);
  int getSizeR() const;
  int getSizeC() const;
  int Nrows() const;
  int Ncols() const;
  
  double Variance() const;
  double Variance_nonzero() const;
  double Sum() const;
  double Mean_nonzero() const;

  My_Matrix Rows0(int a, int b);
  My_Matrix Columns0(int a, int b);
  My_Matrix SubMatrix0(int ra, int rb,int ca, int cb);

  void setRows0(int a, int b,const My_Matrix & M);
  void setRows0(int a, int b,const My_Matrix_bool_scale & M);
  void setColumns0(int a, int b,const My_Matrix & M);
  void setSubMatrix0(int ra, int rb,int ca, int cb, const My_Matrix & M);
  
  void setRows0(int a, int b,const  double & d);
  void setColumns0(int a, int b,const double & d);
  void setSubMatrix0(int ra, int rb,int ca, int cb, const double & d);

  void setColumn0(int a,const double * d);
  My_Matrix ExtractColumns(vector<int>& which);

  void standardize();

  My_Matrix operator - ( const double & d) ;
  My_Matrix operator - ( My_Matrix & M);
  My_Matrix operator - (const My_Matrix & M) const;
  My_Matrix operator + ( const double & d) ;
  My_Matrix operator + ( My_Matrix & M);
  My_Matrix operator + (const My_Matrix & M) const;
  My_Matrix operator * ( My_Matrix & M);
  My_Matrix operator * (const My_Matrix & M) const;
  My_Matrix operator * (const My_Matrix_bool & M);
  My_Matrix operator * (const My_Matrix_bool_scale & M);
  My_Vector operator * (const My_Vector & V);
  My_Matrix operator * ( const double & d);
  My_Matrix operator % (const My_Matrix & M) const;

  double *operator[](int);
  const double *operator[](int) const;
  const My_Matrix &operator=(const My_Matrix &);
  const My_Matrix &operator=(const double &);

  void RandUniformFill();
  void RandGaussFill();

private:
   
   double **ptr; 
   int sizeR;
   int sizeC;
};

#endif



