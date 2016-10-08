#ifndef MY_MAT_BS_S_H
#define MY_MAT_BS_S_H

#include <iostream>

#include "my_mat.h"
#include "my_vec.h"

using namespace std;

class My_Matrix;
class My_Vector;

class My_Matrix_bool_scale_sparse {

  friend class My_Matrix;
  friend class My_Vector;

  friend ostream &operator<<(ostream &, const My_Matrix_bool_scale_sparse &); 
 
  friend  My_Matrix_bool_scale_sparse Transpose(const My_Matrix_bool_scale_sparse& a);

public:

  My_Matrix_bool_scale_sparse();
  My_Matrix_bool_scale_sparse(const My_Matrix_bool_scale_sparse& M);
  My_Matrix_bool_scale_sparse(int,int);                    
  My_Matrix_bool_scale_sparse(int mSizeR, int mSizeC, double mScale); 
  ~My_Matrix_bool_scale_sparse();                                     
  bool el(int r, int c) const { return ptr[r][c]; };
  // bool & el(int r, int c) { return ptr[r][c]; };
  vector<bool>::reference el(int r, int c) { return ptr[r][c]; };
  void ReDimension(int, int);
  void Enlarge(int, int);
  int getSizeR() const;
  int getSizeC() const;
  int Nrows() const;
  int Ncols() const;

  void SetScale(double s){scale = s;};
  double GetScale() const {return scale;};
  My_Vector multMtVScale(const My_Vector & V, const double& dscale);
  My_Vector WplusMultMV(const My_Vector & V,const My_Vector & W ); 

  
  My_Matrix_bool_scale_sparse Rows0(int a, int b);
  My_Matrix_bool_scale_sparse Columns0(int a, int b);
  My_Matrix_bool_scale_sparse SubMatrix0(int ra, int rb,int ca, int cb);

  void setRows0(int a, int b,const My_Matrix_bool_scale_sparse & M);
  void setColumns0(int a, int b,const My_Matrix_bool_scale_sparse & M);
  void setSubMatrix0(int ra, int rb,int ca, int cb, const My_Matrix_bool_scale_sparse & M);
  
  void setRows0(int a, int b,const  double & d);
  void setColumns0(int a, int b,const double & d);
  void setSubMatrix0(int ra, int rb,int ca, int cb, const double & d);

  My_Matrix_bool_scale_sparse cbind( My_Matrix_bool_scale_sparse &);

  My_Matrix operator * ( My_Matrix_bool_scale_sparse & M);
  My_Matrix operator * (const My_Matrix_bool_scale_sparse & M) const;
  My_Matrix operator * (const My_Matrix & M) const;
  My_Vector operator * (const My_Vector & V);



  const vector<bool>& operator[](int i) const{return ptr[i];};
  vector<bool>& operator[](int i){return ptr[i];};
 

  const My_Matrix_bool_scale_sparse &operator=(const My_Matrix_bool_scale_sparse &);

  My_Matrix_bool_scale_sparse operator | (const My_Matrix_bool_scale_sparse & M);
  
private:
   
  vector< vector< bool> > ptr; 
  int sizeC;
  double scale;
  vector< vector<int> > Index;
};

#endif



