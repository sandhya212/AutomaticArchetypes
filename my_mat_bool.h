#ifndef MY_MAT_B_H
#define MY_MAT_B_H

#include <iostream>

#include "my_mat.h"
#include "my_vec.h"

using namespace std;

class My_Matrix;
class My_Vector;
class My_Matrix_bool {

  friend class My_Matrix;
  friend class My_Vector;

  friend ostream &operator<<(ostream &, const My_Matrix_bool &); 
 
  friend  My_Matrix_bool Transpose(const My_Matrix_bool& a);

public:

  My_Matrix_bool();
  My_Matrix_bool(const My_Matrix_bool& M);
  My_Matrix_bool(int,int);                     
  ~My_Matrix_bool();                                     
  bool el(int r, int c) const { return ptr[r][c]; };
  // bool & el(int r, int c) { return ptr[r][c]; };
  vector<bool>::reference el(int r, int c) { return ptr[r][c]; };
  void ReDimension(int, int);
  void Enlarge(int, int);
  int getSizeR() const;
  int getSizeC() const;
  int Nrows() const;
  int Ncols() const;

  
  
  My_Matrix_bool Rows0(int a, int b);
  My_Matrix_bool Columns0(int a, int b);
  My_Matrix_bool SubMatrix0(int ra, int rb,int ca, int cb);

  void setRows0(int a, int b,const My_Matrix_bool & M);
  void setColumns0(int a, int b,const My_Matrix_bool & M);
  void setSubMatrix0(int ra, int rb,int ca, int cb, const My_Matrix_bool & M);
  
  void setRows0(int a, int b,const  double & d);
  void setColumns0(int a, int b,const double & d);
  void setSubMatrix0(int ra, int rb,int ca, int cb, const double & d);

  My_Matrix_bool cbind( My_Matrix_bool &);

  My_Matrix operator * ( My_Matrix_bool & M);
  My_Matrix operator * (const My_Matrix_bool & M) const;
  My_Matrix operator * (const My_Matrix & M) const;
  My_Vector operator * (const My_Vector & V);



  const vector<bool>& operator[](int i) const{return ptr[i];};
  vector<bool>& operator[](int i){return ptr[i];};
 

  const My_Matrix_bool &operator=(const My_Matrix_bool &);

  My_Matrix_bool operator | (const My_Matrix_bool & M);
  
private:
   
   vector< vector< bool> > ptr; 
   int sizeC;
};

#endif



