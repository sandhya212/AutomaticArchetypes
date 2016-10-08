#ifndef MY_VEC1_H
#define MY_VEC1_H

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "my_mat.h"
using namespace std;

class  My_Matrix_bool_scale;
class My_Matrix;
class My_Vector {

  friend ostream& operator<<(ostream &, const My_Vector &); 
  friend istream& operator>>(istream &, My_Vector &);

  friend double Norm(const My_Vector& V){return V.Norm();};

  friend double Norm(const My_Vector& V, const int& a, const int& b);


  friend double Norm_infty(const My_Vector& V){return V.Norm_infty();};

  friend double  Norm_infty(const My_Vector& V, const int& a, const int& b);

  friend double Norm_L1(const My_Vector& V){return V.Norm_L1();};

  friend double  Norm_L1(const My_Vector& V, const int& a, const int& b);

  friend double pNorm(const My_Vector& V, const double& p){return V.pNorm(p);};

  friend double pNorm(const My_Vector& V, const int& a, const int& b, const double& p);
  
  friend double gammaNorm(const My_Vector& V, const double& gamma){return V.gammaNorm(gamma);};
  
  friend double gammaNorm(const My_Vector& V, const int& a, const int& b, const double& gamma);



  friend My_Vector operator*(const double&, const My_Vector &);
  friend My_Vector operator*(const int&, const My_Vector &);

  friend class My_Matrix_bool_scale_sparse;
  friend class My_Matrix_bool_scale;
  friend class My_Matrix;

 friend double multi_interact_5_multTvecNorm( My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3, My_Matrix_bool_scale & M_4,My_Matrix_bool_scale & M_5, const My_Vector& V,  const double& exp);

 friend double multi_interact_4_multTvecNorm(My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2, My_Matrix_bool_scale & M_3,  My_Matrix_bool_scale & M_4, const My_Vector& V,  const double& exp);

 friend double  multi_interact_3_multTvecNorm(My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2,My_Matrix_bool_scale & M_3,const My_Vector& V,  const double& exp);

 friend double  multi_interact_2_multTvecNorm( My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2, const My_Vector& V,  const double& exp);


friend double multi_interact_5_multTvecNorm_L1(My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2, My_Matrix_bool_scale & M_3,My_Matrix_bool_scale & M_4,My_Matrix_bool_scale & M_5, const My_Vector& V,  const double& exp);

  friend double multi_interact_4_multTvecNorm_L1(My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3,  My_Matrix_bool_scale & M_4, const My_Vector& V,  const double& exp);

 friend double  multi_interact_3_multTvecNorm_L1( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2, My_Matrix_bool_scale & M_3,const My_Vector& V,  const double& exp);

 friend double  multi_interact_2_multTvecNorm_L1( My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2, const My_Vector& V,  const double& exp);


 friend   My_Vector multi_interact_5multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3,  const My_Matrix_bool_scale & M_4, const My_Matrix_bool_scale & M_5, const My_Vector& V, const double& exp);

 friend   My_Vector multi_interact_4multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3,  const My_Matrix_bool_scale & M_4, const My_Vector& V, const double& exp);

 friend   My_Vector multi_interact_3multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3,  const My_Vector& V, const double& exp);

 friend   My_Vector multi_interact_2multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Vector& V, const double& exp);

 friend  My_Vector A_shifted_times_x(const My_Matrix & A, const My_Vector & x); /* computes Ax on a "shifted" matrix A (multi-task learning)*/

 friend  My_Vector A_shifted_T_times_x(const My_Matrix & A,  const My_Vector & x); // computes A^t x on a "shifted" matrix A (multi-task learning)

 friend  My_Vector  ABS_SORT(const My_Vector & SV);
public:
  My_Vector();
  My_Vector(const My_Vector& M);
  My_Vector(int);                     
  ~My_Vector();                                     
  double el(int r) const { return ptr[r]; };
  double & el(int r) { return ptr[r]; };
  
  void ReDimension(int);
  void Enlarge(int);
  int getLength() const {return size;};
  int Length() const {return size;};
  
  //int Max() const {return max;}; //added by sandhya - 10th Feb 12
//  int Max(int);

  double* getPTR() { return ptr; };
  void VfromPTR(double* c, const int& n, bool negate);
 
 
  double Norm() const;
  double Norm_L1() const;
  double Norm_infty() const;
  double pNorm(const double& p) const;
  double gammaNorm(const double& gamma) const;
 
  My_Vector SubVector0(int a, int b);
  void setSubVector0(int a, int b, const My_Vector & SV);
  void setSubVector0(int a, int b, const double & d);
  void setSubVector0_fromOtherSubvector_Scale(int a, int b, const My_Vector & SV, const double& scale);

  void setSubVector0_fromOtherSubvector_ABS_SORT(int a, int b, const My_Vector & SV);
 
   


  double operator * (My_Vector & V);
  double operator * (const My_Vector & V)  const;
  double operator * (const vector<double> & V);
  My_Vector operator * (const double & d);
  My_Vector operator * (const int & d);
  My_Vector operator / (const double & d);
  My_Vector operator + (const double & d);
  My_Vector operator - (const double & d);

  My_Vector operator+(My_Vector & V);
  My_Vector operator+( const My_Vector & V) const;
  My_Vector operator-(My_Vector & V);
  My_Vector operator-( const My_Vector & V) const;
  

  My_Vector operator - () const;

  double & operator[](int);
  const double & operator[](int index) const { return ptr[index]; } ;
  const My_Vector &operator=(const My_Vector &);
  const My_Vector &operator=(const double &);

  void Add_A_shifted_times_x(const My_Matrix & A,  const My_Vector & x); // adds A^t x ("shifted" matrix A);

  void RandGaussFill();

  private:
  double *ptr; 
  int size;
 
 
};

#endif



