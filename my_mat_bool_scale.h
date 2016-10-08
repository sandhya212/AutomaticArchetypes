#ifndef MY_MAT_BS_H
#define MY_MAT_BS_H

#include <iostream>

#include "my_mat.h"
#include "my_vec.h"

using namespace std;

class My_Matrix;
class My_Vector;

class My_Matrix_bool_scale {

  friend class My_Matrix;
  friend class My_Vector;

  friend ostream &operator<<(ostream &, const My_Matrix_bool_scale &); 
 
  friend  My_Matrix_bool_scale Transpose(const My_Matrix_bool_scale& a);


  friend My_Matrix_bool_scale multi_interact_5(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3, const My_Matrix_bool_scale & M_4, const My_Matrix_bool_scale & M_5);

  friend My_Matrix_bool_scale multi_interact_4(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3, const My_Matrix_bool_scale & M_4);

  friend My_Matrix_bool_scale multi_interact_3(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3);
  
  friend My_Matrix_bool_scale multi_interact_2(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2);


  friend double multi_interact_5_multTvecNorm(My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2, My_Matrix_bool_scale & M_3,My_Matrix_bool_scale & M_4,My_Matrix_bool_scale & M_5, const My_Vector& V,  const double& exp);

  friend double multi_interact_4_multTvecNorm(My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3,  My_Matrix_bool_scale & M_4, const My_Vector& V,  const double& exp);

 friend double  multi_interact_3_multTvecNorm( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2, My_Matrix_bool_scale & M_3,const My_Vector& V,  const double& exp);

 friend double  multi_interact_2_multTvecNorm( My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2, const My_Vector& V,  const double& exp);
 
friend double multi_interact_5_multTvecNorm_L1(My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2, My_Matrix_bool_scale & M_3,My_Matrix_bool_scale & M_4,My_Matrix_bool_scale & M_5, const My_Vector& V,  const double& exp);

  friend double multi_interact_4_multTvecNorm_L1(My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3,  My_Matrix_bool_scale & M_4, const My_Vector& V,  const double& exp);

 friend double  multi_interact_3_multTvecNorm_L1( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2, My_Matrix_bool_scale & M_3,const My_Vector& V,  const double& exp);

 friend double  multi_interact_2_multTvecNorm_L1( My_Matrix_bool_scale & M_1, My_Matrix_bool_scale & M_2, const My_Vector& V,  const double& exp);

 friend  My_Vector multi_interact_5multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3,  const My_Matrix_bool_scale & M_4, const My_Matrix_bool_scale & M_5, const My_Vector& V, const double& exp);

 friend  My_Vector multi_interact_4multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3,  const My_Matrix_bool_scale & M_4, const My_Vector& V, const double& exp);

 friend  My_Vector multi_interact_3multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3,  const My_Vector& V, const double& exp);

 friend My_Vector  multi_interact_2multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Vector& V, const double& exp);
public:

  My_Matrix_bool_scale();
  My_Matrix_bool_scale(const My_Matrix_bool_scale& M);
  My_Matrix_bool_scale(int,int);                    
  My_Matrix_bool_scale(int mSizeR, int mSizeC, double mScale); 
  ~My_Matrix_bool_scale();                                     
  bool el(int r, int c) const { return ptr[r][c]; };
  // bool & el(int r, int c) { return ptr[r][c]; };
  vector<bool>::reference el(int r, int c) { return ptr[r][c]; };
  void ReDimension(int, int);
  void Enlarge(int, int);
  int getSizeR() const;
  int getSizeC() const;
  int Nrows() const;
  int Ncols() const;
  int NNZcols() const;

  void SetScale(double s){scale = s;};
  double GetScale() const {return scale;};

  My_Vector multMtV(const My_Vector & V);
  double multMtVNorm(const My_Vector & V);
  double multMtVNorm_L1(const My_Vector & V);
  My_Vector multMtVScale(const My_Vector & V, const double& dscale);
  My_Vector WplusMultMV(const My_Vector & V,const My_Vector & W ); 

  
  My_Matrix_bool_scale Rows0(int a, int b);
  My_Matrix_bool_scale Columns0(int a, int b);
  My_Matrix_bool_scale SubMatrix0(int ra, int rb,int ca, int cb);

  void setRows0(int a, int b,const My_Matrix_bool_scale & M);
  void setColumns0(int a, int b,const My_Matrix_bool_scale & M);
  void setSubMatrix0(int ra, int rb,int ca, int cb, const My_Matrix_bool_scale & M);
  
  void setRows0(int a, int b,const  double & d);
  void setColumns0(int a, int b,const double & d);
  void setSubMatrix0(int ra, int rb,int ca, int cb, const double & d);

  My_Matrix_bool_scale cbind( My_Matrix_bool_scale &);

  My_Matrix operator * ( My_Matrix_bool_scale & M);
  My_Matrix operator * (const My_Matrix_bool_scale & M) const;
  My_Matrix operator * (const My_Matrix & M) const;
  My_Vector operator * (const My_Vector & V);



  const vector<bool>& operator[](int i) const{return ptr[i];};
  vector<bool>& operator[](int i){return ptr[i];};
 

  const My_Matrix_bool_scale &operator=(const My_Matrix_bool_scale &);

  My_Matrix_bool_scale operator | (const My_Matrix_bool_scale & M);
  
  void update_RowIndices();


private:
   
  vector< vector< bool> > ptr; 
  int sizeC;
  double scale;
  bool is_up_to_date;
  vector< vector <int> > rowInd;
};

#endif



