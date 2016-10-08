#ifndef MY_SUBVEC1_H
#define MY_SUBVEC1_H

#include <iostream>

using namespace std;

class My_Vector;

class My_SubVector{

 
public:
  My_SubVector(){};
  My_SubVector(const My_Vector& M);      
  ~My_Vector(){};                                     
 
  double & operator[](int)  { return ptr[index]; };
  const double & operator[](int index) const { return ptr[index]; } ;
  const My_Vector &operator=(const My_Vector &);
  const My_Vector &operator=(const double &);
  
private:
   
  double *ptr; 
  int size;
  int start;
  int stop;
};

#endif



