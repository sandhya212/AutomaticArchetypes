
#include "my_vec.h"
#include <cmath>

int  compareDecrease(const void * a, const void * b)
{
   double* arg1 = (double*) a;
   double* arg2 = (double*) b;

   if( *arg1 < *arg2 ) return 1;
     else if( *arg1 == *arg2 ) return 0;
     else return -1;

}


My_Vector::My_Vector(){
  int i;
  size = 1;
 
  ptr = new double[size];

 
  for (i = 0; i < size; i++)
      ptr[i] = 0.;            // initialize array
}

My_Vector::My_Vector(const My_Vector& M){

   int i;
   size = M.size;
  
   ptr = new double[size];
  
 
   for (i = 0; i < size; i++)
     ptr[i] = M.ptr[i];      
}


My_Vector::My_Vector(int mSize)
{
   int i; 
   size = mSize;
  
   ptr = new double[size];
 
 
   for (i = 0; i < size; i++)
     ptr[i] = 0.;            // initialize array
}

My_Vector::~My_Vector()
{  
  delete [] ptr;
}

void My_Vector::Enlarge(int r){

    int i,j,k;
    double* ptr_buff;
    ptr_buff = new double[size];
   

    for (i = 0; i < size; i++)
      ptr_buff[i] = ptr[i];
    
    delete [] ptr;
    
    
    ptr = new double[r];
  
    
    int min_r;
    if(size > r)
	min_r = r;
    else
	min_r = size;
  

    for (i = 0; i < min_r; i++)
      ptr[i] = ptr_buff[i];  
    
    for (i = min_r; i < r; i++)
      ptr[i] = 0.0;
    
   
    delete [] ptr_buff;
    
    size = r;
   
    
}
void My_Vector::ReDimension(int r){
    int i,j;
       
    if(r != size){

      delete [] ptr; 
      size = r;          
      ptr = new double[size];
    }     
  
}



void My_Vector::VfromPTR(double* c, const int& n, bool negate = false){
  if(size != n)
    cout << "wrong size in void VfromPTR(double* c, const int& n, const double& scale)" << endl;
  if(negate){
    for (int i = 0; i < size; i++){
      ptr[i] = -c[i];
    }
  }
  else{
    for (int i = 0; i < size; i++){
      ptr[i] = c[i];
    }
  }
  
  
}

void My_Vector::setSubVector0(int a, int b, const My_Vector & SV){
   int i;
   if(SV.size != b-a+1)
     cout <<"Wrong Size in void setSubVector0(int a, int b, My_Vector & SV)" << endl;
   

   for(i = a; i <= b; i++)
     ptr[i] = SV.ptr[i-a];

}
void My_Vector::setSubVector0_fromOtherSubvector_Scale(int a, int b, const My_Vector & SV, const double& scale){
   
  int i;  

   for(i = a; i <= b; i++)
     ptr[i] = SV.ptr[i]*scale;

}



void My_Vector::setSubVector0_fromOtherSubvector_ABS_SORT(int a, int b, const My_Vector & SV){

  int group_size = b-a+1;

  double* abs_v = new double[group_size];
  for(int i = a; i <= b; i++){
    abs_v[i-a] = fabs(SV.ptr[i]);
    //  cout <<  abs_v[i-a]<< ' ';
  }
  //cout << endl;

  //cout << " -------------" << endl;
  qsort(abs_v,group_size, sizeof(double), compareDecrease);
  for(int i = a; i <= b; i++){
    ptr[i] = abs_v[i-a];
    // cout << ptr[i]<< ' ';
  }
  // cout << endl;
  delete [] abs_v;
 
}

void My_Vector::setSubVector0(int a, int b, const double & d){
   int i;
   

   for(i = a; i <= b; i++)
     ptr[i] = d;

}
My_Vector My_Vector::SubVector0(int a, int b){
    int i;
    int dim = b-a+1;
  

    My_Vector sub(dim); 

    for(i = a; i <= b; i++)
      sub[i-a]= ptr[i];

    return sub;
   
}

double My_Vector::Norm() const 
{
  double sum =0.;
  for (int k = 0; k < size; k++)
    sum += ptr[k]* ptr[k];
  return sqrt(sum);
}

double My_Vector::Norm_infty() const 
{
  double sum = -1.0;
  for (int k = 0; k < size; k++){
    if(sum < fabs( ptr[k]))
      sum =  fabs( ptr[k]);
  }
  return sum;
}

double My_Vector::Norm_L1() const 
{
  double sum =0.;
  for (int k = 0; k < size; k++)
    sum += fabs(ptr[k]);
  return sum;
}

double My_Vector::pNorm( const double& p ) const 
{
  double sum =0.;
  for (int k = 0; k < size; k++)
    sum += pow(fabs(ptr[k]),p);
  //  cout << "p=" <<p <<endl;
  return pow(sum,(1/p));
  
}

 double My_Vector::gammaNorm( const double& gamma) const
{
  double sum =0.;
// cout << "gamma" <<gamma<<"---" <<endl; exit (1);
  for (int k = 0; k < size; k++)
     sum += pow(fabs(ptr[k]),gamma);
  return pow(sum,1/gamma);
 
}




/*
double My_Vector::operator VplusMultMV (const My_Matrix_bool_scale & M, const My_Vector & V) const { 

  int i, j;
  double sum = 0;
  int sizeR = ptr.size();
  int len = V.Length() ;
  if(sizeC!= len)
    cout << "Wrong Size in my_vector = multMtVScale(my_vector, dscale)" << endl;
 
  My_Vector prod(sizeR);
 
   
  for(i = 0; i < sizeR; i++){
    sum = 0.0;
    for(j = 0; j < len; j++){ 
      if(ptr[i][j])
	sum += V[j];
    }
    prod[i] = scale*sum*dscale; 
  }
  return prod; 
}
*/ 

// Overloaded subscript operator
double& My_Vector::operator[](int subscript)
{
   return ptr[subscript];    
}



My_Vector My_Vector::operator * (const double & d){
  
  My_Vector V(size);
  for (int k = 0; k < size; k++)
    V[k]  = d* ptr[k];
  return V;
}
My_Vector My_Vector::operator * (const int & d){
  
  My_Vector V(size);
  for (int k = 0; k < size; k++)
    V[k]  = d* ptr[k];
  return V;
}

My_Vector My_Vector::operator / (const double & d){
  
  My_Vector V(size);
  for (int k = 0; k < size; k++)
    V[k]  = ptr[k]/d;
  return V;
}

My_Vector My_Vector::operator + (const double & d){
  
  My_Vector V(size);
  for (int k = 0; k < size; k++)
    V[k]  = d + ptr[k];
  return V;
}
My_Vector My_Vector::operator - (const double & d){
  
  My_Vector V(size);
  for (int k = 0; k < size; k++)
    V[k]  = ptr[k]-d;
  return V;
}

My_Vector My_Vector::operator+(My_Vector & V){
  
  if(size != V.size)
    cout <<"Wrong Size in My_Vector My_Vector::operator + (const  My_Vector & V)" << endl;
  My_Vector W(size);
  for (int k = 0; k < size; k++)
    W[k]  = ptr[k]+V.ptr[k];
  return W;
}

My_Vector My_Vector::operator-(My_Vector & V){
  
  if(size != V.size)
    cout <<"Wrong Size in My_Vector My_Vector::operator - (const  My_Vector & V)" << endl;
  My_Vector W(size);
  for (int k = 0; k < size; k++)
    W[k]  = ptr[k]-V.ptr[k];
  return W;
}

My_Vector My_Vector::operator+(const My_Vector & V) const {
  
  if(size != V.size)
    cout <<"Wrong Size in My_Vector My_Vector::operator + (const  My_Vector & V)" << endl;
  My_Vector W(size);
  for (int k = 0; k < size; k++)
    W[k]  = ptr[k]+V.ptr[k];
  return W;
}

My_Vector My_Vector::operator-(const My_Vector & V) const {
  
  if(size != V.size)
    cout <<"Wrong Size in My_Vector My_Vector::operator - (const  My_Vector & V)" << endl;
  My_Vector W(size);
  for (int k = 0; k < size; k++)
    W[k]  = ptr[k]-V.ptr[k];
  return W;
}
 

My_Vector My_Vector::operator - () const{
   
  My_Vector V(size);
  for (int k = 0; k < size; k++)
    V[k]  = -ptr[k];
  return V;
}
const My_Vector &My_Vector::operator=(const My_Vector &right)
{
   int i,j,k,l,m;

   if (&right != this) {    // check for self-assignment
 
 
          
     if(size != right.size){
       delete [] ptr; 
       
       size = right.size;    // resize this object
     
       
       ptr = new double[size];
     }     
     
     for (k = 0; k < size; k++){
       
       ptr[k] = right.ptr[k];  // copy My_Vector into object
     }
   }
   
   return *this;   // enables x = y = z;
}

const My_Vector &My_Vector::operator=(const double &right)
{
   int k;  
     
   for (k = 0; k < size; k++) 
	 ptr[k] = right; 
   
   return *this;  
}

void My_Vector::RandGaussFill(){
 double x1, x2, w, y1, y2;
 
 for(int i = 0; i < size; i++){
     do {
       x1 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
       x2 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
       w = x1 * x1 + x2 * x2;
     } while ( w >= 1.0 );
     
     w = sqrt( (-2.0 * log( w ) ) / w );
     y1 = x1 * w;
     // y2 = x2 * w;
     ptr[i] = y1;
   }
}
 

void My_Vector::Add_A_shifted_times_x(const My_Matrix & A,  const My_Vector & x){ // adds A x ("shifted" matrix A);

  int K = A.sizeC;
  int N = A.sizeR;
  int idx =0; 
  
  for(int k =0; k<K;k++){
    for(int n =0; n< N; n++){
      ptr[idx] += A.ptr[n][k]*x.ptr[k];
      idx++;
    }
  }
 
}
ostream &operator<<(ostream &output, const My_Vector &m)
{
 int i;
 int j;
   for (i = 0; i < m.size; i++) 
	 output << m.ptr[i] << ' ';
   output << endl;


   return output;   // enables cout << x << y;

}

istream &operator>>(istream &input, My_Vector &m)
{
 int i;
 int j;
   for (i = 0; i < m.size; i++) 
	 input >> m.ptr[i];
   
   
   return input;   // enables cout << x << y;
   
}
 

// double = my_vector * my_vector
double My_Vector::operator * ( My_Vector & V){ 

  int i;
  double sum = 0;
  if(size!= V.size )
    cout << "Wrong Size in my_vector * my_vector" << endl;
 
 
  
   
  for(i = 0; i < size; i++)
    sum += ptr[i] * V.ptr[i];
	      

  return sum; 
}

double My_Vector::operator * (const My_Vector & V) const { 

  int i;
  double sum = 0;
  if(size!= V.size )
    cout << "Wrong Size in my_vector * my_vector" << endl;
 
 
  
   
  for(i = 0; i < size; i++)
    sum += ptr[i] * V.ptr[i];
	      

  return sum; 
}

// double = my_vector * STL_vector
double My_Vector::operator * (const vector<double> & V){ 

  int i;
  double sum = 0;
  if(size!= V.size() )
    cout << "Wrong Size in  my_vector * STL_vector" << endl;
 
 
  
   
  for(i = 0; i < size; i++)
    sum += ptr[i] * V[i];
	      

  return sum; 
}




My_Vector operator*(const double& d, const My_Vector & V){
  My_Vector W(V.size);
  for (int k = 0; k < V.size; k++)
    W[k]  = d* V.ptr[k];
  return W;
}
My_Vector operator*(const int& d, const My_Vector & V){
  My_Vector W(V.size);
  for (int k = 0; k < V.size; k++)
    W[k]  = d* V.ptr[k];
  return W;
}
 double Norm(const My_Vector & V, const int & a, const int &b){ 

   double sum =0.;
   for (int k = a; k <= b; k++)
     sum += V[k]* V[k];
   return sqrt(sum);
 }


double Norm_infty(const My_Vector & V, const int & a, const int &b){ 

  double sum = -1.0;
  for (int k = a; k <= b; k++){
    if(sum < fabs( V[k]))
      sum =  fabs( V[k]);
  }
  return sum;
}


double Norm_L1(const My_Vector & V, const int & a, const int &b){ 

  double sum = 0.0;
  for (int k = a; k <= b; k++)
    sum += fabs(V[k]);
  return sum;
}


double pNorm(const My_Vector & V, const int & a, const int &b,const double& p){ 
    double sum =0.;
   for (int k = a; k <= b; k++)
     sum += pow(fabs(V[k]),p);
     return pow(sum,(1/p));
 }



 double gammaNorm(const My_Vector & V, const int & a, const int &b,const double& gamma){ 
 double sum =0.;
   for (int k = a; k <= b; k++)
     sum += pow(fabs(V[k]),gamma);
   return pow(sum,(1/gamma));
 }





My_Vector  ABS_SORT(const My_Vector & SV){

  int s = SV.size;
  My_Vector res(s);
  for(int i = 0; i < s; i++){
    res[i] =  fabs(SV.ptr[i]);
  }
  qsort(res.ptr,s, sizeof(double), compareDecrease);
  return res;
}
