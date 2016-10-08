
#include "my_mat.h"
#include <cmath>


My_Matrix::My_Matrix(){
  int i; int j; int k;
  sizeR = 1;
  sizeC = 1;
  
  ptr = new double*[sizeR];
 
  for ( k = 0; k < sizeR; k++){  
    ptr[k] = new double[sizeC];
   
  }
  for (i = 0; i < sizeR; i++)
    for (j = 0; j < sizeC; j++)
      ptr[i][j] = 0.;            // initialize array
}

My_Matrix::My_Matrix(const My_Matrix& M){

   int i; int j; int k;
   sizeR = M.sizeR;
   sizeC = M.sizeC;
   
   ptr = new double*[sizeR];
 
   for ( k = 0; k < sizeR; k++){  
       ptr[k] = new double[sizeC];
     
   }
   for (i = 0; i < sizeR; i++)
     for (j = 0; j < sizeC; j++)
	 ptr[i][j] = M.ptr[i][j];      
}
My_Matrix::My_Matrix(int mSizeR, int mSizeC)
{
   int i; int j; int k;
   sizeR = mSizeR;
   sizeC = mSizeC;
   
   ptr = new double*[sizeR];
 
   for ( k = 0; k < sizeR; k++){  
       ptr[k] = new double[sizeC];
      
   }
   for (i = 0; i < sizeR; i++)
     for (j = 0; j < sizeC; j++)
	 ptr[i][j] = 0.;            // initialize array
}

My_Matrix::~My_Matrix()
{  
  int i;
  for( i = 0; i < sizeR; i++)
    delete [] ptr[i];
  delete [] ptr;
}


void My_Matrix::standardize(){
  double mean, var;
  for(int j =0; j<sizeC; j++){
    mean =0.0;
    for(int i =0; i<sizeR; i++)
      mean += ptr[i][j];
    mean /= sizeR;
    var =0.0;
    for(int i =0; i<sizeR; i++)
      var += pow(ptr[i][j]-mean,2);
    var /= sizeR;
   
    for(int i =0; i<sizeR; i++){
      ptr[i][j] =  (ptr[i][j]-mean) /sqrt(var+1e-6);
    }
  }
}


void My_Matrix::Enlarge(int r, int c){

    int i,j,k;
    double** ptr_buff;
    ptr_buff = new double*[sizeR];
    for ( k = 0; k < sizeR; k++){  
	 ptr_buff[k] = new double[sizeC];
     }

    for (i = 0; i < sizeR; i++)
	for (j = 0; j < sizeC; j++)
	    ptr_buff[i][j] = ptr[i][j];
    
    for( i = 0; i < sizeR; i++)
	delete [] ptr[i];
    delete [] ptr;
    
    
    ptr = new double*[r];
  
    for ( k = 0; k < r; k++){  
	ptr[k] = new double[c];

    }
    
    int min_r, min_c;
    if(sizeR > r)
	min_r = r;
    else
	min_r = sizeR;
    if(sizeC > c)
	min_c = c;
    else
	min_c = sizeC;

    for (i = 0; i < min_r; i++){
	for (j = 0; j < min_c; j++)
	    ptr[i][j] = ptr_buff[i][j];  
	for (j = min_c; j < c; j++)
	    ptr[i][j] = 0.0;
    }
    for (i = min_r; i < r; i++)
	for (j = 0; j < c; j++)
	    ptr[i][j] = 0.0;
    
    for( i = 0; i < sizeR; i++)
	delete [] ptr_buff[i];
    delete [] ptr_buff;
    
    sizeR = r;
    sizeC = c;
    
}
void My_Matrix::ReDimension(int r, int c){
    int i,j;
    
    for(i = 0; i < sizeR; i++)
	delete [] ptr[i];        
    delete [] ptr; 
    
    sizeR = r;   
    sizeC = c;
       
    ptr = new double*[sizeR];     
    for (j = 0; j < sizeR; j++){
	ptr[j] = new double[sizeC]; 
    }
}

My_Matrix My_Matrix::Rows0(int a, int b){
    int i,j;
    int rdim = b-a+1;

    My_Matrix sub(rdim,sizeC); 

    for(i = a; i <= b; i++)
      { 
     	for(j = 0; j < sizeC; j++)
	  {
          	sub[i-a][j] = ptr[i][j]; 
          }
      }

     return sub;
   
}

My_Matrix My_Matrix::Columns0(int a, int b){
    int i,j;
    int cdim = b-a+1;

    My_Matrix sub(sizeR,cdim); 

    for(i = 0; i < sizeR; i++)
      { 
     	for(j = a; j <= b; j++)
	  {
          	sub[i][j-a] = ptr[i][j]; 
          }
      }

     return sub;
   
}

My_Matrix My_Matrix::SubMatrix0(int ra, int rb,int ca, int cb){
    int i,j;
    int rdim = rb-ra+1;
    int cdim = cb-ca+1;

    My_Matrix sub(rdim,cdim); 

    for(i = ra; i <= rb; i++)
      { 
     	for(j = ca; j <= cb; j++)
	  {
          	sub[i-ra][j-ca] = ptr[i][j]; 
          }
      }

     return sub;
   
}

void My_Matrix::setRows0(int a, int b,const My_Matrix & M){

    int i,j;
    
    for(i = a; i <= b; i++)
      for(j = 0; j < sizeC; j++)
	ptr[i][j] = M.ptr[i-a][j]; 
        
}
void My_Matrix::setRows0(int a, int b,const My_Matrix_bool_scale & M){
  int i,j;
    
  for(i = a; i <= b; i++)
    for(j = 0; j < sizeC; j++)
      ptr[i][j] = M.scale * double(M.ptr[i-a][j]);
}
void My_Matrix::setRows0(int a, int b,const double & d){

    int i,j;
   

    for(i = a; i <= b; i++)
      for(j = 0; j < sizeC; j++)
	ptr[i][j] = d; 
        
}

void My_Matrix::setColumns0(int a, int b,const My_Matrix & M){

    int i,j;
   
    for(i = 0; i < sizeR; i++)
      for(j = a; j <= b; j++)
	ptr[i][j] = M.ptr[i][j-a]; 
        
}

void My_Matrix::setColumns0(int a, int b,const double & d){

    int i,j;
   
    for(i = 0; i < sizeR; i++)
      for(j = a; j <= b; j++)
	ptr[i][j] = d; 
        
}

void My_Matrix::setColumn0(int a,const double * d){

    int i;
   
    for(i = 0; i < sizeR; i++)
      ptr[i][a] = d[i]; 
        
}

My_Matrix My_Matrix::ExtractColumns(vector<int>& which){

  int cdim = which.size();
  My_Matrix sub(sizeR,cdim); 
  for( int i =0; i<sizeR; i++){
    for(int j =0; j<cdim;j++)
      sub[i][j] = ptr[i][which[j]];
   
  }
  return sub;
}

void My_Matrix::setSubMatrix0(int ra, int rb, int ca, int cb, const My_Matrix & M){

    int i,j;
   
    for(i = ra; i <= rb; i++)
      for(j = ca; j <= cb; j++)
	ptr[i][j] = M.ptr[i-ra][j-ca]; 
        
}

void My_Matrix::setSubMatrix0(int ra, int rb, int ca, int cb, const double & d){

    int i,j;
   
    for(i = ra; i <= rb; i++)
      for(j = ca; j <= cb; j++)
	ptr[i][j] = d; 
        
}

// Overloaded subscript operator
double *My_Matrix::operator[](int subscript)
{
   return ptr[subscript];    
}

const double * My_Matrix::operator[](int subscript) const {
   return ptr[subscript];  
}

const My_Matrix &My_Matrix::operator=(const double &d){
   for(int i = 0; i < sizeR; i++)
     for (int l = 0; l < sizeC; l++)
       ptr[i][l] = d;
   return *this; 
}

const My_Matrix &My_Matrix::operator=(const My_Matrix &right)
{
   int i,j,k,l,m;

   if (&right != this) {    // check for self-assignment
 
 
       for(i = 0; i < sizeR; i++)
	 delete [] ptr[i];        
   
       delete [] ptr; 
      
       sizeR = right.sizeR;    // resize this object
       sizeC = right.sizeC;
       
       ptr = new double*[sizeR];     
       for (j = 0; j < sizeR; j++){
	   ptr[j] = new double[sizeC];  // create space for My_Matrix copy
	 
	 }
       for (k = 0; k < sizeR; k++){
	   for (l = 0; l < sizeC; l++)
	       ptr[k][l] = right.ptr[k][l];  // copy My_Matrix into object
	 }
     }
   
   return *this;   // enables x = y = z;
}

ostream &operator<<(ostream &output, const My_Matrix &m)
{
 int i;
 int j;
   for (i = 0; i < m.sizeR; i++) 
     {
       for (j = 0; j < m.sizeC; j++)
	 output << m.ptr[i][j] << ' ';
       output << endl;
     }

   return output;   // enables cout << x << y;

}

istream &operator>>(istream &input, My_Matrix &m)
{
 int i;
 int j;
   for (i = 0; i < m.sizeR; i++) 
       for (j = 0; j < m.sizeC; j++)
	 input >> m.ptr[i][j] ;
       
     
   return input;   // enables cout << x << y;

}
// my_matrix = my_matrix - double
My_Matrix My_Matrix::operator - ( const double & d){ 

  int i, j;
 
  
  My_Matrix sum(sizeR,sizeC);
  
  for(i = 0; i < sizeR; i++)
      for(j = 0; j < sizeC; j++)
	sum[i][j] = ptr[i][j] - d; 
	
  return sum; 
}

// my_matrix = my_matrix + double
My_Matrix My_Matrix::operator + ( const double & d){ 

  int i, j;
 
  
  My_Matrix sum(sizeR,sizeC);
  
  for(i = 0; i < sizeR; i++)
      for(j = 0; j < sizeC; j++)
	sum[i][j] = ptr[i][j] + d; 
	
  return sum; 
}


// my_matrix = my_matrix * double
My_Matrix My_Matrix::operator * ( const double & d){ 

  int i, j;
 
 
  My_Matrix sum(sizeR,sizeC);
  
  for(i = 0; i < sizeR; i++)
      for(j = 0; j < sizeC; j++)
	sum[i][j] = ptr[i][j] * d; 
	
  return sum; 
}


// my_matrix = my_matrix - my_matrix
My_Matrix My_Matrix::operator - ( My_Matrix & M){ 

  int i, j;
 
  if(sizeC!= M.sizeC || sizeR!= M.sizeR )
    cout << "Wrong Size in my_matrix = my_matrix - my_matrix" << endl;

  My_Matrix sum(sizeR,sizeC);
  
  for(i = 0; i < sizeR; i++)
      for(j = 0; j < sizeC; j++)
	sum[i][j] = ptr[i][j] - M.ptr[i][j]; 
	
  return sum; 
}


My_Matrix My_Matrix::operator % (const My_Matrix & M) const{ 

  int i, j;
 
  if(sizeC!= M.sizeC || sizeR!= M.sizeR )
    cout << "Wrong Size in my_matrix = my_matrix ** my_matrix" << endl;

  My_Matrix prod(sizeR,sizeC);
  
  for(i = 0; i < sizeR; i++)
      for(j = 0; j < sizeC; j++)
	prod[i][j] = ptr[i][j] * M.ptr[i][j]; 
	
  return prod; 
}


My_Matrix My_Matrix::operator - (const  My_Matrix & M) const { 

  int i, j;
 
  if(sizeC!= M.sizeC || sizeR!= M.sizeR )
    cout << "Wrong Size in my_matrix = my_matrix + my_matrix" << endl;

  My_Matrix sum(sizeR,sizeC);
  
  for(i = 0; i < sizeR; i++)
      for(j = 0; j < sizeC; j++)
	sum[i][j] = ptr[i][j] - M.ptr[i][j]; 
	
  return sum; 
}

// my_matrix = my_matrix + my_matrix
My_Matrix My_Matrix::operator + ( My_Matrix & M){ 

  int i, j;
 
  if(sizeC!= M.sizeC || sizeR!= M.sizeR )
    cout << "Wrong Size in my_matrix = my_matrix + my_matrix" << endl;

  My_Matrix sum(sizeR,sizeC);
  
  for(i = 0; i < sizeR; i++)
      for(j = 0; j < sizeC; j++)
	sum[i][j] = ptr[i][j] + M.ptr[i][j]; 
	
  return sum; 
}

My_Matrix My_Matrix::operator + (const  My_Matrix & M) const { 

  int i, j;
 
  if(sizeC!= M.sizeC || sizeR!= M.sizeR )
    cout << "Wrong Size in my_matrix = my_matrix + my_matrix" << endl;

  My_Matrix sum(sizeR,sizeC);
  
  for(i = 0; i < sizeR; i++)
      for(j = 0; j < sizeC; j++)
	sum[i][j] = ptr[i][j] + M.ptr[i][j]; 
	
  return sum; 
}

// my_matrix = my_matrix * my_matrix
My_Matrix My_Matrix::operator * ( My_Matrix & M){ 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.sizeR )
    cout << "Wrong Size in  my_matrix * my_matrix" << endl;
  int ra = sizeR; 
  int cb = M.sizeC; 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  
   
  for(i = 0; i < ra; i++){
      for(j = 0; j < cb; j++){
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    sum += ptr[i][k] * M.ptr[k][j];
	      
	  prod[i][j] = sum; 
	}
    }
  return prod; 
}

My_Matrix My_Matrix::operator * (const My_Matrix & M) const { 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.sizeR )
    cout << "Wrong Size in  my_matrix * my_matrix" << endl;
  int ra = sizeR; 
  int cb = M.sizeC; 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  
   
  for(i = 0; i < ra; i++){
      for(j = 0; j < cb; j++){
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    sum += ptr[i][k] * M.ptr[k][j];
	      
	  prod[i][j] = sum; 
	}
    }
  return prod; 
}

// my_matrix = my_matrix * my_matrix_bool
My_Matrix My_Matrix::operator * (const My_Matrix_bool & M){ 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.ptr.size() )
    cout << "Wrong Size in my_matrix = my_matrix * my_matrix_bool" << endl;
  int ra = sizeR; 
  int cb = M.sizeC; 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  
   
  for(i = 0; i < ra; i++){
      for(j = 0; j < cb; j++){
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    if(M[k][j])
	      sum +=  ptr[i][k];
	  prod[i][j] = sum; 
	}
    }
  return prod; 
}

// my_matrix = my_matrix * my_matrix_bool_scale
My_Matrix My_Matrix::operator * (const My_Matrix_bool_scale & M){ 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.ptr.size() )
    cout << "Wrong Size in my_matrix = my_matrix * my_matrix_bool_scale" << endl;
  int ra = sizeR; 
  int cb = M.sizeC; 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  double scale = M.scale;
   
  for(i = 0; i < ra; i++){
      for(j = 0; j < cb; j++){
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    if(M[k][j])
	      sum +=  ptr[i][k];
	  prod[i][j] = scale*sum; 
	}
    }
  return prod; 
}

// my_vector = my_matrix * my_vector
My_Vector My_Matrix::operator * (const My_Vector & V){ 

  int i, j, k;
  double sum = 0;
  int len = V.size;
  if(sizeC!= len)
    cout << "Wrong Size in my_vector = my_matrix * my_vector" << endl;
 
  My_Vector prod(sizeR);
  
   
  for(i = 0; i < sizeR; i++){
    sum = 0;
    for(j = 0; j < len; j++){ 
      	sum +=ptr[i][j]* V.ptr[j];
    }
    prod[i] = sum; 
  }
  return prod; 
}



int My_Matrix::getSizeR() const { return sizeR; }

int My_Matrix::getSizeC() const { return sizeC; }

int My_Matrix::Nrows() const { return sizeR; }

int My_Matrix::Ncols() const { return sizeC; }

void My_Matrix::RandUniformFill(){

  for(int i = 0; i < sizeR; i++)
    for(int j = 0; j < sizeC; j++)
      ptr[i][j] = (double)rand()/(double)RAND_MAX;
}


double My_Matrix::Variance() const {
  double mean = 0.0;
  for(int i = 0; i < sizeR; i++)
    for(int j = 0; j < sizeC; j++)
      mean +=ptr[i][j];
  mean /= ( sizeR* sizeC);

  double var = 0.0;
  for(int i = 0; i < sizeR; i++)
    for(int j = 0; j < sizeC; j++)
      var += (ptr[i][j]-mean)*(ptr[i][j]-mean);

  return var/( sizeR* sizeC);
  
}

double My_Matrix::Variance_nonzero() const {
  double mean = 0.0;
  int counter =0;

  for(int i = 0; i < sizeR; i++)
    for(int j = 0; j < sizeC; j++){
      if(ptr[i][j] != 0){
	mean +=ptr[i][j];
	counter++;
      }
    }
  mean /= double(counter);

  double var = 0.0;
  for(int i = 0; i < sizeR; i++)
    for(int j = 0; j < sizeC; j++)
      if(ptr[i][j] != 0){
	var += (ptr[i][j]-mean)*(ptr[i][j]-mean);
      }
  return var/double(counter);
  
}

double My_Matrix::Sum() const {
  double sum = 0.0;
  for(int i = 0; i < sizeR; i++)
    for(int j = 0; j < sizeC; j++)
      sum +=ptr[i][j];
 
  return sum;
  
}

double My_Matrix::Mean_nonzero() const {
  double sum = 0.0;
  int counter =0;
  for(int i = 0; i < sizeR; i++)
    for(int j = 0; j < sizeC; j++){
      if(ptr[i][j] != 0){
	sum +=ptr[i][j];
	counter++;
      }
    }
 
  return sum/double(counter);
  
}

void My_Matrix::RandGaussFill(){
 double x1, x2, w, y1, y2;
 
 for(int i = 0; i < sizeR; i++)
   for(int j = 0; j < sizeC; j++){
     do {
       x1 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
       x2 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
       w = x1 * x1 + x2 * x2;
     } while ( w >= 1.0 );
     
     w = sqrt( (-2.0 * log( w ) ) / w );
     y1 = x1 * w;
     // y2 = x2 * w;
     ptr[i][j] = y1;
   }
}
 


My_Matrix Transpose(const My_Matrix& a){
    int i, j;
    int nr = a.sizeR; 
    int nc = a.sizeC;
    My_Matrix tra(nc,nr); 

    for(i = 0; i < nr; i++)
      { 
     	for(j = 0; j < nc; j++)
	  {
          	tra[j][i] = a[i][j]; 
          }
      }

     return tra;
     
}


My_Vector A_shifted_times_x(const My_Matrix &A, const My_Vector& x){ // computes Ax on a "shifted" matrix A (multi-task learning)

  int K = A.sizeC;
  int N = A.sizeR;
  int idx =0; 
  My_Vector res(K*N);
  for(int k =0; k<K;k++){
    for(int n =0; n< N; n++){
      res.ptr[idx] = A.ptr[n][k]*x.ptr[k];
      idx++;
    }
  }
  return res;
}
My_Vector A_shifted_T_times_x(const My_Matrix &A,  const My_Vector& x){ // computes A^t x on a "shifted" matrix A (multi-task learning)
 int K = A.sizeC;
 int N = A.sizeR;
  double sum;
  My_Vector res(K);
  int idx =0; 
  for(int k =0; k<K;k++){
    sum =0.0;
    for(int n =0; n< N; n++){
      sum +=  A.ptr[n][k]*x.ptr[idx];
      idx++;
    }
    res[k] = sum;
  }
  return res;
}
