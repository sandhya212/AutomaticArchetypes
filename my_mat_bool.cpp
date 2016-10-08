
#include "my_mat_bool.h"
#include <cmath>

My_Matrix_bool::My_Matrix_bool(){
 

  vector<bool> v(1,false);
  ptr.push_back(v);
  sizeC = 1;
}

My_Matrix_bool::My_Matrix_bool(const My_Matrix_bool& M){

   int k;

   ptr = M.ptr;
 
   sizeC = M.sizeC;

  
}


My_Matrix_bool::My_Matrix_bool(int mSizeR, int mSizeC)
{

  for (int  k = 0; k < mSizeR; k++){  
    vector<bool> v(mSizeC,false);
    ptr.push_back(v);
  }
  sizeC = mSizeC;
}

My_Matrix_bool::~My_Matrix_bool()
{  
  ptr.clear();
}



void My_Matrix_bool::Enlarge(int r, int c){

  
    vector < vector < bool> > ptr_buff = ptr;

    ptr.clear();

    for (int i = 0; i < ptr_buff.size(); i++){
      vector<bool> v = ptr[i];
      for(int j =sizeC; j< c; j++)
	v.push_back(false);
      ptr.push_back(v);
    }
    
    ptr_buff.clear();

    for (int i = ptr_buff.size(); i < r; i++){
      vector<bool> v(c,false); 
      ptr.push_back(v);
    }

    sizeC = c;

    
}
void My_Matrix_bool::ReDimension(int r, int c){

  ptr.clear();
  for (int i = 0; i < r; i++){
    vector<bool> v(c,false); 
    ptr.push_back(v);
  }
  sizeC = c;
}

My_Matrix_bool My_Matrix_bool::cbind( My_Matrix_bool & M){

  if(ptr.size() != M.ptr.size()){
    cout << "wrong number of rows in cbind" <<endl;
  }
  My_Matrix_bool ret(ptr.size(), sizeC+M.sizeC);
   
  for(int i = 0; i < ptr.size(); i++){
    for(int j = 0; j < sizeC; j++)
      ret.ptr[i][j] = ptr[i][j]; 
    for(int j = sizeC; j < sizeC+M.sizeC; j++)
      ret.ptr[i][j] = M.ptr[i][j-sizeC]; 
  }
  return ret;
  
}

My_Matrix_bool My_Matrix_bool::Rows0(int a, int b){
    int i;
    int rdim = b-a+1;

    My_Matrix_bool sub(rdim,sizeC); 

    for(i = a; i <= b; i++)
      { 
	sub[i-a] = ptr[i]; 
      }

     return sub;
   
}

My_Matrix_bool My_Matrix_bool::Columns0(int a, int b){
    int i,j;
    int cdim = b-a+1;
    int sizeR = ptr.size();

    My_Matrix_bool sub(sizeR,cdim); 

    for(i = 0; i < sizeR; i++)
      { 
     	for(j = a; j <= b; j++)
	  {
          	sub[i][j-a] = ptr[i][j]; 
          }
      }

     return sub;
   
}

My_Matrix_bool My_Matrix_bool::SubMatrix0(int ra, int rb,int ca, int cb){
    int i,j;
    int rdim = rb-ra+1;
    int cdim = cb-ca+1;

    My_Matrix_bool sub(rdim,cdim); 

    for(i = ra; i <= rb; i++)
      { 
     	for(j = ca; j <= cb; j++)
	  {
          	sub[i-ra][j-ca] = ptr[i][j]; 
          }
      }

     return sub;
   
}


void My_Matrix_bool::setRows0(int a, int b,const My_Matrix_bool & M){

    int i,j;
    
    for(i = a; i <= b; i++)
      ptr[i] = M.ptr[i-a]; 
        
}
void My_Matrix_bool::setRows0(int a, int b,const double & d){

    int i,j;
   

    for(i = a; i <= b; i++)
      ptr[i].assign(sizeC,d);
    
        
}

void My_Matrix_bool::setColumns0(int a, int b,const My_Matrix_bool & M){

    int i,j;
   
    for(i = 0; i < ptr.size(); i++)
      for(j = a; j <= b; j++)
	ptr[i][j] = M.ptr[i][j-a]; 
        
}

void My_Matrix_bool::setColumns0(int a, int b,const double & d){

    int i,j;
   
    for(i = 0; i < ptr.size(); i++)
      for(j = a; j <= b; j++)
	ptr[i][j] = d; 
        
}

void My_Matrix_bool::setSubMatrix0(int ra, int rb, int ca, int cb, const My_Matrix_bool & M){

    int i,j;
   
    for(i = ra; i <= rb; i++)
      for(j = ca; j <= cb; j++)
	ptr[i][j] = M.ptr[i-ra][j-ca]; 
        
}

void My_Matrix_bool::setSubMatrix0(int ra, int rb, int ca, int cb, const double & d){

    int i,j;
   
    for(i = ra; i <= rb; i++)
      for(j = ca; j <= cb; j++)
	ptr[i][j] = d; 
        
}




const My_Matrix_bool &My_Matrix_bool::operator=(const My_Matrix_bool &M)
{
   int k;

   if (&M != this) {    // check for self-assignment
 
 

     ptr = M.ptr;      
     
     sizeC = M.sizeC;

     
   }
   
   return *this;   // enables x = y = z;
}

// my_matrix = my_matrix_bool * my_matrix_bool
My_Matrix My_Matrix_bool::operator * ( My_Matrix_bool & M){ 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.ptr.size() )
    cout << "Wrong Size in  my_matrix = my_matrix_bool * my_matrix_bool" << endl;
  int ra = ptr.size(); 
  int cb = M.sizeC; 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  
   
  for(i = 0; i < ra; i++)
    {
      for(j = 0; j < cb; j++)
	{
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    if(ptr[i][k] && M.ptr[k][j])
	      sum += 1;
	  prod[i][j] = sum; 
	}
    }
  return prod; 
}

My_Matrix My_Matrix_bool::operator * (const My_Matrix_bool & M) const { 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.ptr.size() )
    cout << "Wrong Size in  my_matrix = my_matrix_bool * my_matrix_bool" << endl;
  int ra = ptr.size(); 
  int cb = M.sizeC; 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  
   
  for(i = 0; i < ra; i++)
    {
      for(j = 0; j < cb; j++)
	{
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    if(ptr[i][k] && M.ptr[k][j])
	      sum += 1;
	  prod[i][j] = sum; 
	}
    }
  return prod; 
}

// my_matrix = my_matrix_bool * my_matrix
My_Matrix My_Matrix_bool::operator * (const My_Matrix & M) const { 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.Nrows() )
    cout << "Wrong Size in my_matrix = my_matrix_bool * my_matrix" << endl;
  int ra = ptr.size(); 
  int cb = M.Ncols(); 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  
   
  for(i = 0; i < ra; i++)
    {
      for(j = 0; j < cb; j++)
	{
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    if(ptr[i][k])
	      sum +=  M.ptr[k][j];
	  prod[i][j] = sum; 
	}
    }
  return prod; 
}


// factor interaction
My_Matrix_bool My_Matrix_bool::operator | (const My_Matrix_bool& M){
  int sizeR = ptr.size(), newC = sizeC*M.sizeC;
  My_Matrix_bool inter(sizeR,newC);
  for(int i = 0; i < sizeR; i++){
 
    for(int j = 0; j < sizeC; j++){
      if(ptr[i][j]){
	for(int k=0; k< M.sizeC; k++)
	  if(M.ptr[i][k]){
	     inter[i][j*M.sizeC+k] = true;
	     // cout <<"| "<< i<<' '<<j<<' '<<k << endl;
	  }
      }
    }
    
  } 
  bool exclude_zero_columns = true;
  if(exclude_zero_columns){
    bool flag;
    for(int j = 0; j < newC; j++){
      flag = false;
      for(int i = 0; i < sizeR; i++){
	if( inter[i][j]){
	  flag = true;
	  break;
	}
      }
      if(!flag){
	for(int i = 0; i < sizeR; i++)
	  inter.ptr[i].erase(inter.ptr[i].begin()+j);
	inter.sizeC--;
      }
    }
  }
  
  return inter;
}
// my_vector = my_matrix_bool * my_vector
My_Vector My_Matrix_bool::operator * (const My_Vector & V){ 

  int i, j, k;
  double sum = 0;
  int sizeR = ptr.size();
  int len = V.Length() ;
  if(sizeC!= len)
    cout << "Wrong Size in my_vector = my_matrix_bool * my_vector" << endl;
 
  My_Vector prod(sizeR);
  
   
  for(i = 0; i < sizeR; i++){
    sum = 0;
    for(j = 0; j < len; j++){ 
      if(ptr[i][j])
	sum += V[j];
    }
    prod[i] = sum; 
  }
  return prod; 
}




ostream &operator<<(ostream &output, const My_Matrix_bool &m)
{
 int i;
 int j;
   for (i = 0; i < m.ptr.size(); i++) 
     {
       for (j = 0; j < m.sizeC; j++)
	 output << bool(m.ptr[i][j]) << ' ';
       output << endl;
     }

   return output;   // enables cout << x << y;

}


 

int My_Matrix_bool::getSizeR() const { return ptr.size(); }

int My_Matrix_bool::getSizeC() const { return sizeC; }

int My_Matrix_bool::Nrows() const { return  ptr.size(); }

int My_Matrix_bool::Ncols() const { return sizeC; }



















My_Matrix_bool Transpose(const My_Matrix_bool& a)
{
    int i, j;
    int nr = a.ptr.size(); 
    int nc = a.sizeC;
    My_Matrix_bool tra(nc,nr); 

    for(i = 0; i < nr; i++)
      { 
     	for(j = 0; j < nc; j++)
	  {
          	tra[j][i] = a[i][j]; 
          }
      }

     return tra;
     
}

/*

main(){
  My_Matrix_bool A(3,3);
  A[0][0] = 1;
  A[1][0] = 1;
  A.el(1,1) = true;
  cout << A << endl;
  My_Matrix_bool B(3,2);
  B[0][0] = 1;
  B[0][1] = 1;
  B[1][0] = true;
  cout <<B;
  My_Matrix_bool C =  A | B;
  cout << C;
  My_Matrix_bool D = A.cbind(C);

  cout << D;
}
*/
