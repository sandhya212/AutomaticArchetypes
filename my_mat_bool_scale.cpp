
#include "my_mat_bool_scale.h"
#include <cmath>
#include <algorithm>


My_Matrix_bool_scale::My_Matrix_bool_scale(){
 

  vector<bool> v(1,false);
  ptr.push_back(v);
  sizeC = 1;
  scale = 1.0;
  is_up_to_date = false;
}

My_Matrix_bool_scale::My_Matrix_bool_scale(const My_Matrix_bool_scale& M){

   int k;

   ptr = M.ptr;
 
   sizeC = M.sizeC;
   scale = M.scale;
   is_up_to_date = M.is_up_to_date;
   rowInd = M.rowInd;
  
}


My_Matrix_bool_scale::My_Matrix_bool_scale(int mSizeR, int mSizeC)
{

  for (int  k = 0; k < mSizeR; k++){  
    vector<bool> v(mSizeC,false);
    ptr.push_back(v);
  }
  sizeC = mSizeC;
  scale = 1.0;
  is_up_to_date = false;
}

My_Matrix_bool_scale::My_Matrix_bool_scale(int mSizeR, int mSizeC, double mScale)
{

  for (int  k = 0; k < mSizeR; k++){  
    vector<bool> v(mSizeC,false);
    ptr.push_back(v);
  }
  sizeC = mSizeC;
  scale = mScale;
  is_up_to_date = false;
}

My_Matrix_bool_scale::~My_Matrix_bool_scale()
{  
  ptr.clear();
}



void My_Matrix_bool_scale::Enlarge(int r, int c){

  
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
void My_Matrix_bool_scale::ReDimension(int r, int c){

  ptr.clear();
  for (int i = 0; i < r; i++){
    vector<bool> v(c,false); 
    ptr.push_back(v);
  }
  sizeC = c;
  scale = 1.0;
}

My_Matrix_bool_scale My_Matrix_bool_scale::cbind( My_Matrix_bool_scale & M){

  if(ptr.size() != M.ptr.size()){
    cout << "wrong number of rows in cbind" <<endl;
  }
   if(scale != M.scale){
    cout << "wrong scales in cbind, keeping orig scale!" <<endl;
  }
  My_Matrix_bool_scale ret(ptr.size(), sizeC+M.sizeC,scale);
   
  for(int i = 0; i < ptr.size(); i++){
    for(int j = 0; j < sizeC; j++)
      ret.ptr[i][j] = ptr[i][j]; 
    for(int j = sizeC; j < sizeC+M.sizeC; j++)
      ret.ptr[i][j] = M.ptr[i][j-sizeC]; 
  }
  return ret;
  
}

My_Matrix_bool_scale My_Matrix_bool_scale::Rows0(int a, int b){
    int i;
    int rdim = b-a+1;

    My_Matrix_bool_scale sub(rdim,sizeC,scale); 
    
    for(i = a; i <= b; i++)
      { 
	sub[i-a] = ptr[i]; 
      }
    
    is_up_to_date = false;
    return sub;
   
}

My_Matrix_bool_scale My_Matrix_bool_scale::Columns0(int a, int b){
    int i,j;
    int cdim = b-a+1;
    int sizeR = ptr.size();

    My_Matrix_bool_scale sub(sizeR,cdim,scale); 

    for(i = 0; i < sizeR; i++)
      { 
     	for(j = a; j <= b; j++)
	  {
          	sub[i][j-a] = ptr[i][j]; 
          }
      }
    is_up_to_date = false;
    return sub;
   
}

My_Matrix_bool_scale My_Matrix_bool_scale::SubMatrix0(int ra, int rb,int ca, int cb){
    int i,j;
    int rdim = rb-ra+1;
    int cdim = cb-ca+1;

    My_Matrix_bool_scale sub(rdim,cdim,scale); 

    for(i = ra; i <= rb; i++)
      { 
     	for(j = ca; j <= cb; j++)
	  {
          	sub[i-ra][j-ca] = ptr[i][j]; 
          }
      }
    is_up_to_date = false;
    return sub;
    
}


void My_Matrix_bool_scale::setRows0(int a, int b,const My_Matrix_bool_scale & M){

    int i,j;
    
    for(i = a; i <= b; i++)
      ptr[i] = M.ptr[i-a]; 

    is_up_to_date = false;
        
}
void My_Matrix_bool_scale::setRows0(int a, int b,const double & d){

    int i,j;
   

    for(i = a; i <= b; i++)
      ptr[i].assign(sizeC,d);

    is_up_to_date = false;
    
        
}

void My_Matrix_bool_scale::setColumns0(int a, int b,const My_Matrix_bool_scale & M){

    int i,j;
   
    for(i = 0; i < ptr.size(); i++)
      for(j = a; j <= b; j++)
	ptr[i][j] = M.ptr[i][j-a]; 

    is_up_to_date = false;
        
}

void My_Matrix_bool_scale::setColumns0(int a, int b,const double & d){

    int i,j;
   
    for(i = 0; i < ptr.size(); i++)
      for(j = a; j <= b; j++)
	ptr[i][j] = d; 

    is_up_to_date = false;
        
}

void My_Matrix_bool_scale::setSubMatrix0(int ra, int rb, int ca, int cb, const My_Matrix_bool_scale & M){

    int i,j;
   
    for(i = ra; i <= rb; i++)
      for(j = ca; j <= cb; j++)
	ptr[i][j] = M.ptr[i-ra][j-ca]; 
        
    is_up_to_date = false;
}

void My_Matrix_bool_scale::setSubMatrix0(int ra, int rb, int ca, int cb, const double & d){

    int i,j;
   
    for(i = ra; i <= rb; i++)
      for(j = ca; j <= cb; j++)
	ptr[i][j] = d; 
    is_up_to_date = false;
}

My_Vector My_Matrix_bool_scale::multMtV(const My_Vector & V){ 

  int i, j;
  double sum = 0;
  int sizeR = ptr.size();
  int len = V.size;
  if(sizeR!= len)
    cout << "Wrong Size in my_vector = multMtV(my_vector)" << endl;
 
  My_Vector prod(sizeC);
 
   
  for(i = 0; i < sizeC; i++){
    sum = 0.0;
    for(j = 0; j < len; j++){ 
      if(ptr[j][i])
	sum += V.ptr[j];
    }
    prod[i] = scale*sum; 
  }
  return prod; 
}

double My_Matrix_bool_scale::multMtVNorm(const My_Vector & V){ 

  int i, j;
  double sum = 0;

  int len = V.size;

  int sizeR = ptr.size();
  if(sizeR!= len)
    cout << "Wrong Size in multMtVNorm(my_vector)" << endl;
 
  double norm =0;
 
  if(!is_up_to_date){
    update_RowIndices();
  }

  for(i = 0; i < sizeC; i++){
    sum = 0.0;
   
    for(j = 0; j < rowInd[i].size(); j++)
      sum += V.ptr[rowInd[i][j]];
    
    norm += sum *sum; 
  }
  return scale*sqrt(norm); 
}
double My_Matrix_bool_scale::multMtVNorm_L1(const My_Vector & V){ 

  int i, j;
  double sum = 0;

  int len = V.size;

  int sizeR = ptr.size();
  if(sizeR!= len)
    cout << "Wrong Size in multMtVNorm(my_vector)" << endl;
 
  double norm =0;
 
  if(!is_up_to_date){
    update_RowIndices();
  }

  for(i = 0; i < sizeC; i++){
    sum = 0.0;
   
    for(j = 0; j < rowInd[i].size(); j++)
      sum += V.ptr[rowInd[i][j]];
    
    norm += fabs(sum); 
  }
  return scale*norm; 
}

My_Vector My_Matrix_bool_scale::multMtVScale(const My_Vector & V, const double& dscale){ 

  int i, j;
  double sum = 0;
  int sizeR = ptr.size();
  int len = V.size;
  if(sizeR!= len)
    cout << "Wrong Size in my_vector = multMtVScale(my_vector, dscale)" << endl;
 
  My_Vector prod(sizeC);
   if(!is_up_to_date){
    update_RowIndices();
  }

   
  for(i = 0; i < sizeC; i++){
    sum = 0.0;
   
    for(j = 0; j < rowInd[i].size(); j++)
      sum += V.ptr[rowInd[i][j]];
    
    prod[i] = scale*sum*dscale; 
  }
  return prod; 
}

My_Vector My_Matrix_bool_scale::WplusMultMV(const My_Vector & V,const My_Vector & W ){ 

  int i, j;
  double sum = 0;
  int sizeR = ptr.size();
  int len = V.size ;
  if(sizeC!= len || sizeR != W.size)
    cout << "Wrong Size in my_vector = WplusMultMV(const My_Vector & V,const My_Vector & W)" << endl;
 
  My_Vector prod(sizeR);
 
   
  for(i = 0; i < sizeR; i++){
    sum = 0.0;
    for(j = 0; j < len; j++){ 
      if(ptr[i][j])
	sum += V.ptr[j];
    }
    prod[i] = scale*sum + W[i]; 
  }
  return prod; 
}




const My_Matrix_bool_scale &My_Matrix_bool_scale::operator=(const My_Matrix_bool_scale &M)
{
   int k;

   if (&M != this) {    // check for self-assignment
 
 

     ptr = M.ptr;      
     
     sizeC = M.sizeC;
     scale = M.scale;
     
   }
   is_up_to_date =  M.is_up_to_date;
   rowInd = M.rowInd;
   return *this;   // enables x = y = z;
}

// my_matrix = my_matrix_bool_scale * my_matrix_bool_scale
My_Matrix My_Matrix_bool_scale::operator * ( My_Matrix_bool_scale & M){ 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.ptr.size() )
    cout << "Wrong Size in  my_matrix = my_matrix_bool_scale * my_matrix_bool_scale" << endl;
  int ra = ptr.size(); 
  int cb = M.sizeC; 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  
  double sp = scale*M.scale; 
  for(i = 0; i < ra; i++)
    {
      for(j = 0; j < cb; j++)
	{
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    if(ptr[i][k] && M.ptr[k][j])
	      sum += sp;
	  prod[i][j] = sum; 
	}
    }
  return prod; 
}

My_Matrix My_Matrix_bool_scale::operator * (const My_Matrix_bool_scale & M) const { 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.ptr.size() )
    cout << "Wrong Size in  my_matrix = my_matrix_bool_scale * my_matrix_bool_scale" << endl;
  int ra = ptr.size(); 
  int cb = M.sizeC; 
  int ca = sizeC;
  My_Matrix prod(ra,cb);
  double sp = scale*M.scale; 
   
  for(i = 0; i < ra; i++)
    {
      for(j = 0; j < cb; j++)
	{
	  sum = 0;
	  for(k = 0; k < ca; k++) 
	    if(ptr[i][k] && M.ptr[k][j])
	      sum += sp;
	  prod[i][j] = sum; 
	}
    }
  return prod; 
}

// my_matrix = my_matrix_bool_scale * my_matrix
My_Matrix My_Matrix_bool_scale::operator * (const My_Matrix & M) const { 

  int i, j, k;
  double sum = 0;
  if(sizeC!= M.Nrows() )
    cout << "Wrong Size in my_matrix = my_matrix_bool_scale * my_matrix" << endl;
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
	  prod[i][j] = scale * sum; 
	}
    }
  return prod; 
}


// factor interaction
My_Matrix_bool_scale My_Matrix_bool_scale::operator | (const My_Matrix_bool_scale& M){
  int sizeR = ptr.size(), newC = sizeC*M.sizeC;
  My_Matrix_bool_scale inter(sizeR,newC,scale);
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
  /* 
  bool exclude_zero_columns = false;
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
  */
  
  return inter;
}



// my_vector = my_matrix_bool_scale * my_vector
My_Vector My_Matrix_bool_scale::operator * (const My_Vector & V){ 

  int i, j, k;
  double sum = 0;
  int sizeR = ptr.size();
  int len = V.size ;
  if(sizeC!= len)
    cout << "Wrong Size in my_vector = my_matrix_bool_scale * my_vector" << endl;
 
  My_Vector prod(sizeR);
 
   
  for(i = 0; i < sizeR; i++){
    sum = 0.0;
    for(j = 0; j < len; j++){ 
      if(ptr[i][j])
	sum += V.ptr[j];
    }
    prod[i] = scale*sum; 
  }
  return prod; 
}




ostream &operator<<(ostream &output, const My_Matrix_bool_scale &m)
{
 int i;
 int j;
   for (i = 0; i < m.ptr.size(); i++) 
     {
       for (j = 0; j < m.sizeC; j++)
	 output << m.scale*double(m.ptr[i][j]) << ' ';
       output << endl;
     }

   return output;   // enables cout << x << y;

}

// multiple factor interaction multV
My_Vector multi_interact_2multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Vector& V, const double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC;
  double scale = 1.0/pow(double(newC),exp);  

  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);

  int* jsj = new int[M_1.sizeC];
 
  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i* M_2.sizeC;
  }
 
  double sum;
  My_Vector P(sizeR);

  for(int i = 0; i < sizeR; i++){
    
    sum =0.;
    for(int j = 0; j < M_1.sizeC; j++){
	if(M_1.ptr[i][j]){
	  
	  for(int k=0; k< M_2.sizeC; k++)
	    if(M_2.ptr[i][k]){
	      
	      sum += V.ptr[jsj[j]  + k];
	   
	      // cout <<"| "<< i<<' '<<j<<' '<<k << endl;
		  
	      
	    }
	}
    }
    P[i] = sum * scale;
  }
  delete [] jsj;
  return P;
}

// multiple factor interaction multV
My_Vector multi_interact_3multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3, const My_Vector& V, const double& exp){
  
 
  int sj =  M_2.sizeC *  M_3.sizeC;
  int newC =  M_1.sizeC* sj;
  double scale = 1.0/pow(double(newC),exp);  
  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);

  int* jsj = new int[M_1.sizeC];
  int* ksk = new int[M_2.sizeC];

  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i*sj;
  }
  for(int i = 0; i < M_2.sizeC; i++){
    ksk[i] = i* M_3.sizeC;
  }
 

  double sum;
  My_Vector P(sizeR);
  for(int i = 0; i < sizeR; i++){
    sum =0.;
    for(int j = 0; j < M_1.sizeC; j++){
	if(M_1.ptr[i][j]){
	  
	  for(int k=0; k< M_2.sizeC; k++)
	    if(M_2.ptr[i][k]){
	      
	      for(int l=0; l< M_3.sizeC; l++)
		if(M_3.ptr[i][l]){
		  
	
		   sum += V.ptr[jsj[j] + ksk[k]  + l];
		  // cout <<"| "<< i<<' '<<j<<' '<<k << endl;
		
		}
	    }
	}
    }
    P[i] = sum * scale;
  }
  delete [] jsj;
  delete [] ksk;
  
  return P;
}

My_Vector multi_interact_4multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3,  const My_Matrix_bool_scale & M_4, const My_Vector& V, const double& exp){

  
 
  int sk =  M_3.sizeC* M_4.sizeC;
  int sj =  M_2.sizeC * sk;
  int newC =  M_1.sizeC* sj;
  double scale = 1.0/pow(double(newC),exp);  
  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);

  int* jsj = new int[M_1.sizeC];
  int* ksk = new int[M_2.sizeC];
  int* lsl = new int[M_3.sizeC];
  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i*sj;
  }
  for(int i = 0; i < M_2.sizeC; i++){
    ksk[i] = i*sk;
  }
  for(int i = 0; i < M_3.sizeC; i++){
    lsl[i] = i*M_4.sizeC;
   }
  double sum;
  My_Vector P(sizeR);
  for(int i = 0; i < sizeR; i++){
    sum =0.;
    for(int j = 0; j < M_1.sizeC; j++){
	if(M_1.ptr[i][j]){
	  
	  for(int k=0; k< M_2.sizeC; k++)
	    if(M_2.ptr[i][k]){
	      
	      for(int l=0; l< M_3.sizeC; l++)
		if(M_3.ptr[i][l]){
		  
		  for(int m=0; m< M_4.sizeC; m++)
		    if(M_4.ptr[i][m]){
		       sum += V.ptr[jsj[j] + ksk[k] + lsl[l] + m];
	     
		  }
	      }
	  }
	}
    }
    P[i] = sum * scale;
  }
  delete [] jsj;
  delete [] ksk;
  delete [] lsl;
  return P;
}



My_Vector multi_interact_5multV(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2,  const My_Matrix_bool_scale & M_3,  const My_Matrix_bool_scale & M_4, const My_Matrix_bool_scale & M_5, const My_Vector& V, const double& exp){

  
 
  int sl =  M_4.sizeC* M_5.sizeC;
  int sk =  M_3.sizeC* sl;
  int sj =  M_2.sizeC * sk;
  int newC =  M_1.sizeC* sj;
  double scale = 1.0/pow(double(newC),exp);  
  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);

  int* jsj = new int[M_1.sizeC];
  int* ksk = new int[M_2.sizeC];
  int* lsl = new int[M_3.sizeC];
  int* msm = new int[M_4.sizeC];
  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i*sj;
  }
  for(int i = 0; i < M_2.sizeC; i++){
    ksk[i] = i*sk;
  }
  for(int i = 0; i < M_3.sizeC; i++){
    lsl[i] = i*sl;
   }
 for(int i = 0; i < M_4.sizeC; i++){
    msm[i] = i*M_5.sizeC;
   }

 double sum;
 My_Vector P(sizeR);
  for(int i = 0; i < sizeR; i++){
    sum =0.;
    for(int j = 0; j < M_1.sizeC; j++){
	if(M_1.ptr[i][j]){
	  
	  for(int k=0; k< M_2.sizeC; k++)
	    if(M_2.ptr[i][k]){
	      
	      for(int l=0; l< M_3.sizeC; l++)
		if(M_3.ptr[i][l]){
		  
		  for(int m=0; m< M_4.sizeC; m++)
		    if(M_4.ptr[i][m]){

		      for(int n=0; n< M_5.sizeC; n++)
			if(M_5.ptr[i][n]){
			   sum += V.ptr[jsj[j] + ksk[k] + lsl[l] + msm[m] +n];
	     
			}
		    }
		}
	    }
	}
    }
    P[i] = sum * scale;
  }
  delete [] jsj;
  delete [] ksk;
  delete [] lsl;
  delete [] msm;
  return P;
}




// multiple factor interaction
My_Matrix_bool_scale multi_interact_2(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC;
  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);

  int* jsj = new int[M_1.sizeC];
 
  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i* M_2.sizeC;
  }
 

  for(int i = 0; i < sizeR; i++){
    
    for(int j = 0; j < M_1.sizeC; j++){
	if(M_1.ptr[i][j]){
	  
	  for(int k=0; k< M_2.sizeC; k++)
	    if(M_2.ptr[i][k]){
	      
	    
	      inter[i][jsj[j]  + k] = true;
	      // cout <<"| "<< i<<' '<<j<<' '<<k << endl;
		  
	      
	    }
	}
    }
  }
  delete [] jsj;
  return inter;
}



// multiple factor interaction
My_Matrix_bool_scale multi_interact_3(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3){
  
 
 
  int sj =  M_2.sizeC *  M_3.sizeC;
  int newC =  M_1.sizeC* sj;
  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);

  int* jsj = new int[M_1.sizeC];
  int* ksk = new int[M_2.sizeC];

  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i*sj;
  }
  for(int i = 0; i < M_2.sizeC; i++){
    ksk[i] = i* M_3.sizeC;
  }
 
  for(int i = 0; i < sizeR; i++){
    
    for(int j = 0; j < M_1.sizeC; j++){
	if(M_1.ptr[i][j]){
	  
	  for(int k=0; k< M_2.sizeC; k++)
	    if(M_2.ptr[i][k]){
	      
	      for(int l=0; l< M_3.sizeC; l++)
		if(M_3.ptr[i][l]){
		  
	
		  inter[i][jsj[j] + ksk[k]  + l] = true;
		  // cout <<"| "<< i<<' '<<j<<' '<<k << endl;
		
		}
	    }
	}
    }
  }
  delete [] jsj;
  delete [] ksk;
  
  return inter;
}



// multiple factor interaction
My_Matrix_bool_scale multi_interact_4(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3, const My_Matrix_bool_scale & M_4){
  
 
  int sk =  M_3.sizeC* M_4.sizeC;
  int sj =  M_2.sizeC * sk;
  int newC =  M_1.sizeC* sj;
  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);

  int* jsj = new int[M_1.sizeC];
  int* ksk = new int[M_2.sizeC];
  int* lsl = new int[M_3.sizeC];
  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i*sj;
  }
  for(int i = 0; i < M_2.sizeC; i++){
    ksk[i] = i*sk;
  }
  for(int i = 0; i < M_3.sizeC; i++){
    lsl[i] = i*M_4.sizeC;
   }

  for(int i = 0; i < sizeR; i++){
    
    for(int j = 0; j < M_1.sizeC; j++){
	if(M_1.ptr[i][j]){
	  
	  for(int k=0; k< M_2.sizeC; k++)
	    if(M_2.ptr[i][k]){
	      
	      for(int l=0; l< M_3.sizeC; l++)
		if(M_3.ptr[i][l]){
		  
		  for(int m=0; m< M_4.sizeC; m++)
		    if(M_4.ptr[i][m]){
		      inter[i][jsj[j] + ksk[k] + lsl[l] + m] = true;
	     
		  }
	      }
	  }
	}
    }
  }
  delete [] jsj;
  delete [] ksk;
  delete [] lsl;
  return inter;
}




// multiple factor interaction
My_Matrix_bool_scale multi_interact_5(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3, const My_Matrix_bool_scale & M_4, const My_Matrix_bool_scale & M_5){
  
 
  int sl =  M_4.sizeC* M_5.sizeC;
  int sk =  M_3.sizeC* sl;
  int sj =  M_2.sizeC * sk;
  int newC =  M_1.sizeC* sj;
  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);

  int* jsj = new int[M_1.sizeC];
  int* ksk = new int[M_2.sizeC];
  int* lsl = new int[M_3.sizeC];
  int* msm = new int[M_4.sizeC];
  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i*sj;
  }
  for(int i = 0; i < M_2.sizeC; i++){
    ksk[i] = i*sk;
  }
  for(int i = 0; i < M_3.sizeC; i++){
    lsl[i] = i*sl;
   }
 for(int i = 0; i < M_4.sizeC; i++){
    msm[i] = i*M_5.sizeC;
   }

  for(int i = 0; i < sizeR; i++){
    
    for(int j = 0; j < M_1.sizeC; j++){
	if(M_1.ptr[i][j]){
	  
	  for(int k=0; k< M_2.sizeC; k++)
	    if(M_2.ptr[i][k]){
	      
	      for(int l=0; l< M_3.sizeC; l++)
		if(M_3.ptr[i][l]){
		  
		  for(int m=0; m< M_4.sizeC; m++)
		    if(M_4.ptr[i][m]){

		      for(int n=0; n< M_5.sizeC; n++)
			if(M_5.ptr[i][n]){
			  inter[i][jsj[j] + ksk[k] + lsl[l] + msm[m] +n] = true;
	     
			}
		    }
		}
	    }
	}
    }
  }
  delete [] jsj;
  delete [] ksk;
  delete [] lsl;
  delete [] msm;
  return inter;
}





// multiple factor interaction multTvecNorm
double multi_interact_2_multTvecNorm( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2,const  My_Vector& V,const  double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC;
  int sizeR =  V.size;
  double scale = 1.0/pow(double(newC),exp);  
  double res =0.;
  double sum;    


  if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   

  

  for(int j = 0; j < M_1.sizeC; j++){
  
    for(int k=0; k< M_2.sizeC; k++){
	      
      vector<int> setIntersection;
      set_intersection(M_1.rowInd[j].begin(),M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
		       back_inserter(setIntersection));
      
      sum =0.;

      for(int i = 0; i < setIntersection.size(); i++){
	sum += V.ptr[setIntersection[i]];
      }
      res += sum*sum;  
    }
  }
  
  

  return scale*sqrt(res);
}


// multiple factor interaction multTvecNorm
double multi_interact_3_multTvecNorm( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3, const My_Vector& V,  const double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC * M_3.sizeC;
  int sizeR =  V.size;
  double scale = 1.0/pow(double(newC),exp);  
  double res =0.;
  double sum;    


  if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   if(!M_3.is_up_to_date){
     M_3.update_RowIndices();
   }
   

 
  for(int j = 0; j < M_1.sizeC; j++){
  
    for(int k=0; k< M_2.sizeC; k++){

      vector<int> setIntersection;
      set_intersection(M_1.rowInd[j].begin(), M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
		       back_inserter(setIntersection));
      
	      
	      
      for(int l=0; l< M_3.sizeC; l++){

	vector<int> setIntersection2;
	set_intersection( setIntersection.begin(),  setIntersection.end(), M_3.rowInd[l].begin(),  M_3.rowInd[l].end(),
			 back_inserter(setIntersection2));
		  
	sum =0.;

	for(int i = 0; i < setIntersection2.size(); i++){
	    sum += V.ptr[setIntersection2[i]];
	  }
	res += sum*sum;
	
      }
    }
  }
  

  return scale*sqrt(res);
}


double multi_interact_4_multTvecNorm( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3,  My_Matrix_bool_scale & M_4, const My_Vector& V,  const double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC * M_3.sizeC * M_4.sizeC;
  int sizeR =  V.size;
  double scale = 1.0/pow(double(newC),exp);  
  double sum;    
  double res=0.;

 


  if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   if(!M_3.is_up_to_date){
     M_3.update_RowIndices();
   }
   if(!M_4.is_up_to_date){
     M_4.update_RowIndices();
   }
   



  

 
 
    
  for(int j = 0; j < M_1.sizeC; j++){
  
    for(int k=0; k< M_2.sizeC; k++){


      vector<int> setIntersection;
      set_intersection(M_1.rowInd[j].begin(), M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
		       back_inserter(setIntersection));

	      
      for(int l=0; l< M_3.sizeC; l++){

	vector<int> setIntersection2;
	set_intersection( setIntersection.begin(),  setIntersection.end(), M_3.rowInd[l].begin(), M_3.rowInd[l].end(),
			 back_inserter(setIntersection2));
		  
	for(int m=0; m< M_4.sizeC; m++){

	  vector<int> setIntersection3;
	  set_intersection( setIntersection2.begin(),  setIntersection2.end(), M_4.rowInd[m].begin(),  M_4.rowInd[m].end(),
			    back_inserter(setIntersection3));

	  sum =0.;

	
	  for(int i = 0; i < setIntersection3.size(); i++){
	    sum += V.ptr[setIntersection3[i]];
	  }
	  
	  res += sum*sum;
	  
	}
      }
    }
  }
  

  return scale*sqrt(res);
}



double multi_interact_5_multTvecNorm( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3,  My_Matrix_bool_scale & M_4,  My_Matrix_bool_scale & M_5, const My_Vector& V,  const double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC * M_3.sizeC * M_4.sizeC * M_5.sizeC;
  int sizeR =  V.size;
  double scale = 1.0/pow(double(newC),exp);  
  double sum;    
  double res=0.;

 

 if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   if(!M_3.is_up_to_date){
     M_3.update_RowIndices();
   }
   if(!M_4.is_up_to_date){
     M_4.update_RowIndices();
   }
   if(!M_5.is_up_to_date){
     M_5.update_RowIndices();
   }
   
 
    
  for(int j = 0; j < M_1.sizeC; j++){
  
    for(int k=0; k< M_2.sizeC; k++){


      vector<int> setIntersection;
      set_intersection(M_1.rowInd[j].begin(), M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
		       back_inserter(setIntersection));

	      
      for(int l=0; l< M_3.sizeC; l++){

	vector<int> setIntersection2;
	set_intersection( setIntersection.begin(),  setIntersection.end(), M_3.rowInd[l].begin(),  M_3.rowInd[l].end(),
			 back_inserter(setIntersection2));
		  
	for(int m=0; m< M_4.sizeC; m++){

	  vector<int> setIntersection3;
	  set_intersection( setIntersection2.begin(),  setIntersection2.end(), M_4.rowInd[m].begin(),  M_4.rowInd[m].end(),
			    back_inserter(setIntersection3));


	  for(int n=0; n< M_5.sizeC; n++){
	    
	    vector<int> setIntersection4;
	    set_intersection( setIntersection3.begin(),  setIntersection3.end(), M_5.rowInd[n].begin(),  M_5.rowInd[n].end(),
			      back_inserter(setIntersection4));
	    
	    
	    sum =0.;
	  
	  
	    for(int i = 0; i < setIntersection4.size(); i++){
	      sum += V.ptr[setIntersection4[i]];
	    }
	  
	    res += sum*sum;
	  }
	  
	}
      }
    }
  }
  

  return scale*sqrt(res);
}




// multiple factor interaction multTvecNorm
double multi_interact_2_multTvecNorm_L1( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2,const  My_Vector& V,const  double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC;
  int sizeR =  V.size;
  double scale = 1.0/pow(double(newC),exp);  
  double res =0.;
  double sum;    


  if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   

  

  for(int j = 0; j < M_1.sizeC; j++){
  
    for(int k=0; k< M_2.sizeC; k++){
	      
      vector<int> setIntersection;
      set_intersection(M_1.rowInd[j].begin(),M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
		       back_inserter(setIntersection));
      
      sum =0.;

      for(int i = 0; i < setIntersection.size(); i++){
	sum += V.ptr[setIntersection[i]];
      }
      res += fabs(sum);  
    }
  }
  
  

  return scale*res;
}


// multiple factor interaction multTvecNorm
double multi_interact_3_multTvecNorm_L1( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3, const My_Vector& V,  const double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC * M_3.sizeC;
  int sizeR =  V.size;
  double scale = 1.0/pow(double(newC),exp);  
  double res =0.;
  double sum;    


  if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   if(!M_3.is_up_to_date){
     M_3.update_RowIndices();
   }
   

 
  for(int j = 0; j < M_1.sizeC; j++){
  
    for(int k=0; k< M_2.sizeC; k++){

      vector<int> setIntersection;
      set_intersection(M_1.rowInd[j].begin(), M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
		       back_inserter(setIntersection));
      
	      
	      
      for(int l=0; l< M_3.sizeC; l++){

	vector<int> setIntersection2;
	set_intersection( setIntersection.begin(),  setIntersection.end(), M_3.rowInd[l].begin(),  M_3.rowInd[l].end(),
			 back_inserter(setIntersection2));
		  
	sum =0.;

	for(int i = 0; i < setIntersection2.size(); i++){
	    sum += V.ptr[setIntersection2[i]];
	  }
	res += fabs(sum);
	
      }
    }
  }
  

  return scale*res;
}


double multi_interact_4_multTvecNorm_L1( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3,  My_Matrix_bool_scale & M_4, const My_Vector& V,  const double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC * M_3.sizeC * M_4.sizeC;
  int sizeR =  V.size;
  double scale = 1.0/pow(double(newC),exp);  
  double sum;    
  double res=0.;

 


  if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   if(!M_3.is_up_to_date){
     M_3.update_RowIndices();
   }
   if(!M_4.is_up_to_date){
     M_4.update_RowIndices();
   }
   



  

 
 
    
  for(int j = 0; j < M_1.sizeC; j++){
  
    for(int k=0; k< M_2.sizeC; k++){


      vector<int> setIntersection;
      set_intersection(M_1.rowInd[j].begin(), M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
		       back_inserter(setIntersection));

	      
      for(int l=0; l< M_3.sizeC; l++){

	vector<int> setIntersection2;
	set_intersection( setIntersection.begin(),  setIntersection.end(), M_3.rowInd[l].begin(), M_3.rowInd[l].end(),
			 back_inserter(setIntersection2));
		  
	for(int m=0; m< M_4.sizeC; m++){

	  vector<int> setIntersection3;
	  set_intersection( setIntersection2.begin(),  setIntersection2.end(), M_4.rowInd[m].begin(),  M_4.rowInd[m].end(),
			    back_inserter(setIntersection3));

	  sum =0.;

	
	  for(int i = 0; i < setIntersection3.size(); i++){
	    sum += V.ptr[setIntersection3[i]];
	  }
	  
	  res += fabs(sum);
	  
	}
      }
    }
  }
  

  return scale*res;
}



double multi_interact_5_multTvecNorm_L1( My_Matrix_bool_scale & M_1,  My_Matrix_bool_scale & M_2,  My_Matrix_bool_scale & M_3,  My_Matrix_bool_scale & M_4,  My_Matrix_bool_scale & M_5, const My_Vector& V,  const double& exp){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC * M_3.sizeC * M_4.sizeC * M_5.sizeC;
  int sizeR =  V.size;
  double scale = 1.0/pow(double(newC),exp);  
  double sum;    
  double res=0.;

 

 if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   if(!M_3.is_up_to_date){
     M_3.update_RowIndices();
   }
   if(!M_4.is_up_to_date){
     M_4.update_RowIndices();
   }
   if(!M_5.is_up_to_date){
     M_5.update_RowIndices();
   }
   
 
    
  for(int j = 0; j < M_1.sizeC; j++){
  
    for(int k=0; k< M_2.sizeC; k++){


      vector<int> setIntersection;
      set_intersection(M_1.rowInd[j].begin(), M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
		       back_inserter(setIntersection));

	      
      for(int l=0; l< M_3.sizeC; l++){

	vector<int> setIntersection2;
	set_intersection( setIntersection.begin(),  setIntersection.end(), M_3.rowInd[l].begin(),  M_3.rowInd[l].end(),
			 back_inserter(setIntersection2));
		  
	for(int m=0; m< M_4.sizeC; m++){

	  vector<int> setIntersection3;
	  set_intersection( setIntersection2.begin(),  setIntersection2.end(), M_4.rowInd[m].begin(),  M_4.rowInd[m].end(),
			    back_inserter(setIntersection3));


	  for(int n=0; n< M_5.sizeC; n++){
	    
	    vector<int> setIntersection4;
	    set_intersection( setIntersection3.begin(),  setIntersection3.end(), M_5.rowInd[n].begin(),  M_5.rowInd[n].end(),
			      back_inserter(setIntersection4));
	    
	    
	    sum =0.;
	  
	  
	    for(int i = 0; i < setIntersection4.size(); i++){
	      sum += V.ptr[setIntersection4[i]];
	    }
	  
	    res += fabs(sum);
	  }
	  
	}
      }
    }
  }
  

  return scale*res;
}


void  My_Matrix_bool_scale::update_RowIndices(){

  rowInd.clear();
  int sizeR =  ptr.size();
  for(int j = 0; j < sizeC; j++){
    vector<int> r1;
    for(int i = 0; i < sizeR; i++){
      if(ptr[i][j])
	r1.push_back(i);
    }
    rowInd.push_back(r1);
  }
  is_up_to_date = true;
}

int My_Matrix_bool_scale::getSizeR() const { return ptr.size(); }

int My_Matrix_bool_scale::getSizeC() const { return sizeC; }

int My_Matrix_bool_scale::Nrows() const { return  ptr.size(); }

int My_Matrix_bool_scale::Ncols() const { return sizeC; }


int My_Matrix_bool_scale::NNZcols() const {
 int sizeR = ptr.size();
 int nnz =0;
 for (int i = 0; i < sizeC; i++){
   for(int j =0; j< sizeR; j++){
     if(ptr[j][i]){
       nnz++;
       break;
     }
   }
 }
 // cout << sizeC <<' '<<nnz << endl;
 return nnz; 
}



















My_Matrix_bool_scale Transpose(const My_Matrix_bool_scale& a)
{
    int i, j;
    int nr = a.ptr.size(); 
    int nc = a.sizeC;
    My_Matrix_bool_scale tra(nc,nr,a.scale); 

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
  My_Matrix_bool_scale A(3,3);
  A[0][0] = 1;
  A[1][0] = 1;
  A.el(1,1) = true;
  cout << A << endl;
  My_Matrix_bool_scale B(3,2);
  B[0][0] = 1;
  B[0][1] = 1;
  B[1][0] = true;
  cout <<B;
  My_Matrix_bool_scale C =  A | B;
  cout << C;
  My_Matrix_bool_scale D = A.cbind(C);

  cout << D;
}
*/
