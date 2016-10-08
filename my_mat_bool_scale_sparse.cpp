
#include "my_mat_bool_scale.h"
#include <cmath>


My_Matrix_bool_scale::My_Matrix_bool_scale(){
 
  sizeC = 1;
  scale = 1.0;
  is_up_to_date = false;
}

My_Matrix_bool_scale::My_Matrix_bool_scale(const My_Matrix_bool_scale& M){

 
 
   sizeC = M.sizeC;
   scale = M.scale;
   is_up_to_date = M.is_up_to_date;
   rowInd = M.rowInd;
  
}

My_Matrix_bool_scale::My_Matrix_bool_scale(const vector<vector< int> >& v, double s = 1){

 
 
  sizeC = v.size();
  scale = s;
  is_up_to_date = false;
  rowInd = v;
  
}
My_Matrix_bool_scale::My_Matrix_bool_scale(int cdim, double s = 1){

 
 
   sizeC = cdim;
   scale = s;
   is_up_to_date = false;
   rowInd.resize(cdim);
  
}





My_Matrix_bool_scale::~My_Matrix_bool_scale()
{  
  ptr.clear();
}




My_Matrix_bool_scale My_Matrix_bool_scale::cbind( My_Matrix_bool_scale & M){


   if(scale != M.scale){
    cout << "wrong scales in cbind, keeping orig scale!" <<endl;
  }

   vector<vector<int> > v = rowInd ;

 
  for(int j = 0; j < M.sizeC; j++)
    v.push_back(M.rowInd[j]); 

  My_Matrix_bool_scale ret(v,scale);
 
  return ret;
  
}



double My_Matrix_bool_scale::multMtVNorm(const My_Vector & V){ 

  int i, j;
  double sum = 0;

  int len = V.size;

  int sizeR = ptr.size();
  if(sizeR!= len)
    cout << "Wrong Size in multMtVNorm(my_vector)" << endl;
 
  if(! is_up_to_date)
    update_RowIndices();
  double norm =0;
 
   
  for(i = 0; i < sizeC; i++){
    sum = 0.0;
   
    for(j = 0; j < rowInd[i].size(); j++)
      sum += V.ptr[rowInd[i][j]];
    
    norm += sum *sum; 
  }
  return scale*sqrt(norm); 
}

My_Vector My_Matrix_bool_scale::multMtVScale(const My_Vector & V, const double& dscale){ 

  int i, j;
  double sum = 0;
  int sizeR = ptr.size();
  int len = V.size;
  if(sizeR!= len)
    cout << "Wrong Size in my_vector = multMtVScale(my_vector, dscale)" << endl;
 
  My_Vector prod(sizeC);
  if(! is_up_to_date)
    update_RowIndices();
 
   
  for(i = 0; i < sizeC; i++){
    sum = 0.0;
    
    for(j = 0; j < rowInd[i].size(); j++)
      sum += V.ptr[rowInd[i][j]];
    
    prod[i] = scale*sum*dscale; 
  }
  return prod; 
}





const My_Matrix_bool_scale &My_Matrix_bool_scale::operator=(const My_Matrix_bool_scale &M)
{


   if (&M != this) {    // check for self-assignment
 
 

  
     
     sizeC = M.sizeC;
     scale = M.scale;
     
   }
   is_up_to_date =  M.is_up_to_date;
   rowInd = M.rowInd;
   return *this;   // enables x = y = z;
}






// multiple factor interaction
My_Matrix_bool_scale multi_interact_2(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2){
  
 
 
  int newC =  M_1.sizeC*  M_2.sizeC;
  int sizeR =  M_1.ptr.size();
  My_Matrix_bool_scale inter(sizeR,newC,1);




  
  if(!M_1.is_up_to_date){
    M_1.update_RowIndices();
  }
   if(!M_2.is_up_to_date){
    M_2.update_RowIndices();
  }
   
   int* jsj = new int[M_1.sizeC];
 
  for(int i = 0; i < M_1.sizeC; i++){
    jsj[i] = i* M_2.sizeC;
  }
 
  vector<vector<int> > vI(newC, vector<int>(0));

   for(int j = 0; j < M_1.sizeC; j++){
  
     for(int k=0; k< M_2.sizeC; k++){
	      
       
       set_intersection(M_1.rowInd[j].begin(),M_1.rowInd[j].end(), M_2.rowInd[k].begin(),  M_2.rowInd[k].end(),
			back_inserter( vI[jsj[j]  + k]));
       vI[jsj[j]  + k] = setIntersection;
       
     }
   }



 
   My_Matrix_bool_scale inter(vI,1);

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
 
  vector<vector<int> > vI(newC, vector<int>(0));

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


	set_intersection( setIntersection.begin(),  setIntersection.end(), M_3.rowInd[l].begin(),  M_3.rowInd[l].end(),
			 back_inserter(vI[jsj[j] + ksk[k]  + l]));

	
      }
    }
  }
  


  delete [] jsj;
  delete [] ksk;

  My_Matrix_bool_scale inter(vI,1);
    
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
double multi_interact_2_multTvecNorm(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Vector& V, const double& exp){
  
 
 
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
double multi_interact_3_multTvecNorm(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3, const My_Vector& V,  const double& exp){
  
 
 
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


double multi_interact_4_multTvecNorm(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3, const My_Matrix_bool_scale & M_4, const My_Vector& V,  const double& exp){
  
 
 
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



double multi_interact_5_multTvecNorm(const My_Matrix_bool_scale & M_1, const My_Matrix_bool_scale & M_2, const My_Matrix_bool_scale & M_3, const My_Matrix_bool_scale & M_4, const My_Matrix_bool_scale & M_5, const My_Vector& V,  const double& exp){
  
 
 
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
