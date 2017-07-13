#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

//#define P 2 //2 dimensions for now
#define INF 9999999 //no need for c++11
#define BUCKETSIZE 10

//tree
struct Node{
	int index; // which point I am
	Node *left;
	Node *right; 
	Node (int i) { index = i; left = right = NULL; }
};

void print_tree(Node *T){
	 if (T == NULL) return;
  print_tree(T->left);
  std::cout << T->index << "\n";
  print_tree(T->right);
	
}

//for small trees
Node *miniInsert(Node *Tree, double *coords, int index, int d,int n){

  int P = 2;
  
  if(Tree==NULL){
    return new Node(index);
  }
  
  if(coords[index]<=coords[Tree->index]&&d==0){
    Tree->left=miniInsert(Tree->left,coords,index,(d+1)%P,n);
  }
  
  if(coords[index]>coords[Tree->index]&&d==0){ 
    Tree->right=miniInsert(Tree->right,coords,index,(d+1)%P,n);
  }
  
  if(coords[index+n]<=coords[Tree->index+n]&&d==1){
    Tree->left=miniInsert(Tree->left,coords,index,(d+1)%P,n);
  }
  
  if(coords[index+n]>coords[Tree->index+n]&&d==1){ 
    Tree->right=miniInsert(Tree->right,coords,index,(d+1)%P,n);
  }
  
  return Tree;
}

void getNNIndx(int i, int m, int &iNNIndx, int &iNN){
  
  if(i == 0){
    iNNIndx = 0;//this should never be accessed
    iNN = 0;
    return;
  }else if(i < m){
    iNNIndx = static_cast<int>(static_cast<double>(i)/2*(i-1));
    iNN = i;
    return;
  }else{
    iNNIndx = static_cast<int>(static_cast<double>(m)/2*(m-1)+(i-m)*m);
    iNN = m;
    return;
  } 
}

//as long as we have the numbers and know its 2d, just simple formula
double eucdist2d(double &x0, double &y0, double &x1, double &y1){
	return(sqrt(pow(x0-x1,2)+pow(y0-y1,2)));
}

// double eucdist3d(double &x0,double &y0,double &z0, double &x1, double &y1,double &z1){
// 	return(sqrt(pow(x0-x1,2)+pow(y0-y1,2)+pow(z0-z1,2)));
// }

//a dist, b index
void fSort(double *a, int *b, int n){
  
  int j, k, l;
  double v;
  
  for(j = 1; j <= n-1; j++){
    k = j;  
    while(k > 0 && a[k] < a[k-1]) {
      v = a[k]; l = b[k];
      a[k] = a[k-1]; b[k] = b[k-1];
      a[k-1] = v; b[k-1] = l;
      k--;
    }
  }
}

void get_nn(Node *Tree, int index, int d, double *coords,int n, double *nnDist, int *nnIndx, int iNNIndx, int iNN, int check){

  int P = 2;
  
  if(Tree==NULL){
    return;
  }
  
  double disttemp= eucdist2d(coords[index],coords[index+n],coords[Tree->index],coords[Tree->index+n]); 
  
  if(index!=Tree->index && disttemp<nnDist[iNNIndx+iNN-1]){
    nnDist[iNNIndx+iNN-1]=disttemp;
    nnIndx[iNNIndx+iNN-1]=Tree->index;
    fSort(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
  }
  
  Node *temp1=Tree->left;
  Node *temp2=Tree->right;
  
  if(d==0){
    
    if(coords[index]>coords[Tree->index]){
      std::swap(temp1,temp2);
    }
    
    get_nn(temp1,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN, check);
    
    if(fabs(coords[Tree->index]-coords[index])>nnDist[iNNIndx+iNN-1]){
      return;
    }
    
    get_nn(temp2,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN, check);
  }
  
  if(d==1){
    
    if(coords[index+n]>coords[Tree->index+n]){
      std::swap(temp1,temp2);
    }
    
    get_nn(temp1,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN,check);

    if(fabs(coords[Tree->index+n]-coords[index+n])>nnDist[iNNIndx+iNN-1]){
      return;
    }
    
    get_nn(temp2,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN,check);
  }

}

void mkNNIndx(int n, int m, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU){
  
	int i, iNNIndx, iNN;
	double d;
       	int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  
  	for(i = 0; i < nIndx; i++){
	  nnDist[i] = INF;
	}
	
	Node *Tree=NULL;
	int time_through=-1;
	
	for(i=0;i<n;i++){
	  getNNIndx(i, m, iNNIndx, iNN);
	  nnIndxLU[i] = iNNIndx;
	  nnIndxLU[n+i] = iNN;
	  if(time_through==-1){
	    time_through=i;
	  }
	  
	  if(i!=0){
	    for(int j = time_through; j < i; j++){ 
	      getNNIndx(i, m, iNNIndx, iNN);
	      d = eucdist2d(coords[i], coords[i+n], coords[j], coords[n+j]);
	      if(d < nnDist[iNNIndx+iNN-1]){
		nnDist[iNNIndx+iNN-1] = d;
		nnIndx[iNNIndx+iNN-1] = j;
		
		fSort(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN); 	  
	      }
	    }
	    
	    
	    if(i%BUCKETSIZE==0){
	      
              #pragma omp parallel for private(iNNIndx, iNN)
	      for(int j=time_through;j<time_through+BUCKETSIZE;j++){
		
		getNNIndx(j, m, iNNIndx, iNN);
		get_nn(Tree,j,0, coords,n, nnDist,nnIndx,iNNIndx,iNN,i-BUCKETSIZE);
	      }
	      
	      
	      for(int j=time_through;j<time_through+BUCKETSIZE;j++){
		Tree=miniInsert(Tree,coords,j,0, n);
	      }
	      
	      time_through=-1;
	    }
	    if(i==n-1){
	      
              #pragma omp parallel for private(iNNIndx, iNN)
	      for(int j=time_through;j<n;j++){
		getNNIndx(j, m, iNNIndx, iNN);
		get_nn(Tree,j,0, coords,n, nnDist,nnIndx,iNNIndx,iNN,i-BUCKETSIZE);
	      }
	      
	    }
	  } 
	  if(i==0){
	    Tree=miniInsert(Tree,coords,i,0,n);
	    time_through=-1;
	  }
	}
}

//dont put this in r library
int main(){
	int m=30;
	int n=1000000;
	double *coords=new double[n*2];
	int *nnIndxLU=new int[n*2];
	int *nnIndx=new int[n*m];
	double *nnDist=new double[n*m];

	std::ifstream fin("mygenmytest.txt");
	for(int i=0;i<n;i++){
		fin >> coords[i];
		fin >> coords[i+n];
	}
	mkNNIndx(n,m,coords,nnIndx,nnDist,nnIndxLU);
	//zeros coming up at end b/c of the nIndx calculation above
	for(int i=0;i<(n*2);i++){
		std::cout <<"LU:" << nnIndxLU[i] << "\n";
	}
	for(int i=0;i< (n*m);i++){
		std::cout <<"nnIndx:" << nnIndx[i]<<"\n";
		std::cout << "nnDist:" << nnDist[i] << std::endl;
	}
	return 0;
}
