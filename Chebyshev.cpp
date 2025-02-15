#include <iostream>
#include "chebyshev.h"

using namespace std;
  
  //Kernel function
double ker(double x,double y){
             return 1.0/fabs(x-y);
 }
 
 //Main function.
 int main() {
   
   double a1,a2,b1,b2;
   int n1,n2,m1,m2;
   
   //Entering source detail
   std::cout<<"Enter the source interval: a1=";
   std::cin>>a1;
   std::cout<<"   b1=";
   std::cin>>b1;
   std::cout<<"how many source points you have, n1=";
   std::cin>>n1;
   
   double source_point[n1];
   
   for(int i=0;i<n1;i++){
      std::cout<<"Enter source point["<<i<<"]:";
   std::cin>>source_point[i]; 
   }
   std::cout<<"how many source nodes(chebyshev nodes) you want, n2=";
   std::cin>>n2;
   
   
   Chebyshev source(a1,b1,n2);
   //std::cout<<"source nodes are:"<<std::endl;
   //source.printNodes();
    
   double lag1[n1][n2];
   
   for(int i=0;i<n1;i++){
        for(int j=0;j<n2;j++){
        lag1[i][j]= source.lagrange(n2,j,source_point[i]);
    //    std::cout<<lag1[i][j]<<"   ";
        }
      // std::cout<<std::endl;
   }
    
    //Entering target detail
   std::cout<<"Enter the target interval: a2=";
   std::cin>>a2;
   std::cout<<"   b2=";
   std::cin>>b2;
   std::cout<<"how many target points you have, n1="<<std::endl;
   std::cin>>m1;
   
   double target_point[m1];
   
   for(int i=0;i<m1;i++){
      std::cout<<"Enter target point["<<i<<"]:";
   std::cin>>target_point[i]; 
   }
   std::cout<<"how many target nodes(chebyshev nodes) you want, m2="<<std::endl;
   std::cin>>m2;

   
   Chebyshev target(a2,b2,m2);
   //std::cout<<"target nodes are:"<<std::endl;
  // target.printNodes();
   
   double lag2[n1][n2];
   
   for(int i=0;i<m1;i++){
        for(int j=0;j<m2;j++){
        lag2[i][j]= target.lagrange(m2,j,target_point[i]);
     //   std::cout<<lag2[i][j]<<"   ";
        }
        //std::cout<<std::endl;
   }
    std::cout<<"Interaction matrix"<< std::endl;
   double exact[n1][m1];   // interaction between source and target points
   
   for(int i=0;i<n1;i++){
       
       for(int j=0;j<m1;j++){
         
         exact[i][j]= ker(source_point[i],target_point[j]);
         std::cout<< exact[i][j]<<"  ";
         }
       std::cout<< std::endl;
     
}
std::cout<< std::endl;
   double ker_nodes[n2][m2];   // interaction between chebyshev nodes in source and target intervals
  // std::cout<< "Interaction matrix of chebyshev nodes"<<std::endl;
   
   for(int i=0;i<n2;i++){
       
       for(int j=0;j<m2;j++){
         
         ker_nodes[i][j]= ker(source.get(i),target.get(j));
        // std::cout<< ker_nodes[i][j]<<"  ";
         }
     //  std::cout<< std::endl;
     }
     std::cout<<"Interaction matrix obtained by low rank approximation"<< std::endl;
     double approx[n1][m1];   // low rank approximation
     
     for(int i=0;i<n1;i++){
       
          for(int j=0;j<m1;j++){
            double sum=0;
            for(int k=0;k<m2;k++){
              double temp=0;
              for(int l=0;l<n2;l++){
                temp = temp + ker_nodes[l][k]*lag1[i][l];
               }
            sum = sum + temp*lag2[j][k];
            }
            approx[i][j] = sum;
            std::cout<<approx[i][j]<<"  ";
        }
        std::cout<<std::endl;
   }
   
   
  return 0; 
}

    

