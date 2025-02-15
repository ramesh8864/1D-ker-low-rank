#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <iostream>
#include <cmath>

class Chebyshev {
private:
  double* nodes;
public:
    int size;       // Number of nodes
     // Array to store Chebyshev nodes

    // Constructor 1: Chebyshev nodes in [-1,1]
    Chebyshev(int n) {
        size = n;
        nodes = new double[n];
        for (int i = 0; i < n; ++i) {
            nodes[i] = cos(M_PI * (2 * i + 1) / (2 * n));
        }
    }

    // Constructor 2: Chebyshev nodes in [a,b]
    Chebyshev(double a, double b, int n) {
        size = n;
        nodes = new double[n];
        for (int i = 0; i < n; ++i) {
            double x_i = cos(M_PI * (2 * i + 1) / (2 * n));  // Nodes in [-1,1]
            nodes[i] = 0.5 * ((b - a) * x_i + (b + a));      // Map to [a,b]
        }
    }

    // Destructor to free allocated memory
    ~Chebyshev() {
        delete[] nodes;
    }
    
     double lagrange(int n,int i,double z){
         double coeff=1;
          for(int j=0;j<n;j++){
                if(j!=i)
               coeff= coeff*(z-nodes[j])/(nodes[i]-nodes[j]);
              // std::cout<<coeff<<"  ";
          }
          return coeff;
     }
     //Function to get the value  at ith position
     
     double get(int i){
        return nodes[i];
     }

    // Function to print nodes
    void printNodes() const {
        for (int i = 0; i < size; ++i) {
            std::cout << nodes[i] << " ";
        }
        std::cout << std::endl;
    }
};

#endif // CHEBYSHEV_H

