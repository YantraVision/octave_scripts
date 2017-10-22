#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <cmath>
#include <iomanip>


using namespace std;

#define number_of_points 1 //Number of colm.
#define interpolated_point_size 48 // Num of rows

int xv[11] = {0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0} ;
int yv[11] = {0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0} ;


int Input[48] = { 1,  2,  3,  4,  5,  6,  7,  8,
                     9, 10, 11, 12, 13, 14, 15, 16,
                    17, 18, 19, 20, 21, 22, 23, 24,
                    25, 26, 27, 28, 29, 30, 31, 32,
                    33, 34, 35, 36, 37, 38, 39, 40,
                    41, 42, 43, 44, 45, 46, 47, 48
                  };

int ExpOutput[48] = {
   5,    6,    7,         8,       9,    10,        11,       12,
   13,   14,   15,        60,     17,    18,        19,       20,
   21,   22,   23,        24,     25,    26,        550,     -456,
   29,   30,   31,        32,     33,    34,        35,       36,
   37,   63743,  -115774, 5364,   41,    42,        43,       44,
   45,   46,   47,        48,   76081,  -2090380,   191715,  -58512
   };

/*
double Input[48] = {  
                     1.0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0
                   };


double ExpOutput[48] = {
   5.0000000e+00,  -4.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,
   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   4.4000000e+01,  -8.8000000e+01,   4.4000000e+01,   0.0000000e+00,   0.0000000e+00,
   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   5.2800000e+02,  -1.5400000e+03,
   1.4960000e+03,  -4.8400000e+02,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,
   0.0000000e+00,   6.3360000e+03,  -2.4288000e+04,   3.4892000e+04,  -2.2264000e+04,   5.3240000e+03,   0.0000000e+00,   0.0000000e+00,
   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   0.0000000e+00,   7.6032000e+04,  -3.6115200e+05,   6.8587200e+05,  -6.5098000e+05

};

*/


int output_vector[interpolated_point_size];

//double aCoeff[11] = {1.0e-05*0.0448, 0, -1.0e-05*0.2242, 0, 1.0e-05*0.4484, 0, -1.0e-05*0.4484, 0, 1.0e-05*0.2242, 0, -1.0e-05*0.0448};
//double bCoeff[11] = {1, -9.6276, 41.7283, -107.2220, 180.8803, -209.32, 168.3021, -92.8287, 33.6150, -7.2165, 0.6975};    

int aCoeff[11] = {1,2,3,4,5,6,7,8,9,10,11};
int bCoeff[11] = {5,6,7,8,9,10,11,12,13,14,15};

int aCoeffScaled[11];
int bCoeffScaled[11];


//Filtering
void butterWorthFilter(int *input_points, int *xv, int *yv,int gainValue, int *aCoeff, int *bCoeff, int *output_vector){
   for(int u=0; u< interpolated_point_size; u++) {
      xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4]; xv[4]=xv[5]; 
      xv[5] = xv[6]; xv[6] = xv[7]; xv[7] = xv[8]; xv[8] = xv[9]; xv[9]=xv[10];
      xv[10] = (input_points[u] / 1.0);
      yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; yv[4]= yv[5]; 
      yv[5] = yv[6]; yv[6] = yv[7]; yv[7] = yv[8]; yv[8] = yv[9]; yv[9]= yv[10]; 


      yv[10] =  ((bCoeff[0]  * xv[10])
                +(bCoeff[1]  * xv[9])
       	        +(bCoeff[2]  * xv[8])
				+(bCoeff[3]  * xv[7])
				+(bCoeff[4]  * xv[6])
				+(bCoeff[5]  * xv[5])
				+(bCoeff[6]  * xv[4])
				+(bCoeff[7]  * xv[3])
				+(bCoeff[8]  * xv[2])
				+(bCoeff[9]  * xv[1])
				+(bCoeff[10] * xv[0])

                -(aCoeff[1]  * yv[9])
				-(aCoeff[2]  * yv[8])
				-(aCoeff[3]  * yv[7])
				-(aCoeff[4]  * yv[6])
				-(aCoeff[5]  * yv[5])
				-(aCoeff[6]  * yv[4])
				-(aCoeff[7]  * yv[3])
				-(aCoeff[8]  * yv[2])
				-(aCoeff[9]  * yv[1])
				-(aCoeff[10] * yv[0]));

      output_vector[u] = yv[10];		
   }
}


void checkOutput () {
     for (int i=0; i < 48; i++) {
         if ( output_vector[i] <=  0.95*ExpOutput[i] || output_vector[i] >=  1.05*ExpOutput[i] ) {
            cout << "@Index " << i << " output is " <<  output_vector[i] << " expected is " <<  ExpOutput[i] <<  " diff is " <<  (output_vector[i] - ExpOutput[i])/ExpOutput[i] << endl ;
         } else {
            cout << "@Index " << i << " output is " <<  output_vector[i] << " expected is " <<  ExpOutput[i] <<  " diff is " <<  (output_vector[i] - ExpOutput[i])/ExpOutput[i] << endl ;
         }
     }
}


int main(int argc,  char** argv){
    //Filter Parameters
    //int NZEROS = 10;
    //int NPOLES = 10;
    int GAIN  = 1 ; // 2.235408541e+06;
    for (int i =0; i < 11; i++) {
       aCoeffScaled[i] = aCoeff[i]/aCoeff[0]; 
       bCoeffScaled[i] = bCoeff[i]/aCoeff[0]; 
    }
    //butterworthFilter5thOrder(Input, xv, yv, GAIN,	aCoeff, bCoeff, output_vector);
    butterWorthFilter(Input, xv, yv, GAIN, aCoeffScaled, bCoeffScaled, output_vector);
    //butterWorthInternet(Input, xv, yv, GAIN, aCoeffScaled, bCoeffScaled, output_vector);
    checkOutput();
}
