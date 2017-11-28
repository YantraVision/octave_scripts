#include <stdio.h>
#include <math.h>


#define MAX_VALUES 4

extern int dpd_lut[];

int values[] = {0,1,512,789};
int cmplx_values[] = {1,1,512,200,64,78,789,1000};


void main()
{

//Generate the DPD values
int i;
int real;
int cmplx;
int abs_val;

printf("\nDPD\n");

for(i = 0; i < MAX_VALUES ; i++)
{
   real = dpd_lut[2*values[i]];
   cmplx = dpd_lut[2*values[i]+1];

  printf("\nreal %d cmplx %d \n",real,cmplx);
   		
}
#if 1
for(i = 0; i < MAX_VALUES ; i++)
{
   real = cmplx_values[2*values[i]];
   cmplx = cmplx_values[2*values[i]+1];  
	
   abs_val = (int)sqrt(real*real+cmplx*cmplx);
   
 if(abs_val > 1000)
 {
    abs_val = 1000;
 }

   real = dpd_lut[2*abs_val];  //Apply phase correction
   cmplx = dpd_lut[2*abs_val+1];   //Apply phase correction
   		
   printf("\nreal %d cmplx %d",real,cmplx);
}
#endif
}




int dpd_lut[] = {0};
