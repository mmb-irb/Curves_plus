/*********************************************
 ***   fct:   get date                     ***
 ***   call:  nml                          ***
 *********************************************/

# include <time.h>
# include <string.h>

void getdate_(char *d,int x)
     {
       time_t temp;       
       temp=time(&temp);
       strcpy(d,ctime(&temp));
       return;

     }





















