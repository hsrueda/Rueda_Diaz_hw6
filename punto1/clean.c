#include<stdio.h>
 
int main()
{ 
  int i=0;
  for (i=0; i<31; i++){
    remove('poblaciones_'+str(i)+'_20.dat');
    remove('poblaciones_'+str(i)+'_20.pdf');
  }
  
  return 0;
}
