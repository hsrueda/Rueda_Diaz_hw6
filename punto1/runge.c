#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//Funciones 
float dxdt(float t, float x,float y, float A, float B);
float dydt(float t, float x,float y, float C, float D);
float *generatevector(float *vector, int m);
float *RungeKutta(float t_old, float x_old, float y_old, float h, float A, float B, float C, float D);
void  write_txt(float *t, float *x, float *y, float x0, float y0, int m_points);

// Main
int main(int argc, char **argv){

  float *t=NULL, *x=NULL, *y=NULL, *temporal=NULL, m_points_float=0.0;
  float a=0.0, b=3.0, x0 = 0.0, y0 = 0.0, h=0.01, A = 20.0, B = 1.0, C = 30.0, D = 1.0;
  float t_old=0.0, x_old=0.0, y_old=0.0, x_prima = 0.0, y_prima=0.0;
  int m_points=0, i=0, n=3;
  
  x0 = *argv[1];
  y0 = *argv[2];
 
  m_points_float = (b-a)/h;
  m_points = (int) m_points_float;
  
  t = generatevector(t, m_points);
  x = generatevector(x, m_points);
  y = generatevector(y, m_points);
  temporal = generatevector(temporal, n);

  t[0]=a;
  x[0]=x0;
  y[0]=y0;

  for( i=1 ; i<m_points ; i++ ){

    x_prima = dxdt(t[i-1], x[i-1], y[i-1], A, B);
    y_prima = dydt(t[i-1], x[i-1], y[i-1], C, D);
    
    t_old = t[i-1];
    x_old = x[i-1];
    y_old = y[i-1];
    
    temporal = RungeKutta(t_old, x_old, y_old, h, A, B, C, D);
    
    t[i] = temporal[0];
    x[i] = temporal[1];
    y[i] = temporal[2];
  }
  
  write_txt(t, x, y, x0, y0, m_points);
  
  return 0;
}

// Definicion de las funciones
float dxdt(float t, float x,float y, float A, float B){
  return A*x - B*x*y;
}

float dydt(float t, float x,float y, float C, float D){
  return -C*y + D*x*y;
}

float *generatevector(float *vect,int m){
  int i=0;
  if(!(vect = malloc(sizeof(float)*m))){
    fprintf(stderr, "Problem with memory allocation");
  }
  for( i=0 ; i<m ; i++ ){ 
    vect[i] = 0.0;
  }
  return vect;
}

float *RungeKutta(float t_old, float x_old,float y_old, float h, float A, float B, float C, float D){
  
  float *retorno = NULL;
  float k1_x=0.0, k2_x=0.0, k3_x=0.0, k4_x=0.0;
  float k1_y=0.0, k2_y=0.0, k3_y=0.0, k4_y=0.0;
  float t1=0.0, t2=0.0, t3=0.0, t4=0.0;
  float x_1=0.0, x_2=0.0, x_3=0.0, average_x=0.0;
  float y_1=0.0, y_2=0.0, y_3=0.0, average_y=0.0;
  float t_new=0.0, x_new=0.0, y_new=0.0;
 
  retorno = generatevector(retorno, 3);
  
  k1_x = dxdt(t_old, x_old, y_old, A, B);
  k1_y = dydt(t_old, x_old, y_old, C, D);

  //First step
  t1 = t_old + (h/2.0);
  x_1 = x_old + (h/2.0) * k1_x;
  y_1 = y_old + (h/2.0) * k1_y;
  
  k2_x = dxdt(t1, x_1, y_1, A, B);
  k2_y = dydt(t1, x_1, y_1, C, D);
  
  //Second step
  t2 = t_old + (h/2.0);
  x_2 = x_old + (h/2.0) * k2_x;
  y_2 = y_old + (h/2.0) * k2_y;
  
  k3_x = dxdt(t2, x_2, y_2, A, B);
  k3_y = dydt(t2, x_2, y_2, C, D);
  
  //Third step
  t3 = t_old + (h/1.0);
  x_3 = x_old + (h/1.0) * k3_x;
  y_3 = y_old + (h/1.0) * k3_y;
  
  k4_x = dxdt(t3, x_3, y_3, A, B);
  k4_y = dydt(t3, x_3, y_3, C, D);
  
  //Fourth step
  average_x = (1.0/6.0)*(k1_x+2.0*k2_x+2.0*k3_x+k4_x);
  average_y = (1.0/6.0)*(k1_y+2.0*k2_y+2.0*k3_y+k4_y);
  
  t_new = t_old + h;
  x_new = x_old + h*average_x;
  y_new = y_old + h*average_y;
  
  retorno[0] = t_new;
  retorno[1] = x_new;
  retorno[2] = y_new;

 return retorno;
}

void  write_txt(float *t,float *x, float *y, float x0, float y0, int m_points){
  FILE *fileout;
  int i=0;
  char name[20];
  sprintf(name, "poblaciones_%f_%f.dat", x0, y0);
  fileout = fopen(name,"w");
  for(i=0;i<(m_points);i++){
    fprintf(fileout,"%f %f %f\n",t[i],x[i],y[i]);
  }
}
