#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define Rt 6378137
#define m 1.67*(10**-27)
#define c 299792458
#define c2 8.9875518*(10**16)
#define q 1.6*(10**(-19))
#define Bo 3*(10**(-5))
#define PI 3.14159265


int main(int argc, char **argv){
//Definicion de variables
  double Ec0=atof(argv[1]);
  double alpha=atof(argv[2]);
  double Ecj=Ec0*e*1000000;
  double gamma=1+Ecj/(m*c2);
  double vi=c*sqrt(1-(1/pow(gamma,2)));
  double t0=0.0;
  double tf=100.0;
  double npuntos=(tf-t0)/(2*PI*gamma*m/(e*Bo));

}

double funcion1(double x, double y, double z, double x_prima, double y_prima, double z_prima, double gm){

  double r=pow(x,2)+pow(y,2)+pow(z,2);
  double func1=e*(-Bo*(pow(Rt,3)))*(z_prima*(3*y*z)-x_prima*(2*z*(z-x)*(x-y)*y))/(m*gama*pow(r,5.0/2.0));
  return func1; 

}

double funcion2(double x, double y, double z, double x_prima, double y_prima, double z_prima, double gm){

  double r=pow(x,2)+pow(y,2)+pow(z,2);
  double func2=e*(-Bo*(pow(Rt,3)))*(z_prima*(3*x*z)-x_prima*(2*z*(z-x)*(x-y)*y))/(m*gama*pow(r,5.0/2.0));
  return func2;
}

double funcion3(double x, double y, double z, double x_prima, double y_prima, double z_prima, double gm){

  double r=pow(x,2)+pow(y,2)+pow(z,2);
  double func3=e*(-Bo*(pow(Rt,3)))*(z_prima*(3*x*z)-x_prima*(2*z*(z-x)*(x-y)*y))/(m*gama*pow(r,5.0/2.0));
  return func3;
}

"sebastianrosales"



