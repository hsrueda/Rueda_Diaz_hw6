#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define Rt 6378137
#define m 1.67E-27
#define c 299792458
#define c2 8.9875518E16
#define q 1.6E-19
#define Bo 3E-5
#define PI 3.14159265
#define e5 1E-5

double funcion1(double x,double y, double z,double x_prima, double y_prima,double z_prima, double gama);
double funcion2(double x,double y, double z,double x_prima, double y_prima,double z_prima, double gama);
double funcion3(double x,double y, double z,double x_prima, double y_prima,double z_prima, double gama);

int main(int argc, char **argv){
//Definicion de variables
  double Eco=atof(argv[1]);
  double alpha=atof(argv[2]);
  double Ecj=Eco*q*1000000;
  double gamma=1+Ecj/(m*c2);
  double vi=c*sqrt(1-(1/pow(gamma,2)));
  double t0=0.0;
  double tf=100.0;
  double puntos=(tf-t0)/(2.0*PI*gamma*m/(q*Bo));
  double h=(2.0*PI*gamma*m/(q*Bo))/50.0;

//Conversion de unidades

  double radianes=(alpha*PI/180.0);

  //Condiciones iniciales del trabajo
  double xi=2*Rt;
  double yi=0.0;
  double zi=0.0;
  double x_primai=0;
  double y_primai=vi*sin(radianes);
  double z_primai=vi*cos(radianes);

  //Se abre el archivo y se escribe la primera linea del archivo
  FILE *file;
  char nombre[1000];
  sprintf(nombre, "trayectoria_%d_%d.dat", (int)Eco, (int)alpha);
  file=fopen(nombre,"w");
  fprintf(file, "%f %f %f %f \n",0.0,xi,yi,zi);

  //Aplicamos metodo Rongue Kutta para poder esncontrar las posiciones del proton y al final se escribe en el archivo.
  int i;
  for(i=1;i<puntos;i++){

    double Kx1=funcion1(xi,yi,zi,x_primai,y_primai,z_primai,gamma);
    double Ky1=funcion2(xi,yi,zi,x_primai,y_primai,z_primai,gamma);
    double Kz1=funcion3(xi,yi,zi,x_primai,y_primai,z_primai,gamma);
 //Primer Paso
    double x1=xi+h*x_primai/2.0;
    double y1=xi+h*y_primai/2.0;
    double z1=xi+h*z_primai/2.0;
    double x_prima1=x_primai+h*Kx1/2.0;
    double y_prima1=y_primai+h*Ky1/2.0;
    double z_prima1=z_primai+h*Kz1/2.0;

    double Kx2=funcion1(x1,y1,z1,x_prima1,y_prima1,z_prima1,gamma);
    double Ky2=funcion2(x1,y1,z1,x_prima1,y_prima1,z_prima1,gamma);
    double Kz2=funcion3(x1,y1,z1,x_prima1,y_prima1,z_prima1,gamma);

 //Segundo Paso
    double x2=xi+h*x_prima1/2.0;
    double y2=xi+h*y_prima1/2.0;
    double z2=xi+h*z_prima1/2.0;
    double x_prima2=x_primai+h*Kx2/2.0;
    double y_prima2=y_primai+h*Ky2/2.0;
    double z_prima2=z_primai+h*Kz2/2.0;

    double Kx3=funcion1(x2,y2,z2,x_prima2,y_prima2,z_prima2,gamma);
    double Ky3=funcion2(x2,y2,z2,x_prima2,y_prima2,z_prima2,gamma);
    double Kz3=funcion3(x2,y2,z2,x_prima2,y_prima2,z_prima2,gamma);

//Tercer Paso
    double x3=xi+h*x_prima2/2.0;
    double y3=xi+h*y_prima2/2.0;
    double z3=xi+h*z_prima2/2.0;
    double x_prima3=x_primai+h*Kx3/2.0;
    double y_prima3=y_primai+h*Ky3/2.0;
    double z_prima3=z_primai+h*Kz3/2.0;

    double Kx4=funcion1(x3,y3,z3,x_prima3,y_prima3,z_prima3,gamma);
    double Ky4=funcion2(x3,y3,z3,x_prima3,y_prima3,z_prima3,gamma);
    double Kz4=funcion3(x3,y3,z3,x_prima3,y_prima3,z_prima3,gamma);

//Cuarto Paso
    
    double averange_Kx=(1.0/6.0)*(Kx1+2.0*Kx2+2.0*Kx3+Kx4);
    double averange_Ky=(1.0/6.0)*(Ky1+2.0*Ky2+2.0*Ky3+Ky4);
    double averange_Kz=(1.0/6.0)*(Kz1+2.0*Kz2+2.0*Kz3+Kz4);
    
    double x_primaf=x_primai+h*averange_Kx;
    double y_primaf=y_primai+h*averange_Ky;
    double z_primaf=z_primai+h*averange_Kz;
    double xf=xi+h*(x_primaf+x_primai)/2.0;
    double yf=yi+h*(x_primaf+x_primai)/2.0;
    double zf=zi+h*(x_primaf+x_primai)/2.0;

    //Restablecemos valores originales a los actualizados despues del Rongue Kutta, y los imprimimos en el archivo
    x_primai=x_primaf;
    y_primai=y_primaf;
    z_primai=z_primaf;
    xi=xf;
    yi=yf;
    zi=zf;
    double t=i*e5;

    fprintf(file, "%f %f %f %f \n",t,xf,yf,zf);
  }
  
  fclose(file);
  return 0;
}

double funcion1(double x, double y, double z, double x_prima, double y_prima, double z_prima, double gama){

  double r=pow(x,2)+pow(y,2)+pow(z,2);
  double func1=q*(-Bo*(pow(Rt,3)))*(y_prima*(2*z*(z-x)*(x-y)*y)-z_prima*(3*y*z))/(m*gama*pow(r,5.0/2.0));
  return func1; 

}

double funcion2(double x, double y, double z, double x_prima, double y_prima, double z_prima, double gama){

  double r=pow(x,2)+pow(y,2)+pow(z,2);
  double func2=q*(-Bo*(pow(Rt,3)))*(z_prima*(3*x*z)-x_prima*(2*z*(z-x)*(x-y)*y))/(m*gama*pow(r,5.0/2.0));
  return func2;
}

double funcion3(double x, double y, double z, double x_prima, double y_prima, double z_prima, double gama){

  double r=pow(x,2)+pow(y,2)+pow(z,2);
  double func3=q*(-Bo*(pow(Rt,3)))*(x_prima*(3*y*z)-y_prima*(3*x*z))/(m*gama*pow(r,5.0/2.0));
  return func3;
}





