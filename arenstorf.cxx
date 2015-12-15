#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

void f(double* const k, double* const X, const double mu, const double r, const double s);
void k(double* const X, double* const k1, double* const k2, double* const k3, double* const k4, double* const k5, double* const k6, double* const k7, double mu, const double r, const double s, const double dt);
double numerror(const double* const X, const double* const X2);
void Kutta54(double* const X, double * const X2, const double mu, const double dt);

int main(){
 const double mu=0.012277471,eps=1e-10,T_end=17.065216560157,safetyq=0.6, x0=0.994,y0=0,x0dot=0,y0dot=-2.00158510637908; // Initial values
 double dt=1e-3,t=0;
 double X5[4],X4[4],X_anf[4];
  X5[0]=x0;X5[1]=x0dot;X5[2]=y0;X5[3]=y0dot;X4[0]=x0;X4[1]=x0dot;X4[2]=y0;X4[3]=y0dot;  //Arrays

   ofstream out("Data.txt");

  while(t<T_end){
    for(int i=0;i<4;i++) X_anf[i]=X4[i];
   Kutta54(X5,X4,mu,dt);
    if(numerror(X5,X4)>eps){					// Stepsize control. Choose new Stepsize untill error is low enough
     while(numerror(X5,X4)>eps){
      dt*=pow(eps/numerror(X5,X4),0.2)*safetyq;
       for(int i=0;i<4;i++){ X4[i]=X_anf[i]; X5[i]=X_anf[i];}	// return to initial value to calculate again and check for error
      Kutta54(X5,X4,mu,dt);
     }
    } else {							// if error is already low enough, increase stepsize again
      dt*=pow((eps/numerror(X5,X4)),1/5.0)*safetyq;
       for(int i=0;i<4;i++){ X4[i]=X_anf[i]; X5[i]=X_anf[i];}	// return to initial value
      Kutta54(X5,X4,mu,dt);					// calculate again with increased stepsize
      }                    
    for(int i=0;i<4;i++) X5[i]=X4[i];				// continue with the 4th order solution
   out << t << " " << dt << " " << X4[0] << " " << X4[2] << endl;
   t+=dt;
  }

out.close();
return 0;
}

void Kutta54(double* const X, double * const X2, const double mu, const double dt){
 double r=sqrt(pow((X[0]+mu),2)+pow(X[2],2)),s=sqrt(pow((X[0]+mu-1.0),2)+pow(X[2],2)); //,r2=sqrt(pow(X2[0]+mu,2)+pow(X2[2],2)),s2=sqrt(pow(X2[0]+mu-1,2)+pow(X2[2],2));
 double k1[4],k2[4],k3[4],k4[4],k5[4],k6[4],k7[4];
  k(X,k1,k2,k3,k4,k5,k6,k7,mu,r,s,dt);
  //if(!X){ for(int i=0;i<4;i++) X2[i]+=dt*(5179/57600.0*k1[i]+7571/16695.0*k3[i]+393/640.0*k4[i]-92097339200.0*k5[i]+187/2100.0*k6[i]+1/40.0*k7[i]);}
  //else{ 
   for(int i=0;i<4;i++){
     X[i]+=dt*((35.0/384.0)*k1[i]+(500.0/1113.0)*k3[i]+(125.0/192.0)*k4[i]+(-2187.0/6784.0)*k5[i]+(11.0/84.0)*k6[i]);
     X2[i]+=dt*((5179.0/57600.0)*k1[i]+(7571.0/16695.0)*k3[i]+(393.0/640.0)*k4[i]+(-92097.0/339200.0)*k5[i]+(187.0/2100.0)*k6[i]+(1.0/40.0)*k7[i]);
    }
   //}
}

void f(double* const k, double* const X, const double mu, const double r, const double s){

  k[0]=X[1];
  k[2]=X[3];
  k[1]=X[0]+2*X[3]-(1.0-mu)*(X[0]+mu)/pow(r,3)-mu*(X[0]-1.0+mu)/pow(s,3);
  k[3]=X[2]-2*X[1]-(1.0-mu)*X[2]/pow(r,3)-mu*X[2]/pow(s,3);
}

void k(double* const X, double* const k1, double* const k2, double* const k3, double* const k4, double* const k5, double* const k6, double* const k7, double mu, const double r, const double s, const double dt){

 double a[4];
  f(k1,X,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt/5.0*k1[i];
  f(k2,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+(3.0*dt/40.0)*(k1[i]+3*k2[i]);
  f(k3,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*((44.0/45.0)*k1[i]+(-56.0/15.0)*k2[i]+(32.0/9.0)*k3[i]);
  f(k4,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*((19372.0/6561.0)*k1[i]+(-25360.0/2187.0)*k2[i]+(64448.0/6561.0)*k3[i]+(-212.0/729.0)*k4[i]);
  f(k5,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*((9017.0/3168.0)*k1[i]+(-355.0/33.0)*k2[i]+(46732.0/5247.0)*k3[i]+(49.0/176.0)*k4[i]+(-5103.0/18656.0)*k5[i]);
  f(k6,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*((35.0/384.0)*k1[i]+(500.0/1113.0)*k3[i]+(125.0/192.0)*k4[i]+(-2187.0/6784.0)*k5[i]+(11.0/84.0)*k6[i]);
  f(k7,a,mu,r,s);

}

double numerror(const double* const X, const double* const X2){
 double a=abs(X[0]-X2[0]);
  for(int i=1;i<4;i++){
   if(abs(X[i]-X2[i])>=a) a=abs(X[i]-X2[i]);
  }
  return a;
}
