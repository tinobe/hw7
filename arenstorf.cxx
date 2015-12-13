#include <fstream>
#include <cmath>

using namespace std;

void f(double* const k, double* const X, double mu, const double r, const double s);
void k(double* const X, double* const k1, double* const k2, double* const k3, double* const k4, double* const k5, double* const k6, double* const k7, double mu, const double r, const double s, const double dt);

int main(){
 double dt=0.001, x0=0.994,y0=0,x0dot=0,y0dot=-2.00158510637908,eps,mu=0.012277471,tol=1e-10,T_end=50,t=0,r,s;
 double X[4],X2[4],k1[4],k2[4],k3[4],k4[4],k5[4],k6[4],k7[4];
  X[0]=x0;X[1]=x0dot;X[2]=y0;X[3]=y0dot;X2[0]=x0;X2[1]=x0dot;X2[2]=y0;X2[3]=y0dot;

   ofstream out("Data.txt");

  while(t<T_end){
   r=sqrt(pow(X[0]+mu,2)+pow(X[2],2));
   s=sqrt(pow(X[0]+mu-1,2)+pow(X[2],2));
   k(X,k1,k2,k3,k4,k5,k6,k7,mu,r,s,dt);
    for(int i=0;i<4;i++){
     X[i]+=dt*(35/384.0*k1[i]+500/1113.0*k3[i]+125/192.0*k4[i]-2187/6784.0*k5[i]+11/84.0*k6[i]);
     X2[i]+=dt*(5179/57600.0*k1[i]+7571/16695.0*k3[i]+393/640.0*k4[i]-92097339200.0*k5[i]+187/2100.0*k6[i]+1/40.0*k7[i]);
    }
   out << t << " " << X[0] << " " << X[2] << endl;
   t+=dt;
  }

out.close();
return 0;
}

void f(double* const k, double* const X, double mu, const double r, const double s){

  k[0]=X[1];
  k[2]=X[3];
  k[1]=X[0]+2*X[3]-(1-mu)*(1+mu)/pow(r,3)-mu*(X[0]-1+mu)/pow(s,3);
  k[3]=X[2]-2*X[1]-(1-mu)*X[2]/pow(r,3)-mu*X[2]/pow(s,3);
}

void k(double* const X, double* const k1, double* const k2, double* const k3, double* const k4, double* const k5, double* const k6, double* const k7, double mu, const double r, const double s, const double dt){

 double a[4];
  f(k1,X,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt/5.0*k1[i];
  f(k2,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*(3/40.0*k1[i]+9/40.0*k2[i]);
  f(k3,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*(44/45.0*k1[i]-56/15.0*k2[i]+32/9.0*k3[i]);
  f(k4,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*(19372/6561.0*k1[i]-25360/2187.0*k2[i]+64448/6561.0*k3[i]-212/729.0*k4[i]);
  f(k5,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*(9017/3168.0*k1[i]-355/33.0*k2[i]+46732/5247.0*k3[i]+49/176.0*k4[i]-5103/18656.0*k5[i]);
  f(k6,a,mu,r,s);

  for(int i=0;i<4;i++) a[i]=X[i]+dt*(35/384.0*k1[i]+500/1113.0*k3[i]+125/192.0*k4[i]-2187/6784.0*k5[i]+11/84.0*k6[i]);
  f(k7,a,mu,r,s);

}
