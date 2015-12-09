// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx,double* k1, double* k2, double* k3, double* k4);
void poly(const double ytemp,double* k1, double* k2, double* k3, double* k4, const double dx, double theta, double& ytheta);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
	ofstream out2("blabla");
  const int dim = 2;
  
	double dx = 0.1,x=0;
	const double L = 10000;
	double p0 = 3.0;
	
	for (p0 = 0.1; p0 < 5; p0 += 0.1){
	   x = 0;
	    double y0[dim] = {p0 , 0.0};
	    double yn[dim];
	    double ytemp = 0.0;
	    double k1[dim], k2[dim], k3[dim], k4[dim];

      out2 << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	    while(x<=L)
	    {	
		    x += dx;
		    RKstep(yn, y0, x, dx,k1,k2,k3,k4);
		    
		    
		    if(ytemp>0.0 && y0[1]<0.0) break;
		    
		    
		    
		    ytemp=y0[1];
		    
	    for(int i=0; i<dim; i++) y0[i] = yn[i];
			out2 << x << "\t" << y0[0] << "\t" << y0[1] << endl;
		}
	// 	cout << t1 << "\t" << t2 << endl;
		
	    double theta = 0.5;
	    double thL = 0;
	    double thR = 1;
	    double ytheta = ytemp;
	    
	    while (abs(ytheta) > 1e-8){
		poly(ytemp, k1, k2,k3, k4, dx, theta, ytheta);
		if (ytheta > 0)
		  thL = theta;
		if (ytheta <= 0)
		  thR = theta;
		theta = (thL+thR)/2;
	// 	cout << theta << endl;
	    }
	out << p0 << "\t" << x + theta*dx << endl;
	}
	out.close();
	out2.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx,double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;
	

  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
  f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// Lorenz model
void f(double* const y0, const double x)
{
	double y[2] = { y0[0], y0[1] };

	y0[0] = y[1];
	y0[1] = -(y[0] / (sqrt(1+y[0]*y[0])));
	
}

void poly(const double ytemp,double* k1, double* k2, double* k3, double* k4, const double dx, double theta, double& ytheta){
  double b[4];
  b[0] = theta - (3*theta*theta)/2.0 + (2.0*pow(theta,3))/3.0;
  b[1] = b[2] = theta*theta - (2.0*pow(theta,3))/3.0;
  b[3] = -(theta*theta)/2.0 + (2.0*pow(theta,3))/3.0;
  
  ytheta = ytemp + dx*(b[0]*k1[1] + b[1]*k2[1] + b[2]*k3[1] + b[3]*k4[1]);
  
}
