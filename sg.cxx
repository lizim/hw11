#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);


void step(cmplx* const psi1, cmplx* const psi0, const double hq, 
	  const double dt, const double m, const double dx, 
	  const double xmin, const double omega, const int Nx);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
	const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt = 0.1*dx;
	double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double m = 1;
	const double hq = 1;
	const double lambda = 10;
	const double omega = 0.2;
	const double alpha = sqrt((m*omega)/sqrt(hq));

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0,alpha,lambda,dx,dt,Nx,xmin);

	writeToFile(psi0,"psi_0",dx,Nx,xmin,alpha,lambda,omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
		    
		    step(psi1,psi0,hq,dt,m,dx,xmin,omega,Nx);
		    
		    h    = psi0;
		    psi0 = psi1;
		    psi1 = h;
		    
		    t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
	cout << "t = " << t << endl;

	delete[] psi0;
	delete[] psi1;
	
	return 0;
}
//-----------------------------------
void step(cmplx* const psi1, cmplx* const psi0, const double hq, 
	  const double dt, const double m, const double dx, 
	  const double xmin, const double omega, const int Nx)
{
  double x;
  double V;
  cmplx* H  = new cmplx[Nx]; // Zwischenberechnung von A*PsiN
  cmplx* d  = new cmplx[Nx];
  cmplx* dk = new cmplx[Nx]; // Komplex konjugiertes
  cmplx a  = cmplx(0.0,-(hq*dt)/(4*m*dx*dx));
  cmplx ak = cmplx(0.0,(hq*dt)/(4*m*dx*dx));
  
  for(int i=0;i<Nx;i++){
	x = xmin + i * dx;
	V = 0.5*m * pow(omega*x,2);
	d[i]  = cmplx(1.0,(hq*dt)/2*m*dx*dx) + dt/(2*hq)*V;
	dk[i] = cmplx(1.0,-(hq*dt)/2*m*dx*dx) - dt/(2*hq)*V;
  }
  
  H[0] = dk[0]*psi0[0] + ak*psi0[1];
  for(int i=1;i<(Nx-1);i++){
    H[i] = dk[i]*psi0[i] + ak*(psi0[i-1]+psi0[i+1]);
  }
  H[Nx-1] = ak*psi0[Nx-2] + dk[Nx-1]*psi0[Nx-1];
  
  // forward substitution
  
  for(int i=1;i<Nx;i++){
  
    d[i] -= a*a/d[i-1];
    H[i] -= H[i-1]*a/d[i-1];
  }
    
  // backward substitution
  
  psi1[Nx-1] = H[Nx-1]/d[Nx-1];
  
  for(int i=(Nx-2);i>=0;i--){
    
    psi1[i] = (H[i] - a*psi1[i+1]) / d[i];
  }
  
  delete[] H;
  delete[] d;
  delete[] dk;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
		xi = alpha * x;
		xil = alpha * lambda;
		h1 = -0.5 * pow(xi - xil*cos(omega*t),2);
		h2 = 0.5 * omega*t + xi * xil * sin(omega*t);
		h3 = -0.25 * xil*xil* sin(2*omega*t);
		ana = cmplx(h1,h2+h3);
		ana = sqrt(alpha/sqrt(M_PI)) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
		    << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	// const double x0 = dx*Nx * 0.5; // wofuer ist das?
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(-0.5 * pow(alpha*(x-lambda),2));
	}
}
