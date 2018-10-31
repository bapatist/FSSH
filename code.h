#include <ctime>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixC;

class FSSH
{
	public:
		//Basic Variables
		double x; //distance from origin
		double K; //momentum 
		double F[2]; //Force
		int state;
		double t;
		bool hopped;
		const double dt=1;
		const double M=2000; //(mass of molecule in atomic units)i
		const double hbar=1.0; 

		//boundary variables
		double xmax=10.01, xmin=-10.01; //boundaries of x
		
		double iK;
		int Right2=0;
		int Right1=0;
		int Left2=0;
		int Left1=0;
		
		double PR1, PR2, PL1, PL2;

		//Potential Energy, Wavefunction and Energies of States, dVij = derivative of Vij 
		Matrix Vij, dVij, Dij;
		MatrixC C, Cnew, A, B;
		double initE, finalE, force, Knew, xnew;
	       	int statenew;	
		double D12;
		double TE, E[2], dE[2]; //Total Energy; Energy Array for 2 states; derivative of energy wrt x
	FSSH();	
	void Build();
	void StateE();
	void Velocity();
	void CouplingD();
	void RK4();
	void Position();
	void Initializer();
	void Hopping();
};
