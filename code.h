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

typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> MatrixC;

class FSSH
{
	public:
		//Basic Variables
		double x; //distance from origin
		double K; //momentum 
		double F[2]; //Force
		int state;
		complex<double> c1, c2;
		double t;
		double dt=1;
		double M=2000; //(mass of molecule in atomic units)

		//boundary variables
		double xmax=10.1, xmin=-10.1; //boundaries of x
		
		//Potential Energy, Wavefunction and Energies of States, dVij = derivative of Vij 
		Matrix Vij, dVij;
		MatrixC C;
		double initE, finalE, force, Knew, xnew; 
		double D12;
		double TE, E[2], dE[2]; //Total Energy; Energy Array for 2 states; derivative of energy wrt x
	FSSH();	
	void Build();
	void StateE();
	void Velocity();
	double CouplingD();
	void RK4();
	void Position();
	void Initializer();
};
