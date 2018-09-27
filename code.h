#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
using namespace std;
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

class FSSH
{
	public:
		//Basic Variables
		double x; //distance from origin
		double K; //momentum 
		int state;
		double dt=0.1;

		//boundary variables
		double xmax=10.0, xmin=-10.0; //boundaries of x
		
		//Potential Energy, Wavefunction and Energies of States, dVij = derivative of Vij 
		Matrix Vij, dVij;
		double initE, finalE, force, initK, finalK, initX, finalX; 
		double D12;
		double E[2];
	FSSH();	
	void Build();
	void StateE();
	double Velocity();
	double CouplingD();
	double Position();
};
