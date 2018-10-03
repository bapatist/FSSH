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
		double dt=1;

		//boundary variables
		double xmax=10.1, xmin=-10.1; //boundaries of x
		
		//Potential Energy, Wavefunction and Energies of States, dVij = derivative of Vij 
		Matrix Vij, dVij;
		double initE, finalE, force, initK, finalK, initX, finalX; 
		double D12;
		double E[2], dE[2];
	FSSH();	
	void Build();
	void StateE();
	void Velocity();
	double CouplingD();
	void Position();
};
