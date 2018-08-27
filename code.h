#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
using namespace std;
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

class FSSH
{
	public:
		//Basic Variables
		double x; //distance from origin
		double v; //velocity 
		int state;

		//boundary variables
		double xmax=10.001, xmin=-10.001; //boundaries of x
		
		//Potential Energy, Wavefunction and Energies of States, dVij = derivative of Vij 
		Matrix Vij, WF, E, dVij, phi1, phi2, D12;
	FSSH();	
	void Build(double x);
};
