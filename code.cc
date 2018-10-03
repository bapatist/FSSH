#include "code.h"

void print_mat(const Matrix& a);

FSSH::FSSH(){
       	x=-10.0;
	state = 0;
	K = 1;
	while(x>xmin and x<xmax){
	
		//for (double n=-10; n<10; n++){
			
			Build();
			StateE();
			Position();
			cout << "Position is at " << x << endl;
			cout << "K value is " << K << '\n' << endl;
			//cout << "Energy of States at x= " << x << endl << E[0] << "\t" << E[1] << endl;
			//cout << "Derivative Energy of States at x= " << x << endl << dE[0] << "\t" << dE[1] << endl;

		}	
	//CouplingD();
	//Velocity(); //calculates initial and final E and changes momentum 
	//Position();
	//cout << "Energy of state is= " << E(state) << "value of x is = " << x << endl; 

	
}

void FSSH::Build(){
	//Creating PE Matrix from Tully's Paper for Simple Avoided Crossing
	Vij = Matrix(2,2);
	if (x>0) Vij(0,0) = 0.01*(1-exp(-(1.6*x)));
	else     Vij(0,0) = -0.01*(1-exp(1.6*x));
	Vij(1,1) = -Vij(0,0);
	Vij(1,0) = Vij(0,1) = 0.005*exp(-(pow(x,2)));

	dVij = Matrix(2,2);
	if(x>0) dVij(0,0) = 0.016*exp(-1.6*x);
	else	dVij(0,0) = 0.016*exp(1.6*x);
	dVij(1,1) = -dVij(0,0);
	dVij(1,0) = dVij(0,1) = -0.01*x*exp(-pow(x,2));

}

double FSSH::CouplingD(){

	//Failed attempt to calculate coupling vector 
	
////////Eigen::SelfAdjointEigenSolver<Matrix> solver(Vij);
////////WF = Matrix(2,2);
////////WF = solver.eigenvectors();
////////E = Matrix(2,1);
////////E = solver.eigenvalues();
////////
////////phi1 = Matrix(2,1);
////////phi2 = Matrix(2,1);
////////for(int i=0; i<2; i++){
////////       phi1(i,0) = WF(i,0);
////////       phi2(i,0) = WF(i,1);
////////}
//////////print_mat(WF);
//////////print_mat(phi1);
//////////print_mat(phi2);
	//print_mat(dVij);
//	D12 = Matrix(1,1);
//	double deltaE = E(0) - E(1);
//	D12 = -1 * phi1.transpose() * dVij * phi2 / deltaE;
	//cout << deltaE << endl;
	//cout << endl << "For x=" << x << "  D12 is = " << D12 << endl;

	//CHECKING WITH PRASHANT'S EXPRESSION

	double n = (Vij(0,0) - Vij(1,1)) / 2*Vij(0,1);
	double D12 = (1/(4*(Vij(0,1)*Vij(0,1))*(1+n*n))) * ((Vij(0,1)*(dVij(0,0)-dVij(1,1))) - (dVij(0,1)*(Vij(0,0)-Vij(1,1))));
	return D12; 
}

void FSSH::StateE(){
	E[1] = ((Vij(0,0) + Vij(1,1)) + pow( pow(Vij(0,0)+Vij(1,1),2) - 4*(Vij(0,0)*Vij(1,1)- pow(Vij(0,1), 2)) , 0.5))/2;

	E[0] = ((Vij(0,0) + Vij(1,1)) - pow( pow(Vij(0,0)+Vij(1,1),2) - 4*(Vij(0,0)*Vij(1,1)- pow(Vij(0,1), 2)) , 0.5))/2;

	dE[1] = 0.5*(dVij(0,0) + dVij(1,1) + 0.5*(pow(pow(Vij(0,0)+Vij(1,1), 2) - 4*(Vij(0,0)*Vij(1,1)-pow(Vij(0,1),2)), -0.5))*( 2*(Vij(0,0)+Vij(1,1))*(dVij(0,0) + dVij(1,1)) - 4*(Vij(0,0)*dVij(1,1) + Vij(1,1)*dVij(0,0)) + 8*Vij(0,1)*dVij(0,1)));
	 
	dE[0] = 0.5*(dVij(0,0) + dVij(1,1) - 0.5*(pow(pow(Vij(0,0)+Vij(1,1), 2) - 4*(Vij(0,0)*Vij(1,1)-pow(Vij(0,1),2)), -0.5))*( 2*(Vij(0,0)+Vij(1,1))*(dVij(0,0) + dVij(1,1)) - 4*(Vij(0,0)*dVij(1,1) + Vij(1,1)*dVij(0,0)) + 8*Vij(0,1)*dVij(0,1)));
}

void FSSH::Velocity(){
	K = K + dE[0]*dt;
}
void FSSH::Position(){
	x = x + K*dt - 0.5*dE[0]*pow(dt,2);
}

//print matrix function
void print_mat(const Matrix& a)
{
      FILE* out = stdout;
      int m = a.rows();
      int n = a.cols();
      int ii,jj,kk,nn,ll;
      int i,j;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=n;
      if (nn > kk) nn=kk;
      ll = 2*(nn-ii+1)+1;
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < m; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",a(i,j));
            }
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
}
	
	
	
	
	
