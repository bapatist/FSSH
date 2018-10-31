#include "code.h"

void print_mat(const Matrix& a);

FSSH::FSSH(){
	for(iK=10.0; iK<31.0; iK=iK+5){
		Right1 = Right2 = Left1 = Left2 = 0;
		for(int traj=0; traj<2000; traj++){
			Initializer();
			while(x>xmin and x<xmax){
				Build();     // Calculates Vij and dVij matrices
				StateE();    // Calculates State energies and total energy using Vij
				CouplingD(); // Calculates coupling vector Dij
//				std::ofstream logfile(
//					"logfile.txt", std::ios_base::out | std::ios_base::app
//				);
//				std::ofstream TotalE(	
//					"TotalE.txt", std::ios_base::out | std::ios_base::app
//				);
//				std::ofstream State(
//					"State.txt", std::ios_base::out | std::ios_base::app
//				);
//				logfile << std::setprecision(4) << std::fixed << t << '\t' << x << endl;
//				TotalE << std::setprecision(4) << std::fixed << t << '\t' << TE << endl;
//				State << std::setprecision(4) << std::fixed << x << '\t' << state << endl;
							
				Hopping();  //Updates State
				RK4();      //Updates C vector
				Position(); //Updates x
				Velocity(); //Updates K
				C = Cnew;
				x = xnew;
				K = Knew;
				t=t+dt;
				}
			if(x>0 and state==0) Right1 += 1;
			if(x<0 and state==0) Left1 += 1;
			if(x>0 and state==1) Right2 += 1;
			if(x<0 and state==1) Left2 += 1;	
		}
	
		PR1 = (double)Right1/2000.0;
		PR2 = (double)Right2/2000.0;
		PL1 = (double)Left1/2000.0;
		cout << "For K=" << iK << endl;
		cout << Right1 << endl << Right2 << endl << Left1 << endl << Left2 << endl;		
		cout << PR1 << endl;
		cout << PR2 << endl;
		cout << PL1 << endl << endl;
	}

//		std::ofstream P1(
//			"P1.txt", std::ios_base::out | std::ios_base::app
//		);
//		P1 << std::setprecision(4) << std::fixed << K << '\t' << PR1 << endl;
//		std::ofstream P2(
//			"P2.txt", std::ios_base::out | std::ios_base::app
//		);
//		P2 << std::setprecision(4) << std::fixed << K << '\t' << PR2 << endl;
	//}	
}
void FSSH::Initializer(){
	x=-10.0;
	state = 0;
	K = iK;
	t=0.0;
	C = MatrixC (2,1);
	C(0) = {1,0};
	C(1) = {0,0};
	srand(time(0));
}
void FSSH::Hopping(){
	if(2*M*(E[state]-E[(state+1)%2]) + K*K > 0){
		MatrixC A(2,2), B(2,2);
		for(int l=0; l<2; l++){
			for(int m=0; m<2; m++){
				A(l,m) = C(l)*conj(C(m));
				B(l,m) = (2/hbar)*(imag(conj(A(l,m))*Vij(l,m))) - 2*(real(conj(A(l,m))*x*Dij(l,m)));
				//cout << A(l,m) << endl;
			}
		}
		double rando = rand()/(double)RAND_MAX;
		hopped = 0;
		if(state==0){
		//	cout << dt*B(1,0)/A(0,0) << endl;
			if(real(dt*B(1,0)/A(0,0)) > rando){
				statenew = 1;
				hopped = 1;
			}
		}
		if(state==1){
			if(real(dt*B(0,1)/A(1,1)) > rando){
				statenew = 0;
				hopped = 1;
			}
		}
		if(hopped){
			K = pow((2*M*(E[state]-E[statenew]) + K*K),0.5);
			state = statenew;
		}
	}
}


void FSSH::RK4(){
	MatrixC A(2,2), B(2,2);
	MatrixC k1(2,1), k2(2,1), k3(2,1), k4(2,1);
	MatrixC Mmat(2,2);
	complex<double> iota = {0,1};
	Mmat = -((iota*(Vij.cast<std::complex<double>>()))/hbar) - ((K/M)*(Dij.cast<std::complex<double>>()));
	k1 = dt*Mmat*C;
	k2 = dt*Mmat*(C+0.5*k1);
	k3 = dt*Mmat*(C+0.5*k2);
	k4 = dt*Mmat*(C+k3);
	Cnew = C + (k1 + 2*k2 + 2*k3 + k4)/6;	
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

void FSSH::CouplingD(){
	double n = (Vij(0,0) - Vij(1,1)) / (2*Vij(0,1));
	double D12 = (1/(4*(Vij(0,1)*Vij(0,1))*(1+n*n))) * ((Vij(0,1)*(dVij(0,0)-dVij(1,1))) - (dVij(0,1)*(Vij(0,0)-Vij(1,1))));

	Dij = Matrix(2,2);
	Dij(0,0) = 0; 
	Dij(0,1) = D12; 
	Dij(1,0) = -D12; 
	Dij(1,1) = 0;
}

void FSSH::StateE(){
	E[1] = ((Vij(0,0) + Vij(1,1)) + pow( pow(Vij(0,0)+Vij(1,1),2) - 4*(Vij(0,0)*Vij(1,1)- pow(Vij(0,1), 2)) , 0.5))/2;

	E[0] = ((Vij(0,0) + Vij(1,1)) - pow( pow(Vij(0,0)+Vij(1,1),2) - 4*(Vij(0,0)*Vij(1,1)- pow(Vij(0,1), 2)) , 0.5))/2;

	dE[1] = 0.5*(dVij(0,0) + dVij(1,1) + 0.5*(pow(pow(Vij(0,0)+Vij(1,1), 2) - 4*(Vij(0,0)*Vij(1,1)-pow(Vij(0,1),2)), -0.5))*( 2*(Vij(0,0)+Vij(1,1))*(dVij(0,0) + dVij(1,1)) - 4*(Vij(0,0)*dVij(1,1) + Vij(1,1)*dVij(0,0)) + 8*Vij(0,1)*dVij(0,1)));
	 
	dE[0] = 0.5*(dVij(0,0) + dVij(1,1) - 0.5*(pow(pow(Vij(0,0)+Vij(1,1), 2) - 4*(Vij(0,0)*Vij(1,1)-pow(Vij(0,1),2)), -0.5))*( 2*(Vij(0,0)+Vij(1,1))*(dVij(0,0) + dVij(1,1)) - 4*(Vij(0,0)*dVij(1,1) + Vij(1,1)*dVij(0,0)) + 8*Vij(0,1)*dVij(0,1)));

	F[0] = -dE[0];
	F[1] = -dE[1];
	
	TE = K*K/(2*M) + E[state];
}

void FSSH::Velocity(){
	Knew = K + F[state]*dt;
}
void FSSH::Position(){
	xnew = x + (K*dt + 0.5*F[state]*pow(dt,2))/M;
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
