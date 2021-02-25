#include "TMath.h"

void evaluate_min_dist_to_jet(Jet_struct jt, TVector3 first, TVector3 last, TVector3 *PCA,  double *min_dist){

	long int N_ITER=100;
	TVector3 dir_trk,dir_jet,A,B;
	dir_trk=(last-first).Unit();
	dir_jet=jt.j.Vect().Unit();
	/*double t_min=0;
	double t_max=(last.X()-first.X())/dir_trk.X(); //N.B. può essere anche negativo*/
	double dist_A, dist_B;
	
	A=first;
	B=last;
	dist_A = A.Cross(dir_jet).Mag();
	dist_B = B.Cross(dir_jet).Mag();
	
	for(int n=0;n<N_ITER;n++){
		
		/*cout<<"Interazione "<<n<<"-esima"<<endl;
		cout<<"dist_A: "<<dist_A<<endl;
		cout<<"dist_B: "<<dist_B<<endl;*/
		
		if(dist_A<=dist_B){
			B.SetXYZ((A.X()+B.X())/2.,(A.Y()+B.Y())/2.,(A.Z()+B.Z())/2.);
			dist_B = B.Cross(dir_jet).Mag();
			}
		else if(dist_A>dist_B){
			A.SetXYZ((A.X()+B.X())/2.,(A.Y()+B.Y())/2.,(A.Z()+B.Z())/2.);
			dist_A = A.Cross(dir_jet).Mag();
			}
		}
	if(dist_A<=dist_B){
		*min_dist=dist_A;
		*PCA=A;
		}
	else{
		*min_dist=dist_B;
		*PCA=B;
		}
	
	}
void evaluate_min_dist_to_jet2(Jet_struct jt, TVector3 first, TVector3 last, TVector3 *PCA,  double *min_dist){
	TVector3 dir_trk,dir_jet;
	dir_trk=(last-first).Unit();
	dir_jet=jt.j.Vect().Unit();
	TMatrixD coeff(3,3);
	double data[9]={first.X(),first.Y(),first.Z(),dir_trk.X(), dir_trk.Y(),dir_trk.Z(),dir_jet.X(),dir_jet.Y(),dir_jet.Z()};
	coeff.SetMatrixArray(data,"T");
	if(abs(coeff.Determinant())<1e-10){
		cout<<"determinante: "<<coeff.Determinant()<<endl;
	}
	TMatrixD A(2,2);
	TMatrixD C(2,1);
	TMatrixD Ainv(2,2);
	TMatrixD Y(2,1);
	
	double data_A[4];
	data_A[0]=dir_jet.Dot(dir_trk);
	data_A[1]=-dir_trk.Mag2();
	data_A[2]=dir_jet.Mag2();
	data_A[3]=-dir_jet.Dot(dir_trk);
	A.SetMatrixArray(data_A,"T");
	
	/*cout<<"00: "<<data_A[0]<<" -> "<< TMatrixDRow(A,0)(0)<<endl;
	cout<<"01: "<<data_A[1]<<" -> "<< TMatrixDRow(A,0)(1)<<endl;
	cout<<"10: "<<data_A[2]<<" -> "<< TMatrixDRow(A,1)(0)<<endl;
	cout<<"11: "<<data_A[3]<<" -> "<< TMatrixDRow(A,1)(1)<<endl;*/
	double data_C[2];
	data_C[0]=first.Dot(dir_trk);
	data_C[1]=first.Dot(dir_jet);
	C.SetMatrixArray(data_C,"T");
	
	Ainv=A.Invert();
	Y=Ainv*C;
	double u = TMatrixDRow(Y,0)(0);
	double t = TMatrixDRow(Y,1)(0);
	
	/*cout<<"u: "<<u<<endl;
	cout<<"t: "<<t<<endl;*/
	
	TVector3 PCA_jet,PCA_trk;
	PCA_jet.SetXYZ(dir_jet.X()*u,dir_jet.Y()*u,dir_jet.Z()*u);
	PCA_trk.SetXYZ(dir_trk.X()*t+first.X(),dir_trk.Y()*t+first.Y(),dir_trk.Z()*t+first.Z());
	
	*min_dist=(PCA_jet-PCA_trk).Mag();
	*PCA=PCA_trk;
	if(*min_dist>50.){
	cout<<"min dist: "<<*min_dist<<endl;}
	
	}
