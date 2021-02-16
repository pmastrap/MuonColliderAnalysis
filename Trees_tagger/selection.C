#include "mystruct.h"
#include "TMath.h"

bool trk_selection(Track_struct trk){
	
	//cout<<"trk.t.Pt: "<<trk.t.Pt()<<endl;
	if(trk.t.Pt()>1.){
	
		if(trk.chi2/trk.ndf<5.){
			if(trk.sub_det_hits[0]>0 || trk.sub_det_hits[2]>0){
				//cout<<"traccia selezionata"<<endl;
				return true;
				}
			}
		}
	return false;
		
	}
	
bool vertex_selection(Vertex_struct vtx){
	
	if(vtx.fd_sig_2 > 1.5){
		return true;
		}
	return false;
		
	}	
	
bool jet_selection(Track_struct* trk ,Jet_struct jt){

	int i=0;
	int counter=0;
	if(abs(jt.j.Eta())<2.5){
		if(jt.j.Pt()>10.){
			for(i=0;i<jt.num_trk;i++){
				if(trk[jt.trk_indices[i]].selected == true){
					counter++;
					}
				}
			if(counter>1){
				return true;
				}
			}
		}					
		
	return false;
	
	}
	
int find_first_sv(int dim, Jet_struct jt, Vertex_struct *vtx){
	int index[100];
	double vect[100];
	int k=0;
	if(dim==2){
		//sort sv by uncertainty on 2D flight distance
  		for(k=0;k<jt.num_SV;k++){
  			vect[k]=1000000.;
  			if(vtx[jt.SV_indices[k]].selected==true){
  				vect[k]=vtx[jt.SV_indices[k]].fd_2/vtx[jt.SV_indices[k]].fd_sig_2;
  				}
  			
  			}
  		TMath::Sort(jt.num_SV,vect,index,false);
  		}
  		
  	else if(dim==3){
		//sort sv by uncertainty on 3D flight distance
  		for(k=0;k<jt.num_SV;k++){
  			vect[k]=1000000.;
  			if(vtx[jt.SV_indices[k]].selected==true){
  				vect[k]=vtx[jt.SV_indices[k]].fd_3/vtx[jt.SV_indices[k]].fd_sig_3;
  				}
  			}
  		TMath::Sort(jt.num_SV,vect,index,false);
  		}
  	/*if(jt.num_SV>1){
  	cout<<"---------------------------------------------------------"<<endl;
  	cout<<"----------------QUI ----------------------QUI---------------------------"<<endl;
  	cout<<"indeces"<<endl;
  	for(k=0;k<jt.num_SV;k++){
  		cout<<"uncert "<< dim <<"d "<<vect[k]<<endl;
  		cout<<"index: "<<index[k]<<endl;
  		}
  	}*/
  	return jt.SV_indices[index[0]];
  	}

