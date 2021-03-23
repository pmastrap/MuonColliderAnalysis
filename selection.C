#include "mystruct.h"
#include "TMath.h"


bool trk_selection(Track_struct trk){
	
	//cout<<"trk.t.Pt: "<<trk.t.Pt()<<endl;
	if(trk.t.Pt()>1.){
	
		if(trk.chi2/trk.ndf<5.){
			if((trk.sub_det_hits[0]+trk.sub_det_hits[2])>4){
				//cout<<"traccia selezionata"<<endl;
				return true;
				}
			}
		}
	return false;
		
	}
	
bool lep_selection(Track_struct trk){
	
	//cout<<"trk.t.Pt: "<<trk.t.Pt()<<endl;
	if(trk.t.Pt()>2.){
		if((trk.sub_det_hits[0]+trk.sub_det_hits[2])>4){
			//cout<<"traccia selezionata"<<endl;
			return true;
			}
		}
	return false;
		
	}
	
bool vertex_selection(Vertex_struct vtx){
	if(vtx.id==104){
	if(vtx.fd_sig_2 > 1.5){
		return true;
		}}
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

void order_tracks(Jet_struct jt, Track_struct *trk, int *index){
	double vect[1000];
	int k=0;
	//sort trk by 2D impact-param significance (decreasing)
	for(k=0;k<jt.num_trk;k++){
	  	vect[k]=-9999;
	  	if(trk[jt.trk_indices[k]].selected==true){
	  		vect[k]=trk[jt.trk_indices[k]].sd0/trk[jt.trk_indices[k]].sd0_err;
	  		}		
	  	}
	 TMath::Sort(jt.num_trk,vect,index,true);
	  }
	
void find_first_trk(Jet_struct jt, Track_struct *trk,int *one,int *two, bool lepyn){
	int index[1000];
	double vect[1000];
	int k=0;
	if(lepyn==false){
		//sort trk by 2D impact-param significance (decreasing)
	  	for(k=0;k<jt.num_trk;k++){
	  		vect[k]=-9999;
	  		if(trk[jt.trk_indices[k]].selected==true){
	  			vect[k]=trk[jt.trk_indices[k]].sd0/trk[jt.trk_indices[k]].sd0_err;
	  			}
	  			
	  		}
	  	TMath::Sort(jt.num_trk,vect,index,true);
	  	*one=jt.trk_indices[index[0]];
	  	*two=jt.trk_indices[index[1]];
	  }
	if(lepyn==true){
		//sort trk by 2D impact-param significance (decreasing)
	  	for(k=0;k<jt.num_lep;k++){
	  		vect[k]=-9999;
	  		if(trk[jt.lep_indices[k]].selected==true){
	  			vect[k]=trk[jt.lep_indices[k]].sd0/trk[jt.lep_indices[k]].sd0_err;
	  			}
	  		else{cout<<"ops"<<endl;}
	  			
	  		}
	  	TMath::Sort(jt.num_lep,vect,index,true);
	  	*one=jt.lep_indices[index[0]];
	  	*two=jt.lep_indices[index[1]];
	  }
	  	
  	/*if(jt.num_trk>1){
  	cout<<"---------------------------------------------------------"<<endl;
  	cout<<"----------------QUI ----------------------QUI---------------------------"<<endl;
  	cout<<"indeces"<<endl;
  	cout<<"one: "<<*one<<" two: "<<*two<<endl;
  	cout<<"one: "<<jt.trk_indices[index[0]]<<" two: "<<jt.trk_indices[index[1]]<<endl;
  	for(k=0;k<jt.num_trk;k++){
  		cout<<"d0 significance "<<vect[k]<<endl;
  		cout<<"index: "<<index[k]<<endl;
  		}
  	}*/
  	
  	return;
	}
	
int find_first_sv(int dim, Jet_struct jt, Vertex_struct *vtx){
	int index[100];
	double vect[100];
	int k=0;
	if(dim==2){
		//sort sv by uncertainty on 2D flight distance (increasing)
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

