//variables for tagger
#include "heavyhadron_selection.C"
#include "selection.C"

#define dR_Jet 1.5
#define DEFAULT_VALUE -999

void tagger_variables(){
	
	long int ientry,n,k,i;
	Int_t nmcp,nrec,r2mnrel,nvt,nj,ntrk;
	
	//mc particles
	int num=2000000;
	Float_t *mcmox,*mcmoy,*mcmoz;
	Int_t *mcpdg,*mcpa0,*mcda0,*mcda1,*mcda2,*mcda3,*mcda4;
	mcpdg = (int*) malloc(sizeof(int)*num);
	mcpa0 = (int*) malloc(sizeof(int)*num);
	mcda0 = (int*) malloc(sizeof(int)*num);
	mcda1 = (int*) malloc(sizeof(int)*num);
	mcda2 = (int*) malloc(sizeof(int)*num);
	mcda3 = (int*) malloc(sizeof(int)*num);
	mcda4 = (int*) malloc(sizeof(int)*num);
	mcmox = (float*) malloc(sizeof(float)*num);
	mcmoy = (float*) malloc(sizeof(float)*num);
	mcmoz = (float*) malloc(sizeof(float)*num);
	
	//relations
	int num2=100000;
	Int_t *r2mt,*r2mf,*rcftr,*rcvts;
	Float_t *r2mw;
	r2mt = (int*) malloc(sizeof(int)*num2);
	r2mf = (int*) malloc(sizeof(int)*num2);
	r2mw = (float*) malloc(sizeof(float)*num2);
	rcftr = (int*) malloc(sizeof(int)*num2);
	rcvts = (int*) malloc(sizeof(int)*num2);
	
	//reco particles
	Float_t *rcmox,*rcmoy,*rcmoz,*rcene;
	Int_t *rccha,*rctyp;
	rcmox = (float*) malloc(sizeof(float)*num2);
	rcmoy = (float*) malloc(sizeof(float)*num2);
	rcmoz = (float*) malloc(sizeof(float)*num2);
	rcene = (float*) malloc(sizeof(float)*num2);
	rccha = (int*) malloc(sizeof(int)*num2);
	rctyp = (int*) malloc(sizeof(int)*num2);
	
	//tracks
	Float_t *trk_z0,*trk_d0,*trch2;
	Int_t *trk_atIP, *trndf;
	trk_z0 = (float*) malloc(sizeof(float)*num2);
	trk_d0 = (float*) malloc(sizeof(float)*num2);
	trch2 = (float*) malloc(sizeof(float)*num2);
	trk_atIP = (int*) malloc(sizeof(int)*num2);
	trndf = (int*) malloc(sizeof(int)*num2);
	Int_t trshn[10000][12];
	Float_t cov[10000][15];
	
	//vertices
	Float_t x[1000],y[1000],z[1000],chi[1000],vtcov[1000][6];
	Int_t id[100];
	
	//jets
	Float_t jmox[100],jmoy[100],jmoz[100],jene[100];
	
	Jet_struct Jet[100];
	Vertex_struct Vertex[100];
	Track_struct Track[10000];
	
	
	////TREES FOR MVA//////////////////////////////////////
	int vtx_cat; //0 for novertex and 1 for recovertex
	int n_sv; //number of SV per jet in RecoVertex category
	int n_trk_from_sv;
	double fd_sig_xy_sv; 
	double fd_sig_xyz_sv;
	double corrected_mass_sv;
	double energy_ratio_sv;
	double dR_sv_jet;
	double MassEnergyFraction_sv;
	double Boost_sv;
	
	
   	TTree *tree = new TTree("bb_tree", "bb_tree");
   	tree->Branch("vertex_cat", &vtx_cat, "vertex_cat/I");
        tree->Branch("n_SV", &n_sv, "n_SV/I");
        tree->Branch("n_trk_from_SV", &n_trk_from_sv, "n_trk_from_SV/I");
        tree->Branch("fd_sig_xy_SV", &fd_sig_xy_sv, "fd_sig_xy_SV/D");
        tree->Branch("fd_sig_xyz_SV", &fd_sig_xyz_sv, "fd_sig_xyz_SV/D");
        tree->Branch("corrected_mass_SV", &corrected_mass_sv, "corrected_mass_SV/D");
        tree->Branch("energy_ratio_SV", &energy_ratio_sv, "energy_ratio_SV/D");
        tree->Branch("dR_SV_JET", &dR_sv_jet, "dR_SV_JET/D");
        tree->Branch("MassEnergyFraction_SV", &MassEnergyFraction_sv, "MassEnergyFraction_SV/D");
        tree->Branch("Boost_SV", &Boost_sv, "Boost_SV/D");
        ////////////////////////////////////////////////////////
	
	
	//WHEN YOU CHANGE THE PATH CHANGE ALSO THIS FLAG
	char choice[]="b";

	TChain* fChain = new TChain("fChain");
   	fChain->Add("/home/paola/Scrivania/xml_aggiornati/3/ntuple_bb_1000_mod_3.root/MyLCTuple");
	
	fChain->SetBranchAddress("nmcp", &nmcp);
	fChain->SetBranchAddress("mcpdg", mcpdg);
	fChain->SetBranchAddress("mcpa0", mcpa0);
	fChain->SetBranchAddress("mcda0", mcda0);
	fChain->SetBranchAddress("mcda1", mcda1);
	fChain->SetBranchAddress("mcda2", mcda2);
	fChain->SetBranchAddress("mcda3", mcda3);
	fChain->SetBranchAddress("mcda4", mcda4);
	fChain->SetBranchAddress("mcmox", mcmox);
	fChain->SetBranchAddress("mcmoy", mcmoy);
	fChain->SetBranchAddress("mcmoz", mcmoz);
	
	fChain->SetBranchAddress("r2mnrel", &r2mnrel);
	fChain->SetBranchAddress("r2mt", r2mt);
	fChain->SetBranchAddress("r2mf", r2mf);	
        fChain->SetBranchAddress("r2mw", r2mw);
        
        fChain->SetBranchAddress("nrec", &nrec);
        fChain->SetBranchAddress("rctyp", rctyp);   
	fChain->SetBranchAddress("rcftr", rcftr);
	fChain->SetBranchAddress("rcvts", rcvts);
	fChain->SetBranchAddress("rcmox", rcmox);
  	fChain->SetBranchAddress("rcmoy", rcmoy);
  	fChain->SetBranchAddress("rcmoz", rcmoz);
	fChain->SetBranchAddress("rcene", rcene);
	fChain->SetBranchAddress("rccha", rccha);
	
	fChain->SetBranchAddress("ntrk",&ntrk);
	fChain->SetBranchAddress("tszze",trk_z0);
	fChain->SetBranchAddress("trsip", trk_atIP);
        fChain->SetBranchAddress("tsdze", trk_d0);
        fChain->SetBranchAddress("tscov", cov);
        fChain->SetBranchAddress("trch2", trch2);
	fChain->SetBranchAddress("trndf", trndf);
        fChain->SetBranchAddress("trshn",trshn);
        
        fChain->SetBranchAddress("nvt", &nvt);
        fChain->SetBranchAddress("vtori", id);
	fChain->SetBranchAddress("vtxxx", x);
	fChain->SetBranchAddress("vtyyy", y);
	fChain->SetBranchAddress("vtzzz", z);
	fChain->SetBranchAddress("vtchi", chi);
	fChain->SetBranchAddress("vtcov", vtcov);
	
	fChain->SetBranchAddress("nj", &nj);
   	fChain->SetBranchAddress("jmox", jmox);
	fChain->SetBranchAddress("jmoy", jmoy);
	fChain->SetBranchAddress("jmoz", jmoz);
	fChain->SetBranchAddress("jene", jene);
	
	TLorentzVector reco,vtx_momentum,all_trk_in_jet;
	TVector3  B_hadron, C_hadron;
	TVector3  jet_mo,jet_axis,PV,err_PV,SV,reco_v,flight_dir;
	int  index_of_jet,index_of_vertex,di,counter_sv;
	double delta_R,delta_R_min;
	double mom,theta;
	double X_val,Y_val,pt_v=0,energy_v=0,energy_j=0;;
	double xy_flight_dist,xyz_flight_dist,err_xy_fd,err_xyz_fd,termX,termY,termZ;
	int j_b[10]={-1}, j_c[10]={-1}, j_other[10]={-1}, j_light[10]={-1}, j_index[10]={-1};
	int nj_b=0,nj_other=0,nj_c=0,nj_light=0,is_equal=0,nj_set=0;
	
	int is_matched_with_vertex[100]={0},is_matched_with_jet[100]={-1},vtx_jet_matrix[100][100]={0};
	int vtx_trk_vector[100]={0},vtx_trk_matrix[100][100]={0};
	int jet_to_trk_matrix[100][500]={0},jet_to_trk_vector[100]={0},trk_is_matched_with_jet[100]={-1};
	
	//entry cycle
	//fChain->GetEntries();
	for(ientry=0; ientry< fChain->GetEntries(); ++ientry){
  		fChain->GetEntry(ientry);
  		
  		//initialization
  		for(k=0;k<10;k++){ j_b[k]=-1;j_c[k]=-1;j_other[k]=-1;j_light[k]=-1;}
  		
  		for(k=0;k<100;k++){ 
  			vtx_trk_vector[k]=0;	is_matched_with_jet[k]=-1;	 is_matched_with_vertex[k]=0;
  			jet_to_trk_vector[k]=0; 	trk_is_matched_with_jet[k]=-1;
  			for(n=0;n<100;n++){
  				vtx_trk_matrix[k][n]=-1;	vtx_jet_matrix[k][n]=-1;	jet_to_trk_matrix[k][n]=-1;
  				}
  			}
  		nj_b=0;		nj_c=0;		nj_other=0;	nj_light=0;
  		for(k=0;k<10;k++){ 
  			j_b[k]=-1;	j_c[k]=-1;	j_other[k]=-1;		j_light[k]=-1;		j_index[k]=-1;
  			}
  			
  		//jet
  		for(k=0;k<100;k++){
  			Jet[k].j.SetPxPyPzE(0.,0.,0.,0.);
  			Jet[k].selected="FALSE";
  			sprintf(Jet[k].flavour_tag,"o");
  			sprintf(Jet[k].vertex_category,"NoVertex");
  			sprintf(Jet[k].lepton_category,"NoSoftLepton");
  			Jet[k].num_SV=0;
  			for(i=0;i<50;i++){
  				Jet[k].SV_indices[i]=-1;
  				}
  			Jet[k].num_trk=0;
  			for(i=0;i<500;i++){
  				Jet[k].trk_indices[i]=-1;
  				}
  			}
  		//vertex
  		for(k=0;k<100;k++){
  			Vertex[k].v.SetXYZ(0.,0.,0.);
  			Vertex[k].errv.SetXYZ(0.,0.,0.);
  			Vertex[k].fd_2=0;
  			Vertex[k].fd_3=0;
  			Vertex[k].fd_sig_2=0;
  			Vertex[k].fd_sig_3=0;
  			Vertex[k].selected="FALSE";
  			Vertex[k].id=0;
  			Vertex[k].chi2=0;
  			Vertex[k].jet_index=-1;
  			Vertex[k].num_trk=0;
  			for(i=0;i<500;i++){
  				Vertex[k].trk_indices[i]=-1;
  				}
  			}
  			
  		//tracks
  		for(k=0;k<10000;k++){
  			Track[k].selected="FALSE";
  			Track[k].t.SetPxPyPzE(0.,0.,0.,0.);
  			Track[k].d0=0;
  			Track[k].z0=0;
  			Track[k].chi2=0;
  			Track[k].ndf=0;
  			for(i=0;i<15;i++){
  				Track[k].cov[i]=0;
  				}
  			for(i=0;i<12;i++){
  				Track[k].sub_det_hits[i]=0;
  				}
  			Track[k].jet_index=-1;
  			Track[k].SV_index=-1;
  			}
  			
  		
  		
  		
  		////////////////////////////////////////////////////////////////////////////////////////////
  		//////////////////////////// TRUE FLAVOUR TAG //////////////////////////////////////////////
  		////////////////////////////////////////////////////////////////////////////////////////////
  		// cycle on mc particle to find a B-hadron which has not a B-hadron daughter
  		for(k=0;k<nmcp;k++){
  			if(selection_bb(mcpdg[k])==0 && selection_bb(mcpdg[mcda0[k]])!=0 && selection_bb(mcpdg[mcda1[k]])!=0 && 			selection_bb(mcpdg[mcda2[k]])!=0 && selection_bb(mcpdg[mcda3[k]])!=0 && selection_bb(mcpdg[mcda4[k]])!=0){
  				//cout<<"B"<<endl;
  				B_hadron.SetXYZ(mcmox[k],mcmoy[k],mcmoz[k]);
  				delta_R_min=10000;
  				//find nearest jet
  				for(i=0;i<nj;i++){
  					jet_mo.SetXYZ(jmox[i],jmoy[i],jmoz[i]);
					jet_axis=(jet_mo);
					
					delta_R=jet_axis.DeltaR(B_hadron);
					//cout<<" index jet: "<<i<<" deltaR "<<delta_R<<endl;
					
					if(delta_R<delta_R_min){
						delta_R_min=delta_R;
						index_of_jet=i;
						}
					}
				if(delta_R_min < dR_Jet){
					is_equal=0;
					for(n=0;n<nj_b;n++){
						if(j_b[n]==index_of_jet){is_equal=1;}
						}
					if(is_equal==0){
						//jet is classified as a b-jet 
						j_b[nj_b]=index_of_jet;
						nj_b++;
						}
					}
				}
			}
		//cycle to identify jets not associated to B-hadrons and fill j_other
		for(i=0;i<nj;i++){
			is_equal=0;
			for(n=0;n<nj_b;n++){
				if(j_b[n]==i){is_equal=1;}
				}
			if(is_equal==0){
				j_other[nj_other]=i;
				nj_other++;
				}
			}
		// cycle on mc particle to find a C-hadron which has not a C-hadron daughter
		for(k=0;k<nmcp;k++){
			if(selection_cc(mcpdg[k])==0 && selection_cc(mcpdg[mcda0[k]])!=0 && selection_cc(mcpdg[mcda1[k]])!=0 && 			selection_cc(mcpdg[mcda2[k]])!=0 && selection_cc(mcpdg[mcda3[k]])!=0 && selection_cc(mcpdg[mcda4[k]])!=0){
				//cout<<"C"<<endl;
				C_hadron.SetXYZ(mcmox[k],mcmoy[k],mcmoz[k]);
  				delta_R_min=10000;
  				//find nearest jet between not b-jets
  				for(i=0;i<nj_other;i++){
  					jet_mo.SetXYZ(jmox[j_other[i]],jmoy[j_other[i]],jmoz[j_other[i]]);
					jet_axis=(jet_mo);
					
					delta_R=jet_axis.DeltaR(C_hadron);
					//cout<<" index jet: "<<j_other[i]<<" deltaR "<<delta_R<<endl;
					
					if(delta_R<delta_R_min){
						delta_R_min=delta_R;
						index_of_jet=j_other[i];
						}
					}
				if(delta_R_min < dR_Jet){
					is_equal=0;
					for(n=0;n<nj_c;n++){
						if(j_c[n]==index_of_jet){is_equal=1;}
						}
					if(is_equal==0){
						//jet is classified as a b-jet 
						j_c[nj_c]=index_of_jet;
						nj_c++;
						}
					}
				}
			}
			
		//cycle to identify jets not associated to B/C-hadrons and fill j_light
		for(i=0;i<nj_other;i++){
			jet_mo.SetXYZ(jmox[j_other[i]],jmoy[j_other[i]],jmoz[j_other[i]]);
			is_equal=0;
			for(n=0;n<nj_c;n++){
				if(j_c[n]==j_other[i]){is_equal=1;}
				}
			if(is_equal==0){
				j_light[nj_light]=j_other[i];
				nj_light++;
				}
			}
			
	/*cout<<"--------------------------------------------"<<endl;
	cout<<" In this event there are "<< nj<<" jets"<<endl;
	cout<<" In this event there are "<< nj_b<<" b-jets"<<endl;
	cout<< "Their indeces are: "<<endl;
	for(i=0;i<nj_b;i++){
		cout<<"index: "<<j_b[i]<<endl;
		}
	cout<<" In this event there are "<< nj_c<<" c-jets"<<endl;
	cout<< "Their indeces are: "<<endl;
	for(i=0;i<nj_c;i++){
		cout<<"index: "<<j_c[i]<<endl;
		}
	cout<<" In this event there are "<< nj_light<<" light-flavour-jets"<<endl;
	cout<< "Their indeces are: "<<endl;
	for(i=0;i<nj_light;i++){
		cout<<"index: "<<j_light[i]<<endl;
		}
		
	cout<<"END EVENT"<<endl;
	
	cout<<"--------------------------------------------"<<endl;*/
	
		//if the sample is bb take just b jets
		//if it is cc just the c jets
	
		if(choice[0]=='b'){
			nj_set=nj_b;
			for(i=0;i<nj_b;i++){
				j_index[i]=j_b[i];
				}
			}
			
		if(choice[0]=='c'){
			nj_set=nj_c;
			for(i=0;i<nj_c;i++){
				j_index[i]=j_c[i];
				}
			}
			
		////////////////////////////////////////////////////////////////////////////////////////////
  		////////////////////////////////////////////////////////////////////////////////////////////
  		//cycle on vertices to find the PV (id=103) of the event
  		 for (k=0;k<nvt; k++){
  		 	if(id[k]==103){
  		 		PV.SetXYZ(x[k],y[k],z[k]);
  		 		err_PV.SetXYZ(sqrt(vtcov[k][0]),sqrt(vtcov[k][2]),sqrt(vtcov[k][5]));
  		 		break;
  		 		}			 	
  		 	}
  		 
  		 //cycle on vertices to count SVs(id=104)	
  		 for (k=0;k<nvt; k++){
  		 	if(id[k]==104){
  		 		//SV flight direction
	  	 		SV.SetXYZ(x[k],y[k],z[k]);
	  	 		flight_dir=SV-PV;
	  	 		delta_R_min=10000;
  		 		for(i=0;i<nj;i++){
					//jet axis
					jet_mo.SetXYZ(jmox[i],jmoy[i],jmoz[i]);
					jet_axis=(jet_mo);
					delta_R=jet_axis.DeltaR(flight_dir);
					if(delta_R<delta_R_min){
						delta_R_min=delta_R;
						index_of_jet=i;
						}
					}
				if(delta_R_min < dR_Jet){
					is_matched_with_jet[k]=index_of_jet;
					//it contains the indeces of the verteces related to the jets
					vtx_jet_matrix[index_of_jet][is_matched_with_vertex[index_of_jet]]=k;
					is_matched_with_vertex[index_of_jet]++;
					}
				}
			}
			
		//cycle on reco-particles
  		for(i=0;i<nrec;i++){
  			//if it has a vertex
  			if(rcvts[i]!=-1){
  				//row index is the id number of the vertex, the content is the id number of the tracks related to that vertex
  				vtx_trk_matrix[rcvts[i]][vtx_trk_vector[rcvts[i]]]=i;
  				//array index is the id number of the vertex, the content is the number of tracks related to that vertex 
  				vtx_trk_vector[rcvts[i]]++;
  				}
  			}
  			
  		for(i=0;i<nrec;i++){
  			if(rcftr[i]!=-1){
  				reco.SetPxPyPzE(rcmox[i],rcmoy[i],rcmoz[i],rcene[i]);
  				reco_v=reco.Vect();
  				delta_R_min=10000;
  				for(k=0; k<nj; k++){
  					jet_mo.SetXYZ(jmox[k],jmoy[k],jmoz[k]);
					jet_axis=(jet_mo);
  					delta_R=jet_axis.DeltaR(reco_v);
					if(delta_R<delta_R_min){
						delta_R_min=delta_R;
						index_of_jet=k;
						}
					
					}
				if(delta_R_min < dR_Jet){
					jet_to_trk_matrix[index_of_jet][jet_to_trk_vector[index_of_jet]]=i;
					jet_to_trk_vector[index_of_jet]++;
					trk_is_matched_with_jet[i]=index_of_jet;
					}
				}
			}
			
		//fill Track Structure
  		for(i=0;i<nrec;i++){
  			if(rcftr[i]!=-1){
  				Track[i].t.SetPxPyPzE(rcmox[i],rcmoy[i],rcmoz[i],rcene[i]);
  				Track[i].d0=trk_d0[trk_atIP[rcftr[i]]];
  				Track[i].z0=trk_z0[trk_atIP[rcftr[i]]];
  				Track[i].chi2=trch2[rcftr[i]];
  				Track[i].ndf=trndf[rcftr[i]];
  				for(n=0;n<15;n++){
  					Track[i].cov[n]=cov[trk_atIP[rcftr[i]]][n];
  					}
  				for(n=0;n<12;n++){
  					Track[i].sub_det_hits[n]=trshn[rcftr[i]][n];
  					}
  				Track[i].SV_index=rcvts[i];
  				Track[i].jet_index=trk_is_matched_with_jet[i];
  				Track[i].selected=trk_selection(Track[i]);
  				}
  			}
	
  		
  		//fill Vertex Structure
  		for (k=0;k<nvt; k++){
  			Vertex[k].v.SetXYZ(x[k],y[k],z[k]);
  			Vertex[k].errv.SetXYZ(sqrt(vtcov[k][0]),sqrt(vtcov[k][2]),sqrt(vtcov[k][5]));
  			
  			xy_flight_dist = sqrt(pow(Vertex[k].v.X()-PV.X(),2)+pow(Vertex[k].v.Y()-PV.Y(),2));
  			termX=pow(Vertex[k].v.X()-PV.X(),2)*pow(Vertex[k].errv.X(),2);
  			termY=pow(Vertex[k].v.Y()-PV.Y(),2)*pow(Vertex[k].errv.Y(),2);
  			err_xy_fd = sqrt(termX+termY)/xy_flight_dist;
  			Vertex[k].fd_2=xy_flight_dist;
  			Vertex[k].fd_sig_2=xy_flight_dist/err_xy_fd;
  			
  			xyz_flight_dist = (Vertex[k].v-PV).Mag();
  			termZ=pow(Vertex[k].v.Z()-PV.Z(),2)*pow(Vertex[k].errv.Z(),2);
  			err_xyz_fd = sqrt(termX+termY+termZ)/xyz_flight_dist;
  			Vertex[k].fd_3=xyz_flight_dist;
  			Vertex[k].fd_sig_3=xyz_flight_dist/err_xyz_fd;
  			
  			Vertex[k].id=id[k];
  			Vertex[k].chi2=chi[k];
  			Vertex[k].jet_index=is_matched_with_jet[k];
  			Vertex[k].num_trk=vtx_trk_vector[k];
  			
  			for(i=0;i<vtx_trk_vector[k];i++){
  				Vertex[k].trk_indices[i]=vtx_trk_matrix[k][i];
  				}
  				
  			Vertex[k].selected=vertex_selection(Vertex[k]);
  			}
  			
  		// fill Jet Structure
  		for(n=0;n<nj;n++){
  			//4-mom vector
  			Jet[n].j.SetPxPyPzE(jmox[n],jmoy[n],jmoz[n],jene[n]);
  			//true flavour tag
  			for(i=0;i<nj_set;i++){
  				if(j_index[i]==n){
  					sprintf(Jet[n].flavour_tag,"%s",choice);
  					}
  				}
  				
  			//num trk
  			Jet[n].num_trk=jet_to_trk_vector[n];
  			//trk_indices
  			for(k=0;k<jet_to_trk_vector[n];k++){
  				Jet[n].trk_indices[k]=jet_to_trk_matrix[n][k];
  				}	
  			Jet[n].selected=jet_selection(Track,Jet[n]);
  			
  			
  			//num SV
  			Jet[n].num_SV=is_matched_with_vertex[n];
  			
  			counter_sv=0;
  			//SV indices
  			for(k=0;k<is_matched_with_vertex[n];k++){
  				
  				Jet[n].SV_indices[k]=vtx_jet_matrix[n][k];
  				
  				if(Vertex[vtx_jet_matrix[n][k]].selected==true){
  					counter_sv++;
  					}
  				}
  			if(counter_sv>0){
  				sprintf(Jet[n].vertex_category,"RecoVertex");
  				}
  			for(k=0;k<Jet[n].num_trk;k++){
  				if(abs(rctyp[Jet[n].trk_indices[k]])==13&&Track[Jet[n].trk_indices[k]].selected==true){
  					sprintf(Jet[n].lepton_category,"SoftMuon");
  					break;
  					}
  				if(abs(rctyp[Jet[n].trk_indices[k]])==11&&Track[Jet[n].trk_indices[k]].selected==true){
  					sprintf(Jet[n].lepton_category,"SoftElectron");
  					}
  				}
  			}
  		if(ientry<5){
  			cout<<" "<<endl;
  			cout<<"---------------Jet struct----------------------"<<endl;
  			//jet
  			for(k=0;k<nj;k++){
  				cout<<"INDICE JET: "<<k<<endl;
  				cout<<"4-mom: E Px Py Pz  "<<Jet[k].j.E()<<" "<<Jet[k].j.Px()<<" "<<Jet[k].j.Py()<<" " <<Jet[k].j.E()<<endl;
  				cout<<"pt: "<<Jet[k].j.Pt()<<endl;
  				cout<<"Eta: "<<Jet[k].j.Eta()<<endl;
  				cout<<"selected: "<<Jet[k].selected<<endl;
  				cout<<"Falvour tag: "<<Jet[k].flavour_tag<<endl;
  				cout<<"Vertex cat: "<<Jet[k].vertex_category<<endl;
  				cout<<"lep cat: "<<Jet[k].lepton_category<<endl;
  				cout<<"Num SV: "<<Jet[k].num_SV<<endl;
  				for(i=0;i<Jet[k].num_SV;i++){
  					cout<<"SV indices: "<<Jet[k].SV_indices[i]<<endl;
  					}
  				cout<<"num trk: "<<Jet[k].num_trk<<endl;
  				for(i=0;i<Jet[k].num_trk;i++){
  					cout<<"trk indices: "<<Jet[k].trk_indices[i]<<endl;
  					}
  				}
  			cout<<" "<<endl;
  			cout<<"---------------Vertex struct----------------------"<<endl;
  			//vertex
  			for(k=0;k<nvt;k++){
  				cout<<"INDICE VERTEX: "<<k<<endl;
  				cout<<"coor:  X Y Z  "<<Vertex[k].v.X()<<" "<<Vertex[k].v.Y()<<" "<<Vertex[k].v.Z()<<endl;
  				cout<<"err coor:  X Y Z  "<<Vertex[k].errv.X()<<" "<<Vertex[k].errv.Y()<<" "<<Vertex[k].errv.Z()<<endl;
  				cout<<"2d flight dist sig: "<<Vertex[k].fd_sig_2<<endl;
  				cout<<"3d flight dist sig: "<<Vertex[k].fd_sig_3<<endl;
  				cout<<"selected: "<<Vertex[k].selected<<endl;
  				cout<<"Vtx id: "<<Vertex[k].id<<endl;
  				cout<<"chi2: "<<Vertex[k].chi2<<endl;
  				cout<<"Jet index: "<<Vertex[k].jet_index<<endl;
  				cout<<"Num trk: "<<Vertex[k].num_trk<<endl;
  				for(i=0;i<Vertex[k].num_trk;i++){
  					cout<<"trk indices: "<<Vertex[k].trk_indices[i]<<endl;
  					}
  				}
  			cout<<" "<<endl;
  			/*cout<<"---------------Track struct----------------------"<<endl;
  			//trk
  			for(k=0;k<nrec;k++){
  			if(rcftr[k]!=-1){
  				cout<<"INDICE TRACK: "<<k<<endl;
  				cout<<"pt: "<<Track[k].t.Pt()<<endl;
  				cout<<"4-mom: E Px Py Pz "<<Track[k].t.E()<<" "<<Track[k].t.Px()<<" "<<Track[k].t.Py()<<" "<<Track[k].t.Pz()<<" "<<endl;
  				cout<<"trk d0: "<<Track[k].d0<<endl;
  				cout<<"trk z0: "<<Track[k].z0<<endl;
  				cout<<"trk chi2/ndf: "<<Track[k].chi2/Track[k].ndf<<endl;
  				for(n=0;n<12;n++){
  					cout<<"sub det hits: "<<n<<" : "<<Track[k].sub_det_hits[n]<<endl;
  					}
  				cout<<"SV index: "<<Track[k].SV_index<<endl;
  				cout<<"jet index: "<<Track[k].jet_index<<endl;
  				cout<<"selected: "<<Track[k].selected<<endl;
  				}
  			}	*/
  			}
  			
  		for(k=0;k<nj;k++){
  			if(Jet[k].selected==true && strcmp(Jet[k].flavour_tag,choice)==0){
  			
  				if(strcmp(Jet[k].vertex_category,"RecoVertex")==0){
  				
  					 vtx_cat=1;
  					 counter_sv=0;
  					 for(n=0;n<Jet[k].num_SV;n++){
	  				 	if(Vertex[Jet[k].SV_indices[n]].selected==true){
							counter_sv++;	
							}
						}
	  				 n_sv=counter_sv;
	  				 di=3;
	  				 index_of_vertex=find_first_sv(di,Jet[k],Vertex);
	  				 fd_sig_xy_sv=Vertex[index_of_vertex].fd_sig_2;
	  				 
	  				 fd_sig_xyz_sv=Vertex[index_of_vertex].fd_sig_3;
	  				 
	  				 n_trk_from_sv=Vertex[index_of_vertex].num_trk;
	  				 
	  				 //cycle on number of tracks associated to the vertex
	  				 vtx_momentum.SetPxPyPzE(0.,0.,0.,0.);
	  				 all_trk_in_jet.SetPxPyPzE(0.,0.,0.,0.);
	  				 
  					for(n=0;n<Vertex[index_of_vertex].num_trk;n++){
  						
  						vtx_momentum=vtx_momentum+Track[Vertex[index_of_vertex].trk_indices[n]].t;
  						}
  					flight_dir=Vertex[index_of_vertex].v-PV;
	  	 			mom=sqrt(pow(vtx_momentum.Px(),2)+pow(vtx_momentum.Py(),2)+pow(vtx_momentum.Pz(),2));
	  	 			theta=flight_dir.Angle(vtx_momentum.Vect());
	  	 			corrected_mass_sv=sqrt(pow(vtx_momentum.M(),2)+pow(mom*sin(theta),2))+mom*sin(theta);
	  	 			
	  	 			for(n=0;n<Jet[k].num_trk;n++){
	  	 				if(Track[Jet[k].trk_indices[n]].selected==true){
	  	 					all_trk_in_jet=all_trk_in_jet+Track[Jet[k].trk_indices[n]].t;
	  	 					}
	  	 				}
	  				energy_ratio_sv=vtx_momentum.E()/all_trk_in_jet.E();
	  				dR_sv_jet=flight_dir.DeltaR(Jet[k].j.Vect());
	  				pt_v=0;
	  				energy_v=0;
	  				energy_j=0;
	  				for(n=0;n<Vertex[index_of_vertex].num_trk;n++){
  						pt_v=pt_v+abs(Track[Vertex[index_of_vertex].trk_indices[n]].t.Pt());
  						energy_v=energy_v+Track[Vertex[index_of_vertex].trk_indices[n]].t.E();
  						}
	  				for(n=0;n<Jet[k].num_trk;n++){
	  	 				if(Track[Jet[k].trk_indices[n]].selected==true){
	  	 					energy_j=energy_j+Track[Jet[k].trk_indices[n]].t.E();
	  	 					}
	  	 				}
	  				X_val=corrected_mass_sv*(energy_v/energy_j)/5.2794;
	  				MassEnergyFraction_sv=X_val/(X_val+0.04);
	  				Y_val=sqrt(5.2794)*pt_v/(corrected_mass_sv*sqrt(Jet[k].j.Pt()));
	  				Boost_sv=pow(Y_val,2.)/(pow(Y_val,2.)+10.);
	  				
	  				}
	  			else if(strcmp(Jet[k].vertex_category,"NoVertex")==0){
	  				 vtx_cat=0;
	  				 n_sv=DEFAULT_VALUE;
	  				 fd_sig_xy_sv=DEFAULT_VALUE;
	  				 fd_sig_xyz_sv=DEFAULT_VALUE;
	  				 n_trk_from_sv=DEFAULT_VALUE;
	  				 corrected_mass_sv=DEFAULT_VALUE;
	  				 energy_ratio_sv=DEFAULT_VALUE;
	  				 dR_sv_jet=DEFAULT_VALUE;
	  				 MassEnergyFraction_sv=DEFAULT_VALUE;
	  				 Boost_sv=DEFAULT_VALUE;
	  				 }
	  			else{
	  				cout<<"stem a sbaggjò qualche cos"<<endl;
	  				}
	  			
	  			tree->Fill();
  				}
  			
  			}
  		}
  		
  	tree->Print();
  		
  	//FILL TREES FOR MVA	
	TFile *outfile = new TFile("./bb_tree_forMVA.root", "recreate");
      	outfile->cd("");	
      	tree->Write();
      	outfile->Close();
	
	}
