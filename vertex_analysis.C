////////DESCRIPTION///////////////////////////////////////////////////////////////////////////////////////////////////////
//Con bb sample:													//
//1. Vertici secondari ricostruiti per evento										//
//2. Vertici secondari per jet												//
//3. Distanza PV-SV													//
//4. deltaR tra vertice secondario e jet più vicino									//
//5. massa corretta SV (+ singole variabili della formula) vedi pag 10 paper:						//
//	"Identification of heavy-flavour jets with the CMS detector in pp collisions at 13 TeV"	 			//
//6. Numero di tracce con cui è fittato primario/secondario								//
//7. chi2 vertice secondario/primario VS numero di tracce con cui è ricostruito						//
//8. quante tracce che provengono da PV sono associate erroneamente a SV per evento					//
//9. quante tracce che provengono dal decadimento di B/C sono associate erroneamente al PV per evento			//
//10. numero di SV ricostruiti con particelle che provengono da B/C/resto per evento					//
//11. Per 5 eventi stampo info dettagliate su quali sono le particelle a MC associate a quali vertici ricostruiti/veri  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TLorentzVector.h"

//SELECTION FOR CC
int selection_cc(int mcid){

	if(abs(mcid)==411 || abs(mcid)== 421 || abs(mcid)==10411 || abs(mcid)==10421 || abs(mcid)==413 ||
		   abs(mcid)==423 || abs(mcid)==10413 || abs(mcid)==10423 || abs(mcid)==20413 || abs(mcid)==20423 ||
     		   abs(mcid)==415 || abs(mcid)==425 || abs(mcid)==431 || abs(mcid)==10431 || abs(mcid)==433 || 			   abs(mcid)==10433 ||abs(mcid)==20433 || abs(mcid)==435 || abs(mcid)==4122 || abs(mcid)==4222 || abs(mcid)==4212 || abs(mcid)==4112 || abs(mcid)==4224 || abs(mcid)==4214 || abs(mcid)==4114 || abs(mcid)==4232 || abs(mcid)==4132 || abs(mcid)==4322|| abs(mcid)==4312 || abs(mcid)==4324 || abs(mcid)==4314 || abs(mcid)==4332 || abs(mcid)==4334 || abs(mcid)==4412 || abs(mcid)==4422 || abs(mcid)==4414 || abs(mcid)==4424 || abs(mcid)==4432 || abs(mcid)==4434 || abs(mcid)==4444 || abs(mcid)==441 || abs(mcid)==10441 || abs(mcid)==100441 || abs(mcid)==443 || abs(mcid)==10443 || abs(mcid)==20443 || abs(mcid)==100443 || abs(mcid)==30443 || abs(mcid)==9000443 || abs(mcid)==9010443 || abs(mcid)==9020443 || abs(mcid)==445 || abs(mcid)==100445 ){ return 0;}

	else return -1;

}


//SELECTION FOR BB
int selection_bb(int mcid){

	if(abs(mcid)==511 || abs(mcid)==521 || abs(mcid)==10511 || abs(mcid)==10521 || abs(mcid)==513 || abs(mcid)==523 || abs(mcid)==10513 || abs(mcid)==10523 || abs(mcid)==20513 || abs(mcid)==20523 || abs(mcid)==515 || abs(mcid)==525 || abs(mcid)==531 || abs(mcid)==10531 || abs(mcid)==533 || abs(mcid)==10533 || abs(mcid)==20533 || abs(mcid)==535 || abs(mcid)==541 || abs(mcid)==10541 || abs(mcid)==543 || abs(mcid)==10543 || abs(mcid)==20543 || abs(mcid)==545 || abs(mcid)==551 || abs(mcid)==10551 || abs(mcid)==100551 || abs(mcid)==110551 || abs(mcid)==200551 || abs(mcid)==210551 || abs(mcid)==553 || abs(mcid)==10553 || abs(mcid)==20553 || abs(mcid)==30553 || abs(mcid)==100553 || abs(mcid)==110553 || abs(mcid)==120553 || abs(mcid)==130553 || abs(mcid)==200553 || abs(mcid)==210553 || abs(mcid)==220553 || abs(mcid)==300553 || abs(mcid)==9000553 || abs(mcid)==9010553 || abs(mcid)==555 || abs(mcid)==10555 || abs(mcid)==20555 || abs(mcid)==100555 || abs(mcid)==110555 || abs(mcid)==120555 || abs(mcid)==200555 || abs(mcid)==557 || abs(mcid)==100557 || abs(mcid)==5122 || abs(mcid)==5112 || abs(mcid)==5212 || abs(mcid)==5222 || abs(mcid)==5114 || abs(mcid)==5214 || abs(mcid)==5224 || abs(mcid)==5132 || abs(mcid)==5232 || abs(mcid)==5312 || abs(mcid)==5322 || abs(mcid)==5314 || abs(mcid)==5324 || abs(mcid)==5332 || abs(mcid)==5334 || abs(mcid)==5142 || abs(mcid)==5242 || abs(mcid)==5412 || abs(mcid)==5422 || abs(mcid)==5414 || abs(mcid)==5424 || abs(mcid)==5342 || abs(mcid)==5432 || abs(mcid)==5434 || abs(mcid)==5442 || abs(mcid)==5444 || abs(mcid)==5512 || abs(mcid)==5522 || abs(mcid)==5514 || abs(mcid)==5524 || abs(mcid)==5532 || abs(mcid)==5534 || abs(mcid)==5542 || abs(mcid)==5544 || abs(mcid)==5554){ return 0;}

	else return -1;
}


void vertex_analysis(){
	
	gStyle->SetOptStat("nemrou");
	
	Int_t nmcp,nrec,r2mnrel,nvt,nj,ntrk;
	
	//mc particles
	int num=2000000;
	Float_t *mcvtx,*mcvty,*mcvtz;
	Int_t *mcpdg,*mcpa0;
	mcpdg = (int*) malloc(sizeof(int)*num);
	mcpa0 = (int*) malloc(sizeof(int)*num);
	mcvtx = (float*) malloc(sizeof(float)*num);
	mcvty = (float*) malloc(sizeof(float)*num);
	mcvtz = (float*) malloc(sizeof(float)*num);
	
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
	rcmox = (float*) malloc(sizeof(float)*num2);
	rcmoy = (float*) malloc(sizeof(float)*num2);
	rcmoz = (float*) malloc(sizeof(float)*num2);
	rcene = (float*) malloc(sizeof(float)*num2);
	
	//tracks
	Float_t *trk_z0,*trk_d0;
	Int_t *trk_atIP;
	trk_z0 = (float*) malloc(sizeof(float)*num2);
	trk_d0 = (float*) malloc(sizeof(float)*num2);
	trk_atIP = (int*) malloc(sizeof(int)*num2);
	Float_t cov[100000][15];
	
	//vertices
	Float_t x[100],y[100],z[100],chi[100];
	Int_t id[100];
	
	//jets
	Float_t jmox[100],jmoy[100],jmoz[100],jene[100];
	
	
	TChain* fChain = new TChain("fChain");
   	fChain->Add("/home/paola/Scrivania/xml_aggiornati/ntuple_bb_1000_newgeom.root/MyLCTuple");
	
	fChain->SetBranchAddress("nmcp", &nmcp);
	fChain->SetBranchAddress("mcpdg", mcpdg);
	fChain->SetBranchAddress("mcpa0", mcpa0);
	fChain->SetBranchAddress("mcvtx", mcvtx);
	fChain->SetBranchAddress("mcvty", mcvty);
	fChain->SetBranchAddress("mcvtz", mcvtz);
	
	fChain->SetBranchAddress("r2mnrel", &r2mnrel);
	fChain->SetBranchAddress("r2mt", r2mt);
	fChain->SetBranchAddress("r2mf", r2mf);	
        fChain->SetBranchAddress("r2mw", r2mw);
        
        fChain->SetBranchAddress("nrec", &nrec);   
	fChain->SetBranchAddress("rcftr", rcftr);
	fChain->SetBranchAddress("rcvts", rcvts);
	fChain->SetBranchAddress("rcmox", rcmox);
  	fChain->SetBranchAddress("rcmoy", rcmoy);
  	fChain->SetBranchAddress("rcmoz", rcmoz);
	fChain->SetBranchAddress("rcene", rcene);
	
	fChain->SetBranchAddress("ntrk",&ntrk);
	fChain->SetBranchAddress("tszze",trk_z0);
	fChain->SetBranchAddress("trsip", trk_atIP);
        fChain->SetBranchAddress("tsdze", trk_d0);
        fChain->SetBranchAddress("tscov", cov);
        
        fChain->SetBranchAddress("nvt", &nvt);
        fChain->SetBranchAddress("vtori", id);
	fChain->SetBranchAddress("vtxxx", x);
	fChain->SetBranchAddress("vtyyy", y);
	fChain->SetBranchAddress("vtzzz", z);
	fChain->SetBranchAddress("vtchi", chi);
	
	fChain->SetBranchAddress("nj", &nj);
   	fChain->SetBranchAddress("jmox", jmox);
	fChain->SetBranchAddress("jmoy", jmoy);
	fChain->SetBranchAddress("jmoz", jmoz);
	fChain->SetBranchAddress("jene", jene);
	
	/// PLOT 1 //////////////////////////////
	TH1F* SV_reco_perevent=new TH1F("SV_reco_perevent","SV_reco_perevent",12,-0.5,11.5);
	SV_reco_perevent->SetTitle("Number of SVs reconstructed per event");
   	SV_reco_perevent->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_perevent->SetLineWidth(2);
	/// PLOT 2 //////////////////////////////
	TH1F* SV_reco_perjet=new TH1F("SV_reco_perjet","SV_reco_perjet",12,-0.5,11.5);
	SV_reco_perjet->SetTitle("Number of SVs reconstructed per jet");
   	SV_reco_perjet->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_perjet->SetLineWidth(2);
	/// PLOT 3 //////////////////////////////
	TH1F* SV_to_PV_distance=new TH1F("SV_to_PV_distance","SV_to_PV_distance",150,0.,1.5);
	SV_to_PV_distance->SetTitle("SV flight distance");
   	SV_to_PV_distance->GetXaxis()->SetTitle("distance[mm]");
	SV_to_PV_distance->SetLineWidth(2);
	/// PLOT 4 //////////////////////////////
	TH1F* deltaR_SV_nearestJet=new TH1F("deltaR_SV_nearestJet","deltaR_SV_nearestJet",100,0.,10);
	deltaR_SV_nearestJet->SetTitle("#DeltaR between SV flight direction and its nearest jet axis");
   	deltaR_SV_nearestJet->GetXaxis()->SetTitle("#DeltaR");
	deltaR_SV_nearestJet->SetLineWidth(2);
	/// PLOT 5 //////////////////////////////
	TH1F* vtx_invariant_mass=new TH1F("vtx_invariant_mass","vtx_invariant_mass",100,0.,12);
	vtx_invariant_mass->SetTitle("Invariant mass of tracks associated to SV");
   	vtx_invariant_mass->GetXaxis()->SetTitle("M_{INV}[GeV]");
	vtx_invariant_mass->SetLineWidth(2);
	
	TH1F* vtx_corrected_mass=new TH1F("vtx_corrected_mass","vtx_corrected_mass",100,0.,40);
	vtx_corrected_mass->SetTitle("Corrected SV mass");
   	vtx_corrected_mass->GetXaxis()->SetTitle("M[GeV]");
	vtx_corrected_mass->SetLineWidth(2);
	
	TH1F* vtx_mom=new TH1F("vtx_mom","vtx_mom",100,0.,50);
	vtx_mom->SetTitle("SV momentum");
   	vtx_mom->GetXaxis()->SetTitle("P[GeV/c^{2}]");
	vtx_mom->SetLineWidth(2);
	
	TH1F* vtx_theta=new TH1F("vtx_theta","vtx_theta",50,0.,3.3);
	vtx_theta->SetTitle("Angle between SV momentum and SV flight direction");
   	vtx_theta->GetXaxis()->SetTitle("#theta[rad]");
	vtx_theta->SetLineWidth(2);
	
	// PLOT 6 ///////////////////////////////
	TH1F* num_trk_forPV=new TH1F("num_trk_forPV","num_trk_forPV",30,-0.5,29.5);
	num_trk_forPV->SetTitle("Number of tracks used to fit PV");
   	num_trk_forPV->GetXaxis()->SetTitle("N_{trk}");
	num_trk_forPV->SetLineWidth(2);
	
	TH1F* num_trk_forSV=new TH1F("num_trk_forSV","num_trk_forSV",17,-0.5,16.5);
	num_trk_forSV->SetTitle("Number of tracks used to fit SV");
   	num_trk_forSV->GetXaxis()->SetTitle("N_{trk}");
	num_trk_forSV->SetLineWidth(2);
	// PLOT 7 ///////////////////////////////
	TH2F* chi2_VS_num_trk_PV=new TH2F("chi2_VS_num_trk_PV","chi2_VS_num_trk_PV",30,-0.5,29.5,55,0,55.);
	chi2_VS_num_trk_PV->SetTitle("PV #chi^{2} VS number of tracks used to fit PV");
   	chi2_VS_num_trk_PV->GetXaxis()->SetTitle("N_{trk}");
   	chi2_VS_num_trk_PV->GetYaxis()->SetTitle("#chi^{2}");
	
	TH2F* chi2_VS_num_trk_SV=new TH2F("chi2_VS_num_trk_SV","chi2_VS_num_trk_SV",17,-0.5,16.5,80,0,20.);
	chi2_VS_num_trk_SV->SetTitle("SV #chi^{2} VS number of tracks used to fit SV");
   	chi2_VS_num_trk_SV->GetXaxis()->SetTitle("N_{trk}");
   	chi2_VS_num_trk_SV->GetYaxis()->SetTitle("#chi^{2}");
	// PLOT 8 ///////////////////////////////
	TH1F* from_hh_wrong_to_PV=new TH1F("from_hh_wrong_to_PV","from_hh_wrong_to_PV",10,-0.5,9.5);
	from_hh_wrong_to_PV->SetTitle("Number of tracks from B/C decay wrongly associated to PV per event");
   	from_hh_wrong_to_PV->GetXaxis()->SetTitle("N_{trk}");
	from_hh_wrong_to_PV->SetLineWidth(2);
	// PLOT 9 ///////////////////////////////
	TH1F* from_IP_wrong_to_SV=new TH1F("from_IP_wrong_to_SV","from_IP_wrong_to_SV",10,-0.5,9.5);
	from_IP_wrong_to_SV->SetTitle("Number of tracks from hadronization wrongly associated to SV per event");
   	from_IP_wrong_to_SV->GetXaxis()->SetTitle("N_{trk}");
	from_IP_wrong_to_SV->SetLineWidth(2);
	
	// PLOT 10 ///////////////////////////////
	TH1F* SV_reco_B=new TH1F("SV_reco_B","SV_reco_B",10,-0.5,9.5);
	SV_reco_B->SetTitle("Number of SV (per event) recostructed with at least a track from B decay");
   	SV_reco_B->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_B->SetLineWidth(2);
	
	TH1F* SV_reco_C=new TH1F("SV_reco_C","SV_reco_C",10,-0.5,9.5);
	SV_reco_C->SetTitle("Number of SV (per event) recostructed with at least a track from C decay");
   	SV_reco_C->GetXaxis()->SetTitle("N_{SV}");
   	SV_reco_C->SetLineColor(kRed);
	SV_reco_C->SetLineWidth(2);
	
	TH1F* SV_reco_BC=new TH1F("SV_reco_BC","SV_reco_BC",10,-0.5,9.5);
	SV_reco_BC->SetTitle("Number of SV (per event) recostructed with at least a track both from B and C decay");
   	SV_reco_BC->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_BC->SetLineWidth(2);
	SV_reco_BC->SetLineColor(kGreen+3);
	
	TH1F* SV_reco_noBC=new TH1F("SV_reco_noBC","SV_reco_noBC",10,-0.5,9.5);
	SV_reco_noBC->SetTitle("Number of SV (per event) recostructed with no tracks from B/C decay");
   	SV_reco_noBC->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_noBC->SetLineWidth(2);
	SV_reco_noBC->SetLineColor(kMagenta+1);
	
	long int ientry,k,i,n,m;
	int num_SV_perevent=0,num_SV_perjet=0,index_of_jet=-1,vtx_trk_vector[100]={0},vtx_trk_matrix[100][100]={0},wrong_PV=0,wrong_SV=0;
	int num_SV_B=0,num_SV_C=0,num_SV_BC=0,num_SV_noBC=0,num_trk_B=0,num_trk_C=0,num_trk_BC=0;
	double delta_R=0,delta_R_min=10000,corrected_mass=0,mom=0,theta=0;
	TLorentzVector vtx_momentum,reco;
	TVector3 PV,SV,flight_dir,jet_mo,jet_axis;
	
	//entry cycle
	for(ientry=0; ientry< fChain->GetEntries(); ++ientry){
  		fChain->GetEntry(ientry);
  		
  		//cycle on vertices to find the PV (id=103) of the event and count SVs(id=104)
  		 for (k=0;k<nvt; k++){
  		 	if(id[k]==103){
  		 		PV.SetX(x[k]);
  		 		PV.SetY(y[k]);
  		 		PV.SetZ(z[k]);
  		 		//cout<<"PV =  X: "<<x[k]<<" Y: "<<y[k]<<" Z: "<<z[k]<<endl;
  		 		break;
  		 		}			 	
  		 	}
  		 //cycle on vertices to count SVs(id=104), evaluate flight distance, find nearest jet axis and plot deltaR	
  		 for (k=0;k<nvt; k++){
  		 	if(id[k]==104){
  		 		num_SV_perevent++;
  		 		//SV flight direction
	  	 		SV.SetX(x[k]);
	  	 		SV.SetY(y[k]);
	  	 		SV.SetZ(z[k]);
	  	 		flight_dir=SV-PV;
  		 		//fill histo with SV flight distance//////////////////////////////////////////////
  		 		SV_to_PV_distance->Fill(flight_dir.Mag());///////////// PLOT 3 ///////////////////
  		 		//////////////////////////////////////////////////////////////////////////////////
  		 		
  		 		//find the nearest jet axis and plot deltaR
  		 		for(i=0;i<nj;i++){
					//jet axis
					jet_mo.SetX(jmox[i]);
					jet_mo.SetY(jmoy[i]);
					jet_mo.SetZ(jmoz[i]);
					jet_axis=(jet_mo-PV);
					delta_R=sqrt(pow(jet_axis.Eta()-flight_dir.Eta(),2)+pow(jet_axis.Phi()-flight_dir.Phi(),2));
					if(delta_R<delta_R_min){
						delta_R_min=delta_R;
						index_of_jet=i;
						}
					}
				//fill histo with deltaR min//////////////////////////////////////////////////////
  		 		deltaR_SV_nearestJet->Fill(delta_R_min);///////////// PLOT 4//////////////////////
  		 		//////////////////////////////////////////////////////////////////////////////////
				delta_R_min=10000;
  		 		}
  		 	}//vertex cycle
  		 
  		//fill histo with number of SV per event//////////////////////////////////////////////////////////	
  		SV_reco_perevent->Fill(num_SV_perevent);//////////////////// PLOT 1 //////////////////////////////
  		//////////////////////////////////////////////////////////////////////////////////////////////////
  		
  		num_SV_perevent=0;
  		
  		//cyle on jets
		for(i=0;i<nj;i++){
			//jet axis
			jet_mo.SetX(jmox[i]);
			jet_mo.SetY(jmoy[i]);
			jet_mo.SetZ(jmoz[i]);
			jet_axis=(jet_mo-PV);
			for (n=0;n<nvt; n++){ 
  	 			//SV
  	 			if(id[n]==104){
  	 	 			//SV flight direction
	  	 			SV.SetX(x[n]);
	  	 			SV.SetY(y[n]);
	  	 			SV.SetZ(z[n]);
	  	 			flight_dir=SV-PV;
	  	 			//deltaR between jet_axis and flight_dir
	  	 			delta_R=sqrt(pow(jet_axis.Eta()-flight_dir.Eta(),2)+pow(jet_axis.Phi()-flight_dir.Phi(),2));
	  	 			if(delta_R<1.){
	  	 				num_SV_perjet++;
		  	 			}
  	 				}
  	 			}
			//fill histo with number of SV per jet//////////////////////////////////////////////////////////
  	 		SV_reco_perjet->Fill(num_SV_perjet);//////////////////// PLOT 2 //////////////////////////////
  			//////////////////////////////////////////////////////////////////////////////////////////////////

  	 		num_SV_perjet=0;
  	 		}//jet cycle	
  		
  		//cycle on reco-particles
  		for(i=0;i<nrec;i++){
  			//if it has a vertex
  			if(rcvts[i]!=-1){
  				//row index is the id number of the vertex, the content is the id number of the tracks related to that vertex
  				vtx_trk_matrix[rcvts[i]][vtx_trk_vector[rcvts[i]]]=i;
  				//array index is the id number of the vertex, the content is the number of tracks related to that vertex 
  				vtx_trk_vector[rcvts[i]]++;
  				//cycle on relations to find the related MC
  				for(n=0;n<r2mnrel;n++){
  					if(r2mf[n]==i){
  						break;
  						}
  					}
  				if(r2mw[n]>0.9){
  					//tracks from B/C wrongly associated to PV
  					if((selection_cc(mcpdg[mcpa0[r2mt[n]]])==0||selection_bb(mcpdg[mcpa0[r2mt[n]]])==0)&&id[rcvts[i]]==103){
  						wrong_PV++;
  						}
  					//tracks from hadronization wrongly associated to SV
  					if((abs(mcpdg[mcpa0[r2mt[n]]])==5||abs(mcpdg[mcpa0[r2mt[n]]])==4 ||abs(mcpdg[mcpa0[r2mt[n]]])==3 || abs(mcpdg[mcpa0[r2mt[n]]])==2 || abs(mcpdg[mcpa0[r2mt[n]]])==1 || abs(mcpdg[mcpa0[r2mt[n]]])==21) &&id[rcvts[i]]==104){
  						wrong_SV++;
  						}
  					}
  				}
  		 	}//reco-p cycle
  		////////////////////////////////////////////////////////////////////////////////////////
  		from_hh_wrong_to_PV->Fill(wrong_PV);////////////// PLOT 9///////////////////////////////
  		////////////////////////////////////////////////////////////////////////////////////////
  		from_IP_wrong_to_SV->Fill(wrong_SV);////////////// PLOT 8 //////////////////////////////
  		////////////////////////////////////////////////////////////////////////////////////////
  		
  		wrong_PV=0;
  		wrong_SV=0;
  		
  		for(k=0;k<nvt;k++){
  			if(id[k]==103){
  				////////////////////////////////////////////////////////////////////////////////////////
  				num_trk_forPV->Fill(vtx_trk_vector[k]); ////////////// PLOT 6 //////////////////////////
  				////////////////////////////////////////////////////////////////////////////////////////
  				chi2_VS_num_trk_PV->Fill(vtx_trk_vector[k],chi[k]); ////////////// PLOT 7 //////////////
  				////////////////////////////////////////////////////////////////////////////////////////
  				}
  			else if(id[k]==104){
  				num_trk_forSV->Fill(vtx_trk_vector[k]);
  				chi2_VS_num_trk_SV->Fill(vtx_trk_vector[k],chi[k]);
  				//cycle on number of tracks associated to the vertex
  				for(n=0;n<vtx_trk_vector[k];n++){
  					reco.SetPxPyPzE(rcmox[vtx_trk_matrix[k][n]],rcmoy[vtx_trk_matrix[k][n]],rcmoz[vtx_trk_matrix[k][n]],rcene[vtx_trk_matrix[k][n]]);
  					vtx_momentum=vtx_momentum+reco;
  					}
  				//SV flight direction
	  	 		SV.SetX(x[k]);
	  	 		SV.SetY(y[k]);
	  	 		SV.SetZ(z[k]);
	  	 		flight_dir=SV-PV;
	  	 		mom=sqrt(pow(vtx_momentum.Px(),2)+pow(vtx_momentum.Py(),2)+pow(vtx_momentum.Pz(),2));
	  	 		theta=flight_dir.Angle(vtx_momentum.Vect());
	  	 		corrected_mass=sqrt(pow(vtx_momentum.M(),2)+pow(mom*sin(theta),2))+mom*sin(theta);
	  	 		////////////////////////////////////////////////////////////////////////////////////////
  				vtx_invariant_mass->Fill(vtx_momentum.M());/////////////// PLOT 5 //////////////////////
  				vtx_mom->Fill(mom); ////////////////////////////////////////////////////////////////////
  				vtx_theta->Fill(theta);  ///////////////////////////////////////////////////////////////
  				vtx_corrected_mass->Fill(corrected_mass);///////////////////////////////////////////////
  				////////////////////////////////////////////////////////////////////////////////////////
  				
  				vtx_momentum.SetPxPyPzE(0.,0.,0.,0.);  
  				
  				//cycle on number of tracks associated to the vertex
  				for(n=0;n<vtx_trk_vector[k];n++){
  					//cycle on relations to find the related MC
  					for(m=0;m<r2mnrel;m++){
  						if(r2mf[m]==vtx_trk_matrix[k][n]){
  							break;
  							}
  						}
  					//if the parent is a B/C hadron
  					if(selection_cc(mcpdg[mcpa0[r2mt[m]]])==0){
  						num_trk_C++;
  						}
  					else if(selection_bb(mcpdg[mcpa0[r2mt[m]]])==0){
  						num_trk_B++;
  						}
  					}
  				
  				if(num_trk_B>0){
  					num_SV_B++;
  					}
  				if(num_trk_C>0){
  					num_SV_C++;
  					}
  				if(num_trk_C>0 && num_trk_B>0){
  					num_SV_BC++;
  					}
  				if(num_trk_C==0 && num_trk_B==0){
  					num_SV_noBC++;
  					}
  				num_trk_C=0;
  				num_trk_B=0;
  				}
  			}//vertex cycle
  		////////////////////////////////////////////////////////////////////////////////////
  		SV_reco_B->Fill(num_SV_B); /////////////////////////////////////////////////////////
  		SV_reco_C->Fill(num_SV_C); /////////////////////// PLOT 10 /////////////////////////
  		SV_reco_BC->Fill(num_SV_BC); ///////////////////////////////////////////////////////
  		SV_reco_noBC->Fill(num_SV_noBC);  //////////////////////////////////////////////////
  		////////////////////////////////////////////////////////////////////////////////////
  		num_SV_C=0;num_SV_B=0;num_SV_BC=0;num_SV_noBC=0;
  		
  		if(ientry<3){
  			cout<<"VERTICI SECONDARI RICOSTRUITI NELL'EVENTO: "<<ientry<<endl;
  			for(k=0;k<nvt;k++){
  				if(id[k]==104){
  				cout<<""<<endl;
  				cout<<"Vertice secondario con id number "<<k<<" di coordinate: X= "<<x[k]<<" Y: "<<y[k]<<" Z: "<<z[k]<<endl;
  				cout <<"Viene ricostruito con "<<vtx_trk_vector[k]<<" tracce: "<<endl;
  					for(n=0;n<vtx_trk_vector[k];n++){
  						i=0;
						while( r2mf[i]!=vtx_trk_matrix[k][n] ){i++;}
  						cout<< "Particella: "<<mcpdg[r2mt[i]]<<"(id number "<<r2mf[i]<<") che ha come parent: "<<mcpdg[mcpa0[r2mt[i]]]<<endl;
  						}
  					}
  				}
  			//TEST//
  			for(k=0; k<nrec;k++){
				//se la particella ha un vertice secondario associato
  				if(rcvts[k]!=-1&&id[rcvts[k]]==104){
  					cout<< "La particella con id number "<<k<<"proviene dal vertice (id number "<<rcvts[k]<<")  di tipo "<<id[rcvts[k]]<<"e di coordinate X: "<<x[rcvts[k]]<<" Y: "<<y[rcvts[k]]<<" Z: "<< z[rcvts[k]] <<endl;
  					cout<<" "<<endl;
  					}
  				}
  					
  			cout<<"FINE EVENTO"<<endl;
  			}
  			
  		for(k=0;k<100;k++){ vtx_trk_vector[k]=0; for(n=0;n<100;n++){vtx_trk_matrix[k][n]=0;}}
  			
  		}//entry cycle
  		
  	TCanvas *c1=new TCanvas();
  	SV_reco_perevent->Scale(1./SV_reco_perevent->Integral());
  	SV_reco_perevent->Draw("s hist");	
  	
  	TCanvas *c2=new TCanvas();
  	SV_reco_perjet->Scale(1./SV_reco_perjet->Integral());
  	SV_reco_perjet->Draw("s hist");	
  	
  	TCanvas *c3=new TCanvas();
  	SV_to_PV_distance->Scale(1./SV_to_PV_distance->Integral());
  	SV_to_PV_distance->Draw("s hist");	
  	
  	TCanvas *c4=new TCanvas();
  	deltaR_SV_nearestJet->Scale(1./deltaR_SV_nearestJet->Integral());
  	deltaR_SV_nearestJet->Draw("s hist");
  	
  	TCanvas *c5=new TCanvas();
  	num_trk_forPV->Scale(1./num_trk_forPV->Integral());
  	num_trk_forPV->Draw("s hist");
  	
  	TCanvas *c6=new TCanvas();
  	num_trk_forSV->Scale(1./num_trk_forSV->Integral());
  	num_trk_forSV->Draw("s hist");
  	
  	TCanvas *c7=new TCanvas();
  	chi2_VS_num_trk_PV->Draw("colz");
  	
  	TCanvas *c8=new TCanvas();
  	chi2_VS_num_trk_SV->Draw("colz");
  	
  	TCanvas *c9=new TCanvas();
  	vtx_invariant_mass->Scale(1./vtx_invariant_mass->Integral());
  	vtx_invariant_mass->Draw("s hist");
  	
  	TCanvas *c10=new TCanvas();
  	vtx_corrected_mass->Scale(1./vtx_corrected_mass->Integral());
  	vtx_corrected_mass->Draw("s hist");
  	
  	TCanvas *c11=new TCanvas();
  	vtx_theta->Scale(1./vtx_theta->Integral());
  	vtx_theta->Draw("s hist");
  	
  	TCanvas *c12=new TCanvas();
  	vtx_mom->Scale(1./vtx_mom->Integral());
  	vtx_mom->Draw("s hist");
  	
  	TCanvas *c13=new TCanvas();
  	from_hh_wrong_to_PV->Scale(1./from_hh_wrong_to_PV->Integral());
  	from_hh_wrong_to_PV->Draw("s hist");
  	
  	TCanvas *c14=new TCanvas();
  	from_IP_wrong_to_SV->Scale(1./from_IP_wrong_to_SV->Integral());
  	from_IP_wrong_to_SV->Draw("s hist");
  	
  	TCanvas *c15=new TCanvas();
  	
  	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
	legend2->AddEntry(SV_reco_B, "#splitline{SV per event with at least}{a track from B decay}", "l");
	legend2->AddEntry(SV_reco_C, "#splitline{SV per event with at least}{a track from C decay}" , "l");
	//legend2->AddEntry(SV_reco_BC, "#splitline{Number of SV per event with}{at least a track from both B and C decay}" , "l");
	legend2->AddEntry(SV_reco_noBC, "#splitline{SV per event with}{no tracks from B/C decay}" , "l");
	legend2->SetTextSize(0.04);
	
  	SV_reco_B->Scale(1./SV_reco_B->Integral());
  	SV_reco_B->SetTitle("Number of different types of SV per event");
  	SV_reco_B->Draw("s hist");
  	SV_reco_C->Scale(1./SV_reco_C->Integral());
  	SV_reco_C->Draw("sames hist");
  	//SV_reco_BC->Scale(1./SV_reco_BC->Integral());
  	//SV_reco_BC->Draw("sames hist");
  	SV_reco_noBC->Scale(1./SV_reco_noBC->Integral());
  	SV_reco_noBC->Draw("sames hist");
  	legend2->Draw();
  	
  		
	free(mcvtx);free(mcvty);free(mcvtz);free(mcpdg);free(mcpa0);free(r2mt);free(r2mf);free(rcftr);free(r2mw);free(rcvts);
	free(rcmox);free(rcmoy);free(rcmoz);free(rcene);free(trk_z0);free(trk_d0);free(trk_atIP);	
  		
	}
