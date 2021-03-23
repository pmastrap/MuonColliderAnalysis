////////DESCRIPTION///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																    //
// macro to compare reco_tracks that at gen level comes from quarks/gluons to reco_tracks that at gen level comes from heavy hadrons//
// plot of D0,Z0,D0err,Z0err,pT,radius of innermost hit										    //
// these plots are important to choose the values of the parameters needed to fit primary and secondary verteces		    //
//																    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TLorentzVector.h"
#include "My_Beautiful_Canvas_TrkVar.C"

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


void compare_tracks(){

	//gStyle->SetOptStat("nemrou");
	
	TString extraText = "Simulation"; 
	TString lumi_1_5TeV = "#sqrt{s}=1.5 TeV, 500 fb^{-1}";
	TString softwareText="v02-06-MC";

	long int ientry,n,h,l;
	Float_t sign_z0,sign_d0,trk_theta,trk_eta,distance,vtx_higgs=0,vty_higgs=0,vtz_higgs=0,delta_R=0.,Lxy=0.,d0_true;
	Int_t nrec,nrel,ntrk,t2mnrel,nmcp;

	TLorentzVector reco,mc;
	
	//mc particle
	int num=2000000;
	Float_t *mcvtx,*mcvty,*mcvtz,*mcmox,*mcmoy,*mcmoz,*mcene,*mcepx,*mcepy,*mcepz;
	Int_t *mcpdg,*mcpa,*mcda0,*mcda1;
	mcpdg = (int*) malloc(sizeof(int)*num);
	mcpa = (int*) malloc(sizeof(int)*num);
	mcda0 = (int*) malloc(sizeof(int)*num);
	mcda1 = (int*) malloc(sizeof(int)*num);
	mcvtx = (float*) malloc(sizeof(float)*num);
	mcvty = (float*) malloc(sizeof(float)*num);
	mcvtz = (float*) malloc(sizeof(float)*num);
	mcepx = (float*) malloc(sizeof(float)*num);
	mcepy = (float*) malloc(sizeof(float)*num);
	mcepz = (float*) malloc(sizeof(float)*num);
	mcmox = (float*) malloc(sizeof(float)*num);
	mcmoy = (float*) malloc(sizeof(float)*num);
	mcmoz = (float*) malloc(sizeof(float)*num);
	mcene = (float*) malloc(sizeof(float)*num);
	
	//relation
	int num2=100000;
	Int_t *r2mt,*r2mf,*rcftr;
	Float_t *r2mw;
	r2mt = (int*) malloc(sizeof(int)*num2);
	r2mf = (int*) malloc(sizeof(int)*num2);
	r2mw = (float*) malloc(sizeof(float)*num2);
	rcftr = (int*) malloc(sizeof(int)*num2);
	
	//reco particle
	Float_t *rcmox,*rcmoy,*rcmoz,*ene;
	Int_t *rctyp;
	rctyp=(int*) malloc(sizeof(int)*num2);
	rcmox = (float*) malloc(sizeof(float)*num2);
	rcmoy = (float*) malloc(sizeof(float)*num2);
	rcmoz = (float*) malloc(sizeof(float)*num2);
	ene = (float*) malloc(sizeof(float)*num2);
	
	//tracks
	Float_t *trk_z0,*trk_d0,*trrih;
	Int_t *trk_atIP,*trthn,*thori;
	trk_z0 = (float*) malloc(sizeof(float)*num2);
	trk_d0 = (float*) malloc(sizeof(float)*num2);
	trrih = (float*) malloc(sizeof(float)*num2);
	trk_atIP = (int*) malloc(sizeof(int)*num2);
	trthn = (int*) malloc(sizeof(int)*num2);
	Int_t trshn[1000][12];
	int num3=5000;
	thori=(int*) malloc(sizeof(int)*num3);
	Float_t cov[100000][15];
	Int_t trthi[5000][100];

	TChain* fChain = new TChain("fChain");
   	fChain->Add("/home/paola/Scrivania/xml_aggiornati/Fixed_version/10k/ntuple_bb_10000_fix.root/MyLCTuple");
	
	fChain->SetBranchAddress("nmcp", &nmcp);
	fChain->SetBranchAddress("mcpdg", mcpdg);
	fChain->SetBranchAddress("mcpa0", mcpa);
	fChain->SetBranchAddress("mcda0", mcda0);
	fChain->SetBranchAddress("mcda1", mcda1);
	fChain->SetBranchAddress("mcvtx", mcvtx);
	fChain->SetBranchAddress("mcvty", mcvty);
	fChain->SetBranchAddress("mcvtz", mcvtz);
	fChain->SetBranchAddress("mcepx", mcepx);
	fChain->SetBranchAddress("mcepy", mcepy);
	fChain->SetBranchAddress("mcepz", mcepz);
	fChain->SetBranchAddress("mcmox", mcmox);
  	fChain->SetBranchAddress("mcmoy", mcmoy);
  	fChain->SetBranchAddress("mcmoz", mcmoz);
	fChain->SetBranchAddress("mcene", mcene);
	
	fChain->SetBranchAddress("nrec", &nrec);  
	fChain->SetBranchAddress("rctyp", rctyp);  
	fChain->SetBranchAddress("r2mnrel", &nrel);
	fChain->SetBranchAddress("r2mt", r2mt);
	fChain->SetBranchAddress("r2mf", r2mf);	
        fChain->SetBranchAddress("r2mw", r2mw);
	fChain->SetBranchAddress("rcftr", rcftr);
	fChain->SetBranchAddress("rcmox", rcmox);
  	fChain->SetBranchAddress("rcmoy", rcmoy);
  	fChain->SetBranchAddress("rcmoz", rcmoz);
	fChain->SetBranchAddress("rcene", ene);
	
	fChain->SetBranchAddress("ntrk",&ntrk);
	fChain->SetBranchAddress("tszze",trk_z0);
	fChain->SetBranchAddress("trsip", trk_atIP);
        fChain->SetBranchAddress("tsdze", trk_d0);
        fChain->SetBranchAddress("tscov", cov);
        fChain->SetBranchAddress("trrih",trrih);
        fChain->SetBranchAddress("trthn",trthn);
        fChain->SetBranchAddress("trshn",trshn);
        fChain->SetBranchAddress("trthi",trthi);
        
        fChain->SetBranchAddress("thori",thori);
       
       
        
        TH1F* radiusInnerHit=new TH1F("radiusInnerHit_from_hh","radiusInnerHit",102,29,80);
	radiusInnerHit->SetTitle("Radius of the innermost hit");
	radiusInnerHit->SetLineWidth(2);
    	radiusInnerHit->GetXaxis()->SetTitle("r [mm]");
    	
    	TH1F* radiusInnerHit_o=new TH1F("radiusInnerHit_from_q/g","radiusInnerHit_o",102,29,80);
	radiusInnerHit_o->SetTitle("Radius of the innermost hit");
	radiusInnerHit_o->SetLineWidth(2);
	radiusInnerHit_o->SetLineColor(kRed);
    	radiusInnerHit_o->GetXaxis()->SetTitle("r [mm]");

	TH1F* pt=new TH1F("pt_from_hh","pt",80,0.,20.);
	pt->SetLineWidth(2);
        pt->SetTitle("pT");
    	pt->GetXaxis()->SetTitle("pt [GeV]");
	
	TH1F* pt_o=new TH1F("pt_from_q/g","pt_o",80,0.,20.);
	pt_o->SetLineWidth(2);
	pt_o->SetLineColor(kRed);
	pt_o->SetTitle("pT");
    	pt_o->GetXaxis()->SetTitle("pt [GeV]");
        
        TH1F* reco_Z0=new TH1F("reco_Z0_from_hh","reco_Z0",80,-2,2);
	reco_Z0->SetTitle("Longitudinal impact parameter");
	reco_Z0->SetLineWidth(2);
    	reco_Z0->GetXaxis()->SetTitle("Z0 [mm]");
    	
    	TH1F* reco_D0=new TH1F("reco_D0_from_hh","reco_D0",80,-2,2);
	reco_D0->SetTitle("Transverse impact parameter");
	reco_D0->SetLineWidth(2);
    	reco_D0->GetXaxis()->SetTitle("D0 [mm]");

	TH1F* reco_Z0_o=new TH1F("reco_Z0_from_q/g","reco_Z0_o",80,-2,2);
	reco_Z0_o->SetTitle("Longitudinal impact parameter");
	reco_Z0_o->SetLineWidth(2);
	reco_Z0_o->SetLineColor(kRed);
    	reco_Z0_o->GetXaxis()->SetTitle("Z0 [mm]");
    	
    	TH1F* reco_D0_o=new TH1F("reco_D0_from_q/g","reco_D0_o",80,-2,2);
	reco_D0_o->SetTitle("Transverse impact parameter");
	reco_D0_o->SetLineWidth(2);
	reco_D0_o->SetLineColor(kRed);
    	reco_D0_o->GetXaxis()->SetTitle("D0 [mm]");
    	
    	TH1F* err_Z0=new TH1F("err_Z0_from_hh","err_Z0",50,0,0.5);
	err_Z0->SetTitle("Z0 error");
	err_Z0->SetLineWidth(2);
    	err_Z0->GetXaxis()->SetTitle("mm");
    	
    	TH1F* err_D0=new TH1F("err_D0_from_hh","err_D0",50,0,0.5);
	err_D0->SetTitle("D0 error");
	err_D0->SetLineWidth(2);
    	err_D0->GetXaxis()->SetTitle("mm");

	TH1F* err_Z0_o=new TH1F("err_Z0_from_q/g","err_Z0_o",50,0,0.5);
	err_Z0_o->SetTitle("Z0 error");
	err_Z0_o->SetLineWidth(2);
	err_Z0_o->SetLineColor(kRed);
    	err_Z0_o->GetXaxis()->SetTitle("mm");
    	
    	TH1F* err_D0_o=new TH1F("err_D0_from_q/g","err_D0_o",50,0,0.5);
	err_D0_o->SetTitle("D0 error");
	err_D0_o->SetLineWidth(2);
	err_D0_o->SetLineColor(kRed);
    	err_D0_o->GetXaxis()->SetTitle("mm");
    	
    	TH1F* reco_D0_sig=new TH1F("reco_D0_sig_from_hh","reco_D0_sig",120,-30,30);
	reco_D0_sig->SetTitle("Transverse impact parameter significance");
	reco_D0_sig->SetLineWidth(2);
    	reco_D0_sig->GetXaxis()->SetTitle("Significance");
    	
    	TH1F* reco_D0_sig_o=new TH1F("reco_D0_sig_from_q/g","reco_D0_sig_o",120,-30,30);
	reco_D0_sig_o->SetTitle("Transverse impact parameter significance");
	reco_D0_sig_o->SetLineWidth(2);
	reco_D0_sig_o->SetLineColor(kRed);
    	reco_D0_sig_o->GetXaxis()->SetTitle("Significance");
    	
    	TH1F* reco_Z0_sig=new TH1F("reco_Z0_sig_from_hh","reco_Z0_sig",120,-30,30);
	reco_Z0_sig->SetTitle("Longitudinal impact parameter significance");
	reco_Z0_sig->SetLineWidth(2);
    	reco_Z0_sig->GetXaxis()->SetTitle("Significance");
    	
    	TH1F* reco_Z0_sig_o=new TH1F("reco_Z0_sig_from_q/g","reco_Z0_sig_o",120,-30,30);
	reco_Z0_sig_o->SetTitle("Longitudinal impact parameter significance");
	reco_Z0_sig_o->SetLineWidth(2);
	reco_Z0_sig_o->SetLineColor(kRed);
    	reco_Z0_sig_o->GetXaxis()->SetTitle("Significance");
    	
    	TH1F* NHits=new TH1F("NHits_from_hh","NHits",15,-0.5,14.5);
    	NHits->SetTitle("Number of hits (VTX+FTD) associated to the track");
    	NHits->SetLineWidth(2);
    	
    	TH1F* NHits_o=new TH1F("NHits_from_q/g","NHits_o",15,-0.5,14.5);
    	NHits_o->SetTitle("Number of hits (VTX+FTD) associated to the track");
    	NHits_o->SetLineWidth(2);
	NHits_o->SetLineColor(kRed);
    	
    	TH1F* Z0D0_sig=new TH1F("Z0D0_sig_from_hh","Z0D0_sig",120,0,60);
	Z0D0_sig->SetTitle("D0-Z0 significance");
	Z0D0_sig->SetLineWidth(2);
	Z0D0_sig->GetYaxis()->SetTitle("Events / 0.25");
    	Z0D0_sig->GetXaxis()->SetTitle("Significance");
    	
    	TH1F* Z0D0_sig_o=new TH1F("Z0D0_sig_from_q/g","Z0D0_sig_o",120,0,60);
	Z0D0_sig_o->SetTitle("D0-Z0 significance");
	Z0D0_sig_o->SetLineWidth(2);
	Z0D0_sig_o->SetLineColor(kRed);
	Z0D0_sig_o->GetYaxis()->SetTitle("Events / 0.25");
    	Z0D0_sig_o->GetXaxis()->SetTitle("Significance");
    	
	
    	int nh=0;
    	
    	//counter per tagli
    	double tot_trk_from_hh=0, tot_trk_from_qg=0;
    	double hh_excluded_from_SV=0, qg_excluded_from_SV=0;
    	double hh_excluded_from_PV=0, qg_excluded_from_PV=0;
  	
    	double zd_sig=0;
	
    	TVector3 END,START;

    	
    	//ciclo su entries
    	//fChain->GetEntries()
        for(ientry=0; ientry<fChain->GetEntries();++ientry){
  		fChain->GetEntry(ientry);	
		
		//ciclo su numero di relazioni tra reco e mc
  		for(n=0;n<nrel;n++){
  		
			reco.SetPxPyPzE(rcmox[r2mf[n]],rcmoy[r2mf[n]],rcmoz[r2mf[n]],ene[r2mf[n]]);
			mc.SetPxPyPzE(mcmox[r2mt[n]],mcmoy[r2mt[n]],mcmoz[r2mt[n]],mcene[r2mt[n]]);
			
			//se il peso della relazione reco-mc è maggiore di 0.9 e la reco particle ha una traccia associata
			if(rcftr[r2mf[n]]!=-1 && r2mw[n]>0.9){
				
				START.SetXYZ(mcvtx[r2mt[n]],mcvty[r2mt[n]],mcvtz[r2mt[n]]);
  				END.SetXYZ(mcvtx[mcpa[r2mt[n]]],mcvty[mcpa[r2mt[n]]],mcvtz[mcpa[r2mt[n]]]);
  				
				//se la particella proviene da adroni B o C
				if((selection_cc(mcpdg[mcpa[r2mt[n]]])==0 || selection_bb(mcpdg[mcpa[r2mt[n]]])==0)&&(START-END).Mag()>0){
				
				//se la particella viene dal decadimento (endpoint del padre = start point del figlio)
				if(mcepx[mcpa[r2mt[n]]]==mcvtx[r2mt[n]]&&mcepy[mcpa[r2mt[n]]]==mcvty[r2mt[n]]&& 	mcepz[mcpa[r2mt[n]]]==mcvtz[r2mt[n]]){
					
	  				//fill z0
		   			reco_Z0->Fill(trk_z0[trk_atIP[rcftr[r2mf[n]]]]);
		   			//fill d0
		   			reco_D0->Fill(trk_d0[trk_atIP[rcftr[r2mf[n]]]]);
		   			//err z0
		   			err_Z0->Fill(sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][9]));
		   			//err D0
		   			err_D0->Fill(sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][0]));
		   			
		   			//z0d0- sig
		   			zd_sig=sqrt(pow(trk_d0[trk_atIP[rcftr[r2mf[n]]]]/sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][0]),2)+
		   			pow(trk_z0[trk_atIP[rcftr[r2mf[n]]]]/sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][9]),2));
		   			
		   			Z0D0_sig->Fill(zd_sig);
		   			//pT
					pt->Fill(reco.Perp());
					//Innermost Hit radius
					radiusInnerHit->Fill(trrih[rcftr[r2mf[n]]]);
					//significance d0
					reco_D0_sig->Fill(trk_d0[trk_atIP[rcftr[r2mf[n]]]]/sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][0]));
					//significance z0
					reco_Z0_sig->Fill(trk_z0[trk_atIP[rcftr[r2mf[n]]]]/sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][9]));
					//n hits
					NHits->Fill(trshn[rcftr[r2mf[n]]][0]+trshn[rcftr[r2mf[n]]][2]);
					nh=trshn[rcftr[r2mf[n]]][0]+trshn[rcftr[r2mf[n]]][2];
					
					//tagli
					tot_trk_from_hh++;
					if(nh<4 || reco.Perp()<0.8 || zd_sig<2){
						
						hh_excluded_from_SV++;
						
						}
						
					if(nh<4 || abs(trk_d0[trk_atIP[rcftr[r2mf[n]]]]) > 0.1|| abs(trk_z0[trk_atIP[rcftr[r2mf[n]]]]) > 0.1){
						
						hh_excluded_from_PV++;
						
						}

					
	   				}
	   				}
	   			//se la particella proviene da quark/gluoni
	   			else if(START.Mag()<1e-3){
					//fill z0
		   			reco_Z0_o->Fill(trk_z0[trk_atIP[rcftr[r2mf[n]]]]);
		   			//fill d0
		   			reco_D0_o->Fill(trk_d0[trk_atIP[rcftr[r2mf[n]]]]);
		   			//err z0
					err_Z0_o->Fill(sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][9]));
					//err d0
		   			err_D0_o->Fill(sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][0]));
		   			//pT
					pt_o->Fill(reco.Perp());
					
					//z0d0- sig
		   			zd_sig=sqrt(pow(trk_d0[trk_atIP[rcftr[r2mf[n]]]]/sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][0]),2)+
		   			pow(trk_z0[trk_atIP[rcftr[r2mf[n]]]]/sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][9]),2));
		   			
		   			Z0D0_sig_o->Fill(zd_sig);
					
					//Innermost Hit radius
					radiusInnerHit_o->Fill(trrih[rcftr[r2mf[n]]]);
					//significance d0
					reco_D0_sig_o->Fill(trk_d0[trk_atIP[rcftr[r2mf[n]]]]/sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][0]));
					//significance z0
					reco_Z0_sig_o->Fill(trk_z0[trk_atIP[rcftr[r2mf[n]]]]/sqrt(cov[trk_atIP[rcftr[r2mf[n]]]][9]));
					//n hits
					NHits_o->Fill(trshn[rcftr[r2mf[n]]][0]+trshn[rcftr[r2mf[n]]][2]);
					nh=trshn[rcftr[r2mf[n]]][0]+trshn[rcftr[r2mf[n]]][2];
					
					//tagli
					tot_trk_from_qg++;
					if(nh<4 || reco.Perp()<0.8 || zd_sig<2.){
						
						qg_excluded_from_SV++;
						
						}
						
					if(nh<4 || abs(trk_d0[trk_atIP[rcftr[r2mf[n]]]]) > 0.1 || abs(trk_z0[trk_atIP[rcftr[r2mf[n]]]]) > 0.1){
						
						qg_excluded_from_PV++;
						
						}
					
					
   					}
   				}
   			}
		}
		
	cout<<"Su 1000 eventi --> Totale tracce da hh:  "<< tot_trk_from_hh <<". Totale tracce da qg: "<<tot_trk_from_qg<<endl;
	cout<<" "<<endl;
	cout<<"------------------------------------------"<<endl;
	cout<<" "<<endl;
	
	cout<<"SELEZIONE SU SV"<<endl;
	cout<<" "<<endl;
	cout<<"dal fit di SV sono escluse "<<hh_excluded_from_SV<<" ("<<hh_excluded_from_SV/tot_trk_from_hh*100<<"%) tracce da hh e "<<qg_excluded_from_SV<<" ("<<qg_excluded_from_SV/tot_trk_from_qg*100<<"%) tracce da qg"<<endl;
	
	cout<<" "<<endl;
	cout<<"SELEZIONE SU PV"<<endl;
	cout<<" "<<endl;
	cout<<"dal fit di PV sono escluse "<<hh_excluded_from_PV<<" ("<<hh_excluded_from_PV/tot_trk_from_hh*100<<"%) tracce da hh e "<<qg_excluded_from_PV<<" ("<<qg_excluded_from_PV/tot_trk_from_qg*100<<"%) tracce da qg"<<endl;
	
	
	TString titleX;
	TString titleY;
   	
   	////////////////////////////////////////////////
	titleX = "Radius of the innermost hit [mm]";
	titleY = "Tracks / 0.50 mm";
	radiusInnerHit->SetBinContent(radiusInnerHit->GetNbinsX(), radiusInnerHit->GetBinContent(radiusInnerHit->GetNbinsX()) + radiusInnerHit->GetBinContent(radiusInnerHit->GetNbinsX() + 1));
	radiusInnerHit_o->SetBinContent(radiusInnerHit_o->GetNbinsX(),radiusInnerHit_o->GetBinContent(radiusInnerHit_o->GetNbinsX())+radiusInnerHit_o->GetBinContent(radiusInnerHit_o->GetNbinsX()+1));
   	radiusInnerHit_o->Scale(1./radiusInnerHit_o->Integral());
	radiusInnerHit->Scale(1./radiusInnerHit->Integral());
   	My_Beautiful_Canvas_TrkVar(radiusInnerHit,radiusInnerHit_o,titleX,titleY);
   		
   	////////////////////////////////////////////////
	titleX = "z_{0} [mm]";
	titleY = "Tracks / 0.05 mm";
	reco_Z0->SetBinContent(reco_Z0->GetNbinsX(), reco_Z0->GetBinContent(reco_Z0->GetNbinsX()) + reco_Z0->GetBinContent(reco_Z0->GetNbinsX() + 1));
	reco_Z0->SetBinContent(1, reco_Z0->GetBinContent(1) + reco_Z0->GetBinContent(0));
	reco_Z0_o->SetBinContent(reco_Z0_o->GetNbinsX(),reco_Z0_o->GetBinContent(reco_Z0_o->GetNbinsX())+reco_Z0_o->GetBinContent(reco_Z0_o->GetNbinsX()+1));
	reco_Z0_o->SetBinContent(1, reco_Z0_o->GetBinContent(1) + reco_Z0_o->GetBinContent(0));
   	reco_Z0->Scale(1./reco_Z0->Integral());
	reco_Z0_o->Scale(1./reco_Z0_o->Integral());
	My_Beautiful_Canvas_TrkVar(reco_Z0,reco_Z0_o,titleX,titleY);

	
	////////////////////////////////////////////////
	titleX = "z_{0} significance";
	titleY = "Tracks / 0.50 units";
	reco_Z0_sig->SetBinContent(reco_Z0_sig->GetNbinsX(), reco_Z0_sig->GetBinContent(reco_Z0_sig->GetNbinsX()) + reco_Z0_sig->GetBinContent(reco_Z0_sig->GetNbinsX() + 1));
	reco_Z0_sig->SetBinContent(1, reco_Z0_sig->GetBinContent(1) + reco_Z0_sig->GetBinContent(0));
	reco_Z0_sig_o->SetBinContent(reco_Z0_sig_o->GetNbinsX(),reco_Z0_sig_o->GetBinContent(reco_Z0_sig_o->GetNbinsX())+reco_Z0_sig_o->GetBinContent(reco_Z0_sig_o->GetNbinsX()+1));
	reco_Z0_sig_o->SetBinContent(1, reco_Z0_sig_o->GetBinContent(1) + reco_Z0_sig_o->GetBinContent(0));
   	reco_Z0_sig->Scale(1./reco_Z0_sig->Integral());
	reco_Z0_sig_o->Scale(1./reco_Z0_sig_o->Integral());
	My_Beautiful_Canvas_TrkVar(reco_Z0_sig,reco_Z0_sig_o,titleX,titleY);
   	
   	////////////////////////////////////////////////
	titleX = "d_{0} [mm]";
	titleY = "Tracks / 0.05 mm";
	reco_D0->SetBinContent(reco_D0->GetNbinsX(), reco_D0->GetBinContent(reco_D0->GetNbinsX()) + reco_D0->GetBinContent(reco_D0->GetNbinsX() + 1));
	reco_D0->SetBinContent(1, reco_D0->GetBinContent(1) + reco_D0->GetBinContent(0));
	reco_D0_o->SetBinContent(reco_D0_o->GetNbinsX(),reco_D0_o->GetBinContent(reco_D0_o->GetNbinsX())+reco_D0_o->GetBinContent(reco_D0_o->GetNbinsX()+1));
	reco_D0_o->SetBinContent(1, reco_D0_o->GetBinContent(1) + reco_D0_o->GetBinContent(0));
   	reco_D0->Scale(1./reco_D0->Integral());
	reco_D0_o->Scale(1./reco_D0_o->Integral());
	My_Beautiful_Canvas_TrkVar(reco_D0,reco_D0_o,titleX,titleY);

	
	////////////////////////////////////////////////
	titleX = "d_{0} significance";
	titleY = "Tracks / 0.50 units";
	
	reco_D0_sig->SetBinContent(reco_D0_sig->GetNbinsX(), reco_D0_sig->GetBinContent(reco_D0_sig->GetNbinsX()) + reco_D0_sig->GetBinContent(reco_D0_sig->GetNbinsX() + 1));
	reco_D0_sig->SetBinContent(1, reco_D0_sig->GetBinContent(1) + reco_D0_sig->GetBinContent(0));
	reco_D0_sig_o->SetBinContent(reco_D0_sig_o->GetNbinsX(),reco_D0_sig_o->GetBinContent(reco_D0_sig_o->GetNbinsX())+reco_D0_sig_o->GetBinContent(reco_D0_sig_o->GetNbinsX()+1));
	reco_D0_sig_o->SetBinContent(1, reco_D0_sig_o->GetBinContent(1) + reco_D0_sig_o->GetBinContent(0));
   	reco_D0_sig->Scale(1./reco_D0_sig->Integral());
	reco_D0_sig_o->Scale(1./reco_D0_sig_o->Integral());
	My_Beautiful_Canvas_TrkVar(reco_D0_sig,reco_D0_sig_o,titleX,titleY);

   	/////////////////////////////////////////////////
	titleX = "#sigma_{z_{0}} [mm]";
	titleY = "Tracks / 0.01 mm";
	
	err_Z0->SetBinContent(err_Z0->GetNbinsX(), err_Z0->GetBinContent(err_Z0->GetNbinsX()) + err_Z0->GetBinContent(err_Z0->GetNbinsX() + 1));
	err_Z0_o->SetBinContent(err_Z0_o->GetNbinsX(),err_Z0_o->GetBinContent(err_Z0_o->GetNbinsX())+err_Z0_o->GetBinContent(err_Z0_o->GetNbinsX()+1));
   	err_Z0->Scale(1./err_Z0->Integral());
	err_Z0_o->Scale(1./err_Z0_o->Integral());
	My_Beautiful_Canvas_TrkVar(err_Z0,err_Z0_o,titleX,titleY);
	
   	/////////////////////////////////////////////////
	titleX = "#sigma_{d_{0}} [mm]";
	titleY = "Tracks / 0.01 mm";
	
	err_D0->SetBinContent(err_D0->GetNbinsX(), err_D0->GetBinContent(err_D0->GetNbinsX()) + err_D0->GetBinContent(err_D0->GetNbinsX() + 1));
	err_D0_o->SetBinContent(err_D0_o->GetNbinsX(),err_D0_o->GetBinContent(err_D0_o->GetNbinsX())+err_D0_o->GetBinContent(err_D0_o->GetNbinsX()+1));
   	err_D0->Scale(1./err_D0->Integral());
	err_D0_o->Scale(1./err_D0_o->Integral());
	My_Beautiful_Canvas_TrkVar(err_D0,err_D0_o,titleX,titleY);
 
	/////////////////////////////////////////////////
	titleX = "p_{T} [GeV]";
	titleY = "Tracks / 0.25 GeV";
	
	pt->SetBinContent(pt->GetNbinsX(), pt->GetBinContent(pt->GetNbinsX()) + pt->GetBinContent(pt->GetNbinsX() + 1));
	pt_o->SetBinContent(pt_o->GetNbinsX(),pt_o->GetBinContent(pt_o->GetNbinsX())+pt_o->GetBinContent(pt_o->GetNbinsX()+1));
   	pt->Scale(1./pt->Integral());
	pt_o->Scale(1./pt_o->Integral());
 	My_Beautiful_Canvas_TrkVar(pt,pt_o,titleX,titleY);
   	
   	////////////////////////////////////////////////
	titleX = "Number of hits in vertex detector(B+E)";
	titleY = "Tracks / 1 units";
	
	NHits->SetBinContent(NHits->GetNbinsX(), NHits->GetBinContent(NHits->GetNbinsX()) + NHits->GetBinContent(NHits->GetNbinsX() + 1));
	NHits_o->SetBinContent(NHits_o->GetNbinsX(),NHits_o->GetBinContent(NHits_o->GetNbinsX())+NHits_o->GetBinContent(NHits_o->GetNbinsX()+1));
	NHits->Scale(1./NHits->Integral());
	NHits_o->Scale(1./NHits_o->Integral());
	My_Beautiful_Canvas_TrkVar(NHits,NHits_o,titleX,titleY);
	
	
	///////////////////////////////////////////////
	titleX = "d_{0}-z_{0} significance";
	titleY = "Tracks / 0.50 units";
	
	Z0D0_sig->SetBinContent(Z0D0_sig->GetNbinsX(), Z0D0_sig->GetBinContent(Z0D0_sig->GetNbinsX()) + Z0D0_sig->GetBinContent(Z0D0_sig->GetNbinsX() + 1));
	Z0D0_sig_o->SetBinContent(Z0D0_sig_o->GetNbinsX(), Z0D0_sig_o->GetBinContent(Z0D0_sig_o->GetNbinsX()) + Z0D0_sig_o->GetBinContent(Z0D0_sig_o->GetNbinsX() + 1));
	Z0D0_sig->Scale(1./Z0D0_sig->Integral());
	Z0D0_sig_o->Scale(1./Z0D0_sig_o->Integral());
	My_Beautiful_Canvas_TrkVar(Z0D0_sig,Z0D0_sig_o,titleX,titleY);
	
	free(mcvtx);free(mcvty);free(mcvtz);free(mcpdg);free(mcpa);free(r2mt);free(r2mf);free(rcftr);free(r2mw);
	free(rcmox);free(rcmoy);free(rcmoz);free(trk_z0);free(ene);free(trk_d0);free(trrih);free(trk_atIP);
	
	}

