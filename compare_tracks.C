////////DESCRIPTION///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																    //
// macro to compare reco_tracks that at gen level comes from quarks/gluons to reco_tracks that at gen level comes from heavy hadrons//
// plot of D0,Z0,D0err,Z0err,pT,radius of innermost hit, distance of the particle vertex from Higgs decay vertex		    //
// additional plot of distance of the Higgs decay vertex from IP (gen level)							    //
// these plots are important to choose the values of the parameters needed to fit primary and secondary verteces		    //
//																    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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


void compare_tracks(){

	gStyle->SetOptStat("nemrou");
	
	long int ientry,n;
	Float_t sign_z0,sign_d0,delta_R,trk_theta,trk_eta,distance,vtx_higgs=0,vty_higgs=0,vtz_higgs=0;
	Int_t nrec,nrel,ntrk,t2mnrel,nmcp;

	TLorentzVector reco;
	
	//mc particle
	int num=2000000;
	Float_t *mcvtx,*mcvty,*mcvtz;
	Int_t *mcpdg,*mcpa;
	mcpdg = (int*) malloc(sizeof(int)*num);
	mcpa = (int*) malloc(sizeof(int)*num);
	mcvtx = (float*) malloc(sizeof(float)*num);
	mcvty = (float*) malloc(sizeof(float)*num);
	mcvtz = (float*) malloc(sizeof(float)*num);
	
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
	rcmox = (float*) malloc(sizeof(float)*num2);
	rcmoy = (float*) malloc(sizeof(float)*num2);
	rcmoz = (float*) malloc(sizeof(float)*num2);
	ene = (float*) malloc(sizeof(float)*num2);
	
	//tracks
	Float_t *trk_z0,*trk_d0,*trrih;
	Int_t *trk_atIP;
	trk_z0 = (float*) malloc(sizeof(float)*num2);
	trk_d0 = (float*) malloc(sizeof(float)*num2);
	trrih = (float*) malloc(sizeof(float)*num2);
	trk_atIP = (int*) malloc(sizeof(int)*num2);
	Float_t cov[100000][15];

	TChain* fChain = new TChain("fChain");
   	fChain->Add("/home/paola/Scrivania/xml_aggiornati/ntuple_bb_1000_newgeom.root/MyLCTuple");
	
	fChain->SetBranchAddress("nmcp", &nmcp);
	fChain->SetBranchAddress("mcpdg", mcpdg);
	fChain->SetBranchAddress("mcpa0", mcpa);
	fChain->SetBranchAddress("mcvtx", mcvtx);
	fChain->SetBranchAddress("mcvty", mcvty);
	fChain->SetBranchAddress("mcvtz", mcvtz);
	
	fChain->SetBranchAddress("nrec", &nrec);   
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
        
        TH1F* radiusInnerHit=new TH1F("radiusInnerHit_from_hh","radiusInnerHit",100,29,35);
	radiusInnerHit->SetTitle("Radius of the innermost hit");
	radiusInnerHit->SetLineWidth(2);
    	radiusInnerHit->GetXaxis()->SetTitle("r [mm]");
    	
    	TH1F* radiusInnerHit_o=new TH1F("radiusInnerHit_from_q/g","radiusInnerHit_o",100,29,35);
	radiusInnerHit_o->SetTitle("Radius of the innermost hit");
	radiusInnerHit_o->SetLineWidth(2);
	radiusInnerHit_o->SetLineColor(kRed);
    	radiusInnerHit_o->GetXaxis()->SetTitle("r [mm]");

	TH1F* pt=new TH1F("pt_from_hh","pt",100,0.,50.);
	pt->SetLineWidth(2);
        pt->SetTitle("pT");
    	pt->GetXaxis()->SetTitle("pt [GeV]");
	
	TH1F* pt_o=new TH1F("pt_from_q/g","pt_o",100,0.,50.);
	pt_o->SetLineWidth(2);
	pt_o->SetLineColor(kRed);
	pt_o->SetTitle("pT");
    	pt_o->GetXaxis()->SetTitle("pt [GeV]");
        
        TH1F* reco_Z0=new TH1F("reco_Z0_from_hh","reco_Z0",100,-0.2,0.2);
	reco_Z0->SetTitle("Longitudinal impact parameter");
	reco_Z0->SetLineWidth(2);
    	reco_Z0->GetXaxis()->SetTitle("Z0 [mm]");
    	
    	TH1F* reco_D0=new TH1F("reco_D0_from_hh","reco_D0",100,-0.2,0.2);
	reco_D0->SetTitle("Transverse impact parameter");
	reco_D0->SetLineWidth(2);
    	reco_D0->GetXaxis()->SetTitle("D0 [mm]");

	TH1F* reco_Z0_o=new TH1F("reco_Z0_from_q/g","reco_Z0_o",100,-0.2,0.2);
	reco_Z0_o->SetTitle("Longitudinal impact parameter");
	reco_Z0_o->SetLineWidth(2);
	reco_Z0_o->SetLineColor(kRed);
    	reco_Z0_o->GetXaxis()->SetTitle("Z0 [mm]");
    	
    	TH1F* reco_D0_o=new TH1F("reco_D0_from_q/g","reco_D0_o",100,-0.2,0.2);
	reco_D0_o->SetTitle("Transverse impact parameter");
	reco_D0_o->SetLineWidth(2);
	reco_D0_o->SetLineColor(kRed);
    	reco_D0_o->GetXaxis()->SetTitle("D0 [mm]");
    	
    	TH1F* err_Z0=new TH1F("err_Z0_from_hh","err_Z0",50,0,0.2);
	err_Z0->SetTitle("Z0 error");
	err_Z0->SetLineWidth(2);
    	err_Z0->GetXaxis()->SetTitle("mm");
    	
    	TH1F* err_D0=new TH1F("err_D0_from_hh","err_D0",50,0,0.2);
	err_D0->SetTitle("D0 error");
	err_D0->SetLineWidth(2);
    	err_D0->GetXaxis()->SetTitle("mm");

	TH1F* err_Z0_o=new TH1F("err_Z0_from_q/g","err_Z0_o",50,0,0.2);
	err_Z0_o->SetTitle("Z0 error");
	err_Z0_o->SetLineWidth(2);
	err_Z0_o->SetLineColor(kRed);
    	err_Z0_o->GetXaxis()->SetTitle("mm");
    	
    	TH1F* err_D0_o=new TH1F("err_D0_from_q/g","err_D0_o",50,0,0.2);
	err_D0_o->SetTitle("D0 error");
	err_D0_o->SetLineWidth(2);
	err_D0_o->SetLineColor(kRed);
    	err_D0_o->GetXaxis()->SetTitle("mm");
    	
    	TH2F* colz=new TH2F("colz","colz",100,-0.5,0.5,100,-0.3,0.3);  	
    	
    	TH1F* dist=new TH1F("dist_hh","dist",100,0.,0.1);
	dist->SetTitle("particle vertex distance from IP");
	dist->SetLineWidth(2);
    	dist->GetXaxis()->SetTitle("dist[mm]");
    	
    	TH1F* dist_o=new TH1F("dist_q/g","dist_o",100,0.,0.1);
	dist_o->SetTitle("particle vertex distance from IP");
	dist_o->SetLineWidth(2);
	dist_o->SetLineColor(kRed);
    	dist_o->GetXaxis()->SetTitle("dist[mm]");
    	
    	TH1F* higgs_decay=new TH1F("higgs_decay","higgs_decay",100,0.,10e-9);
    	
    	//ciclo su entries
        for(ientry=0; ientry<fChain->GetEntries() ;++ientry){
  		fChain->GetEntry(ientry);
		//trova higgs
		n=0;
		while( abs(mcpdg[mcpa[n]])!=25 && n<nmcp ){n++;}
		if(n<nmcp){
			cout<<"Vertice decadimento Higgs: X="<<mcvtx[n]<<" Y="<<mcvty[n]<<" Z="<<mcvtz[n]<<endl;
			//cout<<"figlio Higgs: "<<mcpdg[n]<<endl;
			//cout<<"Higgs: "<<mcpdg[mcpa[n]]<<endl;
			vtx_higgs=mcvtx[n];
			vty_higgs=mcvty[n];
			vtz_higgs=mcvtz[n];
			higgs_decay->Fill(sqrt(vtx_higgs*vtx_higgs+vty_higgs*vty_higgs+vtz_higgs*vtz_higgs));
			}		
		cout<<"FINE EVENTO"<<endl;
		//ciclo su numero di relazioni tra reco e mc
  		for(n=0;n<nrel;n++){
			reco.SetPx(rcmox[r2mf[n]]);
			reco.SetPy(rcmoy[r2mf[n]]);
			reco.SetPz(rcmoz[r2mf[n]]);
			reco.SetE(ene[r2mf[n]]);
			//se il peso della relazione reco-mc è maggiore di 0.9 e la reco particle ha una traccia associata
			if(rcftr[r2mf[n]]!=-1 && r2mw[n]>0.9 && reco.M()>0){
				distance=sqrt((mcvtx[r2mt[n]]-vtx_higgs)*(mcvtx[r2mt[n]]-vtx_higgs)+(mcvty[r2mt[n]]-vty_higgs)*(mcvty[r2mt[n]]-vty_higgs)+(mcvtz[r2mt[n]]-vtz_higgs)*(mcvtz[r2mt[n]]-vtz_higgs));
				//se la particella proviene da adroni B o C
				if(selection_cc(mcpdg[mcpa[r2mt[n]]])==0 || selection_bb(mcpdg[mcpa[r2mt[n]]])==0){
					//fill dist
					dist->Fill(distance);
	  				//fill z0
		   			reco_Z0->Fill(trk_z0[trk_atIP[rcftr[r2mf[n]]]]);
		   			//fill d0
		   			reco_D0->Fill(trk_d0[trk_atIP[rcftr[r2mf[n]]]]);
		   			colz->Fill(trk_z0[trk_atIP[rcftr[r2mf[n]]]],trk_d0[trk_atIP[rcftr[r2mf[n]]]]);
		   			err_Z0->Fill(cov[trk_atIP[rcftr[r2mf[n]]]][9]);
		   			err_D0->Fill(cov[trk_atIP[rcftr[r2mf[n]]]][0]);
					pt->Fill(reco.Perp());
					radiusInnerHit->Fill(trrih[rcftr[r2mf[n]]]);
	   				}
	   			//se la particella proviene da quark/gluoni
				else if(abs(mcpdg[mcpa[r2mt[n]]])==5 || abs(mcpdg[mcpa[r2mt[n]]])==4 || abs(mcpdg[mcpa[r2mt[n]]])==21 || abs(mcpdg[mcpa[r2mt[n]]])==3 || abs(mcpdg[mcpa[r2mt[n]]])==2 || abs(mcpdg[mcpa[r2mt[n]]])==1){
				//else if(abs(mcpdg[mcpa[t2mt[n]]])==211 || abs(mcpdg[mcpa[t2mt[n]]])==310 || abs(mcpdg[mcpa[t2mt[n]]])==130 ||abs(mcpdg[mcpa[t2mt[n]]])==321){
					dist_o->Fill(distance);
					//fill z0
		   			reco_Z0_o->Fill(trk_z0[trk_atIP[rcftr[r2mf[n]]]]);
		   			//fill d0
		   			reco_D0_o->Fill(trk_d0[trk_atIP[rcftr[r2mf[n]]]]);
					err_Z0_o->Fill(cov[trk_atIP[rcftr[r2mf[n]]]][9]);
		   			err_D0_o->Fill(cov[trk_atIP[rcftr[r2mf[n]]]][0]);
					pt_o->Fill(reco.Perp());
					radiusInnerHit_o->Fill(trrih[rcftr[r2mf[n]]]);
   					}
   				}
   			
   			}
		}
	
	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
	legend2->AddEntry(reco_Z0, "from heavy hadron decay", "l");
	legend2->AddEntry(reco_Z0_o, "from hadronization" , "l");
	legend2->SetTextSize(0.04);
   	
   	TCanvas *c100= new TCanvas("c100","c100",900,600);
   	radiusInnerHit_o->Scale(1./radiusInnerHit_o->Integral());
	radiusInnerHit->Scale(1./radiusInnerHit->Integral());
   	radiusInnerHit->Draw("s hist");
	radiusInnerHit_o->Draw("sames hist");
	legend2->Draw();
   		
   	TCanvas *c1= new TCanvas("c1","c1",900,600);
   	reco_Z0->Scale(1./reco_Z0->Integral());
	reco_Z0_o->Scale(1./reco_Z0_o->Integral());
   	reco_Z0->Draw("s hist");
	reco_Z0_o->Draw("sames hist");
	legend2->Draw();
   	
   	TCanvas *c2= new TCanvas("c2","c2",900,600);
   	reco_D0->Scale(1./reco_D0->Integral());
   	reco_D0->Draw("s hist");
	reco_D0_o->Scale(1./reco_D0_o->Integral());
	reco_D0_o->Draw("sames hist");
	legend2->Draw();
   	
   	TCanvas *c3= new TCanvas("c3","c3",900,600);
   	//colz->Scale(1./colz->Integral());
   	colz->Draw("colz");
   	
   	TCanvas *c10= new TCanvas("c10","c10",900,600);
   	err_Z0->Scale(1./err_Z0->Integral());
	err_Z0_o->Scale(1./err_Z0_o->Integral());
   	err_Z0->Draw("s hist");
	err_Z0_o->Draw("sames hist");
	legend2->Draw();
   	
   	TCanvas *c20= new TCanvas("c20","c20",900,600);
   	err_D0->Scale(1./err_D0->Integral());
	err_D0_o->Scale(1./err_D0_o->Integral());
   	err_D0->Draw("s hist");
  	err_D0_o->Draw("sames hist");	
  	legend2->Draw();

	TCanvas *c30= new TCanvas("c30","c30",900,600);
   	pt->Scale(1./pt->Integral());
	pt_o->Scale(1./pt_o->Integral());
   	pt_o->Draw("s hist");
        pt->Draw("sames hist");
   	legend2->Draw();
   	
   	TCanvas *c200= new TCanvas("c200","c200",900,600);
   	dist->Scale(1./dist->Integral());
	dist_o->Scale(1./dist_o->Integral());
  	dist_o->Draw("s hist");	
  	dist->Draw("sames hist");
  	legend2->Draw();
        
	TCanvas *c500= new TCanvas("c500","c500",900,600);
	higgs_decay->Draw();
	
	free(mcvtx);free(mcvty);free(mcvtz);free(mcpdg);free(mcpa);free(r2mt);free(r2mf);free(rcftr);free(r2mw);
	free(rcmox);free(rcmoy);free(rcmoz);free(trk_z0);free(ene);free(trk_d0);free(trrih);free(trk_atIP);
	
	}

