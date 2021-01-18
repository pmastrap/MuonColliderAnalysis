////////DESCRIPTION///////////////////////////////////
//						    //
// evaluate resolution of reco-p variables	    //
// plot of energy, eta, px, py, pz, d0, z0	    //
//////////////////////////////////////////////////////

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

void resolution(){

	gStyle->SetOptStat("nemrou");
	
	long int ientry,n,l;
	Int_t nmc,nrel;
	
	int num=1000000;
	Int_t *mcpdg,*mcpa,*mcda0;
	Float_t *mcmox,*mcmoy,*mcmoz,*mcene;
	Float_t *mcvtx,*mcvty,*mcvtz;
	Float_t *mcepx,*mcepy,*mcepz;
	
	//mc particles
	mcpdg = (int*) malloc(sizeof(int)*num);
	mcpa = (int*) malloc(sizeof(int)*num);
	mcmox = (float*) malloc(sizeof(float)*num);
	mcmoy = (float*) malloc(sizeof(float)*num);
	mcmoz = (float*) malloc(sizeof(float)*num);
	mcene = (float*) malloc(sizeof(float)*num);

	mcvtx = (float*) malloc(sizeof(float)*num);
	mcvty = (float*) malloc(sizeof(float)*num);
	mcvtz = (float*) malloc(sizeof(float)*num);

	mcepx = (float*) malloc(sizeof(float)*num);
	mcepy = (float*) malloc(sizeof(float)*num);
	mcepz = (float*) malloc(sizeof(float)*num);
	
	int num2=1000000;
	Int_t *trk_atIP,*r2mf,*r2mt,*rcftr;
	Float_t *r2mw,*trk_d0,*trk_z0,*trk_omega;
	
	//relations
	r2mf = (int*) malloc(sizeof(int)*num2);
	r2mt = (int*) malloc(sizeof(int)*num2);
	rcftr = (int*) malloc(sizeof(int)*num2);
	r2mw = (float*) malloc(sizeof(float)*num2);
	
	//tracks
	trk_atIP = (int*) malloc(sizeof(int)*num2);
	trk_d0 = (float*) malloc(sizeof(float)*num2);
	trk_z0 = (float*) malloc(sizeof(float)*num2);
	trk_omega = (float*) malloc(sizeof(float)*num2);
	
	//reco particles
	Float_t *rcmox,*rcmoy,*rcmoz,*ene;
	rcmox = (float*) malloc(sizeof(float)*num2);
	rcmoy = (float*) malloc(sizeof(float)*num2);
	rcmoz = (float*) malloc(sizeof(float)*num2);
	ene = (float*) malloc(sizeof(float)*num2);

	TChain* fChain = new TChain("fChain");
   	fChain->Add("/home/paola/Scrivania/xml_aggiornati/ntuple_bb_1000_newgeom.root/MyLCTuple");
   	
   	fChain->SetBranchAddress("nmcp", &nmc);
	fChain->SetBranchAddress("mcpdg", mcpdg);
	fChain->SetBranchAddress("mcpa0", mcpa);
	fChain->SetBranchAddress("mcmox", mcmox);
	fChain->SetBranchAddress("mcmoy", mcmoy);
	fChain->SetBranchAddress("mcmoz", mcmoz);
	fChain->SetBranchAddress("mcene", mcene);
	
	//end point
	fChain->SetBranchAddress("mcepx", mcepx);
	fChain->SetBranchAddress("mcepy", mcepy);
	fChain->SetBranchAddress("mcepz", mcepz);
	
	//production point
	fChain->SetBranchAddress("mcvtx", mcvtx);
	fChain->SetBranchAddress("mcvty", mcvty);
	fChain->SetBranchAddress("mcvtz", mcvtz);
	
	fChain->SetBranchAddress("r2mnrel", &nrel);
	fChain->SetBranchAddress("r2mt", r2mt);
	fChain->SetBranchAddress("r2mf", r2mf);	
        fChain->SetBranchAddress("r2mw", r2mw);
	fChain->SetBranchAddress("rcftr", rcftr);
	
        fChain->SetBranchAddress("tszze",trk_z0);
	fChain->SetBranchAddress("trsip", trk_atIP);
        fChain->SetBranchAddress("tsdze", trk_d0);
        fChain->SetBranchAddress("tsome", trk_omega);
        
        fChain->SetBranchAddress("rcmox", rcmox);
  	fChain->SetBranchAddress("rcmoy", rcmoy);
  	fChain->SetBranchAddress("rcmoz", rcmoz);
	fChain->SetBranchAddress("rcene", ene);
	
	/*TH1F* pt_gen=new TH1F("pt_gen","pt_gen",10,0.,50.);
	pt_gen->SetLineWidth(2);
        pt_gen->SetTitle("pT");
    	pt_gen->GetXaxis()->SetTitle("pt [GeV]");
	
	TH1F* pt_gen_matched=new TH1F("pt_gen_matched","pt_gen_matched",10,0.,50.);
	pt_gen_matched->SetLineWidth(2);
	pt_gen_matched->SetLineColor(kRed);
	pt_gen_matched->SetTitle("pT");
    	pt_gen_matched->GetXaxis()->SetTitle("pt [GeV]");
    	
    	TH1F* dist=new TH1F("dist","dist",100,0.,0.5);*/
	
    	TProfile* resol_px=new TProfile("resol_px","resol_px",20,-2.66,2.66,-200,200.);
    	TProfile* resol_py=new TProfile("resol_py","resol_py",20,-2.66,2.66,-200,200.);
    	TProfile* resol_pz=new TProfile("resol_pz","resol_pz",20,-2.66,2.66,-200,200.);
    	TProfile* resol_ene=new TProfile("resol_ene","resol_ene",20,-2.66,2.66,-200,200.);
    	TProfile* resol_p=new TProfile("resol_p","resol_p",20,-2.66,2.66,-200,200.);
    	
    	TProfile* resol_eta=new TProfile("resol_eta","resol_eta",20,-2.66,2.66,-200,200.);
    	TProfile* resol_pt=new TProfile("resol_pt","resol_pt",20,-2.66,2.66,-200,200.);
    	
    	TProfile* resol_d0=new TProfile("resol_d0","resol_d0",20,-2.66,2.66,-200,1000.);
    	TProfile* resol_z0=new TProfile("resol_z0","resol_z0",20,-2.66,2.66,-200,1000.);
    	
    	TH1F* boost=new TH1F("boost","boost",100,0.,20);
    	
    	TLorentzVector mc,reco,mc_pa;
    	TVector3 start_point,end_point;
    	float pt_track,reco_mom,mc_mom,mc_pa_mom,d0,z0,Lxy;
	int count=0,pdg[20]={-1},length=0;
	
	for(ientry=0; ientry<fChain->GetEntries() ;++ientry){
  		fChain->GetEntry(ientry);
  		cout<<"nrel: "<<nrel<<endl;	
  		for(n=0;n<nrel;n++){
			if(r2mw[n]>0.9){
			
  				mc.SetPxPyPzE(mcmox[r2mt[n]],mcmoy[r2mt[n]],mcmoz[r2mt[n]],mcene[r2mt[n]]);
  				reco.SetPxPyPzE(rcmox[r2mf[n]],rcmoy[r2mf[n]],rcmoz[r2mf[n]],ene[r2mf[n]]);
  				
  				//pt_track=3*10e-5*abs(3.57/trk_omega[trk_atIP[r2mf[n]]]);
  				reco_mom=sqrt(reco.Px()*reco.Px()+reco.Py()*reco.Py()+reco.Pz()*reco.Pz());
				mc_mom=sqrt(mc.Px()*mc.Px()+mc.Py()*mc.Py()+mc.Pz()*mc.Pz());
				
				resol_eta->Fill(mc.Eta(),abs((mc.Eta()-reco.Eta())/mc.Eta()));
  				resol_pt->Fill(mc.Eta(),abs((mc.Pt()-reco.Pt())/mc.Pt()));
  				resol_ene->Fill(mc.Eta(),abs((mc.E()-reco.E())/mc.E()));
  				resol_px->Fill(mc.Eta(),abs((mc.Px()-reco.Px())/mc.Px()));
  				resol_py->Fill(mc.Eta(),abs((mc.Py()-reco.Py())/mc.Py()));
  				resol_pz->Fill(mc.Eta(),abs((mc.Pz()-reco.Pz())/mc.Pz()));
  				resol_p->Fill(mc.Eta(),abs((mc_mom-reco_mom)/mc_mom));
  				
  				//se alla reco-particle è associata una track
  				if(rcftr[r2mf[n]]!=-1){
  					//se il parent della particella al gen è un adrone pesante
  					if(selection_cc(mcpdg[mcpa[r2mt[n]]])==0 || selection_bb(mcpdg[mcpa[r2mt[n]]])==0){
  						
  						//confronto d0 misurato con quello atteso
  						//distanza tra B/C production point e decay point nel piano xy
  						Lxy=sqrt(pow(mcvtx[mcpa[r2mt[n]]]-mcepx[mcpa[r2mt[n]]],2)+pow(mcvty[mcpa[r2mt[n]]]-mcepy[mcpa[r2mt[n]]],2));
  						d0=Lxy*sin(mc.Pt()/mc_mom);
  						z0=mcvtz[mcpa[r2mt[n]]]-mcepz[mcpa[r2mt[n]]];
  						resol_d0->Fill(mc.Eta(),abs((d0-trk_d0[trk_atIP[rcftr[r2mf[n]]]])/d0));
  						resol_z0->Fill(mc.Eta(),abs((z0-trk_z0[trk_atIP[rcftr[r2mf[n]]]])/z0));
  						}
  					}
  				}
  			}
  		length=0;
  		count=0;
  		}
  	/*for(ientry=0; ientry<fChain->GetEntries() ;++ientry){
  		fChain->GetEntry(ientry);
	  	for(l=0;l<nmc;l++){
	  		mc.SetPxPyPzE(mcmox[l],mcmoy[l],mcmoz[l],mcene[l]);
	  		mc_mom=sqrt(mc.Px()*mc.Px()+mc.Py()*mc.Py()+mc.Pz()*mc.Pz());
	  		//boost plot degli adroni pesanti
	  		if(selection_cc(mcpdg[n])==0 || selection_bb(mcpdg[n])==0){
	  			if(abs(mcpdg[mcpa[n]])==5 || abs(mcpdg[mcpa[n]])==4 || abs(mcpdg[mcpa[n]])==21 || abs(mcpdg[mcpa[n]])==3 || abs(mcpdg[mcpa[n]])==2 || abs(mcpdg[mcpa[n]])==1){
		  			cout<<mcpdg[n]<<endl;
		  			boost->Fill(mc_mom/mc.Pt());
		  			}
	  			}
	  		}
	  	}*/
  	cout<<count<<endl;
	TCanvas* c1=new TCanvas();
	resol_eta->Draw();
	
	TCanvas* c2=new TCanvas();
	resol_pt->Draw();
	
	TCanvas* c3=new TCanvas();
	resol_ene->Draw();
	
	TCanvas* c4=new TCanvas();
	resol_px->Draw();
	
	TCanvas* c5=new TCanvas();
	resol_py->Draw();
	
	TCanvas* c6=new TCanvas();
	resol_pz->Draw();
	
	TCanvas* c7=new TCanvas();
	resol_p->Draw();
	
	TCanvas* c8=new TCanvas();
	resol_d0->Draw();
	
	TCanvas* c9=new TCanvas();
	resol_z0->Draw();
	
	TCanvas* c10=new TCanvas();
	boost->Scale(1./boost->Integral());
	boost->Draw("s hist");
	
	
  	/*TH1F *h_ratio = (TH1F*)pt_gen_matched->Clone("h_ratio");
  	h_ratio->Sumw2();
  	h_ratio->SetStats(0);
  	h_ratio->Divide(pt_gen);
  	h_ratio->Draw("ep");*/
	free(mcpdg);free(mcpa);free(r2mt);free(r2mf);free(r2mw);free(trk_d0);free(trk_omega);free(rcftr);
	free(rcmox);free(rcmoy);free(rcmoz);free(ene);free(mcmox);free(mcmoy);free(mcmoz);free(mcene);free(trk_z0);free(trk_atIP);
	free(mcepx);free(mcepy);free(mcepz);free(mcvtx);free(mcvty);free(mcvtz);
	}
