#include "TLorentzVector.h"
void mass_reco(){
	gStyle->SetOptStat("nemrou");
	
	long int n,ientry,i;
	
	//relation
	int num2=1000000;
	Int_t *r2mt,*r2mf, nrel;
	Float_t *r2mw;
	r2mt = (int*) malloc(sizeof(int)*num2);
	r2mf = (int*) malloc(sizeof(int)*num2);
	r2mw = (float*) malloc(sizeof(float)*num2);
	
	//reco particle
	Float_t *rcmox,*rcmoy,*rcmoz,*ene,*rcmas;
	Int_t *type;
	rcmox = (float*) malloc(sizeof(float)*num2);
	rcmoy = (float*) malloc(sizeof(float)*num2);
	rcmoz = (float*) malloc(sizeof(float)*num2);
	rcmas = (float*) malloc(sizeof(float)*num2);
	ene = (float*) malloc(sizeof(float)*num2);
	type = (int*) malloc(sizeof(int)*num2);
	
	//mc particle
	Int_t *mcpdg,*mcpa0,*mcpa1,nmc;
	Float_t *mcmox,*mcmoy,*mcmoz,*mcene,*mcmas;
	int num=100000000;
	mcpdg = (int*) malloc(sizeof(int)*num);
	mcpa0 = (int*) malloc(sizeof(int)*num);
	mcpa1 = (int*) malloc(sizeof(int)*num);
	mcmox = (float*) malloc(sizeof(float)*num);
	mcmoy = (float*) malloc(sizeof(float)*num);
	mcmoz = (float*) malloc(sizeof(float)*num);
	mcmas = (float*) malloc(sizeof(float)*num);
	mcene = (float*) malloc(sizeof(float)*num);
	
	TChain* fChain = new TChain("fChain");
   	fChain->Add("/home/paola/Scrivania/xml_aggiornati/ntuple_bb_1000_newgeom.root/MyLCTuple");
   	
	  
	fChain->SetBranchAddress("r2mnrel", &nrel);
	fChain->SetBranchAddress("r2mt", r2mt);
	fChain->SetBranchAddress("r2mf", r2mf);	
        fChain->SetBranchAddress("r2mw", r2mw);
        
        fChain->SetBranchAddress("rcmox", rcmox);
  	fChain->SetBranchAddress("rcmoy", rcmoy);
  	fChain->SetBranchAddress("rcmoz", rcmoz);
	fChain->SetBranchAddress("rcene", ene);
	fChain->SetBranchAddress("rcmas", rcmas);
	fChain->SetBranchAddress("rctyp", type);
	
	fChain->SetBranchAddress("nmcp", &nmc);
	fChain->SetBranchAddress("mcpdg", mcpdg);
	fChain->SetBranchAddress("mcpa0", mcpa0);
	fChain->SetBranchAddress("mcpa1", mcpa1);
	
	fChain->SetBranchAddress("mcmox", mcmox);
	fChain->SetBranchAddress("mcmoy", mcmoy);
	fChain->SetBranchAddress("mcmoz", mcmoz);
	fChain->SetBranchAddress("mcene", mcene);
	fChain->SetBranchAddress("mcmas", mcmas);
	
	TH1F* energy=new TH1F("energy","energy",300,0,10);
	energy->SetLineWidth(1);
	energy->SetLineColor(kRed);
	energy->SetTitle("Energy");
    	energy->GetXaxis()->SetTitle("E [GeV]");
    	
	TH1F* pt=new TH1F("pt","pt",300,0,10);
	pt->SetLineWidth(1);
	pt->SetLineColor(kRed);
	pt->SetTitle("pT");
    	pt->GetXaxis()->SetTitle("pT [GeV]");
    	
	TH1F* eta=new TH1F("eta","eta",300,-2.66,2.66);
	eta->SetLineWidth(1);
	eta->SetLineColor(kRed);
	eta->SetTitle("#eta");
    	eta->GetXaxis()->SetTitle("#eta");
	
	TH1F* energy_all=new TH1F("energy_all","energy_all",300,0,10);
	energy_all->SetLineWidth(1);
	energy_all->SetLineColor(kBlue);
	energy_all->SetTitle("Energy");
    	energy_all->GetXaxis()->SetTitle("E [GeV]");
    	
	TH1F* pt_all=new TH1F("pt_all","pt_all",300,0,10);
	pt_all->SetLineWidth(1);
	pt_all->SetLineColor(kBlue);
	pt_all->SetTitle("pT");
    	pt_all->GetXaxis()->SetTitle("pT [GeV]");
    	
	TH1F* eta_all=new TH1F("eta_all","eta_all",300,-2.66,2.66);
	eta_all->SetLineWidth(1);
	eta_all->SetLineColor(kBlue);
	eta_all->SetTitle("#eta");
    	eta_all->GetXaxis()->SetTitle("#eta");
    	
    	TH1F* eta_region=new TH1F("eta_region","eta_region",4,-0.5,3.5);
	eta_region->SetLineWidth(1);
	eta_region->SetLineColor(kBlue);
	eta_region->SetTitle("#eta");
    	eta_region->GetXaxis()->SetTitle("#eta region");
    	
	TH1F* deltaM=new TH1F("deltaM","deltaM",100,-0.5,1);
	deltaM->GetXaxis()->SetTitle("#DeltaM[GeV]");
	TH1F* deltaM_zoom=new TH1F("deltaM_zoom","deltaM_zoom",50,0.9,1);
	deltaM_zoom->GetXaxis()->SetTitle("#DeltaM[GeV]");

	//TH1F* pdg=new TH1F("pdg","pdg",2250,-0.5,2249.5);
	//TH1F* mass=new TH1F("mass","mass",100,-0.08,1);
	
	TLorentzVector reco,mc;
	
	//int id[20][2]={0,0},counter=0,isequal=0;
	
	float reco_mom,mc_mom;

	//int count=0;
	
	//ciclo su entries
        for(ientry=0; ientry<fChain->GetEntries() ;++ientry){
  		fChain->GetEntry(ientry);
		//ciclo su relazioni mc-reco
  		for(n=0;n<nrel;n++){
			//richiedo un buon link mc-reco
			if(r2mw[n]>0.9){
			
  				reco.SetPxPyPzE(rcmox[r2mf[n]],rcmoy[r2mf[n]],rcmoz[r2mf[n]],ene[r2mf[n]]);
				mc.SetPxPyPzE(mcmox[r2mt[n]],mcmoy[r2mt[n]],mcmoz[r2mt[n]],mcene[r2mt[n]]);
			
				eta_all->Fill(mc.Eta());
				pt_all->Fill(mc.Pt());
				energy_all->Fill(mc.E());
				
			   //differenza tra massa calcolata da energia e momento misurati e massa da pdg associata dopo id da Particle Flow
			   //es. se la particella è id come neutrone rcmas sarà 0.939565
				deltaM->Fill(rcmas[r2mf[n]]-reco.M());
				
				//prendo solo particelle problematiche
				if(abs(reco.M()-rcmas[r2mf[n]])>0.5){
				
					if(mc.Eta()>-0.5&&mc.Eta()<0.5){eta_region->Fill(0);}
					if((mc.Eta()>-1.5 &&mc.Eta()<-0.5)||(mc.Eta()>0.5 &&mc.Eta()<1.5)){eta_region->Fill(1);}
					if((mc.Eta()>-2 &&mc.Eta()<-1.5)||(mc.Eta()>1.5 &&mc.Eta()<2.)){eta_region->Fill(2);}
					if(mc.Eta()<-2||mc.Eta()>2){eta_region->Fill(3);}
					//if(mc.Pt()<1){
					
	  				eta->Fill(mc.Eta());
	  				pt->Fill(mc.Pt());
					energy->Fill(mc.E());
					deltaM_zoom->Fill(rcmas[r2mf[n]]-reco.M());
	  				//pdg->Fill(abs(mcpdg[r2mt[n]]));
	  				//if(mcpdg[r2mt[n]]>10000){cout<<mcpdg[r2mt[n]]<<endl;}
	  				//pdg->Fill(abs(type[r2mf[n]]));
	  				/*for(i=0;i<counter;i++){
	  					if(abs(id[i][0])==abs(mcpdg[r2mt[n]])){isequal=1;id[i][1]++;}
	  					}
	  				if(isequal==0){id[counter][0]=abs(mcpdg[r2mt[n]]);id[counter][1]=1;counter++;}
	  				isequal=0;*/
					reco_mom=sqrt(reco.Px()*reco.Px()+reco.Py()*reco.Py()+reco.Pz()*reco.Pz());
					mc_mom=sqrt(mc.Px()*mc.Px()+mc.Py()*mc.Py()+mc.Pz()*mc.Pz());
					
					cout<<"id gen: "<<mcpdg[r2mt[n]]<<" id reco: "<<type[r2mf[n]]<<endl;
					cout<<"Gen eta: "<<mc.Eta()<<"Gen Pt: "<<mc.Pt()<<endl;
					cout<<"RECO Px: "<<reco.Px()<<" Py: "<<reco.Py()<<" Pz: "<<reco.Pz()<<" Energia: "<<reco.E()<<" Mass: "<<reco.M()<<endl;
					cout<<"GEN Px: "<<mc.Px()<<" Py: "<<mc.Py()<<" Pz: "<<mc.Pz()<<" Energia: "<<mc.E()<<" Mass: "<<mc.M()<<endl;
					cout<< " "<<endl;
  					}
				}
  			}
		cout<<"FINE EVENTO"<<endl;
  		}
	

	//for(i=0;i<counter;i++){cout<<"Particella: "<<id[i][0]<<" : "<<id[i][1]<<"volte" <<endl;}
	TCanvas *c0=new TCanvas();
	eta_region->Scale(1./eta_region->Integral());
	eta_region->Draw("s hist");
	
	TCanvas *c00=new TCanvas();
	deltaM->Scale(1./deltaM->Integral());
	deltaM->Draw("s hist");
	
	TCanvas *c000=new TCanvas();
	deltaM_zoom->Draw();
	//pdg->Draw();
	
	TCanvas *c1=new TCanvas();
	energy->Scale(1./energy->Integral());
	energy->Draw("s hist");
	energy_all->Scale(1./energy_all->Integral());
	energy_all->Draw("sames hist");
	
	TCanvas *c2=new TCanvas();
	pt->Scale(1./pt->Integral());
	pt->Draw("s hist");
	pt_all->Scale(1./pt_all->Integral());
	pt_all->Draw("sames hist");
	
	TCanvas *c3=new TCanvas();
	eta->Scale(1./eta->Integral());
	eta->Draw("s hist");
	eta_all->Scale(1./eta_all->Integral());
	eta_all->Draw("sames hist");
	
	free(mcpdg);free(mcpa0);free(r2mt);free(r2mf);free(r2mw);free(mcpa1);
	free(rcmox);free(rcmoy);free(rcmoz);free(ene);free(mcmox);free(mcmoy);free(mcmoz);free(mcene);free(mcmas);free(rcmas);
	}
	
