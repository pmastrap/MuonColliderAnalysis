////////DESCRIPTION///////////////////////////////////////////////////////////////////
//						      				    //
// signed impact parameter of all tracks inside a jet 				    //
// this parameter is not avaible in the collection   				    //
// we need to build it up:   			     				    //
// 1st: find PCA from d0,z0 and phi		      				    //
// 2nd: evaluate the sign taking the scalar product between jet axis and PCA vector //
//										    //
//////////////////////////////////////////////////////////////////////////////////////


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


void SIP(){

	gStyle->SetOptStat("nemrou");
	
	Int_t ntrk,nj,index,t2mnrel;
	
	long int ientry,n,k,l;
	
	int num2=100000,num=2000000;
	Float_t *trk_z0,*trk_d0,*tsphi,*tstnl,*t2mw;
	Int_t *trk_atIP,*t2mf,*t2mt;
	Int_t *mcpdg,*mcpa;
	
	mcpdg = (int*) malloc(sizeof(int)*num);
	mcpa = (int*) malloc(sizeof(int)*num);
	trk_z0 = (float*) malloc(sizeof(float)*num2);
	trk_d0 = (float*) malloc(sizeof(float)*num2);
	tsphi = (float*) malloc(sizeof(float)*num2);
	tstnl = (float*) malloc(sizeof(float)*num2);
	trk_atIP = (int*) malloc(sizeof(int)*num2);
	t2mf = (int*) malloc(sizeof(int)*num2);
	t2mt = (int*) malloc(sizeof(int)*num2);
	t2mw = (float*) malloc(sizeof(float)*num2);
	
	Float_t cov[100000][15];
	
	Float_t jmox[1000], jmoy[1000], jmoz[1000];
	
	Float_t trk_theta,trk_eta,delta_R,dist,min;
	
	TVector3 PCA,jet_mo;
	
	TChain* fChain = new TChain("fChain");
   	fChain->Add("/home/paola/Scrivania/xml_aggiornati/ntuple_cc_1000_newgeom.root/MyLCTuple");
   	fChain->SetBranchAddress("ntrk",&ntrk);
	fChain->SetBranchAddress("tszze",trk_z0);
	fChain->SetBranchAddress("trsip", trk_atIP);
        fChain->SetBranchAddress("tsdze", trk_d0);
        fChain->SetBranchAddress("tsphi", tsphi);
        fChain->SetBranchAddress("tscov", cov);
        fChain->SetBranchAddress("tstnl",tstnl);

	fChain->SetBranchAddress("t2mnrel",&t2mnrel);
	fChain->SetBranchAddress("t2mf",t2mf);
	fChain->SetBranchAddress("t2mt",t2mt);
	fChain->SetBranchAddress("t2mw",t2mw);
        
        fChain->SetBranchAddress("nj", &nj);
   	fChain->SetBranchAddress("jmox", jmox);
	fChain->SetBranchAddress("jmoy", jmoy);
	fChain->SetBranchAddress("jmoz", jmoz);
	
	fChain->SetBranchAddress("mcpdg", mcpdg);
	fChain->SetBranchAddress("mcpa0", mcpa);
	
	TH1F* SIP_plot=new TH1F("SIP_plot_c","SIP_plot",50,-0.5,0.5);
	SIP_plot->SetTitle("Signed 3D impact parameter");
	SIP_plot->SetLineWidth(2);
    	SIP_plot->GetXaxis()->SetTitle("SIP[mm]");
    	SIP_plot->SetLineColor(kRed);
    	
    	TH1F* R_plot=new TH1F("R_plot","R_plot",50,0,10.);
	
	//ciclo su entries
        for(ientry=0; ientry<fChain->GetEntries() ;++ientry){
  		fChain->GetEntry(ientry);
  		//ciclo su tracce
  		for(n=0;n<ntrk;n++){
  		
  			//cerco relazione corrispondente
  			l=0;
			while(t2mf[l]!=n && l<t2mnrel ){l++;}
			//seleziono solo tracce che provengono da adroni pesanti (al gen)
  			//if(t2mw[l]>0.9 && (selection_cc(mcpdg[mcpa[t2mt[l]]])==0 || selection_bb(mcpdg[mcpa[t2mt[l]]])==0)){
  			cout<<"padre: "<<mcpdg[mcpa[t2mt[l]]]<<endl;
  			trk_theta = (TMath::PiOver2()-TMath::ATan(tstnl[trk_atIP[n]]));
			trk_eta=-TMath::Log(TMath::Tan(trk_theta/2));
			min=10000;
			//trovo il jet con l'asse più vicino alla traccia
			for(k=0;k<nj;k++){
				jet_mo.SetX(jmox[k]);
				jet_mo.SetY(jmoy[k]);
				jet_mo.SetZ(jmoz[k]);
  				delta_R=sqrt(pow(jet_mo.Eta()-trk_eta,2)+pow(jet_mo.Phi()-tsphi[trk_atIP[n]],2));
  				if(delta_R<min){
  					min=delta_R;
  					index=k;
  					}
  				}
  			jet_mo.SetX(jmox[index]);
			jet_mo.SetY(jmoy[index]);
			jet_mo.SetZ(jmoz[index]);
			delta_R=sqrt(pow(jet_mo.Eta()-trk_eta,2)+pow(jet_mo.Phi()-tsphi[trk_atIP[n]],2));
			R_plot->Fill(delta_R);
			if(delta_R<1){
			  	PCA.SetX(-(trk_d0[trk_atIP[n]])*TMath::Sin(tsphi[trk_atIP[n]]));
			  	PCA.SetY((trk_d0[trk_atIP[n]])*TMath::Cos(tsphi[trk_atIP[n]]));
				PCA.SetZ(trk_z0[trk_atIP[n]]);
					
				dist=sqrt(PCA.X()*PCA.X()+PCA.Y()*PCA.Y()+PCA.Z()*PCA.Z());
				if(PCA*jet_mo>0){
					SIP_plot->Fill(dist);
					}
				if(PCA*jet_mo<0){
					SIP_plot->Fill(-dist);
					}
				}	
			//}
			}
			
		}
	//bb sample
	TChain* fChainb = new TChain("fChainb");
   	fChainb->Add("/home/paola/Scrivania/xml_aggiornati/ntuple_bb_1000_newgeom.root/MyLCTuple");
   	fChainb->SetBranchAddress("ntrk",&ntrk);
	fChainb->SetBranchAddress("tszze",trk_z0);
	fChainb->SetBranchAddress("trsip", trk_atIP);
        fChainb->SetBranchAddress("tsdze", trk_d0);
        fChainb->SetBranchAddress("tsphi", tsphi);
        fChainb->SetBranchAddress("tscov", cov);
        fChainb->SetBranchAddress("tstnl",tstnl);
        
        fChainb->SetBranchAddress("nj", &nj);
   	fChainb->SetBranchAddress("jmox", jmox);
	fChainb->SetBranchAddress("jmoy", jmoy);
	fChainb->SetBranchAddress("jmoz", jmoz);
	
	TH1F* SIP_plot_b=new TH1F("SIP_plot_b","SIP_plot_b",50,-0.5,0.5);
	SIP_plot_b->SetTitle("Signed 3D impact parameter");
	SIP_plot_b->SetLineWidth(2);
    	SIP_plot_b->GetXaxis()->SetTitle("SIP[mm]");
    	
	
	//ciclo su entries
        for(ientry=0; ientry<fChainb->GetEntries() ;++ientry){
  		fChainb->GetEntry(ientry);
  		//ciclo su tracce
  		for(n=0;n<ntrk;n++){
  			trk_theta = (TMath::PiOver2()-TMath::ATan(tstnl[trk_atIP[n]]));
			trk_eta=-TMath::Log(TMath::Tan(trk_theta/2));
			min=10000;
			for(k=0;k<nj;k++){
				jet_mo.SetX(jmox[k]);
				jet_mo.SetY(jmoy[k]);
				jet_mo.SetZ(jmoz[k]);
  				delta_R=sqrt(pow(jet_mo.Eta()-trk_eta,2)+pow(jet_mo.Phi()-tsphi[trk_atIP[n]],2));
  				if(delta_R<min){
  					min=delta_R;
  					index=k;
  					}
  				}
  			jet_mo.SetX(jmox[index]);
			jet_mo.SetY(jmoy[index]);
			jet_mo.SetZ(jmoz[index]);
			delta_R=sqrt(pow(jet_mo.Eta()-trk_eta,2)+pow(jet_mo.Phi()-tsphi[trk_atIP[n]],2));
			if(delta_R<1){
			  	PCA.SetX(-(trk_d0[trk_atIP[n]])*TMath::Sin(tsphi[trk_atIP[n]]));
			  	PCA.SetY((trk_d0[trk_atIP[n]])*TMath::Cos(tsphi[trk_atIP[n]]));
				PCA.SetZ(trk_z0[trk_atIP[n]]);
					
				dist=sqrt(PCA.X()*PCA.X()+PCA.Y()*PCA.Y()+PCA.Z()*PCA.Z());
				if(PCA*jet_mo>0){
					SIP_plot_b->Fill(dist);
					}
				if(PCA*jet_mo<0){
					SIP_plot_b->Fill(-dist);
					}
				}	
				
			}
			
		}
		
	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
	legend2->AddEntry(SIP_plot, "H #rightarrow c#bar{c}", "l");
	legend2->AddEntry(SIP_plot_b, "H #rightarrow b#bar{b}" , "l");
	legend2->SetTextSize(0.04);
	
	SIP_plot->Scale(1./SIP_plot->Integral());	
	SIP_plot->Draw("s hist");
	
	SIP_plot_b->Scale(1./SIP_plot_b->Integral());	
	SIP_plot_b->Draw("sames hist");
	legend2->Draw();
	
	TCanvas *c1=new TCanvas();
	R_plot->Draw();

	free(trk_z0);free(trk_d0);free(tsphi);free(tstnl);free(trk_atIP);free(t2mf);free(t2mt);free(t2mw);
			
	}
