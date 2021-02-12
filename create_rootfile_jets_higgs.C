#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"



void create_rootfile_jets_higgs(Int_t plot_ID=0){

gStyle->SetOptStat("nemrou");

cout << "Plot IDs are the following" << endl;
cout << "0: number of reconstructed jets" << endl;
cout << "1: pseudorapidity of jets in the event" << endl;
cout << "2: transverse momentum of jets in the event" << endl;
cout << "3: phi of jets in the event" << endl;
cout << "4: energy of jets in the event" << endl;
cout << "5: invariant mass of jet pair associated to Higgs" << endl;

TH1F* n_jets_ev=new TH1F("n_jets_ev","n_jets_ev",10,-0.5,9.5);
	n_jets_ev->SetTitle("Number of jets in the event");
    	n_jets_ev->GetXaxis()->SetTitle("N_{jets}");
TH1F* m_inv=new TH1F("m_inv","m_inv",1000,0,2500);
	m_inv->SetTitle("Invariant mass of jet pair");
   	m_inv->GetXaxis()->SetTitle("M_{H} [GeV]");
TH1F* eta_histo=new TH1F("eta_histo","eta_histo",30,-5,5);
	eta_histo->SetTitle("Jets pseudorapidity");
    	eta_histo->GetXaxis()->SetTitle("#eta");
TH1F* E_histo=new TH1F("E_histo","E_histo",1000,0,2000);
	E_histo->SetTitle("Jets energy");
    	E_histo->GetXaxis()->SetTitle("Energy [GeV]");
TH1F* pt_histo=new TH1F("pt_histo","pt_histo",1000,0,2000);
	pt_histo->SetTitle("Jets transverse momentum");
    	pt_histo->GetXaxis()->SetTitle("p_{T} [GeV]");
TH1F* theta_histo=new TH1F("theta_histo","theta_histo",50,-0.5,4);
	theta_histo->SetTitle("Jets polar angle");
	theta_histo->GetXaxis()->SetTitle("#theta [rad]");
TH1F* phi_histo=new TH1F("phi_histo","phi_histo",100,-4,4);
	phi_histo->SetTitle("Jets azimuthal angle");
    	phi_histo->GetXaxis()->SetTitle("#phi [rad]");
//separazione tra jet
TH1F* Dtheta=new TH1F("Dtheta","Dtheta",50,-0.5,3);
	Dtheta->SetTitle("Jets #Delta#theta");
	Dtheta->GetXaxis()->SetTitle("#Delta#theta [rad]");
TH1F* Dr=new TH1F("Dr","Dr",100,0,10);
	Dr->SetTitle("Jets #DeltaR");
	Dr->GetXaxis()->SetTitle("#DeltaR [rad]");
TH1F* DPhi=new TH1F("DPhi","DPhi",50,-0.5,4);
	DPhi->SetTitle("Jets #Delta#phi");
	DPhi->GetXaxis()->SetTitle("#Delta#phi [rad]");
//Higgs Plots
TH1F* eta_higgs=new TH1F("eta_higgs","eta_higgs",30,-5,5);
	eta_higgs->SetTitle("Higgs pseudorapidity");
    	eta_higgs->GetXaxis()->SetTitle("#eta_{H}");
TH1F* E_higgs=new TH1F("E_higgs","E_higgs",1000,0,2000);
	E_higgs->SetTitle("Higgs energy");
    	E_higgs->GetXaxis()->SetTitle("Energy_{H} [GeV]");
TH1F* pt_higgs=new TH1F("pt_higgs","pt_higgs",1000,0,2000);
	pt_higgs->SetTitle("Higgs transverse momentum");
    	pt_higgs->GetXaxis()->SetTitle("p_{T} [GeV]");
TH1F* theta_higgs=new TH1F("theta_higgs","theta_higgs",50,-0.5,4);
	theta_higgs->SetTitle("Higgs polar angle");
	theta_higgs->GetXaxis()->SetTitle("#theta_{H} [rad]");
TH1F* phi_higgs=new TH1F("phi_higgs","phi_higgs",100,-4,4);
	phi_higgs->SetTitle("Higgs azimuthal angle");
    	phi_higgs->GetXaxis()->SetTitle("#phi_{H} [rad]");
    	
TH1F* h_invmass=new TH1F("h_invmass","h_invmass",1000,0,2500);
TH1F* pt_j1=new TH1F("pt_jet1","pt_jet1",200,0,1500);
pt_j1->SetLineColor(kBlue);
pt_j1->SetLineWidth(2);
TH1F* pt_j2=new TH1F("pt_jet2","pt_jet2",200,0,1500);
pt_j2->SetLineColor(kRed);
pt_j2->SetLineWidth(2);

   Int_t  nj;
   Float_t         jmox[100];   //[njet]
   Float_t         jmoy[100];   //[njet]
   Float_t         jmoz[100];   //[njet]
   Float_t         jene[100];   //[njet]

   TLorentzVector j[100],sum,h,j_selected[100],j_selected_sorted[100]; 

   TChain* fChain = new TChain("fChain");
   //JetHistogramGenJetTuple
   fChain->Add("/home/paola/Scrivania/xml_aggiornati/Irriducible_bkg/ntuple_cc_nohiggs_1000.root/MyLCTuple"); 
   //fChain->Add("/home/paola/Scrivania/xml_aggiornati/3/ntuple_cc_1000_mod_3.root/MyLCTuple");

   fChain->SetBranchAddress("nj", &nj);
   fChain->SetBranchAddress("jmox", jmox);
   fChain->SetBranchAddress("jmoy", jmoy);
   fChain->SetBranchAddress("jmoz", jmoz);
   fChain->SetBranchAddress("jene", jene);

   double inv_mass;

   double delta_min=10000000;
   double delta_m;
   int jet1_num=0;
   int jet2_num=0;
   int n_jet_selected=0;
   int n_jet_tot=0;
   int n_jet_etacut_tot=0,n_jet_etacut=0;
   int n_jet_ptcut_tot=0;
   int n_events_etacut=0;
   int n_events_ptcut=0;
   int n_events_minvcut=0;
   int n_events_tot=fChain->GetEntries();
   
   double pT[100];
   int index[100];

   double delta_R=0,delta_theta=0,delta_phi=0;
   
   double mh=125.;
   unsigned int ientry=0,k=0;

	for(ientry=0; ientry< fChain->GetEntries(); ++ientry){
   		fChain->GetEntry(ientry);
   		
   		n_jet_tot=n_jet_tot+nj;
   		
   		//selezione jet
   		// eta <2.5
   		// pt >10 GeV
   		n_jet_selected=0;
   		n_jet_etacut=0;
   		for (k=0;k<nj; k++){ 
   			j[k].SetPxPyPzE(jmox[k],jmoy[k],jmoz[k],jene[k]);
   			if(abs(j[k].Eta())<2.5){
   				n_jet_etacut_tot++;
   				n_jet_etacut++;
   				if(j[k].Pt()>10.){
   					n_jet_ptcut_tot++;
	   				j_selected[n_jet_selected]=j[k];
	   				n_jet_selected++;
	   				}
   				}
   			}
   		
   		if(n_jet_etacut>1){n_events_etacut++;}
   
  		
  		//sort jets by pT
  		
  		for(k=0;k<n_jet_selected;k++){
  			pT[k]=j_selected[k].Pt();
  			}
  		TMath::Sort(n_jet_selected,pT,index);
  		
  		for(k=0;k<n_jet_selected;k++){
  			j_selected_sorted[k]=j_selected[index[k]];
  			}
  		
  		for(k=0;k<n_jet_selected;k++){	
	  		cout<<pT[k]<<endl;
	  		cout<<index[k]<<endl;
	  		cout<<j_selected_sorted[k].Pt()<<endl;
	  		cout<<"---------------------"<<endl;
	  		}
  		cout<<"FINE EVENTO"<<endl;
  		
  		//selezione eventi
  		if(n_jet_selected>1){ //selection of events with more than one jet selected
  			n_events_ptcut++;
  			h=j_selected_sorted[0]+j_selected_sorted[1];
  			
  			pt_j1->Fill(j_selected_sorted[0].Pt());
  			pt_j2->Fill(j_selected_sorted[1].Pt());
  			
  			if(h.M()<250.){//selection of events with dijet invariant mass less than 250 GeV
	  			n_events_minvcut++;
	  			m_inv->Fill(h.M());
				eta_higgs->Fill(h.PseudoRapidity());
				pt_higgs->Fill(h.Perp());
				theta_higgs->Fill(h.Theta());
				phi_higgs->Fill(h.Phi());
				E_higgs->Fill(h.E());
				
				//histograms of the selected jets of selected events
				n_jets_ev->Fill(n_jet_selected);
		   		eta_histo->Fill(j_selected_sorted[0].PseudoRapidity());
				pt_histo->Fill(j_selected_sorted[0].Perp());
				theta_histo->Fill(j_selected_sorted[0].Theta());
				phi_histo->Fill(j_selected_sorted[0].Phi());
				E_histo->Fill(j_selected_sorted[0].E());
				
				eta_histo->Fill(j_selected_sorted[1].PseudoRapidity());
				pt_histo->Fill(j_selected_sorted[1].Perp());
				theta_histo->Fill(j_selected_sorted[1].Theta());
				phi_histo->Fill(j_selected_sorted[1].Phi());
				E_histo->Fill(j_selected_sorted[1].E());
				
				//separazione tra jet
			
				delta_theta=abs(j_selected_sorted[0].Theta()-j_selected_sorted[1].Theta());
				delta_R=j_selected_sorted[0].DeltaR(j_selected_sorted[1]);
				delta_phi=abs(j_selected_sorted[0].DeltaPhi(j_selected_sorted[1]));
				Dtheta->Fill(delta_theta);
				Dr->Fill(delta_R);
				DPhi->Fill(delta_phi);
				
				}
			
			}
  		
  		//just a check 
  		if(nj>1){
		for(k=0;k<nj;k++){
  			pT[k]=j[k].Pt();
  			}
  		TMath::Sort(nj,pT,index);
  		
  		for(k=0;k<nj;k++){
  			j_selected_sorted[k]=j[index[k]];
  			}
		h=j_selected_sorted[0]+j_selected_sorted[1];
		h_invmass->Fill(h.M());
		}
		for(k=0;k<100;k++){
			j_selected_sorted[k].SetPxPyPzE(0,0,0,0);
			j_selected[k].SetPxPyPzE(0,0,0,0);
			j[k].SetPxPyPzE(0,0,0,0);
			}
		} // end of loop on fChain->GetEntries()
		
		
	h_invmass->Draw();
	
	TCanvas *c1=new TCanvas();
	gStyle->SetOptStat("nemrou");
	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
	legend2->AddEntry(pt_j1, "Fist jet (pT sorted)", "l");
	legend2->AddEntry(pt_j2, "Second jet (pT sorted)", "l");
	legend2->SetTextSize(0.03);
	pt_j2->SetTitle("pT distribution for the 2 highest pT jets");
	pt_j2->GetXaxis()->SetTitle("pT [GeV]");
	pt_j2->Draw("s hist");
	pt_j1->Draw("sames hist");
	legend2->Draw();
	
		
	double eff=0.,err_eff=0;	
	cout<<"JET SELECTION"<<endl;
	cout<<"Numero jet totali (su tutti gli eventi simulati): "<<n_jet_tot<<endl;
	
	cout<<"Numero jet dopo eta-cut: "<<n_jet_etacut_tot<< "+/- "<<sqrt(n_jet_etacut_tot)<<endl;
	eff=double(n_jet_etacut_tot)/double(n_jet_tot);
	err_eff=sqrt(eff*(1-eff)/double(n_jet_tot));
	cout<<"Efficienza assoluta per eta-cut su jet: "<<eff<<" +/- "<<err_eff<<endl;
	
	cout<<"Numero jet dopo pt-cut+eta-cut: "<<n_jet_ptcut_tot<< "+/- "<<sqrt(n_jet_ptcut_tot)<<endl;
	eff=double(n_jet_ptcut_tot)/double(n_jet_tot);
	err_eff=sqrt(eff*(1-eff)/double(n_jet_tot));
	cout<<"Efficienza assoluta per pt+eta-cut su jet: "<<eff<<" +/- "<<err_eff<<endl;
	
	cout<<"Numero eventi con almeno due jet che passano eta-cut: "<<n_events_etacut<<" +/- "<<sqrt(n_events_etacut)<<endl;
	
	cout<<"Numero eventi con almeno due jet che passano pt-cut+eta-cut: "<<n_events_ptcut<<" +/- "<<sqrt(n_events_ptcut)<<endl;
	
	cout<<"EVENT SELECTION"<<endl;
	cout<<"Numero eventi che hanno almeno 2 jet che superano selezioni pt+eta: "<<n_events_ptcut<<endl;
	eff=double(n_events_ptcut)/double(n_events_tot);
	err_eff=sqrt(eff*(1-eff)/double(n_events_tot));
	cout<<"Efficienza per richiesta almeno due jet: "<<eff<<" +/- "<<err_eff<<endl;
	
	cout<<"Numero eventi che hanno almeno 2 jet e massa invariante della coppia a più alto pt <250 GeV: "<<n_events_minvcut<<endl;
	eff=double(n_events_minvcut)/double(n_events_tot);
	err_eff=sqrt(eff*(1-eff)/double(n_events_tot));
	cout<<"efficienza 2jet+minv: "<<eff<<" +/- "<<err_eff<<endl;
	
	
	
	


if (plot_ID==0){
    TCanvas* n_jets = new TCanvas("n_jets","n_jets",900,600);
   // n_jets_ev->Scale(1./n_jets_ev->Integral());
    n_jets_ev->Draw("s hist");
    }


if (plot_ID==1){
    TCanvas* eta_canv = new TCanvas("eta_canv","eta_canv",900,600);
    eta_histo->Draw();
    }

if (plot_ID==2){
    TCanvas* pt_canv = new TCanvas("pt_canv","pt_canv",900,600);
    pt_histo->Draw();
    }

if (plot_ID==3){
    TCanvas* phi_canv = new TCanvas("phi_canv","phi_canv",900,600);
    phi_histo->Draw();
    }

if (plot_ID==4){
    TCanvas* E_canv = new TCanvas("E_canv","E_canv",900,600);
    E_histo->Draw();
    }

if (plot_ID==5){
    TCanvas* minv_canv = new TCanvas("minv_canv","minv_canv",900,600);
    m_inv->Scale(1./m_inv->Integral());
    m_inv->Draw("s hist");
    }


	TFile* f=new TFile("/home/paola/Scrivania/jet_cc.root","RECREATE");
	   if ( f->IsOpen() ) cout << "File opened successfully" << endl;
	n_jets_ev->Write();
	m_inv->Write();
	eta_histo->Write();
	E_histo->Write();
	pt_histo->Write();
	theta_histo->Write();
	phi_histo->Write();
	
	//separazione tra jet
	Dtheta->Write();
	Dr->Write();
	DPhi->Write();

	//higgs
	eta_higgs->Write();
	E_higgs->Write();
	pt_higgs->Write();
	theta_higgs->Write();
	phi_higgs->Write();


}



