#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"



void create_rootfile_jets_higgs(Int_t plot_ID=0)
{
cout << "Plot IDs are the following" << endl;
cout << "0: number of reconstructed jets" << endl;
cout << "1: pseudorapidity of jets in the event" << endl;
cout << "2: transverse momentum of jets in the event" << endl;
cout << "3: phi of jets in the event" << endl;
cout << "4: energy of jets in the event" << endl;
cout << "5: invariant mass of jet pair associated to Higgs with highest pT" << endl;

TH1F* n_jets_ev=new TH1F("n_jets_ev","n_jets_ev",10,-0.5,9.5);
	n_jets_ev->SetTitle("Number of jets in the event");
    	n_jets_ev->GetXaxis()->SetTitle("N_{jets}");
TH1F* m_inv=new TH1F("m_inv","m_inv",25,0,250);
	m_inv->SetTitle("Invariant mass of jet pair");
   	m_inv->GetXaxis()->SetTitle("M_{H} [GeV]");
TH1F* eta_histo=new TH1F("eta_histo","eta_histo",30,-5,5);
	eta_histo->SetTitle("Jets pseudorapidity");
    	eta_histo->GetXaxis()->SetTitle("#eta");
TH1F* E_histo=new TH1F("E_histo","E_histo",100,0,600);
	E_histo->SetTitle("Jets energy");
    	E_histo->GetXaxis()->SetTitle("Energy [GeV]");
TH1F* pt_histo=new TH1F("pt_histo","pt_histo",60,0,300);
	pt_histo->SetTitle("Jets transverse momentum");
    	pt_histo->GetXaxis()->SetTitle("p_{T} [GeV]");
TH1F* theta_histo=new TH1F("theta_histo","theta_histo",45,-0.5,4);
	theta_histo->SetTitle("Jets polar angle");
	theta_histo->GetXaxis()->SetTitle("#theta [rad]");
TH1F* phi_histo=new TH1F("phi_histo","phi_histo",20,-4,4);
	phi_histo->SetTitle("Jets azimuthal angle");
    	phi_histo->GetXaxis()->SetTitle("#phi [rad]");
//separazione tra jet
TH1F* Dtheta=new TH1F("Dtheta","Dtheta",50,-0.5,3);
	Dtheta->SetTitle("Jets #Delta#theta");
	Dtheta->GetXaxis()->SetTitle("#Delta#theta [rad]");
TH1F* Dr=new TH1F("Dr","Dr",50,0,8);
	Dr->SetTitle("Jets #DeltaR");
	Dr->GetXaxis()->SetTitle("#DeltaR [rad]");
TH1F* DPhi=new TH1F("DPhi","DPhi",50,-0.5,3.5);
	DPhi->SetTitle("Jets #Delta#phi");
	DPhi->GetXaxis()->SetTitle("#Delta#phi [rad]");
//Higgs Plots
TH1F* eta_higgs=new TH1F("eta_higgs","eta_higgs",30,-5,5);
	eta_higgs->SetTitle("Higgs pseudorapidity");
    	eta_higgs->GetXaxis()->SetTitle("#eta_{H}");
TH1F* E_higgs=new TH1F("E_higgs","E_higgs",100,0,600);
	E_higgs->SetTitle("Higgs energy");
    	E_higgs->GetXaxis()->SetTitle("Energy_{H} [GeV]");
TH1F* pt_higgs=new TH1F("pt_higgs","pt_higgs",70,0,700);
	pt_higgs->SetTitle("Higgs transverse momentum");
    	pt_higgs->GetXaxis()->SetTitle("p_{T} [GeV]");
TH1F* theta_higgs=new TH1F("theta_higgs","theta_higgs",45,-0.5,4);
	theta_higgs->SetTitle("Higgs polar angle");
	theta_higgs->GetXaxis()->SetTitle("#theta_{H} [rad]");
TH1F* phi_higgs=new TH1F("phi_higgs","phi_higgs",20,-4,4);
	phi_higgs->SetTitle("Higgs azimuthal angle");
    	phi_higgs->GetXaxis()->SetTitle("#phi_{H} [rad]");

double mh=125.;

Int_t  nj;
   Float_t         jmox[100];   //[njet]
   Float_t         jmoy[100];   //[njet]
   Float_t         jmoz[100];   //[njet]
   Float_t         jene[100];   //[njet]

TLorentzVector j[100]; 
TLorentzVector sum;
TLorentzVector h;
TVector3 h_p;
TVector3 j_p;

// List of branches
   TBranch        *b_njet;   
   TBranch        *b_jmox;   
   TBranch        *b_jmoy;   
   TBranch        *b_jmoz;     
   TBranch        *b_jene;   

   TChain* fChain = new TChain("fChain");
   fChain->Add("/home/paola/Scrivania/xml_aggiornati/ntuple_bb_1000_newgeom.root/MyLCTuple");
   

   fChain->SetBranchAddress("nj", &nj, &b_njet);
   fChain->SetBranchAddress("jmox", jmox, &b_jmox);
   fChain->SetBranchAddress("jmoy", jmoy, &b_jmoy);
   fChain->SetBranchAddress("jmoz", jmoz, &b_jmoz);
   fChain->SetBranchAddress("jene", jene, &b_jene);

double inv_mass[1000][1000];

double delta_min=10000000;
double delta_m1;
int jet1_num=0;
int jet2_num=0;

double delta_R=0;
double delta_theta=0;
double delta_phi=0;

for(unsigned int ientry=0; ientry< fChain->GetEntries(); ++ientry)
  {
   fChain->GetEntry(ientry);
   
   	//RIEMPIO IL LORENTZ VECTOR
	
	for (int k=0;k<nj; k++)
  	 { 
  	 	j[k].SetPx(jmox[k]);
  	 	j[k].SetPy(jmoy[k]);
  	 	j[k].SetPz(jmoz[k]);
  	 	j[k].SetE(jene[k]);
  	 	
  	 }

   n_jets_ev->Fill(nj);
   for (int k=0;k<nj; k++)
   { 
	  j_p = j[k].Vect();
          eta_histo->Fill(j_p.PseudoRapidity());
	  pt_histo->Fill(j_p.Perp());
          theta_histo->Fill(j_p.Theta());
          phi_histo->Fill(j_p.Phi());
          E_histo->Fill(j[k].E());
   }

    if(nj>1) //selection of events with more than one jet


{
     for (int i=0; i<nj-1; i++) //Calculation of invariant mass 
     {
        for (int k=i+1; k<nj;k++)
        {
		 sum=j[i]+j[k];
                 inv_mass[i][k]= sum.M();    
        }

     }
    double delta_min=10000000;
      for (int i=0; i<nj-1; i++)
     { 
        for (int k=i+1; k<nj;k++)
        {
		delta_m1=TMath::Abs(mh-inv_mass[i][k]);
		if (delta_m1<delta_min ){delta_min=delta_m1; jet1_num=i;jet2_num=k;}
                
        } //end k
     } //end i

	//una volta individuata la coppia di jet posso ricavare info su higgs e riempire il TLorentzVector h
	h=j[jet1_num]+j[jet2_num];
	h_p = h.Vect();
	m_inv->Fill(inv_mass[jet1_num][jet2_num]);
	eta_higgs->Fill(h_p.PseudoRapidity());
	pt_higgs->Fill(h_p.Perp());
        theta_higgs->Fill(h_p.Theta());
        phi_higgs->Fill(h_p.Phi());
        E_higgs->Fill(h.E());
        
        
        //separazione tra jet
	delta_R=sqrt((j[jet1_num].PseudoRapidity()-j[jet2_num].PseudoRapidity())*(j[jet1_num].PseudoRapidity()-j[jet2_num].PseudoRapidity())+(j[jet1_num].Phi()-j[jet2_num].Phi())*(j[jet1_num].Phi()-j[jet2_num].Phi())); 
	delta_theta=abs(j[jet1_num].Theta()-j[jet2_num].Theta());
	delta_phi=abs(j[jet1_num].Phi()-j[jet2_num].Phi());
	
	Dtheta->Fill(delta_theta);
	Dr->Fill(delta_R);
	if(delta_phi>=TMath::Pi()){
		DPhi->Fill(2*TMath::Pi()-delta_phi);}
	else DPhi->Fill(delta_phi);
	

} //end nj>1

for(int i=0; i<1000; i++) {for(int k=0;k<1000;k++){inv_mass[i][k]=0;}}
delta_min=10000000;
jet1_num=0;
jet2_num=0;


} // end of loop on fChain->GetEntries()


if (plot_ID==0){
    TCanvas* n_jets = new TCanvas("n_jets","n_jets",900,600);
    n_jets_ev->Scale(1./n_jets_ev->Integral());
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


	TFile* f=new TFile("/home/paola/Scrivania/reco_jet_higgs_bb_1000.root","RECREATE");
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



