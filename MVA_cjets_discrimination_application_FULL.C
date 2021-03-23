//MVA application

#define CvsB 0.02
#define CvsL 0.05
#define NUM_EVENTS 10000
#define MAX_JETS 30000


#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/CrossValidation.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/MethodCuts.h"


using namespace TMVA;
//choice is either 'signal' or 'b' 'light'

void MVA_cjets_discrimination_application_FULL(string choice ="signal"){
	gStyle->SetOptStat("nemrou");
	TString weight_path_B ="/home/paola/Scrivania/macro/Trees_tagger/BDT/dataset_10k_3/weights/TMVA_new_BDT.weights.xml";
	TString weight_path_L ="/home/paola/Scrivania/macro/Trees_tagger/BDT/dataset_10k_3_light/weights/TMVA_new_BDT.weights.xml";
	
	// Output file
	
	cout<<"we"<<endl;
	// Get the signal and background trees from TFile source(s);
	TFile *filein;
	TTree *myTree;
	
	if(choice=="signal"){
		// SIGNAL tree	
		filein = new TFile("/home/paola/Scrivania/macro/Trees_tagger/c_tree_forSEL4_10k.root");
		myTree = (TTree*)filein->Get("tree_tagging_var");
		cout<<"===>>  Taking the "<<choice<<" sample"<<endl;
		}
	else if(choice=="b"){
		//BKG from b tree
		filein = new TFile("/home/paola/Scrivania/macro/Trees_tagger/b_tree_forSEL4_10k.root");
		myTree = (TTree*)filein->Get("tree_tagging_var");
		cout<<"===>>  Taking the "<<choice<<" sample"<<endl;
		}
	else if(choice=="light"){
		filein = new TFile("/home/paola/Scrivania/macro/Trees_tagger/o_tree_forSEL4_10k.root");
		myTree = (TTree*)filein->Get("tree_tagging_var");
		cout<<"===>>  Taking the "<<choice<<" sample"<<endl;
		}
		
	else{
		cout<<"Invalid argument; Choose 'signal' or 'background'"<<endl;
		
		}
	
	//for analysis
	double j_px;
	double j_py;
	double j_pz;
	double j_E;
	int j_flavour;
	int event_id;
	////////////////////////////
	
	int vertex_cat;
	int n_SV;
	int n_trk_from_SV;
	double fd_sig_xy_SV;
	double fd_sig_xyz_SV;
	double corrected_mass_SV;
	double MassEnergyFraction_SV;
	double Boost_SV;
	double energy_ratio_SV;
	double dR_SV_JET;
	double SIP_sig_xy_TRK1;
	double SIP_sig_xy_TRK2;
	double SIP_sig_xyz_TRK1;
	double SIP_sig_xyz_TRK2;
	double SIP_sig_xy_TRK_above_trh;
	double SIP_sig_xyz_TRK_above_trh;
	double dR_jet_TRK1;
	int lepton_cat;
	double pt_lj_LEP1;
	double Et_trk_Jet;
	
	
	myTree->SetBranchAddress("j_px", &j_px);
   	myTree->SetBranchAddress("j_py", &j_py);
   	myTree->SetBranchAddress("j_pz", &j_pz);
   	myTree->SetBranchAddress("j_E", &j_E);
   	myTree->SetBranchAddress("j_flavour",&j_flavour);
	myTree->SetBranchAddress("event_id",&event_id);
	
    	myTree->SetBranchAddress("vertex_cat",&vertex_cat);
    	myTree->SetBranchAddress("n_SV",&n_SV);
    	myTree->SetBranchAddress("n_trk_from_SV",&n_trk_from_SV);
    	myTree->SetBranchAddress("fd_sig_xy_SV",&fd_sig_xy_SV);
    	myTree->SetBranchAddress("fd_sig_xyz_SV",&fd_sig_xyz_SV);
    	myTree->SetBranchAddress("corrected_mass_SV",&corrected_mass_SV);
    	myTree->SetBranchAddress("MassEnergyFraction_SV",&MassEnergyFraction_SV);
    	myTree->SetBranchAddress("Boost_SV",&Boost_SV);
    	myTree->SetBranchAddress("energy_ratio_SV",&energy_ratio_SV);
    	myTree->SetBranchAddress("dR_SV_JET",&dR_SV_JET);
    	myTree->SetBranchAddress("SIP_sig_xy_TRK1",&SIP_sig_xy_TRK1);
    	myTree->SetBranchAddress("SIP_sig_xy_TRK2",&SIP_sig_xy_TRK2);
    	myTree->SetBranchAddress("SIP_sig_xyz_TRK1",&SIP_sig_xyz_TRK1);
    	myTree->SetBranchAddress("SIP_sig_xyz_TRK2",&SIP_sig_xyz_TRK2);
    	myTree->SetBranchAddress("SIP_sig_xy_TRK_above_trh",&SIP_sig_xy_TRK_above_trh);
    	myTree->SetBranchAddress("SIP_sig_xyz_TRK_above_trh",&SIP_sig_xyz_TRK_above_trh);
    	myTree->SetBranchAddress("dR_jet_TRK1",&dR_jet_TRK1);
    	myTree->SetBranchAddress("lepton_cat",&lepton_cat);
    	myTree->SetBranchAddress("pt_lj_LEP1",&pt_lj_LEP1);
    	myTree->SetBranchAddress("Et_trk_Jet",& Et_trk_Jet);
    	
    	//TMVA Reader initialization
    	// This loads the library
    	TMVA::Tools::Instance();
  	TMVA::Reader *reader_B = new TMVA::Reader( "!Color:!Silent" );
   	float  forBDTevaluation1; 
   	float  forBDTevaluation2; 
  	float  forBDTevaluation3; 
   	float  forBDTevaluation4; 
   	float  forBDTevaluation5; 
   	float  forBDTevaluation6; 
   	float  forBDTevaluation7; 
   	float  forBDTevaluation8; 
   	float  forBDTevaluation9; 
   	float  forBDTevaluation10;
   	float  forBDTevaluation11;
   	float  forBDTevaluation12;
   	float  forBDTevaluation13;
   	float  forBDTevaluation14;
   	float  forBDTevaluation15;
   	float  forBDTevaluation16;
   	float  forBDTevaluation17;
   	float  forBDTevaluation18;
   	float  forBDTevaluation19;
   	float  forBDTevaluation20;
   	
   	// Variables declaration
   	///////////////////////////////////////////////////////////////////
   	//CvsB								 //
   	///////////////////////////////////////////////////////////////////
   	reader_B->AddVariable("vertex_cat", &forBDTevaluation1 );
   	reader_B->AddVariable("n_SV", &forBDTevaluation2 );
    	reader_B->AddVariable("n_trk_from_SV", &forBDTevaluation3);
    	reader_B->AddVariable("fd_sig_xy_SV", &forBDTevaluation4);
    	reader_B->AddVariable("fd_sig_xyz_SV", &forBDTevaluation5);
    	reader_B->AddVariable("corrected_mass_SV", &forBDTevaluation6);
    	reader_B->AddVariable("MassEnergyFraction_SV", &forBDTevaluation7);
    	reader_B->AddVariable("Boost_SV", &forBDTevaluation8);
    	reader_B->AddVariable("energy_ratio_SV", &forBDTevaluation9);
    	reader_B->AddVariable("dR_SV_JET", &forBDTevaluation10);
    	reader_B->AddVariable("SIP_sig_xy_TRK1", &forBDTevaluation11);
    	reader_B->AddVariable("SIP_sig_xy_TRK2", &forBDTevaluation12);
    	reader_B->AddVariable("SIP_sig_xyz_TRK1", &forBDTevaluation13);
    	reader_B->AddVariable("SIP_sig_xyz_TRK2", &forBDTevaluation14);
    	reader_B->AddVariable("SIP_sig_xy_TRK_above_trh", &forBDTevaluation15);
    	reader_B->AddVariable("SIP_sig_xyz_TRK_above_trh", &forBDTevaluation16);
    	reader_B->AddVariable("dR_jet_TRK1", &forBDTevaluation17);
    	reader_B->AddVariable("lepton_cat", &forBDTevaluation18);
    	reader_B->AddVariable("pt_lj_LEP1", &forBDTevaluation19);
    	reader_B->AddVariable("Et_trk_Jet", &forBDTevaluation20);
    	
    	///////////////////////////////////////////////////////////////////
   	//CvsL							 	 //
   	///////////////////////////////////////////////////////////////////
    	
    	TMVA::Reader *reader_L = new TMVA::Reader( "!Color:!Silent" );
    	reader_L->AddVariable("SIP_sig_xy_TRK1", &forBDTevaluation11);
    	reader_L->AddVariable("SIP_sig_xy_TRK2", &forBDTevaluation12);
    	reader_L->AddVariable("SIP_sig_xyz_TRK1", &forBDTevaluation13);
    	reader_L->AddVariable("SIP_sig_xyz_TRK2", &forBDTevaluation14);
    	reader_L->AddVariable("SIP_sig_xy_TRK_above_trh", &forBDTevaluation15);
    	reader_L->AddVariable("SIP_sig_xyz_TRK_above_trh", &forBDTevaluation16);
    	reader_L->AddVariable("dR_jet_TRK1", &forBDTevaluation17);
    	reader_L->AddVariable("lepton_cat", &forBDTevaluation18);
    	reader_L->AddVariable("pt_lj_LEP1", &forBDTevaluation19);
    	reader_L->AddVariable("Et_trk_Jet", &forBDTevaluation20);
    	
    	
    	
	reader_B->TMVA::Reader::BookMVA("BDT", weight_path_B);
	reader_L->TMVA::Reader::BookMVA("BDT", weight_path_L);
  	cout<<"Using weights in "<<weight_path_B<<" for CvsB" <<endl;
  	cout<<"Using weights in "<<weight_path_L<<" for CvsL" <<endl;
  	
  	long int nentries,i,j,k;
  	nentries=myTree->GetEntries();
  	
  	int selected_B=0, rejected_B=0;
  	int selected_L=0, rejected_L=0;
  	
  	int event_matrix[NUM_EVENTS][10];
  	
  	int event_vector[NUM_EVENTS];
  	
  	double pT[10]={0};
  	int index[10];
  	
  	int jet_tag[MAX_JETS];
  	int jet_selected_tag[10];
  	
  	TLorentzVector jets[MAX_JETS];
  	TLorentzVector jet_selected[10];
  	TLorentzVector jet_selected_sorted[10];
  	TLorentzVector h;
  	
  	int COUNTER_1=0;
  	int COUNTER_2=0;
  	
  	//before BDT selection
  	int event_matrix_BeforeCut[NUM_EVENTS][10];
  	int event_vector_BeforeCut[NUM_EVENTS];
  	TLorentzVector jets_BeforeCut[MAX_JETS];
  	TLorentzVector jet_selected_BeforeCut[10];
  	TLorentzVector jet_selected_sorted_BeforeCut[10];
  	double pT_BeforeCut[10]={0};
  	int index_BeforeCut[10];
  	/////////////////////////////////////////////
  	
  	for(i=0;i<10;i++){
  		pT[i]=0;
  		pT[i]=0;
  		}
  	
  	for(i=0;i<NUM_EVENTS;i++){
  		event_vector[i]=0;
  		event_vector_BeforeCut[i]=0;
  		for(j=0;j<10;j++){
  			event_matrix[i][j]=-1;
  			event_matrix_BeforeCut[i][j]=-1;
  			}
  		}
  	
  	//Start event Loop
  	for (long int jentry=0; jentry<nentries; jentry++) {
  	
		myTree->GetEntry(jentry);
		
		forBDTevaluation1 = vertex_cat;
		forBDTevaluation2 = n_SV;
		forBDTevaluation3 = n_trk_from_SV;
		forBDTevaluation4 = fd_sig_xy_SV;
		forBDTevaluation5 = fd_sig_xyz_SV;
		forBDTevaluation6 = corrected_mass_SV;
		forBDTevaluation7 = MassEnergyFraction_SV;
		forBDTevaluation8 = Boost_SV;
		forBDTevaluation9 = energy_ratio_SV;
		forBDTevaluation10 = dR_SV_JET;
		forBDTevaluation11 = SIP_sig_xy_TRK1;
		forBDTevaluation12 = SIP_sig_xy_TRK2;
		forBDTevaluation13 = SIP_sig_xyz_TRK1;
		forBDTevaluation14 = SIP_sig_xyz_TRK2;
		forBDTevaluation15 = SIP_sig_xy_TRK_above_trh;
		forBDTevaluation16 = SIP_sig_xyz_TRK_above_trh;
		forBDTevaluation17 = dR_jet_TRK1;
		forBDTevaluation18 = lepton_cat;
		forBDTevaluation19 = pt_lj_LEP1;
		forBDTevaluation20 = Et_trk_Jet;
    	
    		double BDT_decision_B = reader_B->EvaluateMVA("BDT");
    		double BDT_decision_L = reader_L->EvaluateMVA("BDT");
    		cout<<"Evaluation done"<<endl;
    		
    		if( BDT_decision_B >= CvsB){//B-selected
           		selected_B++;
           		if(BDT_decision_L >= CvsL){//L-selected
           			selected_L++;
		   		//num of the jet
		   		event_matrix[event_id][event_vector[event_id]]=jentry;
		   		event_vector[event_id]++;
		   		jets[jentry].SetPxPyPzE(j_px,j_py,j_pz,j_E);
		   		jet_tag[jentry]=j_flavour;
		   		//cout<<"E: "<<j_E<<endl;
		   		}
		   	else if(BDT_decision_L < CvsL){//L-rejected
           			rejected_L++;
           			}
           		}
      		else if ( BDT_decision_B < CvsB ){ //B-rejected
           		rejected_B++;
           		}
           		
           	//no BDT selection
           	event_matrix_BeforeCut[event_id][event_vector_BeforeCut[event_id]]=jentry;
           	event_vector_BeforeCut[event_id]++;
           	jets_BeforeCut[jentry].SetPxPyPzE(j_px,j_py,j_pz,j_E);
  		
		}
		
	std::cout << "==> B-Selected jets " << selected_B << std::endl;
    	std::cout << "==> B-Rejected jets " << rejected_B << std::endl;
    	std::cout << "==> L-Selected jets " << selected_L << std::endl;
    	std::cout << "==> L-Rejected jets " << rejected_L << std::endl;
    	std::cout << "==> TMVAApplication is done!" << std::endl;
    	
    	
    	//massa invariante dopo tutte le selezioni
    	TH1F* m_inv=new TH1F("m_inv","m_inv",100,0,250);
	m_inv->SetTitle("Invariant mass of jet pair");
   	m_inv->GetXaxis()->SetTitle("M_{H} [GeV]");
   	m_inv->GetYaxis()->SetTitle("Number of events");
   	
   	//massa invariante dopo tutte le selezioni ECCETTO BDT CUT
    	TH1F* m_inv_BeforeCut=new TH1F("m_inv_BeforeCut","m_inv_BeforeCut",100,0,250);
	m_inv_BeforeCut->SetTitle("Invariant mass of jet pair _ BeforeBDTCut");
   	m_inv_BeforeCut->GetXaxis()->SetTitle("M_{H} [GeV]");
   	m_inv_BeforeCut->GetYaxis()->SetTitle("Number of events");
   	
   	
   	//PLOT DI CONTROLLO
   	//info nel full range
   	TH1F* pt1_histo=new TH1F("pt1_histo","pt1_histo",60,0,300);
	pt1_histo->SetTitle("Transverse momentum of the first jet");
    	pt1_histo->GetXaxis()->SetTitle("p_{T} [GeV]");
    	
    	TH1F* pt2_histo=new TH1F("pt2_histo","pt2_histo",60,0,300);
	pt2_histo->SetTitle("Transverse momentum of the second jet");
    	pt2_histo->GetXaxis()->SetTitle("p_{T} [GeV]");
    	
    	TH1F* m1_histo=new TH1F("m1_histo","m1_histo",60,0,300);
	m1_histo->SetTitle("Mass of the first jet");
    	m1_histo->GetXaxis()->SetTitle("M_{j1} [GeV]");
    	
    	TH1F* m2_histo=new TH1F("m2_histo","m2_histo",60,0,300);
	m2_histo->SetTitle("Mass of the second jet");
    	m2_histo->GetXaxis()->SetTitle("M_{j2} [GeV]");
    	
    	TH1F* E1_histo=new TH1F("E1_histo","E1_histo",80,0,500);
	E1_histo->SetTitle("Energy of the first jet");
    	E1_histo->GetXaxis()->SetTitle("E_{j1} [GeV]");
    	
    	TH1F* E2_histo=new TH1F("E2_histo","E2_histo",80,0,500);
	E2_histo->SetTitle("Energy of the second jet");
    	E2_histo->GetXaxis()->SetTitle("E_{j2} [GeV]");
    	
    	TH1F* flav1_histo=new TH1F("flav1_histo","flav1_histo",2,-0.5,1.5);
	flav1_histo->SetTitle("Flavour of the first jet");
    	flav1_histo->GetXaxis()->SetTitle("");
    	
    	TH1F* flav2_histo=new TH1F("flav2_histo","flav2_histo",2,-0.5,1.5);
	flav2_histo->SetTitle("Flavour of the second jet");
    	flav2_histo->GetXaxis()->SetTitle("");
    	
    	//info nel range di massa invariante <100 GeV
    	TH1F* pt1_histo_r=new TH1F("pt1_histo_r","pt1_histo_r",60,0,300);
	pt1_histo_r->SetTitle("Transverse momentum of the first jet");
    	pt1_histo_r->GetXaxis()->SetTitle("p_{T} [GeV]");
    	
    	TH1F* pt2_histo_r=new TH1F("pt2_histo_r","pt2_histo_r",60,0,300);
	pt2_histo_r->SetTitle("Transverse momentum of the second jet");
    	pt2_histo_r->GetXaxis()->SetTitle("p_{T} [GeV]");
    	
    	TH1F* m1_histo_r=new TH1F("m1_histo_r","m1_histo_r",60,0,300);
	m1_histo_r->SetTitle("Mass of the first jet");
    	m1_histo_r->GetXaxis()->SetTitle("M_{j1} [GeV]");
    	
    	TH1F* m2_histo_r=new TH1F("m2_histo_r","m2_histo_r",60,0,300);
	m2_histo_r->SetTitle("Transverse momentum of the second jet");
    	m2_histo_r->GetXaxis()->SetTitle("M_{j2} [GeV]");
    	
    	TH1F* E1_histo_r=new TH1F("E1_histo_r","E1_histo_r",80,0,500);
	E1_histo_r->SetTitle("Energy of the first jet");
    	E1_histo_r->GetXaxis()->SetTitle("E_{j1} [GeV]");
    	
    	TH1F* E2_histo_r=new TH1F("E2_histo_r","E2_histo_r",80,0,500);
	E2_histo_r->SetTitle("Energy of the second jet");
    	E2_histo_r->GetXaxis()->SetTitle("E_{j2} [GeV]");
    	
    	TH1F* flav1_histo_r=new TH1F("flav1_histo_r","flav1_histo_r",2,-0.5,1.5);
	flav1_histo_r->SetTitle("Flavour of the first jet");
    	flav1_histo_r->GetXaxis()->SetTitle("");
    	
    	TH1F* flav2_histo_r=new TH1F("flav2_histo_r","flav2_histo_r",2,-0.5,1.5);
	flav2_histo_r->SetTitle("Flavour of the second jet");
    	flav2_histo_r->GetXaxis()->SetTitle("");
    	
    	
    	
    	TH2F* pT_vs_flav=new TH2F("colz","colz", 2,-0.5,1.5, 60,0,300);
    	
    	
    	//EFFICIENZA C-TAGGING VS PT
    	double pT_bins[8]={5,20,30,50,70,100,135};
    	TH1F* pt_after=new TH1F("pt_after","pt_after",6,pT_bins);
    	TH1F* pt_before=new TH1F("pt_before","pt_before",6,pT_bins);
    	
    	//No BDT Cut
    	for(k=0;k<NUM_EVENTS;k++){
		if(event_vector_BeforeCut[k]>1){
			//cout<<"------------------------------"<<endl;
			//cout<<"Num event: "<<k<<endl;
			
			for(i=0;i<10;i++){
				index[i]=-1;
				}
			
			for(i=0;i<event_vector_BeforeCut[k];i++){
			
  				jet_selected_BeforeCut[i]=jets_BeforeCut[event_matrix_BeforeCut[k][i]];
  				pT_BeforeCut[i]=jets_BeforeCut[event_matrix_BeforeCut[k][i]].Pt();
  				//cout<<"pt: "<<pT_BeforeCut[i]<<endl;
  				
  				//efficiency//
				pt_before->Fill(pT_BeforeCut[i]);
				/////////////////////////////////
  				}
  			TMath::Sort(event_vector_BeforeCut[k],pT_BeforeCut,index);
  		
  			for(i=0;i<event_vector_BeforeCut[k];i++){
  				jet_selected_sorted_BeforeCut[i]=jet_selected_BeforeCut[index[i]];
  				}
  			//cout<<" "<<endl;
  			//cout<<"first: "<<jet_selected_sorted_BeforeCut[0].Pt()<<endl;
  			//cout<<"second: "<<jet_selected_sorted_BeforeCut[1].Pt()<<endl;
  			h=jet_selected_sorted_BeforeCut[0]+jet_selected_sorted_BeforeCut[1];
  			if(h.M()<250){
  				m_inv_BeforeCut->Fill(h.M());
  				}
  			}
  		}
    	
	for(k=0;k<NUM_EVENTS;k++){
		if(event_vector[k]>1){
			//cout<<"------------------------------"<<endl;
			//cout<<"Num event: "<<k<<endl;
			COUNTER_1++;
			
			for(i=0;i<10;i++){
				index[i]=-1;
				}
			
			for(i=0;i<event_vector[k];i++){
			
  				jet_selected[i]=jets[event_matrix[k][i]];
  				jet_selected_tag[i]=jet_tag[event_matrix[k][i]];
  				pT[i]=jets[event_matrix[k][i]].Pt();
  				//cout<<"pt: "<<pT[i]<<endl;
  				
  				//efficiency//
				pt_after->Fill(pT[i]);
				/////////////////////////////////
  				}
  			TMath::Sort(event_vector[k],pT,index);
  		
  			for(i=0;i<event_vector[k];i++){
  				jet_selected_sorted[i]=jet_selected[index[i]];
  				}
  			/*cout<<" "<<endl;
  			cout<<"first: "<<jet_selected_sorted[0].Pt()<<endl;
  			cout<<"second: "<<jet_selected_sorted[1].Pt()<<endl;*/
  			h=jet_selected_sorted[0]+jet_selected_sorted[1];
  			if(h.M()<250){
  				COUNTER_2++;
  				m_inv->Fill(h.M());
  				pt1_histo->Fill(jet_selected_sorted[0].Pt());
  				pt2_histo->Fill(jet_selected_sorted[1].Pt());
  				E1_histo->Fill(jet_selected_sorted[0].E());
  				E2_histo->Fill(jet_selected_sorted[1].E());
  				
  				flav1_histo->Fill(jet_selected_tag[index[0]]);
  				flav2_histo->Fill(jet_selected_tag[index[1]]);
  				
  				m1_histo->Fill(jet_selected_sorted[0].M());
  				m2_histo->Fill(jet_selected_sorted[1].M());
  				
  				}
  			if(h.M()<100){
  				pt1_histo_r->Fill(jet_selected_sorted[0].Pt());
  				pt2_histo_r->Fill(jet_selected_sorted[1].Pt());
  				
  				E1_histo_r->Fill(jet_selected_sorted[0].E());
  				E2_histo_r->Fill(jet_selected_sorted[1].E());
  				
  				flav1_histo_r->Fill(jet_selected_tag[index[0]]);
  				flav2_histo_r->Fill(jet_selected_tag[index[1]]);
  				
  				m1_histo_r->Fill(jet_selected_sorted[0].M());
  				m2_histo_r->Fill(jet_selected_sorted[1].M());
  				
  				pT_vs_flav->Fill(jet_selected_tag[index[0]],jet_selected_sorted[0].Pt());
  				pT_vs_flav->Fill(jet_selected_tag[index[1]],jet_selected_sorted[1].Pt());
  		
  				}
  			}
  		}
  	
  	//m_inv->Scale(1./m_inv->Integral());
  	//m_inv->Sumw2();
    	m_inv->Draw("s hist");
    	
    	TCanvas *eff_can=new TCanvas();
    	TEfficiency* tag_eff = 0;
	if(TEfficiency::CheckConsistency(*pt_after,*pt_before)){
		tag_eff = new TEfficiency(*pt_after,*pt_before);
		}
	tag_eff->SetLineColor(kBlack);
  	tag_eff->SetMarkerStyle(21);
  	tag_eff->SetTitle("Tagging efficiency VS pT;pT;efficiency");
  	tag_eff->Draw("ap");
    	
    	TCanvas * myc0=new TCanvas();
    	m_inv_BeforeCut->Draw("s hist");
    	
    	TCanvas * myc1=new TCanvas();
    	pt1_histo->Draw("s hist");
    	
    	TCanvas * myc2=new TCanvas();
    	pt2_histo->Draw("s hist");
    	
    	TCanvas * myc3=new TCanvas();
    	flav1_histo->Draw("s hist");
    	
    	TCanvas * myc4=new TCanvas();
    	flav2_histo->Draw("s hist");
    	
    	TCanvas * myc5=new TCanvas();
    	pT_vs_flav->Draw("colz");
    	
    	TFile* f=new TFile(Form("/home/paola/Scrivania/InvMass_%s_FULL.root",choice.c_str()),"RECREATE");
	if ( f->IsOpen() ) cout << "File opened successfully" << endl;
	m_inv->Write();
	m_inv_BeforeCut->Write();
	//efficiency plot
	tag_eff->Write();
	//control plots
	pt1_histo->Write();
	pt2_histo->Write();
	E1_histo->Write();
	E2_histo->Write();
    	flav1_histo->Write();
    	flav2_histo->Write();
    	m1_histo->Write();
    	m2_histo->Write();
    	//restricted region
    	pt1_histo_r->Write();
	pt2_histo_r->Write();
	E1_histo_r->Write();
	E2_histo_r->Write();
    	flav1_histo_r->Write();
    	flav2_histo_r->Write();
    	m1_histo_r->Write();
    	m2_histo_r->Write();
    	
    	cout<<"Numero eventi totali: "<<NUM_EVENTS<<endl;
  	cout<<"Numero eventi con almeno 2 jet che passano le selezioni: "<<COUNTER_1<<endl;
  	cout<<"Numero eventi con almeno mjj<250 : "<<COUNTER_2<<endl;
    	
    	
	
}
