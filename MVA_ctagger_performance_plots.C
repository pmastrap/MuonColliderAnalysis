//MVA application

#define CvsB 0.02
#define CvsL 0.05
#define CvsL_min -0.2
#define CvsL_max 0.4
#define CvsB_min -0.44
#define CvsB_max 0.22
#define step 0.001
#define eff_step 0.02


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

void MVA_ctagger_performance_plots(){
	gStyle->SetOptStat(0);
	TString weight_path_B ="/home/paola/Scrivania/macro/Trees_tagger/BDT/dataset_10k_3/weights/TMVA_new_BDT.weights.xml";
	TString weight_path_L ="/home/paola/Scrivania/macro/Trees_tagger/BDT/dataset_10k_3_light/weights/TMVA_new_BDT.weights.xml";
	
	// Output file
	
	// Get the signal and background trees from TFile source(s);
	TFile *filein_c,*filein_b,*filein_o;
	TTree *myTree_c,*myTree_b,*myTree_o;
	
	// SIGNAL tree	
	filein_c = new TFile("/home/paola/Scrivania/macro/Trees_tagger/c_tree_forMVA3_10k.root");
	myTree_c = (TTree*)filein_c->Get("tree_tagging_var");

	//BKG from b tree
	filein_b = new TFile("/home/paola/Scrivania/macro/Trees_tagger/b_tree_forMVA3_10k.root");
	myTree_b = (TTree*)filein_b->Get("tree_tagging_var");
	
	//BKG from light tree
	filein_o = new TFile("/home/paola/Scrivania/macro/Trees_tagger/o_tree_forMVA3_10k.root");
	myTree_o = (TTree*)filein_o->Get("tree_tagging_var");
	
	long int nentries;
	long int jentry,i,k,m;
	
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
   	
   	//histograms
   	TH1F* h_CvsL_output_c= new TH1F("h_CvsL_output_c","h_CvsL_output_c",(CvsL_max-CvsL_min)/step,CvsL_min,CvsL_max);
   	h_CvsL_output_c->SetLineColor(kGreen+2);
   	h_CvsL_output_c->SetLineWidth(2);
   	TH1F* h_CvsL_output_b= new TH1F("h_CvsL_output_b","h_CvsL_output_b",(CvsL_max-CvsL_min)/step,CvsL_min,CvsL_max);
   	h_CvsL_output_b->SetLineColor(kRed);
   	h_CvsL_output_b->SetLineWidth(2);
   	TH1F* h_CvsL_output_o= new TH1F("h_CvsL_output_o","h_CvsL_output_o",(CvsL_max-CvsL_min)/step,CvsL_min,CvsL_max);
   	h_CvsL_output_o->SetLineColor(kBlue);
   	h_CvsL_output_o->SetLineWidth(2);
   	
   	TH1F* h_CvsB_output_c= new TH1F("h_CvsB_output_c","h_CvsB_output_c",(CvsB_max-CvsB_min)/step,CvsB_min,CvsB_max);
   	h_CvsB_output_c->SetLineColor(kGreen+2);
   	h_CvsB_output_c->SetLineWidth(2);
   	TH1F* h_CvsB_output_b= new TH1F("h_CvsB_output_b","h_CvsB_output_b",(CvsB_max-CvsB_min)/step,CvsB_min,CvsB_max);
   	h_CvsB_output_b->SetLineColor(kRed);
   	h_CvsB_output_b->SetLineWidth(2);
   	TH1F* h_CvsB_output_o= new TH1F("h_CvsB_output_o","h_CvsB_output_o",(CvsB_max-CvsB_min)/step,CvsB_min,CvsB_max);
   	h_CvsB_output_o->SetLineColor(kBlue);
   	h_CvsB_output_o->SetLineWidth(2);
   	
   	
   	
   	TH1F* eff_c_B = new TH1F("eff_c_B","eff_c_B",(CvsB_max-CvsB_min)/step,CvsB_min,CvsB_max);
   	eff_c_B->SetLineColor(kGreen+2);
   	eff_c_B->SetLineWidth(2);
	TH1F* eff_b_B = new TH1F("eff_b_B","eff_b_B",(CvsB_max-CvsB_min)/step,CvsB_min,CvsB_max);
	eff_b_B->SetLineColor(kRed);
   	eff_b_B->SetLineWidth(2);
	TH1F* eff_o_B = new TH1F("eff_o_B","eff_o_B",(CvsB_max-CvsB_min)/step,CvsB_min,CvsB_max);
	eff_o_B->SetLineColor(kBlue);
   	eff_o_B->SetLineWidth(2);
	
	TH1F* eff_c_L = new TH1F("eff_c_L","eff_c_L",(CvsL_max-CvsL_min)/step,CvsL_min,CvsL_max);
	eff_c_L->SetLineColor(kGreen+2);
   	eff_c_L->SetLineWidth(2);
	TH1F* eff_b_L = new TH1F("eff_b_L","eff_b_L",(CvsL_max-CvsL_min)/step,CvsL_min,CvsL_max);
	eff_b_L->SetLineColor(kRed);
   	eff_b_L->SetLineWidth(2);
	TH1F* eff_o_L = new TH1F("eff_o_L","eff_o_L",(CvsL_max-CvsL_min)/step,CvsL_min,CvsL_max);
	eff_o_L->SetLineColor(kBlue);
   	eff_o_L->SetLineWidth(2);
   	
   	
   	//correlation between CvsB and CvsL
   	TGraph* correlation_c = new TGraph();
   	correlation_c->SetMarkerStyle(8);
   	correlation_c->SetMarkerColor(kGreen+2);
   	correlation_c->SetMarkerSize(0.2);
   	TGraph* correlation_b = new TGraph();
   	correlation_b->SetMarkerStyle(8);
   	correlation_b->SetMarkerColor(kRed);
   	correlation_b->SetMarkerSize(0.2);
   	TGraph* correlation_o = new TGraph();
   	correlation_o->SetMarkerStyle(8);
   	correlation_o->SetMarkerColor(kBlue);
   	correlation_o->SetMarkerSize(0.2);
   	
   	
   	//Charm eff contours
   	TH2F* c_cont = new TH2F("c_cont","c_cont",1.1/eff_step,0,1.1,1.1/eff_step,0,1.1);
	
    	//TMVA Reader initialization
    	// This loads the library
    	TMVA::Tools::Instance();
   	
   	// Variables declaration
   	///////////////////////////////////////////////////////////////////
   	//CvsB								 //
   	///////////////////////////////////////////////////////////////////
   	TMVA::Reader *reader_B = new TMVA::Reader( "!Color:!Silent" );
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
  	
  	
  	//C
  	myTree_c->SetBranchAddress("vertex_cat",&vertex_cat);
    	myTree_c->SetBranchAddress("n_SV",&n_SV);
    	myTree_c->SetBranchAddress("n_trk_from_SV",&n_trk_from_SV);
    	myTree_c->SetBranchAddress("fd_sig_xy_SV",&fd_sig_xy_SV);
    	myTree_c->SetBranchAddress("fd_sig_xyz_SV",&fd_sig_xyz_SV);
    	myTree_c->SetBranchAddress("corrected_mass_SV",&corrected_mass_SV);
    	myTree_c->SetBranchAddress("MassEnergyFraction_SV",&MassEnergyFraction_SV);
    	myTree_c->SetBranchAddress("Boost_SV",&Boost_SV);
    	myTree_c->SetBranchAddress("energy_ratio_SV",&energy_ratio_SV);
    	myTree_c->SetBranchAddress("dR_SV_JET",&dR_SV_JET);
    	myTree_c->SetBranchAddress("SIP_sig_xy_TRK1",&SIP_sig_xy_TRK1);
    	myTree_c->SetBranchAddress("SIP_sig_xy_TRK2",&SIP_sig_xy_TRK2);
    	myTree_c->SetBranchAddress("SIP_sig_xyz_TRK1",&SIP_sig_xyz_TRK1);
    	myTree_c->SetBranchAddress("SIP_sig_xyz_TRK2",&SIP_sig_xyz_TRK2);
    	myTree_c->SetBranchAddress("SIP_sig_xy_TRK_above_trh",&SIP_sig_xy_TRK_above_trh);
    	myTree_c->SetBranchAddress("SIP_sig_xyz_TRK_above_trh",&SIP_sig_xyz_TRK_above_trh);
    	myTree_c->SetBranchAddress("dR_jet_TRK1",&dR_jet_TRK1);
    	myTree_c->SetBranchAddress("lepton_cat",&lepton_cat);
    	myTree_c->SetBranchAddress("pt_lj_LEP1",&pt_lj_LEP1);
    	myTree_c->SetBranchAddress("Et_trk_Jet",& Et_trk_Jet);

  	nentries=myTree_c->GetEntries();
  	
  	/////////////////////////////////////////////
  	
  	//Start event Loop
  	for (jentry=0; jentry<nentries; jentry++) {
  	
		myTree_c->GetEntry(jentry);
		
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
    		
    		h_CvsB_output_c->Fill(reader_B->EvaluateMVA("BDT"));
    		h_CvsL_output_c->Fill(reader_L->EvaluateMVA("BDT"));
    		
    		correlation_c->SetPoint(jentry,reader_L->EvaluateMVA("BDT"),reader_B->EvaluateMVA("BDT"));
    		//double BDT_decision_B = reader_B->EvaluateMVA("BDT");
    		//double BDT_decision_L = reader_L->EvaluateMVA("BDT");
    		//cout<<"Evaluation done"<<endl;
		}
		
	//B
  	myTree_b->SetBranchAddress("vertex_cat",&vertex_cat);
    	myTree_b->SetBranchAddress("n_SV",&n_SV);
    	myTree_b->SetBranchAddress("n_trk_from_SV",&n_trk_from_SV);
    	myTree_b->SetBranchAddress("fd_sig_xy_SV",&fd_sig_xy_SV);
    	myTree_b->SetBranchAddress("fd_sig_xyz_SV",&fd_sig_xyz_SV);
    	myTree_b->SetBranchAddress("corrected_mass_SV",&corrected_mass_SV);
    	myTree_b->SetBranchAddress("MassEnergyFraction_SV",&MassEnergyFraction_SV);
    	myTree_b->SetBranchAddress("Boost_SV",&Boost_SV);
    	myTree_b->SetBranchAddress("energy_ratio_SV",&energy_ratio_SV);
    	myTree_b->SetBranchAddress("dR_SV_JET",&dR_SV_JET);
    	myTree_b->SetBranchAddress("SIP_sig_xy_TRK1",&SIP_sig_xy_TRK1);
    	myTree_b->SetBranchAddress("SIP_sig_xy_TRK2",&SIP_sig_xy_TRK2);
    	myTree_b->SetBranchAddress("SIP_sig_xyz_TRK1",&SIP_sig_xyz_TRK1);
    	myTree_b->SetBranchAddress("SIP_sig_xyz_TRK2",&SIP_sig_xyz_TRK2);
    	myTree_b->SetBranchAddress("SIP_sig_xy_TRK_above_trh",&SIP_sig_xy_TRK_above_trh);
    	myTree_b->SetBranchAddress("SIP_sig_xyz_TRK_above_trh",&SIP_sig_xyz_TRK_above_trh);
    	myTree_b->SetBranchAddress("dR_jet_TRK1",&dR_jet_TRK1);
    	myTree_b->SetBranchAddress("lepton_cat",&lepton_cat);
    	myTree_b->SetBranchAddress("pt_lj_LEP1",&pt_lj_LEP1);
    	myTree_b->SetBranchAddress("Et_trk_Jet",& Et_trk_Jet);
	
  	nentries=myTree_b->GetEntries();
  	
  	/////////////////////////////////////////////
  	
  	//Start event Loop
  	for (jentry=0; jentry<nentries; jentry++) {
  	
		myTree_b->GetEntry(jentry);
		
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

    		h_CvsB_output_b->Fill(reader_B->EvaluateMVA("BDT"));
    		h_CvsL_output_b->Fill(reader_L->EvaluateMVA("BDT"));
    		
    		correlation_b->SetPoint(jentry,reader_L->EvaluateMVA("BDT"),reader_B->EvaluateMVA("BDT"));
    		//double BDT_decision_B = reader_B->EvaluateMVA("BDT");
    		//double BDT_decision_L = reader_L->EvaluateMVA("BDT");
    		//cout<<"Evaluation done"<<endl;
		}
		
	//Light
  	myTree_o->SetBranchAddress("vertex_cat",&vertex_cat);
    	myTree_o->SetBranchAddress("n_SV",&n_SV);
    	myTree_o->SetBranchAddress("n_trk_from_SV",&n_trk_from_SV);
    	myTree_o->SetBranchAddress("fd_sig_xy_SV",&fd_sig_xy_SV);
    	myTree_o->SetBranchAddress("fd_sig_xyz_SV",&fd_sig_xyz_SV);
    	myTree_o->SetBranchAddress("corrected_mass_SV",&corrected_mass_SV);
    	myTree_o->SetBranchAddress("MassEnergyFraction_SV",&MassEnergyFraction_SV);
    	myTree_o->SetBranchAddress("Boost_SV",&Boost_SV);
    	myTree_o->SetBranchAddress("energy_ratio_SV",&energy_ratio_SV);
    	myTree_o->SetBranchAddress("dR_SV_JET",&dR_SV_JET);
    	myTree_o->SetBranchAddress("SIP_sig_xy_TRK1",&SIP_sig_xy_TRK1);
    	myTree_o->SetBranchAddress("SIP_sig_xy_TRK2",&SIP_sig_xy_TRK2);
    	myTree_o->SetBranchAddress("SIP_sig_xyz_TRK1",&SIP_sig_xyz_TRK1);
    	myTree_o->SetBranchAddress("SIP_sig_xyz_TRK2",&SIP_sig_xyz_TRK2);
    	myTree_o->SetBranchAddress("SIP_sig_xy_TRK_above_trh",&SIP_sig_xy_TRK_above_trh);
    	myTree_o->SetBranchAddress("SIP_sig_xyz_TRK_above_trh",&SIP_sig_xyz_TRK_above_trh);
    	myTree_o->SetBranchAddress("dR_jet_TRK1",&dR_jet_TRK1);
    	myTree_o->SetBranchAddress("lepton_cat",&lepton_cat);
    	myTree_o->SetBranchAddress("pt_lj_LEP1",&pt_lj_LEP1);
    	myTree_o->SetBranchAddress("Et_trk_Jet",& Et_trk_Jet);
    	
  	nentries=myTree_o->GetEntries();
  	
  	/////////////////////////////////////////////
  	
  	//Start event Loop
  	for (jentry=0; jentry<nentries; jentry++) {
  	
		myTree_o->GetEntry(jentry);
		
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

    		h_CvsB_output_o->Fill(reader_B->EvaluateMVA("BDT"));
    		h_CvsL_output_o->Fill(reader_L->EvaluateMVA("BDT"));
    		
    		correlation_o->SetPoint(jentry,reader_L->EvaluateMVA("BDT"),reader_B->EvaluateMVA("BDT"));
    		//double BDT_decision_B = reader_B->EvaluateMVA("BDT");
    		//double BDT_decision_L = reader_L->EvaluateMVA("BDT");
    		//cout<<"Evaluation done"<<endl;

		}
	int sum1=0;
	int sum2=0;
	int sum3=0;	
	
	double x[100000];
	double y[100000];
	
	//eff plots
	cout<<h_CvsB_output_c->GetXaxis()->GetNbins()<<endl;
	for(i=0;i<h_CvsB_output_c->GetXaxis()->GetNbins();i++){
		sum1=0;
		sum2=0;
		sum3=0;
		for(k=i;k<h_CvsB_output_c->GetXaxis()->GetNbins();k++){
			sum1=sum1+h_CvsB_output_c->GetBinContent(k);
			}
		for(k=i;k<h_CvsB_output_c->GetXaxis()->GetNbins();k++){
			sum2=sum2+h_CvsB_output_b->GetBinContent(k);
			}
		for(k=i;k<h_CvsB_output_c->GetXaxis()->GetNbins();k++){
			sum3=sum3+h_CvsB_output_o->GetBinContent(k);
			}
		x[i]=double(sum1)/double(myTree_c->GetEntries());
		y[i]=double(sum2)/double(myTree_b->GetEntries());
		eff_c_B->SetBinContent(i,double(sum1)/double(myTree_c->GetEntries()));
		eff_b_B->SetBinContent(i,double(sum2)/double(myTree_b->GetEntries()));
		eff_o_B->SetBinContent(i,double(sum3)/double(myTree_o->GetEntries()));
		}
		
	//ROC curve
   	TGraph* ROC_CvsB=new TGraph(h_CvsB_output_c->GetXaxis()->GetNbins(),x,y);
   	ROC_CvsB->SetTitle("CvsB tagger ROC curve");
	ROC_CvsB->GetXaxis()->SetTitle("c jet efficiency");
	ROC_CvsB->GetYaxis()->SetTitle("b jet rejection");
	ROC_CvsB->SetLineWidth(3);
	ROC_CvsB->SetLineColor(kGreen+2);
	ROC_CvsB->SetMarkerColor(kGreen+2);
	
	//eff plots
	cout<<h_CvsL_output_c->GetXaxis()->GetNbins()<<endl;
	for(i=0;i<h_CvsL_output_c->GetXaxis()->GetNbins();i++){
		sum1=0;
		sum2=0;
		sum3=0;
		for(k=i;k<h_CvsL_output_c->GetXaxis()->GetNbins();k++){
			sum1=sum1+h_CvsL_output_c->GetBinContent(k);
			}
		for(k=i;k<h_CvsL_output_c->GetXaxis()->GetNbins();k++){
			sum2=sum2+h_CvsL_output_b->GetBinContent(k);
			}
		for(k=i;k<h_CvsL_output_c->GetXaxis()->GetNbins();k++){
			sum3=sum3+h_CvsL_output_o->GetBinContent(k);
			}
		x[i]=double(sum1)/double(myTree_c->GetEntries());
		y[i]=double(sum3)/double(myTree_o->GetEntries());
		eff_c_L->SetBinContent(i,double(sum1)/double(myTree_c->GetEntries()));
		eff_b_L->SetBinContent(i,double(sum2)/double(myTree_b->GetEntries()));
		eff_o_L->SetBinContent(i,double(sum3)/double(myTree_o->GetEntries()));
		}
	
	//ROC curve
   	TGraph* ROC_CvsL=new TGraph(h_CvsL_output_c->GetXaxis()->GetNbins(),x,y);
   	ROC_CvsL->SetTitle("CvsL tagger ROC curve");
	ROC_CvsL->GetXaxis()->SetTitle("c jet efficiency");
	ROC_CvsL->GetYaxis()->SetTitle("light jet rejection");
	ROC_CvsL->SetLineWidth(3);
	ROC_CvsL->SetLineColor(kGreen+2);
	ROC_CvsL->SetMarkerColor(kGreen+2);
	
	
	//Charm eff countour
	double light_eff=0.;
	double b_eff=0.;
	double c_eff=0.;
	//TGraph2D *g = new TGraph2D();
	TGraph* gr[8];
	for(m=0;m<8;m++){
		gr[m]= new TGraph();
		}
	double epsilon = 0.001;
	double eff[8]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	double n_points[8]={0};
	for(i=0;i<h_CvsL_output_c->GetXaxis()->GetNbins();i++){
		for(k=0;k<h_CvsB_output_c->GetXaxis()->GetNbins();k++){
			cout<<"(CvsL,CvsB) = ("<<(i*step)+CvsL_min<<","<<(k*step)+CvsB_min<<")"<<endl;
			light_eff=eff_o_L->GetBinContent(i)*eff_o_B->GetBinContent(k);
			b_eff=eff_b_L->GetBinContent(i)*eff_b_B->GetBinContent(k);
			c_eff=eff_c_L->GetBinContent(i)*eff_c_B->GetBinContent(k);
			cout<<"b_eff: "<<b_eff<<endl;
			cout<<"c_eff: "<<c_eff<<endl;
			cout<<"light_eff: "<<light_eff<<endl;
			c_cont->SetBinContent(int(b_eff/eff_step)+1,int(light_eff/eff_step)+1,c_eff);
			//fill charm efficiency contour
			for(m=0;m<8;m++){
				if(c_eff>eff[m]-epsilon && c_eff <eff[m]+epsilon){
					gr[m]->SetPoint(n_points[m], b_eff, light_eff);
					n_points[m]++;
					}
				}
			//g->SetPoint(i*h_CvsB_output_c->GetXaxis()->GetNbins()+k, b_eff, light_eff, c_eff);
			}
		}

	auto legend0 = new TLegend(0.67,0.74,0.86,0.90);
	legend0->AddEntry(h_CvsB_output_b, "b-jets(sample H #rightarrow b#bar{b})" , "l");
	legend0->AddEntry(h_CvsB_output_c, "c-jets(sample H #rightarrow c#bar{c})" , "l");
	legend0->AddEntry(h_CvsB_output_o, "light-jets" , "l");
	legend0->SetTextSize(0.030);
	legend0->SetBorderSize(0);
	legend0->SetTextFont(42);
	
	auto legend1 = new TLegend(0.15,0.15,0.45,0.38);
	legend1->SetHeader("  Charm efficiency contours","L");
	legend1->SetNColumns(2);
	for(i=0;i<8;i++){
		legend1->AddEntry(gr[i], Form("#epsilon^{c} = %.1f",eff[i]) ,"p");
		}
	legend1->SetTextSize(0.032);
	legend1->SetBorderSize(0);
	
	//Draw charm efficiency contour
	int W = 800;
  	int H = 600;
	TCanvas* cc = new TCanvas("cc","cc",10,10,W,H);
	//My_Beautiful_Canvas_BDTperform(cc);
	Color_t colors[8]={kBlack,kRed,kMagenta+1,kGreen+1,kBlue+1,kYellow-1,kYellow,kRed+3};
	gPad->SetLogy();
	gPad->SetLogx();
	gPad->SetGrid(1,1);
	gr[0]->GetXaxis()->SetTitle("b-jet misid. probability");
	gr[0]->GetYaxis()->SetTitle("light-jet misid. probability");
	gr[0]->SetLineColor(colors[0]);
	gr[0]->SetLineWidth(1);
	gr[0]->SetMarkerStyle(8);
	gr[0]->SetMarkerColor(colors[0]);
	gr[0]->SetMarkerSize(0.5);
	gr[0]->GetXaxis()->SetLimits(1e-3,2);
	gr[0]->GetYaxis()->SetRangeUser(1e-3,10.);
	gr[0]->Draw("APL");
	cc->Update();
	for(k=1;k<8;k++){
		gr[k]->SetLineColor(colors[k]);
		gr[k]->SetLineWidth(1);
		gr[k]->SetMarkerStyle(8);
		gr[k]->SetMarkerColor(colors[k]);
		gr[k]->SetMarkerSize(0.5);
		gr[k]->Draw("samePL");
		}
	legend1->Draw();
	
	
	
	auto legend00 = new TLegend(0.25,0.65,0.80,0.70);
	legend00->SetNColumns(3);
	legend00->AddEntry(h_CvsB_output_b, "b-jets" , "f");
	legend00->AddEntry(h_CvsB_output_c, "c-jets" , "f");
	legend00->AddEntry(h_CvsB_output_o, "light-jets" , "f");
	legend00->SetTextSize(0.03);
	legend00->SetBorderSize(0);
	TCanvas* c_correlation=new TCanvas();
	correlation_b->GetXaxis()->SetLimits(-0.26,0.4); 
	correlation_b->GetHistogram()->SetMaximum(0.7);
	correlation_b->GetXaxis()->SetTitle("CvsL discriminator");
	correlation_b->GetYaxis()->SetTitle("CvsB discriminator");
	correlation_b->Draw("AP");
	correlation_o->Draw("Psame");
	correlation_c->Draw("Psame");
	legend00->Draw();
	
	auto legend000= new TLegend(0.55,0.35,0.8,0.60);
	legend000->AddEntry(ROC_CvsB, "c tagger CvsB" , "l");
	legend000->SetTextSize(0.05);
	legend000->SetBorderSize(0);
	TCanvas* c_ROC_B=new TCanvas();
	gPad->SetLogy();
	ROC_CvsB->GetXaxis()->SetTitle("c-jet efficiency");
	ROC_CvsB->GetYaxis()->SetTitle("b-jet misid. probability");
	ROC_CvsB->Draw("APL");
	legend000->Draw();
	
	auto legend0000= new TLegend(0.55,0.35,0.8,0.60);
	legend0000->AddEntry(ROC_CvsL, "c tagger CvsL" , "l");
	legend0000->SetTextSize(0.05);
	legend0000->SetBorderSize(0);
	TCanvas* c_ROC_L=new TCanvas();
	gPad->SetLogy();
	ROC_CvsL->GetXaxis()->SetTitle("c-jet efficiency");
	ROC_CvsL->GetYaxis()->SetTitle("b-jet misid. probability");
	ROC_CvsL->Draw("APL");
	legend0000->Draw();
	
	TCanvas* c_eff_CvsB=new TCanvas();
	eff_c_B->SetMaximum(1.4);
	eff_c_B->GetXaxis()->SetTitle("c tagger CvsB efficiency");
	eff_c_B->GetYaxis()->SetTitle("BDT cut");
	eff_c_B->Draw("");
	eff_b_B->Draw("same");
	eff_o_B->Draw("same");
	legend0->Draw();
	
	TCanvas* c_eff_CvsL=new TCanvas();
	eff_c_L->SetMaximum(1.4);
	eff_c_L->GetXaxis()->SetTitle("c tagger CvsL efficiency");
	eff_c_L->GetYaxis()->SetTitle("BDT cut");
	eff_c_L->Draw("");
	eff_b_L->Draw("same");
	eff_o_L->Draw("same");
	legend0->Draw();
	
	TCanvas* c_CvsB = new TCanvas();
	h_CvsB_output_c->Rebin(10);
	h_CvsB_output_b->Rebin(10);
	h_CvsB_output_o->Rebin(10);
	h_CvsB_output_c->GetXaxis()->SetTitle("CvsB discriminator");
	h_CvsB_output_c->GetYaxis()->SetTitle("Jets / 0.01 units");
	h_CvsB_output_c->Scale(1/h_CvsB_output_c->Integral());
	h_CvsB_output_b->Scale(1/h_CvsB_output_b->Integral());
	h_CvsB_output_o->Scale(1/h_CvsB_output_o->Integral());
	h_CvsB_output_c->GetYaxis()->SetRangeUser(0,0.20);
	h_CvsB_output_c->Draw("hist");
	h_CvsB_output_b->Draw("same hist");
	h_CvsB_output_o->Draw("same hist");
	legend0->Draw();
	
	TCanvas* c_CvsL = new TCanvas();
	h_CvsL_output_c->Rebin(10);
	h_CvsL_output_b->Rebin(10);
	h_CvsL_output_o->Rebin(10);
	h_CvsL_output_c->GetXaxis()->SetTitle("CvsL discriminator");
	h_CvsL_output_c->GetYaxis()->SetTitle("Jets / 0.01 units");
	h_CvsL_output_c->Scale(1/h_CvsL_output_c->Integral());
	h_CvsL_output_b->Scale(1/h_CvsL_output_b->Integral());
	h_CvsL_output_o->Scale(1/h_CvsL_output_o->Integral());
	h_CvsL_output_c->GetYaxis()->SetRangeUser(0,0.15);
	h_CvsL_output_c->Draw("hist");
	h_CvsL_output_b->Draw("same hist");
	h_CvsL_output_o->Draw("same hist");
	legend0->Draw();
	
	}
