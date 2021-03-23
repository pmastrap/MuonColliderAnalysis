//MVA for b-c jets discrimination

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
#include "TMVA/TMVAGui.h"

//#include "T3M_common.h"

using namespace TMVA;

void MVA_cjets_discrimination(){

	// Output file
	TFile *fout = new TFile("TMVA_output_10k_3.root", "RECREATE");
	
	// Get the signal and background trees from TFile source(s);
	
        // SIGNAL tree	
	TFile *f_sig = new TFile("/home/paola/Scrivania/macro/Trees_tagger/c_tree_forMVA3_10k.root");
	TTree *signalTree     = (TTree*)f_sig->Get("tree_tagging_var");

	//BKG tree
	TFile *f_bkg = new TFile("/home/paola/Scrivania/macro/Trees_tagger/b_tree_forMVA3_10k.root");
	TTree *bkgTree     = (TTree*)f_bkg->Get("tree_tagging_var");
	
	//weights
	Double_t sigWeight1  = 1.0;
    	Double_t bkgWeight1 = 1.0;
    
    
    	Factory *factory = new Factory("TMVA_new", fout, "");
    	DataLoader *dataloader = new DataLoader("dataset_10k_3");
    	
    	dataloader->AddSignalTree(signalTree,sigWeight1);
   	dataloader->AddBackgroundTree(bkgTree,bkgWeight1);
   	
   	 // Variables declaration
   	dataloader->AddVariable("vertex_cat", 'I');
    	dataloader->AddVariable("n_SV", 'I');
    	dataloader->AddVariable("n_trk_from_SV", 'I');
    	dataloader->AddVariable("fd_sig_xy_SV", 'D');
    	dataloader->AddVariable("fd_sig_xyz_SV", 'D');
    	dataloader->AddVariable("corrected_mass_SV", 'D');
    	dataloader->AddVariable("MassEnergyFraction_SV", 'D');
    	dataloader->AddVariable("Boost_SV", 'D');
    	dataloader->AddVariable("energy_ratio_SV", 'D');
    	dataloader->AddVariable("dR_SV_JET", 'D');
    	dataloader->AddVariable("SIP_sig_xy_TRK1", 'D');
    	dataloader->AddVariable("SIP_sig_xy_TRK2", 'D');
    	dataloader->AddVariable("SIP_sig_xyz_TRK1", 'D');
    	dataloader->AddVariable("SIP_sig_xyz_TRK2", 'D');
    	dataloader->AddVariable("SIP_sig_xy_TRK_above_trh", 'D');
    	dataloader->AddVariable("SIP_sig_xyz_TRK_above_trh", 'D');
    	dataloader->AddVariable("dR_jet_TRK1", 'D');
    	dataloader->AddVariable("lepton_cat", 'I');
    	dataloader->AddVariable("pt_lj_LEP1", 'D');
    	dataloader->AddVariable("Et_trk_Jet", 'D');
    	
    	
    	//cuts
    	TCut cutS="";
    	TCut cutB="";
    	
    	dataloader->PrepareTrainingAndTestTree( cutS, cutB,"SplitMode=Random:NormMode=NumEvents:!V" );
	
	// Booking of MVA methods : BDT
        factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT","!H:!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=50");

        // Training the MVA methods
        factory->TrainAllMethods();
        
        // Testing the MVA methods
        factory->TestAllMethods();
        
        // Evaluating the MVA methods
        factory->EvaluateAllMethods();
		
	// Save the output
    fout->Close();
    
    std::cout << "==> Wrote root file: " << fout->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
    
    delete factory;
    delete dataloader;
    
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()){
        TMVAGui("TMVA_output_10k_3.root");
    	}
	
}
