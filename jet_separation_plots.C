#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

void jet_separation_plots(){

	TFile *f = new TFile("/home/paola/Scrivania/reco_jet_higgs_bb_1000.root");

	TH1F * Dtheta_bb = (TH1F*)f->Get("Dtheta");
	Dtheta_bb->SetLineWidth(2);
	Dtheta_bb->Scale(1./Dtheta_bb->Integral());
	
	TH1F * DR_bb = (TH1F*)f->Get("Dr");
	DR_bb->SetLineWidth(2);
	DR_bb->Scale(1./DR_bb->Integral());

	TH1F * DPhi_bb = (TH1F*)f->Get("DPhi");
	DPhi_bb->SetLineWidth(2);
	DPhi_bb->Scale(1./DPhi_bb->Integral());
	
	TFile *f2 = new TFile("/home/paola/Scrivania/reco_jet_higgs_cc_1000.root");
	
	TH1F * Dtheta_cc = (TH1F*)f2->Get("Dtheta");
	Dtheta_cc->SetLineWidth(2);
	Dtheta_cc->SetLineColor(kRed);
	Dtheta_cc->Scale(1./Dtheta_cc->Integral());
	
	TH1F * DR_cc = (TH1F*)f2->Get("Dr");
	DR_cc->SetLineWidth(2);
	DR_cc->SetLineColor(kRed);
	DR_cc->Scale(1./DR_cc->Integral());

	TH1F * DPhi_cc = (TH1F*)f2->Get("DPhi");
	DPhi_cc->SetLineWidth(2);
	DPhi_cc->SetLineColor(kRed);
	DPhi_cc->Scale(1./DPhi_cc->Integral());
	
	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
     	legend2->AddEntry(Dtheta_cc, "H #rightarrow c#bar{c}", "l");
     	legend2->AddEntry(Dtheta_bb, "H #rightarrow b#bar{b}", "l");
     	legend2->SetTextSize(0.04);
	
	TCanvas *myc1 = new TCanvas("myc1", "canvas", 900,600);
	
	DR_bb->Draw("s hist");
              myc1->Update();
              TPaveStats* stat1 = (TPaveStats*)DR_bb->GetListOfFunctions()->FindObject("stats");
              stat1->SetX1NDC(0.55);
	      stat1->SetX2NDC(0.85);
	      stat1->SetY1NDC(0.54);
	      stat1->SetY2NDC(0.69);
	      stat1->SetTextColor(kBlue);
	      
              DR_cc->Draw("sames hist");
              myc1->Update();
              TPaveStats* stat2 = (TPaveStats*)DR_cc->GetListOfFunctions()->FindObject("stats");
	      stat2->SetX1NDC(0.55);
              stat2->SetX2NDC(0.85);
              stat2->SetY1NDC(0.7);
              stat2->SetY2NDC(0.85);
              stat2->SetTextColor(kRed);
              
        legend2->Draw();
              
        
	TCanvas *myc2 = new TCanvas("myc2", "canvas", 900,600);
	
              Dtheta_cc->Draw("s hist");
              myc2->Update();
              TPaveStats* stat20 = (TPaveStats*)Dtheta_cc->GetListOfFunctions()->FindObject("stats");
	      stat20->SetX1NDC(0.55);
              stat20->SetX2NDC(0.85);
              stat20->SetY1NDC(0.7);
              stat20->SetY2NDC(0.85);
              stat20->SetTextColor(kRed);
              
              Dtheta_bb->Draw("sames hist");
              myc2->Update();
              TPaveStats* stat10 = (TPaveStats*)Dtheta_bb->GetListOfFunctions()->FindObject("stats");
              stat10->SetX1NDC(0.55);
	      stat10->SetX2NDC(0.85);
	      stat10->SetY1NDC(0.54);
	      stat10->SetY2NDC(0.69);
	      stat10->SetTextColor(kBlue);
              
	legend2->Draw();

	TCanvas *myc10 = new TCanvas("myc10", "canvas", 900,600);
	
	DPhi_bb->Draw("s hist");
              myc10->Update();
              TPaveStats* stat11 = (TPaveStats*)DPhi_bb->GetListOfFunctions()->FindObject("stats");
              stat11->SetX1NDC(0.55);
	      stat11->SetX2NDC(0.85);
	      stat11->SetY1NDC(0.54);
	      stat11->SetY2NDC(0.69);
	      stat11->SetTextColor(kBlue);
	      
              DPhi_cc->Draw("sames hist");
              myc10->Update();
              TPaveStats* stat21 = (TPaveStats*)DPhi_cc->GetListOfFunctions()->FindObject("stats");
	      stat21->SetX1NDC(0.55);
              stat21->SetX2NDC(0.85);
              stat21->SetY1NDC(0.7);
              stat21->SetY2NDC(0.85);
              stat21->SetTextColor(kRed);
              
        legend2->Draw();

	
	
}

	
