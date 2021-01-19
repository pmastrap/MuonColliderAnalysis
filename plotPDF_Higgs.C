#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

void plotPDF_Higgs(){

	TFile *f = new TFile("/home/paola/Scrivania/reco_jet_higgs_cc_1000.root");
    	
	TH1F* eta = (TH1F*)f->Get("eta_higgs");
	eta->SetTitle("Higgs pseudorapidity");
    	eta->GetXaxis()->SetTitle("#eta");
    	eta->SetLineColor(kRed);
    	
	TH1F* E = (TH1F*)f->Get("E_higgs");
	E->SetTitle("Higgs energy");
    	E->GetXaxis()->SetTitle("Energy [GeV]");
    	E->SetLineColor(kRed);
    	
	TH1F* pt = (TH1F*)f->Get("pt_higgs");
	pt->SetTitle("Higgs transverse momentum");
    	pt->GetXaxis()->SetTitle("p_{T} [GeV]");
    	pt->SetLineColor(kRed);
    	
    	TH1F* phi= (TH1F*)f->Get("phi_higgs");
	phi->SetTitle("Higgs azimuthal angle");
    	phi->GetXaxis()->SetTitle("#phi [rad]");
    	phi->SetLineColor(kRed);
    	
    	TCanvas *myc1 = new TCanvas("myc1", "canvas", 0,0,1000,1300);
              myc1->Divide(2,2);
              myc1->cd(1);
              //gPad->SetGrid(1,1);
              myc1->SetFillColor(0);
	      
	//HISTO BIB
	
	TFile *f2 = new TFile("/home/paola/Scrivania/reco_jet_higgs_bb_1000.root");
	
   	
	TH1F* eta_BIB = (TH1F*)f2->Get("eta_higgs");
	eta_BIB->SetTitle("Higgs pseudorapidity");
    	eta_BIB->GetXaxis()->SetTitle("#eta");
    	
    	
	TH1F* E_BIB = (TH1F*)f2->Get("E_higgs");
	E_BIB->SetTitle("Higgs energy");
    	E_BIB->GetXaxis()->SetTitle("Energy [GeV]");
    	
    	
	TH1F* pt_BIB = (TH1F*)f2->Get("pt_higgs");
	pt_BIB->SetTitle("Higgs transverse momentum");
    	pt_BIB->GetXaxis()->SetTitle("p_{T} [GeV]");
    	
    	
    	TH1F* phi_BIB= (TH1F*)f2->Get("phi_higgs");
	phi_BIB->SetTitle("Higgs azimuthal angle");
    	phi_BIB->GetXaxis()->SetTitle("#phi [rad]");
    	
    	
    	
    	auto legend2 = new TLegend(0.65,0.55,0.96,0.7);
     	//legend2->SetHeader("Distance");
     	legend2->AddEntry(pt, "H #rightarrow c#bar{c}", "l");
     	legend2->AddEntry(pt_BIB, "H #rightarrow b#bar{b}" , "l");
     	legend2->SetTextSize(0.04);
     	
	      
	      myc1->cd(1);
	      
	      eta_BIB->Scale(1./eta_BIB->Integral());
 	eta_BIB->SetLineWidth(2);
 	eta_BIB->Draw("s hist");
	myc1->Update();
              TPaveStats* stat1 = (TPaveStats*)eta_BIB->GetListOfFunctions()->FindObject("stats");
              stat1->SetX1NDC(0.55);
	      stat1->SetX2NDC(0.85);
	      stat1->SetY1NDC(0.54);
	      stat1->SetY2NDC(0.69);
	      stat1->SetTextColor(kBlue);	
	eta->Scale(1./eta->Integral());
 	eta->SetLineWidth(2); 
 	eta->SetLineColor(kRed);	
	eta->Draw("sames hist");
	myc1->Update();
              TPaveStats* stat2 = (TPaveStats*)eta->GetListOfFunctions()->FindObject("stats");
	      stat2->SetX1NDC(0.55);
              stat2->SetX2NDC(0.85);
              stat2->SetY1NDC(0.7);
              stat2->SetY2NDC(0.85);
	      stat2->SetTextColor(kRed);
	      
	legend2->Draw();
		
	      
	      myc1->cd(2);
	      
	      E_BIB->Scale(1./E_BIB->Integral());
	 	E_BIB->SetLineWidth(2);
	 	E_BIB->Draw("s hist");
		myc1->Update();
              TPaveStats* stat3 = (TPaveStats*)E_BIB->GetListOfFunctions()->FindObject("stats");
              stat3->SetX1NDC(0.55);
	      stat3->SetX2NDC(0.85);
	      stat3->SetY1NDC(0.54);
	      stat3->SetY2NDC(0.69);
	      stat3->SetTextColor(kBlue);	
		E->Scale(1./E->Integral());
	 	E->SetLineWidth(2); 
	 	E->SetLineColor(kRed);	
		E->Draw("sames hist");
		myc1->Update();
              TPaveStats* stat4 = (TPaveStats*)E->GetListOfFunctions()->FindObject("stats");
	      stat4->SetX1NDC(0.55);
              stat4->SetX2NDC(0.85);
              stat4->SetY1NDC(0.7);
              stat4->SetY2NDC(0.85);
	      stat4->SetTextColor(kRed);
	     
		legend2->Draw();
	      
	      
	      myc1->cd(3);
	      pt_BIB->Scale(1./pt_BIB->Integral());
	 	pt_BIB->SetLineWidth(2);
	 	pt_BIB->Draw("s hist");
		myc1->Update();
              TPaveStats* stat5 = (TPaveStats*)pt_BIB->GetListOfFunctions()->FindObject("stats");
              stat5->SetX1NDC(0.55);
	      stat5->SetX2NDC(0.85);
	      stat5->SetY1NDC(0.54);
	      stat5->SetY2NDC(0.69);
	      stat5->SetTextColor(kBlue);	
		pt->Scale(1./pt->Integral());
	 	pt->SetLineWidth(2); 
	 	pt->SetLineColor(kRed);	
		pt->Draw("sames hist");
		myc1->Update();
              TPaveStats* stat6 = (TPaveStats*)pt->GetListOfFunctions()->FindObject("stats");
	      stat6->SetX1NDC(0.55);
              stat6->SetX2NDC(0.85);
              stat6->SetY1NDC(0.7);
              stat6->SetY2NDC(0.85);
	      stat6->SetTextColor(kRed);
	 
	      legend2->Draw();
	      
	      
	      myc1->cd(4);
	      phi_BIB->Scale(1./phi_BIB->Integral());
	 	phi_BIB->SetLineWidth(2);
	 	phi_BIB->Draw("s hist");
		myc1->Update();
              TPaveStats* stat7 = (TPaveStats*)phi_BIB->GetListOfFunctions()->FindObject("stats");
              stat7->SetX1NDC(0.55);
	      stat7->SetX2NDC(0.85);
	      stat7->SetY1NDC(0.54);
	      stat7->SetY2NDC(0.69);
	      stat7->SetTextColor(kBlue);	
		phi->Scale(1./phi->Integral());
	 	phi->SetLineWidth(2); 
	 	phi->SetLineColor(kRed);	
		phi->Draw("sames hist");
		myc1->Update();
              TPaveStats* stat8 = (TPaveStats*)phi->GetListOfFunctions()->FindObject("stats");
	      stat8->SetX1NDC(0.55);
              stat8->SetX2NDC(0.85);
              stat8->SetY1NDC(0.7);
              stat8->SetY2NDC(0.85);
	      stat8->SetTextColor(kRed);

	      legend2->Draw();
	      
	      //myc1->SaveAs("/home/paola/Scrivania/PlotJet.pdf");
	
	}
