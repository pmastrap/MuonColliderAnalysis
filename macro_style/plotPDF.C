#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

void plotPDF(){
	
	double lumi=500.;
	double xsec_sig = 8.9;
	double xsec_irr_bkg = 180.1;
	int N_events=10000;
	//double factor_sig=lumi*xsec_sig/double(N_events);
	//double factor_bkg=lumi*xsec_irr_bkg/double(N_events);
	double factor_sig=1;
	double factor_bkg=1;
	
	
	TFile *f = new TFile("/home/paola/Scrivania/jet_udsg.root");
	
	TH1F * n_jets = (TH1F*)f->Get("n_jets_ev");
	n_jets->SetTitle("Jets per event");
    	n_jets->GetXaxis()->SetTitle("N_{jets}");
    	n_jets->SetLineWidth(2);
	n_jets->Scale(factor_bkg);
	n_jets->Scale(1/n_jets->Integral());
    	
	TH1F * m_inv = (TH1F*)f->Get("m_inv");
	m_inv->SetTitle("Invariant mass of jet pair");
   	m_inv->GetXaxis()->SetTitle("M_{H} [GeV]");
   	m_inv->SetLineWidth(2);
   	m_inv->Rebin(4);
   	m_inv->Scale(factor_bkg);
   	m_inv->Scale(1/m_inv->Integral());

	TH1F* eta = (TH1F*)f->Get("eta_histo");
	eta->SetTitle("Jets pseudorapidity");
    	eta->GetXaxis()->SetTitle("#eta");
    	eta->SetLineWidth(2);
	//eta->GetYaxis()->SetRangeUser(0,0.12);
	eta->Scale(factor_bkg);
	eta->Scale(1/eta->Integral());
    	
	TH1F* E = (TH1F*)f->Get("E_histo");
	E->SetTitle("Jets energy");
    	E->GetXaxis()->SetTitle("Energy [GeV]");
    	E->SetLineWidth(2);
    	E->Rebin(4);
	E->GetXaxis()->SetRangeUser(0,1500.);
  	E->Scale(factor_bkg);
  	E->Scale(1/E->Integral());
    	
	TH1F* pt = (TH1F*)f->Get("pt_histo");
	pt->SetTitle("Jets transverse momentum");
    	pt->GetXaxis()->SetTitle("p_{T} [GeV]");
    	pt->SetLineWidth(2);
    	pt->Rebin(4);
    	pt->GetXaxis()->SetRangeUser(0,1500.);
    	pt->Scale(factor_bkg);
    	pt->Scale(1/pt->Integral());
  
    	TH1F* phi= (TH1F*)f->Get("phi_histo");
	phi->SetTitle("Jets azimuthal angle");
    	phi->GetXaxis()->SetTitle("#phi [rad]");
    	phi->SetLineWidth(2);
    	phi->Scale(factor_bkg);
    	phi->Scale(1/phi->Integral());
    	
    	
    	//HISTO BB
	
	/*TFile *f3 = new TFile("/home/paola/reco_jet_higgs_bb_1000.root");
	TH1F * n_jets_bb = (TH1F*)f3->Get("n_jets_ev");
	n_jets_bb->SetTitle("Jets per event");
    	n_jets_bb->GetXaxis()->SetTitle("N_{jets}");
    	n_jets_bb->SetLineColor(kGreen+3);
    	n_jets_bb->SetLineWidth(2);
	n_jets_bb->Scale(factor/n_jets_bb->Integral());
    	
	TH1F * m_inv_bb = (TH1F*)f3->Get("m_inv");
	m_inv_bb->SetTitle("Invariant mass of jet pair");
   	m_inv_bb->GetXaxis()->SetTitle("M_{H} [GeV]");
   	m_inv_bb->SetLineColor(kGreen+3);
   	m_inv_bb->SetLineWidth(2);
   	m_inv_bb->Scale(factor/m_inv_bb->Integral());

	TH1F* eta_bb = (TH1F*)f3->Get("eta_histo");
	eta_bb->SetTitle("Jets pseudorapidity");
    	eta_bb->GetXaxis()->SetTitle("#eta");
    	eta_bb->SetLineColor(kGreen+3);
    	eta_bb->SetLineWidth(2);
	eta_bb->Scale(factor/eta_bb->Integral());
    	
	TH1F* E_bb = (TH1F*)f3->Get("E_histo");
	E_bb->SetTitle("Jets energy");
    	E_bb->GetXaxis()->SetTitle("Energy [GeV]");
    	E_bb->SetLineColor(kGreen+3);
    	E_bb->SetLineWidth(2);
  	E_bb->Scale(factor/E_bb->Integral());
    	
	TH1F* pt_bb = (TH1F*)f3->Get("pt_histo");
	pt_bb->SetTitle("Jets transverse momentum");
    	pt_bb->GetXaxis()->SetTitle("p_{T} [GeV]");
    	pt_bb->SetLineColor(kGreen+3);
    	pt_bb->SetLineWidth(2);
    	pt_bb->Scale(factor/pt_bb->Integral());
  
    	TH1F* phi_bb= (TH1F*)f3->Get("phi_histo");
	phi_bb->SetTitle("Jets azimuthal angle");
    	phi_bb->GetXaxis()->SetTitle("#phi [rad]");
    	phi_bb->SetLineColor(kGreen+3);
    	phi_bb->SetLineWidth(2);
    	phi_bb->Scale(factor/phi_bb->Integral());*/
    	
    	
  
    	TCanvas *myc1 = new TCanvas("myc1", "canvas", 0,0,1000,1300);
              myc1->Divide(3,2);
              gStyle->SetOptStat("mre");
              myc1->cd(1);
              myc1->SetFillColor(0);
	      
	
	TFile *f2 = new TFile("/home/paola/Scrivania/jet_cc.root");
	
	TH1F * n_jets_BIB = (TH1F*)f2->Get("n_jets_ev");
	n_jets_BIB->SetTitle("Jets per event");
    	n_jets_BIB->GetXaxis()->SetTitle("N_{jets}");
    	n_jets_BIB->SetLineColor(kRed);
    	n_jets_BIB->SetLineWidth(2);
    	n_jets_BIB->Scale(factor_sig);
    	n_jets_BIB->Scale(1/n_jets_BIB->Integral());
	//n_jets_BIB->GetYaxis()->SetRangeUser(0,0.8);
    	
	TH1F * m_inv_BIB = (TH1F*)f2->Get("m_inv");
	m_inv_BIB->SetTitle("Invariant mass of jet pair");
   	m_inv_BIB->GetXaxis()->SetTitle("M_{H} [GeV]");
   	m_inv_BIB->SetLineColor(kRed);
   	m_inv_BIB->SetLineWidth(2);
   	m_inv_BIB->Rebin(4);
	m_inv_BIB->Scale(factor_sig);
	m_inv_BIB->Scale(1/m_inv_BIB->Integral());
	//m_inv_BIB->GetYaxis()->SetRangeUser(0,0.3);
  
   	
	TH1F* eta_BIB = (TH1F*)f2->Get("eta_histo");
	eta_BIB->SetTitle("Jet pseudorapidity");
    	eta_BIB->GetXaxis()->SetTitle("#eta");
    	eta_BIB->SetLineColor(kRed);
    	eta_BIB->SetLineWidth(2);
	//eta_BIB->GetYaxis()->SetRangeUser(0,0.12);
	eta_BIB->Scale(factor_sig);
	eta_BIB->Scale(1/eta_BIB->Integral());
    	
	TH1F* E_BIB = (TH1F*)f2->Get("E_histo");
	E_BIB->SetTitle("Jet energy");
    	E_BIB->GetXaxis()->SetTitle("Energy [GeV]");
    	E_BIB->SetLineColor(kRed);
    	E_BIB->SetLineWidth(2);
    	E_BIB->Rebin(4);
	E_BIB->GetXaxis()->SetRangeUser(0,1500.);
	E_BIB->Scale(factor_sig);
	E_BIB->Scale(1/E_BIB->Integral());
    	
	TH1F* pt_BIB = (TH1F*)f2->Get("pt_histo");
	pt_BIB->SetTitle("Jet transverse momentum");
    	pt_BIB->GetXaxis()->SetTitle("p_{T} [GeV]");
    	pt_BIB->SetLineColor(kRed);
    	pt_BIB->SetLineWidth(2);
    	pt_BIB->Rebin(4);
    	pt_BIB->GetXaxis()->SetRangeUser(0,1500.);
	pt_BIB->Scale(factor_sig);
	pt_BIB->Scale(1/pt_BIB->Integral());
    	
    	TH1F* phi_BIB= (TH1F*)f2->Get("phi_histo");
	phi_BIB->SetTitle("Jet azimuthal angle");
    	phi_BIB->GetXaxis()->SetTitle("#phi [rad]");
    	phi_BIB->SetLineColor(kRed);
    	phi_BIB->SetLineWidth(2);
	phi_BIB->Scale(factor_sig);
	phi_BIB->Scale(1/phi_BIB->Integral());
	//phi_BIB->GetYaxis()->SetRangeUser(0,0.12);
	
    	
    	
    	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
     	//legend2->SetHeader("Distance");
     	legend2->AddEntry(n_jets, "#mu^{+}#mu^{-} #rightarrow #gamma/Z^{0} #rightarrow c#bar{c}", "l");
     	//legend2->AddEntry(n_jets_BIB, "#splitline{H #rightarrow c#bar{c} + BIB(10)}{+correction}" , "l");
     	legend2->AddEntry(n_jets_BIB, "#mu^{+}#mu^{-} #rightarrow H #rightarrow c#bar{c}", "l");
     	legend2->SetTextSize(0.04);
     	
    	
              myc1->cd(1);
              gPad->SetLogy();
              n_jets_BIB->Draw("s hist");
              myc1->Update();
              TPaveStats* stat1 = (TPaveStats*)n_jets_BIB->GetListOfFunctions()->FindObject("stats");
              stat1->SetX1NDC(0.55);
	      stat1->SetX2NDC(0.85);
	      stat1->SetY1NDC(0.54);
	      stat1->SetY2NDC(0.69);
	      stat1->SetTextColor(kRed);
	      
              n_jets->Draw("sames hist");
              myc1->Update();
              TPaveStats* stat2 = (TPaveStats*)n_jets->GetListOfFunctions()->FindObject("stats");
	      stat2->SetX1NDC(0.55);
              stat2->SetX2NDC(0.85);
              stat2->SetY1NDC(0.7);
              stat2->SetY2NDC(0.85);
              stat2->SetTextColor(kBlue);
              
              /*n_jets_bb->Draw("sames hist");
              myc1->Update();
              TPaveStats* stat3 = (TPaveStats*)n_jets_bb->GetListOfFunctions()->FindObject("stats");
	      stat3->SetX1NDC(0.55);
              stat3->SetX2NDC(0.85);
              stat3->SetY1NDC(0.38);
              stat3->SetY2NDC(0.53);
              stat3->SetTextColor(kGreen+3);*/
              
              legend2->Draw();
              
              myc1->cd(2);
              gPad->SetLogx();
	      m_inv_BIB->Draw("s hist");
	      myc1->Update();
              TPaveStats* stat4 = (TPaveStats*)m_inv_BIB->GetListOfFunctions()->FindObject("stats");
              stat4->SetX1NDC(0.57);
	      stat4->SetX2NDC(0.87);
	      stat4->SetY1NDC(0.54);
	      stat4->SetY2NDC(0.69);
	      stat4->SetTextColor(kRed);
	      
	      m_inv->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat5 = (TPaveStats*)m_inv->GetListOfFunctions()->FindObject("stats");
	      stat5->SetX1NDC(0.57);
              stat5->SetX2NDC(0.87);
              stat5->SetY1NDC(0.7);
              stat5->SetY2NDC(0.85);
              stat5->SetTextColor(kBlue);
              
	      /*m_inv_bb->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat6 = (TPaveStats*)m_inv_bb->GetListOfFunctions()->FindObject("stats");
	      stat6->SetX1NDC(0.57);
              stat6->SetX2NDC(0.87);
              stat6->SetY1NDC(0.38);
              stat6->SetY2NDC(0.53);
              stat6->SetTextColor(kGreen+3);*/
              
	      legend2->Draw();
	      
	      myc1->cd(3);
	      gPad->SetLogy();
	      eta_BIB->Draw("s hist");
	      myc1->Update();
              TPaveStats* stat7 = (TPaveStats*)eta_BIB->GetListOfFunctions()->FindObject("stats");
              stat7->SetX1NDC(0.63);
	      stat7->SetX2NDC(0.93);
	      stat7->SetY1NDC(0.54);
	      stat7->SetY2NDC(0.69);
	      stat7->SetTextColor(kRed);
	      
	      eta->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat8 = (TPaveStats*)eta->GetListOfFunctions()->FindObject("stats");
	      stat8->SetX1NDC(0.63);
              stat8->SetX2NDC(0.93);
              stat8->SetY1NDC(0.7);
              stat8->SetY2NDC(0.85);
              stat8->SetTextColor(kBlue);
              
	      /*eta_bb->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat9 = (TPaveStats*)eta_bb->GetListOfFunctions()->FindObject("stats");
	      stat9->SetX1NDC(0.63);
              stat9->SetX2NDC(0.93);
              stat9->SetY1NDC(0.38);
              stat9->SetY2NDC(0.53);
              stat9->SetTextColor(kGreen+3);*/
	      legend2->Draw();
	      
	      myc1->cd(4);
	      gPad->SetLogy();
	      E_BIB->Draw("s hist");
	      myc1->Update();
              TPaveStats* stat10 = (TPaveStats*)E_BIB->GetListOfFunctions()->FindObject("stats");
              stat10->SetX1NDC(0.55);
	      stat10->SetX2NDC(0.85);
	      stat10->SetY1NDC(0.54);
	      stat10->SetY2NDC(0.69);
	      stat10->SetTextColor(kRed);
	      
	      E->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat11 = (TPaveStats*)E->GetListOfFunctions()->FindObject("stats");
	      stat11->SetX1NDC(0.55);
              stat11->SetX2NDC(0.85);
              stat11->SetY1NDC(0.7);
              stat11->SetY2NDC(0.85);
              stat11->SetTextColor(kBlue);
	      /*E_bb->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat12 = (TPaveStats*)E_bb->GetListOfFunctions()->FindObject("stats");
	      stat12->SetX1NDC(0.55);
              stat12->SetX2NDC(0.85);
              stat12->SetY1NDC(0.38);
              stat12->SetY2NDC(0.53);
              stat12->SetTextColor(kGreen+3);*/
	      legend2->Draw();
	      
	      myc1->cd(5);
	      gPad->SetLogy();
	      pt_BIB->Draw("s hist");
	      myc1->Update();
              TPaveStats* stat13 = (TPaveStats*)pt_BIB->GetListOfFunctions()->FindObject("stats");
              stat13->SetX1NDC(0.55);
	      stat13->SetX2NDC(0.85);
	      stat13->SetY1NDC(0.54);
	      stat13->SetY2NDC(0.69);
	      stat13->SetTextColor(kRed);
	      
	      pt->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat14 = (TPaveStats*)pt->GetListOfFunctions()->FindObject("stats");
	      stat14->SetX1NDC(0.55);
              stat14->SetX2NDC(0.85);
              stat14->SetY1NDC(0.7);
              stat14->SetY2NDC(0.85);
              stat14->SetTextColor(kBlue);
              
	      /*pt_bb->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat15 = (TPaveStats*)pt_bb->GetListOfFunctions()->FindObject("stats");
	      stat15->SetX1NDC(0.55);
              stat15->SetX2NDC(0.85);
              stat15->SetY1NDC(0.38);
              stat15->SetY2NDC(0.53);
              stat15->SetTextColor(kGreen+3);*/
	      
	      legend2->Draw();
	      
	      myc1->cd(6);
	      gPad->SetLogy();
	      phi_BIB->Draw("s hist");
	      myc1->Update();
              TPaveStats* stat16 = (TPaveStats*)phi_BIB->GetListOfFunctions()->FindObject("stats");
              stat16->SetX1NDC(0.69);
	      stat16->SetX2NDC(0.99);
	      stat16->SetY1NDC(0.54);
	      stat16->SetY2NDC(0.69);
	      stat16->SetTextColor(kRed);
	      
	      phi->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat17 = (TPaveStats*)phi->GetListOfFunctions()->FindObject("stats");
	      stat17->SetX1NDC(0.69);
              stat17->SetX2NDC(0.99);
              stat17->SetY1NDC(0.7);
              stat17->SetY2NDC(0.85);
              stat17->SetTextColor(kBlue);
              
	      /*phi_bb->Draw("sames hist");
	      myc1->Update();
              TPaveStats* stat18 = (TPaveStats*)phi_bb->GetListOfFunctions()->FindObject("stats");
	      stat18->SetX1NDC(0.69);
              stat18->SetX2NDC(0.99);
              stat18->SetY1NDC(0.38);
              stat18->SetY2NDC(0.53);
              stat18->SetTextColor(kGreen+3);*/
	      
	      legend2->Draw();
	      
	      //myc1->SaveAs("/home/paola/Scrivania/PlotJet.pdf");
	
	}
