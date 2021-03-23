#include "tdrstyle.C"
#include "CMS_lumi_v2.h"



void My_Beautiful_Canvas_TrkVar(TH1F* h1, TH1F* h2, TString Xaxis_title, TString Yaxis_title){

	TString extraText = "Simulation"; 
	TString lumi_1_5TeV = "#sqrt{s}=1.5 TeV, 500 fb^{-1}";
	TString softwareText="v02-06-MC";
	
	setTDRStyle();
	
	int W = 800;
  	int H = 600;
	TCanvas* canv = new TCanvas("cavas","canvas",10,10,W,H);
	int H_ref = 600; 
 	int W_ref = 800; 
  	// references for T, B, L, R
 	float T = 0.08*H_ref;
 	float B = 0.12*H_ref; 
  	float L = 0.12*W_ref;
  	float R = 0.04*W_ref;
  	canv->SetFillColor(0);
  	canv->SetBorderMode(0);
 	canv->SetFrameFillStyle(0);
	canv->SetFrameBorderMode(0);
	canv->SetLeftMargin( L/W );
	canv->SetRightMargin( R/W );
	canv->SetTopMargin( T/H );
	canv->SetBottomMargin( B/H );
	canv->SetTickx(0);
        canv->SetTicky(0);
        
	TH1* fake_hist = new TH1F("fake_hist","fake_hist ",h1->GetNbinsX(),h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
  	//fake_hist ->GetXaxis()->SetNdivisions(6,5,0);
	fake_hist->GetYaxis()->SetTitle(Yaxis_title);
    	fake_hist->GetXaxis()->SetTitle(Xaxis_title); 
  	fake_hist ->GetYaxis()->SetNdivisions(6,5,0);
  	fake_hist ->GetYaxis()->SetTitleOffset(1);

  	fake_hist ->SetMaximum( 1);
  	fake_hist ->Draw();
  	
  	TPad* pad=canv;

	int alignY_=3;
	int alignX_=1;
	int align_= 10*alignX_+ alignY_;
	float HH = pad->GetWh();
	float WW = pad->GetWw();
	float le = pad->GetLeftMargin();
	float t = pad->GetTopMargin();
	float r = pad->GetRightMargin();
	float b = pad->GetBottomMargin();
	float e = 0.025;

	pad->cd();

	TString lumiText=lumi_1_5TeV;

	TLatex latex;
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);

	float extraTextSize = extraOverCmsTextSize*cmsTextSize;

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(lumiTextSize*t);
	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);
	latex.DrawLatex(0.3,1-t+lumiTextOffset*t,softwareText);
	float posX_;
	
	posX_ = le + relPosX*(1-le-r);
	float posY_ = 1-t - relPosY*(1-t-b);

	latex.SetTextFont(cmsTextFont);
	latex.SetTextSize(cmsTextSize*t);
	latex.SetTextAlign(align_);
	latex.DrawLatex(posX_, posY_, cmsText);
	latex.SetTextFont(extraTextFont);
	latex.SetTextAlign(align_);
	latex.SetTextSize(extraTextSize*t);
	latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
	auto legend3 = new TLegend(0.67,1-0.16-t+lumiTextOffset*t,0.86,1-0.02-t+lumiTextOffset*t);
	legend3->AddEntry(h1, "from heavy hadron decay", "l");
	legend3->AddEntry(h2, "from hadronization" , "l");
	legend3->SetTextSize(0.03);
	legend3->SetTextFont(42);
	legend3->SetBorderSize(0);
	legend3->Draw();
	
	h1->Draw("sames hist");
	h2->Draw("sames hist");
	
	canv->SaveAs(Form("/home/paola/Scrivania/macro/macro_style/Plot_belli/Tracks_hh_VS_hadron/%s.png",h1->GetTitle()));
	canv->SaveAs(Form("/home/paola/Scrivania/macro/macro_style/Plot_belli/Tracks_hh_VS_hadron/%s.pdf",h1->GetTitle()));
	}
