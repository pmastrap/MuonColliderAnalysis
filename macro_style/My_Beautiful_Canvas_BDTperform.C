#include "tdrstyle.C"
#include "CMS_lumi_v2.h"



void My_Beautiful_Canvas_BDTperform(TCanvas* canv){

	TString extraText = "Simulation"; 
	TString lumi_1_5TeV = "#sqrt{s}=1.5 TeV, 500 fb^{-1}";
	TString softwareText="v02-06-MC";
	
	setTDRStyle();
	int W = 800;
  	int H = 600;
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

	}
