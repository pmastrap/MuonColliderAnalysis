void overlap(){
	
	double lumi=500.;
	double xsec_sig = 8.9;
	double xsec_bkg = 180.1;
	int N_events=10000;
	double factor_sig=lumi*xsec_sig/double(N_events);
	double factor_bkg=lumi*xsec_bkg/double(N_events);
	
	//////////////////////////////////////////////////////////////////////////
	///////////////// BACKGROUND /////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	
	TFile *f_bkg= new TFile("/home/paola/Scrivania/InvMass_background.root");
	TH1F * m_inv_bkg = (TH1F*)f_bkg->Get("m_inv");
	m_inv_bkg->SetLineWidth(2);
	m_inv_bkg->SetLineColor(kRed);
	m_inv_bkg->Scale(factor_bkg);
	m_inv_bkg->GetYaxis()->SetTitle("Number of events");
	
	TH1F * m_inv_bc_bkg = (TH1F*)f_bkg->Get("m_inv_BeforeCut");
	m_inv_bc_bkg->SetLineWidth(2);
	m_inv_bc_bkg->SetLineColor(kRed);
	m_inv_bc_bkg->Scale(factor_bkg);
	m_inv_bc_bkg->GetYaxis()->SetTitle("Number of events");
	
	TEfficiency* tag_eff_bkg = (TEfficiency*)f_bkg->Get("pt_before_clone");
	tag_eff_bkg->SetTitle("Mis-tagging efficiency VS pT;pT;efficiency");
	
	/////////////////////////////////////////////////////////////////
	////////// full range ///////////////////////////////////////////
	
	TH1F * pt1_bkg = (TH1F*)f_bkg->Get("pt1_histo");
	pt1_bkg->SetLineWidth(2);
	pt1_bkg->SetLineColor(kRed);
	pt1_bkg->Scale(1/pt1_bkg->Integral());
	pt1_bkg->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * pt2_bkg = (TH1F*)f_bkg->Get("pt2_histo");
	pt2_bkg->SetLineWidth(2);
	pt2_bkg->SetLineColor(kRed);
	pt2_bkg->Scale(1/pt2_bkg->Integral());
	pt2_bkg->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * flav1_bkg = (TH1F*)f_bkg->Get("flav1_histo");
	flav1_bkg->SetLineWidth(2);
	flav1_bkg->SetLineColor(kRed);
	flav1_bkg->Scale(1/flav1_bkg->Integral());
	flav1_bkg->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * flav2_bkg = (TH1F*)f_bkg->Get("flav2_histo");
	flav2_bkg->SetLineWidth(2);
	flav2_bkg->SetLineColor(kRed);
	flav2_bkg->Scale(1/flav2_bkg->Integral());
	flav2_bkg->GetYaxis()->SetTitle("Normalized number of events");
	
	///////////////////////////////////////////////////////////
	////////// region <100 GeV ////////////////////////////////
	
	TH1F * pt1_r_bkg = (TH1F*)f_bkg->Get("pt1_histo_r");
	pt1_r_bkg->SetLineWidth(2);
	pt1_r_bkg->SetLineColor(kRed);
	pt1_r_bkg->Scale(1/pt1_r_bkg->Integral());
	pt1_r_bkg->SetTitle("p_{T} of the first jet - Region of M_{jj}<100 GeV");
	pt1_r_bkg->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * pt2_r_bkg = (TH1F*)f_bkg->Get("pt2_histo_r");
	pt2_r_bkg->SetLineWidth(2);
	pt2_r_bkg->SetLineColor(kRed);
	pt2_r_bkg->Scale(1/pt2_r_bkg->Integral());
	pt2_r_bkg->SetTitle("p_{T} of the second jet - Region of M_{jj}<100 GeV");
	pt2_r_bkg->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * flav1_r_bkg = (TH1F*)f_bkg->Get("flav1_histo_r");
	flav1_r_bkg->SetLineWidth(2);
	flav1_r_bkg->SetLineColor(kRed);
	flav1_r_bkg->Scale(1/flav1_r_bkg->Integral());
	flav1_r_bkg->SetTitle("Flavour of the first jet - Region of M_{jj}<100 GeV");
	flav1_r_bkg->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * flav2_r_bkg = (TH1F*)f_bkg->Get("flav2_histo_r");
	flav2_r_bkg->SetLineWidth(2);
	flav2_r_bkg->SetLineColor(kRed);
	flav2_r_bkg->Scale(1/flav2_r_bkg->Integral());
	flav2_r_bkg->SetTitle("Flavour of the second jet - Region of M_{jj}<100 GeV");
	flav2_r_bkg->GetYaxis()->SetTitle("Normalized number of events");
	
	//////////////////////////////////////////////////////////////////////////
	///////////////// SIGNAL //////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	
	TFile *f_sig= new TFile("/home/paola/Scrivania/InvMass_signal.root");
	TH1F * m_inv_sig = (TH1F*)f_sig->Get("m_inv");
	m_inv_sig->SetLineWidth(2);
	m_inv_sig->SetLineColor(kGreen+2);
	m_inv_sig->Scale(factor_sig);
	
	TH1F * m_inv_bc_sig = (TH1F*)f_sig->Get("m_inv_BeforeCut");
	m_inv_bc_sig->SetLineWidth(2);
	m_inv_bc_sig->SetLineColor(kGreen+2);
	m_inv_bc_sig->Scale(factor_sig);
	
	TEfficiency* tag_eff_sig = (TEfficiency*)f_sig->Get("pt_before_clone");
	tag_eff_sig->SetTitle("c-tagging efficiency VS pT;pT;efficiency");
	
	/////////////////////////////////////////////////////////////////
	////////// full range ///////////////////////////////////////////
	
	TH1F * pt1_sig = (TH1F*)f_sig->Get("pt1_histo");
	pt1_sig->SetLineWidth(2);
	pt1_sig->SetLineColor(kGreen+2);
	pt1_sig->Scale(1/pt1_sig->Integral());
	pt1_sig->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * pt2_sig = (TH1F*)f_sig->Get("pt2_histo");
	pt2_sig->SetLineWidth(2);
	pt2_sig->SetLineColor(kGreen+2);
	pt2_sig->Scale(1/pt2_sig->Integral());
	pt2_sig->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * flav1_sig = (TH1F*)f_sig->Get("flav1_histo");
	flav1_sig->SetLineWidth(2);
	flav1_sig->SetLineColor(kGreen+2);
	flav1_sig->Scale(1/flav1_sig->Integral());
	flav1_sig->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * flav2_sig = (TH1F*)f_sig->Get("flav2_histo");
	flav2_sig->SetLineWidth(2);
	flav2_sig->SetLineColor(kGreen+2);
	flav2_sig->Scale(1/flav2_sig->Integral());
	flav2_sig->GetYaxis()->SetTitle("Normalized number of events");
	
	///////////////////////////////////////////////////////////
	////////// region <100 GeV ////////////////////////////////
	
	TH1F * pt1_r_sig = (TH1F*)f_sig->Get("pt1_histo_r");
	pt1_r_sig->SetLineWidth(2);
	pt1_r_sig->SetLineColor(kGreen+2);
	pt1_r_sig->Scale(1/pt1_r_sig->Integral());
	pt1_r_sig->GetYaxis()->SetTitle("Normalized number of events");
	pt1_r_sig->SetTitle("p_{T} of the first jet - Region of M_{jj}<100 GeV");
	
	TH1F * pt2_r_sig = (TH1F*)f_sig->Get("pt2_histo_r");
	pt2_r_sig->SetLineWidth(2);
	pt2_r_sig->SetLineColor(kGreen+2);
	pt2_r_sig->Scale(1/pt2_r_sig->Integral());
	pt2_r_sig->SetTitle("p_{T} of the second jet - Region of M_{jj}<100 GeV");
	pt2_r_sig->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * flav1_r_sig = (TH1F*)f_sig->Get("flav1_histo_r");
	flav1_r_sig->SetLineWidth(2);
	flav1_r_sig->SetLineColor(kGreen+2);
	flav1_r_sig->Scale(1/flav1_r_sig->Integral());
	flav1_r_sig->SetTitle("Flavour of the first jet - Region of M_{jj}<100 GeV");
	flav1_r_sig->GetYaxis()->SetTitle("Normalized number of events");
	
	TH1F * flav2_r_sig = (TH1F*)f_sig->Get("flav2_histo_r");
	flav2_r_sig->SetLineWidth(2);
	flav2_r_sig->SetLineColor(kGreen+2);
	flav2_r_sig->Scale(1/flav2_r_sig->Integral());
	flav2_r_sig->SetTitle("Flavour of the second jet - Region of M_{jj}<100 GeV");
	flav2_r_sig->GetYaxis()->SetTitle("Normalized number of events");
	
	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
     	legend2->AddEntry(m_inv_bkg, "#mu^{+}#mu^{-} #rightarrow H #nu#bar{#nu} #rightarrow b#bar{b}#nu#bar{#nu}" , "l");
     	legend2->AddEntry(m_inv_sig, "#mu^{+}#mu^{-} #rightarrow H #nu#bar{#nu} #rightarrow c#bar{c}#nu#bar{#nu}", "l");
     	legend2->SetTextSize(0.04);
     	
     	TCanvas *myc=new TCanvas();
	      m_inv_bkg->Draw("s ehist");
	      myc->Update();
              TPaveStats* stat4 = (TPaveStats*)m_inv_bkg->GetListOfFunctions()->FindObject("stats");
              stat4->SetX1NDC(0.57);
	      stat4->SetX2NDC(0.87);
	      stat4->SetY1NDC(0.54);
	      stat4->SetY2NDC(0.69);
	      stat4->SetTextColor(kRed);
	      
	      m_inv_sig->Draw("sames ehist");
	      myc->Update();
              TPaveStats* stat5 = (TPaveStats*)m_inv_sig->GetListOfFunctions()->FindObject("stats");
	      stat5->SetX1NDC(0.57);
              stat5->SetX2NDC(0.87);
              stat5->SetY1NDC(0.7);
              stat5->SetY2NDC(0.85);
              stat5->SetTextColor(kGreen+2);
             
             legend2->Draw("");
             
         TCanvas *myc1=new TCanvas();
	      m_inv_bc_bkg->Draw("s ehist");
	      myc1->Update();
              TPaveStats* stat1 = (TPaveStats*)m_inv_bc_bkg->GetListOfFunctions()->FindObject("stats");
              stat1->SetX1NDC(0.57);
	      stat1->SetX2NDC(0.87);
	      stat1->SetY1NDC(0.54);
	      stat1->SetY2NDC(0.69);
	      stat1->SetTextColor(kRed);
	      
	      m_inv_bc_sig->Draw("sames ehist");
	      myc1->Update();
              TPaveStats* stat2 = (TPaveStats*)m_inv_bc_sig->GetListOfFunctions()->FindObject("stats");
	      stat2->SetX1NDC(0.57);
              stat2->SetX2NDC(0.87);
              stat2->SetY1NDC(0.7);
              stat2->SetY2NDC(0.85);
              stat2->SetTextColor(kGreen+2);
          legend2->Draw("");    
             
          TCanvas *myc2=new TCanvas();
          tag_eff_sig->Draw("ap");
          TCanvas *myc3=new TCanvas();
          tag_eff_bkg->Draw("ap");
          
          
          TCanvas *myc4=new TCanvas();
	      pt1_bkg->Draw("s ehist");
	      myc4->Update();
              TPaveStats* stat6 = (TPaveStats*)pt1_bkg->GetListOfFunctions()->FindObject("stats");
              stat6->SetX1NDC(0.57);
	      stat6->SetX2NDC(0.87);
	      stat6->SetY1NDC(0.54);
	      stat6->SetY2NDC(0.69);
	      stat6->SetTextColor(kRed);
	      
	      pt1_sig->Draw("sames ehist");
	      myc4->Update();
              TPaveStats* stat7 = (TPaveStats*)pt1_sig->GetListOfFunctions()->FindObject("stats");
	      stat7->SetX1NDC(0.57);
              stat7->SetX2NDC(0.87);
              stat7->SetY1NDC(0.7);
              stat7->SetY2NDC(0.85);
              stat7->SetTextColor(kGreen+2);
          legend2->Draw("");   
          
          TCanvas *myc5=new TCanvas();
	      pt2_bkg->Draw("s ehist");
	      myc5->Update();
              TPaveStats* stat8 = (TPaveStats*)pt2_bkg->GetListOfFunctions()->FindObject("stats");
              stat8->SetX1NDC(0.57);
	      stat8->SetX2NDC(0.87);
	      stat8->SetY1NDC(0.54);
	      stat8->SetY2NDC(0.69);
	      stat8->SetTextColor(kRed);
	      
	      pt2_sig->Draw("sames ehist");
	      myc5->Update();
              TPaveStats* stat9 = (TPaveStats*)pt2_sig->GetListOfFunctions()->FindObject("stats");
	      stat9->SetX1NDC(0.57);
              stat9->SetX2NDC(0.87);
              stat9->SetY1NDC(0.7);
              stat9->SetY2NDC(0.85);
              stat9->SetTextColor(kGreen+2);
          legend2->Draw("");   
          
          TCanvas *myc6=new TCanvas();
              flav1_bkg->GetYaxis()->SetRangeUser(0.,1);
	      flav1_bkg->Draw("s hist");
	      myc6->Update();
              TPaveStats* stat10 = (TPaveStats*)flav1_bkg->GetListOfFunctions()->FindObject("stats");
              stat10->SetX1NDC(0.57);
	      stat10->SetX2NDC(0.87);
	      stat10->SetY1NDC(0.54);
	      stat10->SetY2NDC(0.69);
	      stat10->SetTextColor(kRed);
	      
	      flav1_sig->Draw("sames hist");
	      myc6->Update();
              TPaveStats* stat11 = (TPaveStats*)flav1_sig->GetListOfFunctions()->FindObject("stats");
	      stat11->SetX1NDC(0.57);
              stat11->SetX2NDC(0.87);
              stat11->SetY1NDC(0.7);
              stat11->SetY2NDC(0.85);
              stat11->SetTextColor(kGreen+2);
          legend2->Draw(""); 
          
          TCanvas *myc7=new TCanvas();
              flav2_bkg->GetYaxis()->SetRangeUser(0.,1);
	      flav2_bkg->Draw("s hist");
	      myc7->Update();
              TPaveStats* stat12 = (TPaveStats*)flav2_bkg->GetListOfFunctions()->FindObject("stats");
              stat12->SetX1NDC(0.57);
	      stat12->SetX2NDC(0.87);
	      stat12->SetY1NDC(0.54);
	      stat12->SetY2NDC(0.69);
	      stat12->SetTextColor(kRed);
	      
	      flav2_sig->Draw("sames hist");
	      myc7->Update();
              TPaveStats* stat13 = (TPaveStats*)flav2_sig->GetListOfFunctions()->FindObject("stats");
	      stat13->SetX1NDC(0.57);
              stat13->SetX2NDC(0.87);
              stat13->SetY1NDC(0.7);
              stat13->SetY2NDC(0.85);
              stat13->SetTextColor(kGreen+2);
          legend2->Draw("");    
          
          ///restricted region
          TCanvas *myc4_r=new TCanvas();
	      pt1_r_bkg->Draw("s ehist");
	      myc4_r->Update();
              TPaveStats* stat6_r = (TPaveStats*)pt1_r_bkg->GetListOfFunctions()->FindObject("stats");
              stat6_r->SetX1NDC(0.57);
	      stat6_r->SetX2NDC(0.87);
	      stat6_r->SetY1NDC(0.54);
	      stat6_r->SetY2NDC(0.69);
	      stat6_r->SetTextColor(kRed);
	      
	      pt1_r_sig->Draw("sames ehist");
	      myc4_r->Update();
              TPaveStats* stat7_r = (TPaveStats*)pt1_r_sig->GetListOfFunctions()->FindObject("stats");
	      stat7_r->SetX1NDC(0.57);
              stat7_r->SetX2NDC(0.87);
              stat7_r->SetY1NDC(0.7);
              stat7_r->SetY2NDC(0.85);
              stat7_r->SetTextColor(kGreen+2);
          legend2->Draw("");   
          
          TCanvas *myc5_r=new TCanvas();
	      pt2_r_bkg->Draw("s ehist");
	      myc5_r->Update();
              TPaveStats* stat8_r = (TPaveStats*)pt2_r_bkg->GetListOfFunctions()->FindObject("stats");
              stat8_r->SetX1NDC(0.57);
	      stat8_r->SetX2NDC(0.87);
	      stat8_r->SetY1NDC(0.54);
	      stat8_r->SetY2NDC(0.69);
	      stat8_r->SetTextColor(kRed);
	      
	      pt2_r_sig->Draw("sames ehist");
	      myc5_r->Update();
              TPaveStats* stat9_r = (TPaveStats*)pt2_r_sig->GetListOfFunctions()->FindObject("stats");
	      stat9_r->SetX1NDC(0.57);
              stat9_r->SetX2NDC(0.87);
              stat9_r->SetY1NDC(0.7);
              stat9_r->SetY2NDC(0.85);
              stat9_r->SetTextColor(kGreen+2);
          legend2->Draw("");   
          
          TCanvas *myc6_r=new TCanvas();
              flav1_r_bkg->GetYaxis()->SetRangeUser(0.,1);
	      flav1_r_bkg->Draw("s hist");
	      myc6_r->Update();
              TPaveStats* stat10_r = (TPaveStats*)flav1_r_bkg->GetListOfFunctions()->FindObject("stats");
              stat10_r->SetX1NDC(0.57);
	      stat10_r->SetX2NDC(0.87);
	      stat10_r->SetY1NDC(0.54);
	      stat10_r->SetY2NDC(0.69);
	      stat10_r->SetTextColor(kRed);
	      
	      flav1_r_sig->Draw("sames hist");
	      myc6_r->Update();
              TPaveStats* stat11_r = (TPaveStats*)flav1_r_sig->GetListOfFunctions()->FindObject("stats");
	      stat11_r->SetX1NDC(0.57);
              stat11_r->SetX2NDC(0.87);
              stat11_r->SetY1NDC(0.7);
              stat11_r->SetY2NDC(0.85);
              stat11_r->SetTextColor(kGreen+2);
          legend2->Draw(""); 
          
          TCanvas *myc7_r=new TCanvas();
              flav2_r_bkg->GetYaxis()->SetRangeUser(0.,1);
	      flav2_r_bkg->Draw("s hist");
	      myc7_r->Update();
              TPaveStats* stat12_r = (TPaveStats*)flav2_r_bkg->GetListOfFunctions()->FindObject("stats");
              stat12_r->SetX1NDC(0.57);
	      stat12_r->SetX2NDC(0.87);
	      stat12_r->SetY1NDC(0.54);
	      stat12_r->SetY2NDC(0.69);
	      stat12_r->SetTextColor(kRed);
	      
	      flav2_r_sig->Draw("sames hist");
	      myc7_r->Update();
              TPaveStats* stat13_r = (TPaveStats*)flav2_r_sig->GetListOfFunctions()->FindObject("stats");
	      stat13_r->SetX1NDC(0.57);
              stat13_r->SetX2NDC(0.87);
              stat13_r->SetY1NDC(0.7);
              stat13_r->SetY2NDC(0.85);
              stat13_r->SetTextColor(kGreen+2);
          legend2->Draw("");    
       
	}
