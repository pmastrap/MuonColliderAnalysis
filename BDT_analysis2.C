
void BDT_analysis2(){
	
	gStyle->SetOptStat("0");
	
	Int_t           classID;
   	Float_t         SIP_sig_xy_TRK1;
   	Float_t         SIP_sig_xy_TRK2;
   	Float_t         SIP_sig_xyz_TRK1;
   	Float_t         SIP_sig_xyz_TRK2;
   	Float_t         SIP_sig_xy_TRK_above_trh;
   	Float_t         SIP_sig_xyz_TRK_above_trh;
  	Float_t         dR_jet_TRK1;
   	Float_t         lepton_cat;
   	Float_t         pt_lj_LEP1;
   	Float_t         Et_trk_Jet;
   	Float_t         BDT;
   	
   	//BDT<0
   	TH1F* SIP_sig_xy_TRK1_h=new TH1F("SIP_sig_xy_TRK1_h","SIP_sig_xy_TRK1_h",100,-10,100);
   	SIP_sig_xy_TRK1_h->SetLineColor(kRed);
	
	TH1F* SIP_sig_xy_TRK2_h=new TH1F("SIP_sig_xy_TRK2_h","SIP_sig_xy_TRK2_h",100,-10,100);
	SIP_sig_xy_TRK2_h->SetLineColor(kRed);

	TH1F* SIP_sig_xyz_TRK1_h=new TH1F("SIP_sig_xyz_TRK1_h","SIP_sig_xyz_TRK1_h",100,-10,100);
	SIP_sig_xyz_TRK1_h->SetLineColor(kRed);
	
	TH1F* SIP_sig_xyz_TRK2_h=new TH1F("SIP_sig_xyz_TRK2_h","SIP_sig_xyz_TRK2_h",100,-10,100);
	SIP_sig_xyz_TRK2_h->SetLineColor(kRed);
	
	TH1F* SIP_sig_xy_TRK_above_trh_h=new TH1F("SIP_sig_xy_TRK_above_trh_h","SIP_sig_xy_TRK_above_trh_h",200,-100,100);
	SIP_sig_xy_TRK_above_trh_h->SetLineColor(kRed);
	
	TH1F* SIP_sig_xyz_TRK_above_trh_h=new TH1F("SIP_sig_xyz_TRK_above_trh_h","SIP_sig_xyz_TRK_above_trh_h",200,-100,100);
	SIP_sig_xyz_TRK_above_trh_h->SetLineColor(kRed);
	
	TH1F* dR_jet_TRK1_h=new TH1F("dR_jet_TRK1_h","dR_jet_TRK1_h",70,0.,1.0);
	dR_jet_TRK1_h->SetLineColor(kRed);
	
	TH1F* lepton_cat_h =new TH1F("lepton_cat_h","lepton_cat_h",3,-0.5,2.5);
	lepton_cat_h->SetLineColor(kRed);
	
	TH1F* pt_lj_LEP1_h =new TH1F("pt_lj_LEP1_h","pt_lj_LEP1_h",50,0.,1.);
	pt_lj_LEP1_h->SetLineColor(kRed);
	
	TH1F* Et_trk_Jet_h=new TH1F("Et_trk_Jet_h","Et_trk_Jet_h",50,0.,1.);
	Et_trk_Jet_h->SetLineColor(kRed);
	
	//BDT>0.1
	TH1F* SIP_sig_xy_TRK1_h2=new TH1F("SIP_sig_xy_TRK1_h2","SIP_sig_xy_TRK1_h2",100,-10,100);
	
	TH1F* SIP_sig_xy_TRK2_h2=new TH1F("SIP_sig_xy_TRK2_h2","SIP_sig_xy_TRK2_h2",100,-10,100);
	
	TH1F* SIP_sig_xyz_TRK1_h2=new TH1F("SIP_sig_xyz_TRK1_h2","SIP_sig_xyz_TRK1_h2",100,-10,100);
	
	TH1F* SIP_sig_xyz_TRK2_h2=new TH1F("SIP_sig_xyz_TRK2_h2","SIP_sig_xyz_TRK2_h2",100,-10,100);
	
	TH1F* SIP_sig_xy_TRK_above_trh_h2=new TH1F("SIP_sig_xy_TRK_above_trh_h2","SIP_sig_xy_TRK_above_trh_h2",200,-100,100);
	
	TH1F* SIP_sig_xyz_TRK_above_trh_h2=new TH1F("SIP_sig_xyz_TRK_above_trh_h2","SIP_sig_xyz_TRK_above_trh_h2",200,-100,100);
	
	TH1F* dR_jet_TRK1_h2=new TH1F("dR_jet_TRK1_h2","dR_jet_TRK1_h2",70,0.,1.0);
	
	TH1F* lepton_cat_h2 =new TH1F("lepton_cat_h2","lepton_cat_h2",3,-0.5,2.5);
	
	TH1F* pt_lj_LEP1_h2 =new TH1F("pt_lj_LEP1_h2","pt_lj_LEP1_h2",50,0.,1.);
	
	TH1F* Et_trk_Jet_h2=new TH1F("Et_trk_Jet_h2","Et_trk_Jet_h2",50,0.,1.);
	
	//BDT<0 BKG
   	TH1F* SIP_sig_xy_TRK1_h3=new TH1F("SIP_sig_xy_TRK1_h3","SIP_sig_xy_TRK1_h3",100,-10,100);
   	SIP_sig_xy_TRK1_h3->SetLineColor(kGreen+2);
	
	TH1F* SIP_sig_xy_TRK2_h3=new TH1F("SIP_sig_xy_TRK2_h3","SIP_sig_xy_TRK2_h3",100,-10,100);
	SIP_sig_xy_TRK2_h3->SetLineColor(kGreen+2);

	TH1F* SIP_sig_xyz_TRK1_h3=new TH1F("SIP_sig_xyz_TRK1_h3","SIP_sig_xyz_TRK1_h3",100,-10,100);
	SIP_sig_xyz_TRK1_h3->SetLineColor(kGreen+2);
	
	TH1F* SIP_sig_xyz_TRK2_h3=new TH1F("SIP_sig_xyz_TRK2_h3","SIP_sig_xyz_TRK2_h3",100,-10,100);
	SIP_sig_xyz_TRK2_h3->SetLineColor(kGreen+2);
	
	TH1F* SIP_sig_xy_TRK_above_trh_h3=new TH1F("SIP_sig_xy_TRK_above_trh_h3","SIP_sig_xy_TRK_above_trh_h3",200,-100,100);
	SIP_sig_xy_TRK_above_trh_h3->SetLineColor(kGreen+2);
	
	TH1F* SIP_sig_xyz_TRK_above_trh_h3=new TH1F("SIP_sig_xyz_TRK_above_trh_h3","SIP_sig_xyz_TRK_above_trh_h3",200,-100,100);
	SIP_sig_xyz_TRK_above_trh_h3->SetLineColor(kGreen+2);
	
	TH1F* dR_jet_TRK1_h3=new TH1F("dR_jet_TRK1_h3","dR_jet_TRK1_h3",70,0.,1.0);
	dR_jet_TRK1_h3->SetLineColor(kGreen+2);
	
	TH1F* lepton_cat_h3 =new TH1F("lepton_cat_h3","lepton_cat_h3",3,-0.5,2.5);
	lepton_cat_h3->SetLineColor(kGreen+2);
	
	TH1F* pt_lj_LEP1_h3 =new TH1F("pt_lj_LEP1_h3","pt_lj_LEP1_h3",50,0.,1.);
	pt_lj_LEP1_h3->SetLineColor(kGreen+2);
	
	TH1F* Et_trk_Jet_h3=new TH1F("Et_trk_Jet_h3","Et_trk_Jet_h3",50,0.,1.);
	Et_trk_Jet_h3->SetLineColor(kGreen+2);
	
	//check on test tree
	TChain* fChain = new TChain("fChain");
   	fChain->Add("/home/paola/Scrivania/macro/Trees_tagger/BDT/TMVA_output_10k_3_light.root/dataset_10k_3_light/TestTree");
   	
   	fChain->SetBranchAddress("classID", &classID);
   	fChain->SetBranchAddress("BDT", &BDT);
   	
  	fChain->SetBranchAddress("SIP_sig_xy_TRK1", &SIP_sig_xy_TRK1);
   	fChain->SetBranchAddress("SIP_sig_xy_TRK2", &SIP_sig_xy_TRK2);
   	fChain->SetBranchAddress("SIP_sig_xyz_TRK1", &SIP_sig_xyz_TRK1);
   	fChain->SetBranchAddress("SIP_sig_xyz_TRK2", &SIP_sig_xyz_TRK2);
   	fChain->SetBranchAddress("SIP_sig_xy_TRK_above_trh", &SIP_sig_xy_TRK_above_trh);
   	fChain->SetBranchAddress("SIP_sig_xyz_TRK_above_trh", &SIP_sig_xyz_TRK_above_trh);
   	fChain->SetBranchAddress("dR_jet_TRK1", &dR_jet_TRK1);
   	fChain->SetBranchAddress("lepton_cat", &lepton_cat);
   	fChain->SetBranchAddress("pt_lj_LEP1", &pt_lj_LEP1);
   	fChain->SetBranchAddress("Et_trk_Jet", &Et_trk_Jet);
   	
   	long int ientry=0;
   	
   	for(ientry=0; ientry< fChain->GetEntries();++ientry){
   		
   		fChain->GetEntry(ientry);	
   		
   		if(classID==0){
   			if(BDT<0.){
   				SIP_sig_xy_TRK1_h->Fill(SIP_sig_xy_TRK1);
   				SIP_sig_xy_TRK2_h->Fill(SIP_sig_xy_TRK2);
			   	SIP_sig_xyz_TRK1_h->Fill(SIP_sig_xyz_TRK1);
			   	SIP_sig_xyz_TRK2_h->Fill(SIP_sig_xyz_TRK2);
			   	SIP_sig_xy_TRK_above_trh_h->Fill(SIP_sig_xy_TRK_above_trh);
			   	SIP_sig_xyz_TRK_above_trh_h->Fill(SIP_sig_xyz_TRK_above_trh);
			  	dR_jet_TRK1_h->Fill(dR_jet_TRK1);
			   	lepton_cat_h->Fill(lepton_cat);
			   	pt_lj_LEP1_h->Fill(pt_lj_LEP1);
			   	Et_trk_Jet_h->Fill(Et_trk_Jet);
			   	}
			else if(BDT>0.1){
				SIP_sig_xy_TRK1_h2->Fill(SIP_sig_xy_TRK1);
   				SIP_sig_xy_TRK2_h2->Fill(SIP_sig_xy_TRK2);
			   	SIP_sig_xyz_TRK1_h2->Fill(SIP_sig_xyz_TRK1);
			   	SIP_sig_xyz_TRK2_h2->Fill(SIP_sig_xyz_TRK2);
			   	SIP_sig_xy_TRK_above_trh_h2->Fill(SIP_sig_xy_TRK_above_trh);
			   	SIP_sig_xyz_TRK_above_trh_h2->Fill(SIP_sig_xyz_TRK_above_trh);
			  	dR_jet_TRK1_h2->Fill(dR_jet_TRK1);
			   	lepton_cat_h2->Fill(lepton_cat);
			   	pt_lj_LEP1_h2->Fill(pt_lj_LEP1);
			   	Et_trk_Jet_h2->Fill(Et_trk_Jet);
				}
			}
		else if(classID==1){
   			if(BDT<0.){
   				SIP_sig_xy_TRK1_h3->Fill(SIP_sig_xy_TRK1);
   				SIP_sig_xy_TRK2_h3->Fill(SIP_sig_xy_TRK2);
			   	SIP_sig_xyz_TRK1_h3->Fill(SIP_sig_xyz_TRK1);
			   	SIP_sig_xyz_TRK2_h3->Fill(SIP_sig_xyz_TRK2);
			   	SIP_sig_xy_TRK_above_trh_h3->Fill(SIP_sig_xy_TRK_above_trh);
			   	SIP_sig_xyz_TRK_above_trh_h3->Fill(SIP_sig_xyz_TRK_above_trh);
			  	dR_jet_TRK1_h3->Fill(dR_jet_TRK1);
			   	lepton_cat_h3->Fill(lepton_cat);
			   	pt_lj_LEP1_h3->Fill(pt_lj_LEP1);
			   	Et_trk_Jet_h3->Fill(Et_trk_Jet);
   				}
   		
   		
   			}
   		}
   		
   	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
	legend2->AddEntry(SIP_sig_xy_TRK1_h, "SIGNAL BDT < 0", "l");
	legend2->AddEntry(SIP_sig_xy_TRK1_h2, "SIGNAL BDT > 0.1" , "l");
	legend2->AddEntry(SIP_sig_xy_TRK1_h3, "BKG BDT < 0" , "l");
	legend2->SetTextSize(0.04);
 
 	TCanvas* c=new TCanvas();
 	SIP_sig_xy_TRK1_h->Scale(1/SIP_sig_xy_TRK1_h->Integral());
 	SIP_sig_xy_TRK1_h->Draw("hist");
 	SIP_sig_xy_TRK1_h2->Scale(1/SIP_sig_xy_TRK1_h2->Integral());
 	SIP_sig_xy_TRK1_h2->Draw("hist same");
 	SIP_sig_xy_TRK1_h3->Scale(1/SIP_sig_xy_TRK1_h3->Integral());
 	SIP_sig_xy_TRK1_h3->Draw("hist same");
 	legend2->Draw();
 	
 	TCanvas* c1=new TCanvas();
 	SIP_sig_xy_TRK2_h->Scale(1/SIP_sig_xy_TRK2_h->Integral());
   	SIP_sig_xy_TRK2_h->Draw("hist");
   	SIP_sig_xy_TRK2_h2->Scale(1/SIP_sig_xy_TRK2_h2->Integral());
   	SIP_sig_xy_TRK2_h2->Draw("hist same");
   	SIP_sig_xy_TRK2_h3->Scale(1/SIP_sig_xy_TRK2_h3->Integral());
   	SIP_sig_xy_TRK2_h3->Draw(" hist same");
   	legend2->Draw();
   	
   	TCanvas* c2=new TCanvas();
   	SIP_sig_xyz_TRK1_h->Scale(1/SIP_sig_xyz_TRK1_h->Integral());
	SIP_sig_xyz_TRK1_h->Draw("hist");
	SIP_sig_xyz_TRK1_h2->Scale(1/SIP_sig_xyz_TRK1_h2->Integral());
	SIP_sig_xyz_TRK1_h2->Draw("hist same");
	SIP_sig_xyz_TRK1_h3->Scale(1/SIP_sig_xyz_TRK1_h3->Integral());
	SIP_sig_xyz_TRK1_h3->Draw("hist same");
	legend2->Draw();
	
	TCanvas* c3=new TCanvas();
	SIP_sig_xyz_TRK2_h->Scale(1/SIP_sig_xyz_TRK2_h->Integral());
	SIP_sig_xyz_TRK2_h->Draw("hist");
	SIP_sig_xyz_TRK2_h2->Scale(1/SIP_sig_xyz_TRK2_h2->Integral());
	SIP_sig_xyz_TRK2_h2->Draw("hist same");
	SIP_sig_xyz_TRK2_h3->Scale(1/SIP_sig_xyz_TRK2_h3->Integral());
	SIP_sig_xyz_TRK2_h3->Draw("hist same");
	legend2->Draw();
	
	TCanvas* c4=new TCanvas();
	SIP_sig_xy_TRK_above_trh_h->Scale(1/SIP_sig_xy_TRK_above_trh_h->Integral());
	SIP_sig_xy_TRK_above_trh_h->Draw("hist");
	SIP_sig_xy_TRK_above_trh_h2->Scale(1/SIP_sig_xy_TRK_above_trh_h2->Integral());
	SIP_sig_xy_TRK_above_trh_h2->Draw("hist same");
	SIP_sig_xy_TRK_above_trh_h3->Scale(1/SIP_sig_xy_TRK_above_trh_h3->Integral());
	SIP_sig_xy_TRK_above_trh_h3->Draw("hist same");
	legend2->Draw();
	
	TCanvas* c5=new TCanvas();
	SIP_sig_xyz_TRK_above_trh_h->Scale(1/SIP_sig_xyz_TRK_above_trh_h->Integral());
	SIP_sig_xyz_TRK_above_trh_h->Draw("hist");
	SIP_sig_xyz_TRK_above_trh_h2->Scale(1/SIP_sig_xyz_TRK_above_trh_h2->Integral());
	SIP_sig_xyz_TRK_above_trh_h2->Draw("hist same");
	SIP_sig_xyz_TRK_above_trh_h3->Scale(1/SIP_sig_xyz_TRK_above_trh_h3->Integral());
	SIP_sig_xyz_TRK_above_trh_h3->Draw("hist same");
	legend2->Draw();
	
	TCanvas* c6=new TCanvas();
	dR_jet_TRK1_h->Scale(1/dR_jet_TRK1_h->Integral());
	dR_jet_TRK1_h->Draw("hist");
	dR_jet_TRK1_h2->Scale(1/dR_jet_TRK1_h2->Integral());	
	dR_jet_TRK1_h2->Draw("hist same");
	dR_jet_TRK1_h3->Scale(1/dR_jet_TRK1_h3->Integral());	
	dR_jet_TRK1_h3->Draw("hist same");
	legend2->Draw();
	
	TCanvas* c7=new TCanvas();
	lepton_cat_h->Scale(1/lepton_cat_h->Integral());
	lepton_cat_h->Draw("hist");
	lepton_cat_h2->Scale(1/lepton_cat_h2->Integral());
	lepton_cat_h2->Draw("hist same");
	lepton_cat_h3->Scale(1/lepton_cat_h3->Integral());
	lepton_cat_h3->Draw("hist same");
	legend2->Draw();
	
	TCanvas* c8=new TCanvas();
	pt_lj_LEP1_h->Scale(1/pt_lj_LEP1_h->Integral());
	pt_lj_LEP1_h->Draw("hist");
	pt_lj_LEP1_h2->Scale(1/pt_lj_LEP1_h2->Integral());
	pt_lj_LEP1_h2->Draw("hist same");
	pt_lj_LEP1_h3->Scale(1/pt_lj_LEP1_h3->Integral());
	pt_lj_LEP1_h3->Draw("hist same");
	legend2->Draw();
	
	TCanvas* c9=new TCanvas();
	Et_trk_Jet_h->Scale(1/Et_trk_Jet_h->Integral());
	Et_trk_Jet_h->Draw("hist");
	Et_trk_Jet_h2->Scale(1/Et_trk_Jet_h2->Integral());
	Et_trk_Jet_h2->Draw("hist same");
	Et_trk_Jet_h3->Scale(1/Et_trk_Jet_h3->Integral());
	Et_trk_Jet_h3->Draw("hist same");
	legend2->Draw();
 
 
 
 
   	
   	}
