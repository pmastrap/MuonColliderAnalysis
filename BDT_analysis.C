// analyze BDT results 
#define NUM_EVT_S 2*8.9*500
#define NUM_EVT_B 2*180.1*500

/*#define NUM_EVT_S 1000
#define NUM_EVT_B 1000*/


void BDT_analysis(){
	//Chain* fChain = new TChain("fChain");
   	//fChain->Add("/home/paola/Scrivania/macro/Trees_tagger/BDT/TMVA_output_10k.root/dataset_10k/");
	
	TFile *f = new TFile("/home/paola/Scrivania/macro/Trees_tagger/BDT/TMVA_output_10k_3.root");
	
	TH1F * bdt_sig = (TH1F*)f->Get("/dataset_10k_3/Method_BDT/BDT/MVA_BDT_S");
	TH1F * bdt_train_sig = (TH1F*)f->Get("/dataset_10k_3/Method_BDT/BDT/MVA_BDT_Train_S");
	
	TH1F * bdt_bkg = (TH1F*)f->Get("/dataset_10k_3/Method_BDT/BDT/MVA_BDT_B");
	TH1F * bdt_train_bkg = (TH1F*)f->Get("/dataset_10k_3/Method_BDT/BDT/MVA_BDT_Train_B");
	
	TH1F * bdt_sig_eff = (TH1F*)f->Get("/dataset_10k_3/Method_BDT/BDT/MVA_BDT_effS");
	TH1F * bdt_bkg_eff = (TH1F*)f->Get("/dataset_10k_3/Method_BDT/BDT/MVA_BDT_effB");
	
	
	cout<<"N bin: "<<bdt_sig_eff->GetNbinsX()<<endl;
	cout<<"Low edge: "<<bdt_sig_eff->GetXaxis()->GetBinLowEdge(1)<<endl;
	cout<<"Up edge: "<<bdt_sig_eff->GetXaxis()->GetBinUpEdge(bdt_sig_eff->GetNbinsX())<<endl;
	
	cout<<"N bin: "<<bdt_bkg_eff->GetNbinsX()<<endl;
	cout<<"Low edge: "<<bdt_bkg_eff->GetXaxis()->GetBinLowEdge(1)<<endl;
	cout<<"Up edge: "<<bdt_bkg_eff->GetXaxis()->GetBinUpEdge(bdt_bkg_eff->GetNbinsX())<<endl;
	
	double sradb[20000];
	double bdt_output[20000];
	double max=0.;
	int N_bins=0;
	long int i=1;
	
	N_bins=bdt_sig_eff->GetNbinsX();
	
	for(i=1;i<N_bins+1;i++){
		bdt_output[i-1]=bdt_sig_eff->GetXaxis()->GetBinCenter(i);
		sradb[i-1]=bdt_sig_eff->GetBinContent(i)*NUM_EVT_S/(sqrt(bdt_sig_eff->GetBinContent(i)*NUM_EVT_S + bdt_bkg_eff->GetBinContent(i)*NUM_EVT_B));
		}
	
	TGraph * s_rad_b= new TGraph(N_bins,bdt_output,sradb);
	max=TMath::MaxElement(s_rad_b->GetN(),s_rad_b->GetY());
	cout<<"Max of S/rad(B): "<<max<<endl;
	
	for(i=0;i<N_bins;i++){
		if(sradb[i]==max){ break;}
		}
	cout<<"BDT output: "<<bdt_output[i]<<endl;
	s_rad_b->SetTitle("Cut efficiency");
	s_rad_b->GetXaxis()->SetTitle("Cut value applied on BDT output");
	s_rad_b->GetYaxis()->SetTitle("S/rad(S+B)");
	s_rad_b->SetLineColor(kBlue);
	s_rad_b->SetMarkerColor(kBlue);
	s_rad_b->SetMarkerStyle(20);
	s_rad_b->SetMarkerSize(0.5);
	s_rad_b->Draw("AP");
	
	
	double KS_sig =0.;
	double KS_bkg =0.;
	double KS_maxdist_sig=0.0;
	double KS_maxdist_bkg=0.0;
	
	KS_sig = bdt_sig->KolmogorovTest(bdt_train_sig, "D");
	KS_bkg = bdt_bkg->KolmogorovTest(bdt_train_bkg, "D");
	
	KS_maxdist_sig = bdt_sig->KolmogorovTest(bdt_train_sig, "MD");
	KS_maxdist_bkg = bdt_bkg->KolmogorovTest(bdt_train_bkg, "MD");
	
	cout<<"KS probability for signal: "<<KS_sig<<endl;
	cout<<"Maximum Kolmogorov distance for signal: "<<KS_maxdist_sig<<endl;
	cout<<"KS probability for bkg: "<<KS_bkg<<endl;
	cout<<"Maximum Kolmogorov distance for bkg: "<<KS_maxdist_bkg<<endl;

	TLatex latex;
  	latex.SetTextSize(0.03);
   	//latex.SetTextAlign(12);  //align at top
   	latex.DrawLatex(-.5,10.,Form("#splitline{For %.0f signal events and %.0f bkg events}{the maximum S/sqrt(S+B) is %.2f when cutting at %.2f}",NUM_EVT_S,NUM_EVT_B,max,bdt_output[i]));
	}
