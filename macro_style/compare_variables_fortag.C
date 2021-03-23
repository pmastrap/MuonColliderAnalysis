#define DIM 49
#include "My_Beautiful_Canvas_TagVar.C"

void compare_variables_fortag(string FILEIN1="/home/paola/Scrivania/macro/Trees_tagger/b_tree_forMVA3_10k.root", string FILEIN2="/home/paola/Scrivania/macro/Trees_tagger/c_tree_forMVA3_10k.root", string FILEIN3 = "/home/paola/Scrivania/macro/Trees_tagger/o_tree_forMVA3_10k.root"){

	string str[DIM]{"vertex_cat","lepton_cat","n_SV","n_trk_from_SV","fd_sig_xy_SV","fd_sig_xyz_SV","corrected_mass_SV",
	"energy_ratio_SV","dR_SV_JET","MassEnergyFraction_SV","Boost_SV","n_trk_in_Jet","SIP_sig_xy_TRK1","SIP_sig_xy_TRK2","SIP_sig_xyz_TRK1",
	"SIP_sig_xyz_TRK2","eta_rel_TRK1","eta_rel_TRK2","pt_rel_TRK1","pt_rel_TRK2","pt_rel_permom_TRK1","pt_rel_permom_TRK2","pl_rel_TRK1",
	"pl_rel_TRK2","pl_rel_permom_TRK1","pl_rel_permom_TRK2","dR_jet_TRK1","dR_jet_TRK2", "dR_jet_tot_TRK", "Et_trk_Jet",
	"SIP_sig_xy_TRK_above_trh","SIP_sig_xyz_TRK_above_trh","dist_trk1_Jet_atPCA","dist_trk2_Jet_atPCA","decay_length_TRK1","decay_length_TRK2",
	"n_LEP","SIP_sig_xyz_LEP1","SIP_sig_xyz_LEP2","eta_rel_LEP1","eta_rel_LEP2","pt_rel_LEP1","pt_rel_LEP2","dR_jet_LEP1","dR_jet_LEP2", 		"pt_lj_LEP1","pt_lj_LEP2","pl_rel_perjetmom_LEP1","pl_rel_perjetmom_LEP2"};
	
	int n_bins[DIM];
	double low_edge[DIM],up_edge[DIM];
	
	TString titleX, titleY;
	
	TFile *f1 = new TFile(Form("%s",FILEIN1.c_str())) ;
  	TTree *tree1 = (TTree*)f1->Get("tree_tagging_var");
  	
  	TFile *f2 = new TFile(Form("%s",FILEIN2.c_str())) ;
  	TTree *tree2 = (TTree*)f2->Get("tree_tagging_var");
  	
  	TH1 *h3[DIM],*h1[DIM],*h2[DIM];
  	
  	TTree *tree3;
  	if(FILEIN3!="nofile"){
  		TFile *f3 = new TFile(Form("%s",FILEIN3.c_str())) ;
  		tree3 = (TTree*)f3->Get("tree_tagging_var");
  		}
	
	//ciclo su tutti branch del tree
	
	for(int i=0;i<DIM;i++){
	//if(str[i]=="SIP_sig_xy_TRK1" || str[i]=="SIP_sig_xy_TRK2"|| str[i]=="SIP_sig_xyz_TRK1"|| str[i]=="SIP_sig_xyz_TRK2"){
	//if(str[i]=="SIP_sig_xy_TRK_above_trh"|| str[i]=="SIP_sig_xyz_TRK_above_trh"){	
	//if(str[i]=="SIP_sig_xy_TRK2"){
	if(str[i]=="SIP_sig_xyz_TRK_above_trh"){
		gEnv->SetValue("Hist.Binning.1D.x",100); 
		gEnv->SetValue("Hist.Binning.2D.x",40); 
		gEnv->SetValue("Hist.Binning.3D.x",20);
		
		n_bins[i]=60;
	  	low_edge[i]=-30;
	  	up_edge[i]=30;
		
		tree1->Draw(Form("%s >> h1_t[%d](%d,%f,%f)",str[i].c_str(),i,n_bins[i],low_edge[i],up_edge[i]),Form("%s >-999",str[i].c_str()),"go off");
	  	//tree1->Draw(Form("%s >> h1_t[%d]",str[i].c_str(),i),Form("%s >-999",str[i].c_str()),"go off");
	  	h1[i] = (TH1*)gDirectory->Get(Form("h1_t[%d]",i));
	  	
	  	n_bins[i]=h1[i]->GetNbinsX();
	  	low_edge[i]=h1[i]->GetXaxis()->GetBinLowEdge(1);
	  	up_edge[i]=h1[i]->GetXaxis()->GetBinUpEdge(n_bins[i]);
	  	
	  	/*cout<<"n bins: "<<n_bins[i]<<endl;
	  	cout<<low_edge[i]<<endl;
	  	cout<<up_edge[i]<<endl;*/
	  	
	  	tree2->Draw(Form("%s >> h2_t[%d](%d,%f,%f)",str[i].c_str(),i,n_bins[i],low_edge[i],up_edge[i]),Form("%s >-999",str[i].c_str()),"go off");
	  	h2[i] = (TH1*)gDirectory->Get(Form("h2_t[%d]",i));
	  	if(FILEIN3!="nofile"){
	  		tree3->Draw(Form("%s >> h3_t[%d](%d,%f,%f)",str[i].c_str(),i,n_bins[i],low_edge[i],up_edge[i]),Form("%s >-999",str[i].c_str()),"go off");
	  		h3[i] =(TH1F*)gDirectory->Get(Form("h3_t[%d]",i));
	  		}
	  	
	  	
	  	auto legend = new TLegend(0.15,0.78,0.5,0.89);
	  	legend->AddEntry(h1[i],"b-jets (sample H #rightarrow b#bar{b})","l");
	  	legend->AddEntry(h2[i],"c-jets (sample H #rightarrow c#bar{c})","l");
	  	if(FILEIN3!="nofile"){
	  		legend->AddEntry(h3[i],"light-jets","l");
	  		}
	  		
	  	h1[i]->SetLineWidth(2);
	  	h1[i]->SetLineColor(kRed);
	  	h1[i]->SetName(Form("%s_b",str[i].c_str()));
	  	h1[i]->SetTitle(Form("%s",str[i].c_str()));
	  	//overflow in the last bin
	  	h1[i]->SetBinContent(h1[i]->GetNbinsX(),h1[i]->GetBinContent(h1[i]->GetNbinsX()) + h1[i]->GetBinContent(h1[i]->GetNbinsX() + 1));
	  	//underflow in the first bin
		h1[i]->SetBinContent(1,h1[i]->GetBinContent(1) + h1[i]->GetBinContent(0));
	  	h1[i]->Scale(1/h1[i]->Integral());

	  	h2[i]->SetLineWidth(2);
	  	h2[i]->SetLineColor(kGreen+2);
	  	h2[i]->SetName(Form("%s_c",str[i].c_str()));
	  	h2[i]->SetTitle(Form("%s",str[i].c_str()));
	  	//overflow in the last bin
	  	h2[i]->SetBinContent(h2[i]->GetNbinsX(),h2[i]->GetBinContent(h2[i]->GetNbinsX()) + h2[i]->GetBinContent(h2[i]->GetNbinsX() + 1));
	  	//underflow in the first bin
		h2[i]->SetBinContent(1,h2[i]->GetBinContent(1) + h2[i]->GetBinContent(0));
	  	h2[i]->Scale(1/h2[i]->Integral());
	  	
	  	
	  	if(FILEIN3!="nofile"){
	  		h3[i]->SetLineWidth(2);
		  	h3[i]->SetLineColor(kBlue);
		  	h3[i]->SetName(Form("%s_light",str[i].c_str()));
		  	h3[i]->SetTitle(Form("%s",str[i].c_str()));
		  	//overflow in the last bin
	  		h3[i]->SetBinContent(h3[i]->GetNbinsX(),h3[i]->GetBinContent(h3[i]->GetNbinsX()) + h3[i]->GetBinContent(h3[i]->GetNbinsX() + 1));
	  		//underflow in the first bin
			h3[i]->SetBinContent(1,h3[i]->GetBinContent(1) + h3[i]->GetBinContent(0));
		  	h3[i]->Scale(1/h3[i]->Integral());
		  	}
		titleX="3D SIP/#sigma of first track above M_{INV}- threhsold";
	  	titleY="Jets / 1 unit";
	  	My_Beautiful_Canvas_TagVar(h1[i],h2[i],h3[i],titleX,titleY);
		
		}
		}
	f1->Close();
	f2->Close();
  	
  	}
