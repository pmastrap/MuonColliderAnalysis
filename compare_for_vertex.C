void compare_for_vertex(){

	gStyle->SetOptStat("nemrou");
	
	//bb
	TFile *f = new TFile("/home/paola/Scrivania/vertex_bb_1000.root");
	
	TH1F * SV_reco_perevent_bb = (TH1F*)f->Get("SV_reco_perevent");
	SV_reco_perevent_bb->SetTitle("Number of SVs reconstructed per event");
   	SV_reco_perevent_bb->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_perevent_bb->SetLineWidth(2);
	
	TH1F * SV_reco_perjet_bb= (TH1F*)f->Get("SV_reco_perjet");
	SV_reco_perjet_bb->SetTitle("Number of SVs reconstructed per jet");
   	SV_reco_perjet_bb->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_perjet_bb->SetLineWidth(2);
	
	TH1F* jet_fraction_vertex_cat_bb=(TH1F*)f->Get("jet_fraction_vertex_cat");
	jet_fraction_vertex_cat_bb->SetTitle("Jet fraction per vertex category");
	jet_fraction_vertex_cat_bb->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_j_bb=(TH1F*)f->Get("SV_to_PV_distance_j");
	SV_to_PV_distance_j_bb->SetTitle("SV flight distance (jet-matched only)");
   	SV_to_PV_distance_j_bb->GetXaxis()->SetTitle("distance[mm]");
   	SV_to_PV_distance_j_bb->SetLineColor(kBlue);
	SV_to_PV_distance_j_bb->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_sig_j_bb=(TH1F*)f->Get("SV_to_PV_distance_sig_j");
	SV_to_PV_distance_sig_j_bb->SetTitle("SV flight distance significance (jet-matched only)");
   	SV_to_PV_distance_sig_j_bb->GetXaxis()->SetTitle("significance");
   	SV_to_PV_distance_sig_j_bb->SetLineColor(kBlue);
	SV_to_PV_distance_sig_j_bb->SetLineWidth(2);
	
	TH1F* deltaR_SV_nearestJet_bb=(TH1F*)f->Get("deltaR_SV_nearestJet");
	deltaR_SV_nearestJet_bb->SetTitle("#DeltaR between SV flight direction and its nearest jet axis");
   	deltaR_SV_nearestJet_bb->GetXaxis()->SetTitle("#DeltaR");
	deltaR_SV_nearestJet_bb->SetLineWidth(2);
	
	TH1F* vtx_corrected_mass_bb=(TH1F*)f->Get("vtx_corrected_mass");
	vtx_corrected_mass_bb->SetTitle("Corrected SV mass (jet-matched only)");
   	vtx_corrected_mass_bb->GetXaxis()->SetTitle("M[GeV]");
   	vtx_corrected_mass_bb->SetLineColor(kBlue);
	vtx_corrected_mass_bb->SetLineWidth(2);
	
	TH1F* chi2_per_num_trk_SV_bb=(TH1F*)f->Get("chi2_per_num_trk_SV");
   	chi2_per_num_trk_SV_bb->SetTitle("SV #chi^{2} per number of tracks used to fit SV (jet-matched only)");
   	chi2_per_num_trk_SV_bb->GetXaxis()->SetTitle("#chi^{2}/N_{trk}");
   	chi2_per_num_trk_SV_bb->GetYaxis()->SetTitle("");
   	chi2_per_num_trk_SV_bb->SetLineWidth(2);
   	
   	
   	TH1F* from_hh_wrong_to_PV_bb=(TH1F*)f->Get("from_hh_wrong_to_PV");
   	from_hh_wrong_to_PV_bb->SetTitle("Number of tracks from B/C decay wrongly associated to PV per event");
   	from_hh_wrong_to_PV_bb->GetXaxis()->SetTitle("N_{trk}");
	from_hh_wrong_to_PV_bb->SetLineWidth(2);
	
  	TH1F* from_IP_wrong_to_SV_bb=(TH1F*)f->Get("from_IP_wrong_to_SV");
  	from_IP_wrong_to_SV_bb->SetTitle("Number of tracks from hadronization wrongly associated to SV per event");
   	from_IP_wrong_to_SV_bb->GetXaxis()->SetTitle("N_{trk}");
	from_IP_wrong_to_SV_bb->SetLineWidth(2);
	
  	TH1F* num_trk_forPV_bb=(TH1F*)f->Get("num_trk_forPV");
  	num_trk_forPV_bb->SetTitle("Number of tracks used to fit PV");
   	num_trk_forPV_bb->GetXaxis()->SetTitle("N_{trk}");
	num_trk_forPV_bb->SetLineWidth(2);
	
  	TH1F* num_trk_forSV_bb=(TH1F*)f->Get("num_trk_forSV");
	num_trk_forSV_bb->SetTitle("Number of tracks used to fit SV");
   	num_trk_forSV_bb->GetXaxis()->SetTitle("N_{trk}");
	num_trk_forSV_bb->SetLineWidth(2);
   	
   	
   	//cc
   	TFile *f2 = new TFile("/home/paola/Scrivania/vertex_cc_1000.root");
   	
   	TH1F * SV_reco_perevent = (TH1F*)f2->Get("SV_reco_perevent");
	SV_reco_perevent->SetTitle("Number of SVs reconstructed per event");
   	SV_reco_perevent->GetXaxis()->SetTitle("N_{SV}");
   	SV_reco_perevent->SetLineColor(kRed);
	SV_reco_perevent->SetLineWidth(2);
	
	TH1F * SV_reco_perjet = (TH1F*)f2->Get("SV_reco_perjet");
	SV_reco_perjet->SetTitle("Number of SVs reconstructed per jet");
   	SV_reco_perjet->GetXaxis()->SetTitle("N_{SV}");
   	SV_reco_perjet->SetLineColor(kRed);
	SV_reco_perjet->SetLineWidth(2);
	
	TH1F* jet_fraction_vertex_cat=(TH1F*)f2->Get("jet_fraction_vertex_cat");
	jet_fraction_vertex_cat->SetTitle("Jet fraction per vertex category");
	jet_fraction_vertex_cat->SetLineColor(kRed);
	jet_fraction_vertex_cat->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_j=(TH1F*)f2->Get("SV_to_PV_distance_j");
	SV_to_PV_distance_j->SetTitle("SV flight distance (jet-matched only)");
   	SV_to_PV_distance_j->GetXaxis()->SetTitle("distance[mm]");
   	SV_to_PV_distance_j->SetLineColor(kRed);
	SV_to_PV_distance_j->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_sig_j=(TH1F*)f2->Get("SV_to_PV_distance_sig_j");
	SV_to_PV_distance_sig_j->SetTitle("SV flight distance significance (jet-matched only)");
   	SV_to_PV_distance_sig_j->GetXaxis()->SetTitle("significance");
   	SV_to_PV_distance_sig_j->SetLineColor(kRed);
	SV_to_PV_distance_sig_j->SetLineWidth(2);
	
	TH1F* deltaR_SV_nearestJet=(TH1F*)f2->Get("deltaR_SV_nearestJet");
	deltaR_SV_nearestJet->SetTitle("#DeltaR between SV flight direction and its nearest jet axis");
   	deltaR_SV_nearestJet->GetXaxis()->SetTitle("#DeltaR");
   	deltaR_SV_nearestJet->SetLineColor(kRed);
	deltaR_SV_nearestJet->SetLineWidth(2);
	
	TH1F* vtx_corrected_mass=(TH1F*)f2->Get("vtx_corrected_mass");
	vtx_corrected_mass->SetTitle("Corrected SV mass (jet-matched only)");
   	vtx_corrected_mass->GetXaxis()->SetTitle("M[GeV]");
   	vtx_corrected_mass->SetLineColor(kRed);
	vtx_corrected_mass->SetLineWidth(2);
	
	TH1F* chi2_per_num_trk_SV=(TH1F*)f2->Get("chi2_per_num_trk_SV");
   	chi2_per_num_trk_SV->SetTitle("SV #chi^{2} per number of tracks used to fit SV (jet-matched only)");
   	chi2_per_num_trk_SV->GetXaxis()->SetTitle("#chi^{2}/N_{trk}");
   	chi2_per_num_trk_SV->GetYaxis()->SetTitle("");
   	chi2_per_num_trk_SV->SetLineColor(kRed);
   	chi2_per_num_trk_SV->SetLineWidth(2);
   	
   	TH1F* from_hh_wrong_to_PV=(TH1F*)f2->Get("from_hh_wrong_to_PV");
   	from_hh_wrong_to_PV->SetTitle("Number of tracks from B/C decay wrongly associated to PV per event");
   	from_hh_wrong_to_PV->GetXaxis()->SetTitle("N_{trk}");
   	from_hh_wrong_to_PV->SetLineColor(kRed);
	from_hh_wrong_to_PV->SetLineWidth(2);
	
  	TH1F* from_IP_wrong_to_SV=(TH1F*)f2->Get("from_IP_wrong_to_SV");
  	from_IP_wrong_to_SV->SetTitle("Number of tracks from hadronization wrongly associated to SV per event");
   	from_IP_wrong_to_SV->GetXaxis()->SetTitle("N_{trk}");
   	from_IP_wrong_to_SV->SetLineColor(kRed);
	from_IP_wrong_to_SV->SetLineWidth(2);
	
  	TH1F* num_trk_forPV=(TH1F*)f2->Get("num_trk_forPV");
  	num_trk_forPV->SetTitle("Number of tracks used to fit PV");
   	num_trk_forPV->GetXaxis()->SetTitle("N_{trk}");
   	num_trk_forPV->SetLineColor(kRed);
	num_trk_forPV->SetLineWidth(2);
	
  	TH1F* num_trk_forSV=(TH1F*)f2->Get("num_trk_forSV");
	num_trk_forSV->SetTitle("Number of tracks used to fit SV");
   	num_trk_forSV->GetXaxis()->SetTitle("N_{trk}");
   	num_trk_forSV->SetLineColor(kRed);
	num_trk_forSV->SetLineWidth(2);
   	
   	//f2->Close();
   	
   	//draw
   	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
     	legend2->AddEntry(SV_reco_perevent, "H #rightarrow c#bar{c}", "l");
     	legend2->AddEntry(SV_reco_perevent_bb, "H #rightarrow b#bar{b}", "l");
     	legend2->SetTextSize(0.04);
   	
   	TCanvas* c1=new TCanvas();
   	SV_reco_perevent->Draw("s hist");
        c1->Update();
        TPaveStats* stat1 = (TPaveStats*)SV_reco_perevent->GetListOfFunctions()->FindObject("stats");
   	stat1->SetTextColor(kRed);
   	SV_reco_perevent_bb->Draw("sames hist");
        c1->Update();
        TPaveStats* stat2 = (TPaveStats*)SV_reco_perevent_bb->GetListOfFunctions()->FindObject("stats");
   	stat2->SetTextColor(kBlue);
   	legend2->Draw();
   	
   	TCanvas* c2=new TCanvas();
   	SV_reco_perjet->Draw("s hist");
        c2->Update();
        TPaveStats* stat3 = (TPaveStats*)SV_reco_perjet->GetListOfFunctions()->FindObject("stats");
   	stat3->SetTextColor(kRed);
   	SV_reco_perjet_bb->Draw("sames hist");
        c2->Update();
        TPaveStats* stat4 = (TPaveStats*)SV_reco_perjet_bb->GetListOfFunctions()->FindObject("stats");
   	stat4->SetTextColor(kBlue);
   	legend2->Draw();
   	
   	TCanvas* c3=new TCanvas();
   	jet_fraction_vertex_cat->Draw("s hist");
        c3->Update();
        /*TPaveStats* stat5 = (TPaveStats*)jet_fraction_vertex_cat->GetListOfFunctions()->FindObject("stats");
   	stat5->SetTextColor(kRed);*/
   	jet_fraction_vertex_cat_bb->Draw("sames hist");
        c3->Update();
        /*TPaveStats* stat6 = (TPaveStats*)jet_fraction_vertex_cat_bb->GetListOfFunctions()->FindObject("stats");
   	stat6->SetTextColor(kBlue);*/
   	legend2->Draw();
   	
   	TCanvas* c4=new TCanvas();
   	SV_to_PV_distance_j->Draw("s hist");
        c4->Update();
        TPaveStats* stat7 = (TPaveStats*)SV_to_PV_distance_j->GetListOfFunctions()->FindObject("stats");
   	stat7->SetTextColor(kRed);
   	SV_to_PV_distance_j_bb->Draw("sames hist");
        c4->Update();
        TPaveStats* stat8 = (TPaveStats*)SV_to_PV_distance_j_bb->GetListOfFunctions()->FindObject("stats");
   	stat8->SetTextColor(kBlue);
   	legend2->Draw();
   	
   	TCanvas* c5=new TCanvas();
   	SV_to_PV_distance_sig_j->Draw("s hist");
        c5->Update();
        TPaveStats* stat9 = (TPaveStats*)SV_to_PV_distance_sig_j->GetListOfFunctions()->FindObject("stats");
   	stat9->SetTextColor(kRed);
   	SV_to_PV_distance_sig_j_bb->Draw("sames hist");
        c5->Update();
        TPaveStats* stat10 = (TPaveStats*)SV_to_PV_distance_sig_j_bb->GetListOfFunctions()->FindObject("stats");
   	stat10->SetTextColor(kBlue);
   	legend2->Draw();
   	
   	
   	TCanvas* c6=new TCanvas();
   	deltaR_SV_nearestJet_bb->Draw("s hist");
        c6->Update();
        TPaveStats* stat11 = (TPaveStats*)deltaR_SV_nearestJet_bb->GetListOfFunctions()->FindObject("stats");
   	stat11->SetTextColor(kBlue);
   	deltaR_SV_nearestJet->Draw("sames hist");
        c6->Update();
        TPaveStats* stat12 = (TPaveStats*)deltaR_SV_nearestJet->GetListOfFunctions()->FindObject("stats");
   	stat12->SetTextColor(kRed);
   	legend2->Draw();
   	
   	TCanvas* c7=new TCanvas();
   	vtx_corrected_mass->Draw("s hist");
        c7->Update();
        TPaveStats* stat300 = (TPaveStats*)vtx_corrected_mass->GetListOfFunctions()->FindObject("stats");
   	stat300->SetTextColor(kRed);
   	vtx_corrected_mass_bb->Draw("sames hist");
        c7->Update();
        TPaveStats* stat30 = (TPaveStats*)vtx_corrected_mass_bb->GetListOfFunctions()->FindObject("stats");
   	stat30->SetTextColor(kBlue);
   	legend2->Draw();

   	TCanvas* c8=new TCanvas();
   	chi2_per_num_trk_SV_bb->Draw("s hist");
        c8->Update();
        TPaveStats* stat13 = (TPaveStats*)chi2_per_num_trk_SV_bb->GetListOfFunctions()->FindObject("stats");
   	stat13->SetTextColor(kBlue);
   	chi2_per_num_trk_SV->Draw("sames hist");
        c8->Update();
        TPaveStats* stat14 = (TPaveStats*)chi2_per_num_trk_SV->GetListOfFunctions()->FindObject("stats");
   	stat14->SetTextColor(kRed);
   	legend2->Draw();
   	
   	TCanvas* c9=new TCanvas();
   	from_hh_wrong_to_PV->Draw("s hist");
        c9->Update();
        TPaveStats* stat15 = (TPaveStats*)from_hh_wrong_to_PV->GetListOfFunctions()->FindObject("stats");
   	stat15->SetTextColor(kBlue);
   	from_hh_wrong_to_PV_bb->Draw("sames hist");
        c9->Update();
        TPaveStats* stat16 = (TPaveStats*)from_hh_wrong_to_PV_bb->GetListOfFunctions()->FindObject("stats");
   	stat16->SetTextColor(kRed);
   	legend2->Draw();
   	
   	
   	TCanvas* c10=new TCanvas();
   	from_IP_wrong_to_SV_bb->Draw("s hist");
        c10->Update();
        TPaveStats* stat17 = (TPaveStats*)from_IP_wrong_to_SV_bb->GetListOfFunctions()->FindObject("stats");
   	stat17->SetTextColor(kBlue);
   	from_IP_wrong_to_SV->Draw("sames hist");
        c10->Update();
        TPaveStats* stat18 = (TPaveStats*)from_IP_wrong_to_SV->GetListOfFunctions()->FindObject("stats");
   	stat18->SetTextColor(kRed);
   	legend2->Draw();
   	
   	TCanvas* c11=new TCanvas();
   	num_trk_forPV->Draw("s hist");
        c11->Update();
        TPaveStats* stat19 = (TPaveStats*)num_trk_forPV->GetListOfFunctions()->FindObject("stats");
   	stat19->SetTextColor(kRed);
   	num_trk_forPV_bb->Draw("sames hist");
        c11->Update();
        TPaveStats* stat20 = (TPaveStats*)num_trk_forPV_bb->GetListOfFunctions()->FindObject("stats");
   	stat20->SetTextColor(kBlue);
   	legend2->Draw();
   	
  	TCanvas* c12=new TCanvas();
   	num_trk_forSV_bb->Draw("s hist");
        c12->Update();
        TPaveStats* stat21 = (TPaveStats*)num_trk_forSV_bb->GetListOfFunctions()->FindObject("stats");
   	stat21->SetTextColor(kBlue);
   	num_trk_forSV->Draw("sames hist");
        c12->Update();
        TPaveStats* stat22 = (TPaveStats*)num_trk_forSV->GetListOfFunctions()->FindObject("stats");
   	stat22->SetTextColor(kRed);
   	legend2->Draw();
  	
   	
   	}
	
	
