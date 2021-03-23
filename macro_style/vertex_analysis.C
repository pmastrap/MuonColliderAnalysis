////////DESCRIPTION///////////////////////////////////////////////////////////////////////////////////////////////////////
//Con bb sample:													//
//1. Vertici secondari ricostruiti per evento										//
//2. Vertici secondari per jet												//
//2b. Jet fraction per vertex category: 0 SV, almeno 1 SV								//
//3. Distanza PV-SV + distanza/errore (flight distance significance)							//
//4. deltaR tra vertice secondario e jet più vicino									//
//5. massa corretta SV (+ singole variabili della formula) vedi pag 10 paper:						//
//	"Identification of heavy-flavour jets with the CMS detector in pp collisions at 13 TeV"	 			//
//5b. stesso ma selezionando solo SV associati a vertici + selezione SV da decadimento B/C/resto			//
//6. Numero di tracce con cui è fittato primario/secondario								//
//7. chi2 vertice secondario/primario VS numero di tracce con cui è ricostruito						//
//8. quante tracce che provengono da PV sono associate erroneamente a SV per evento					//
//9. quante tracce che provengono dal decadimento di B/C sono associate erroneamente al PV per evento			//
//10. numero di SV ricostruiti con particelle che provengono da B/C/resto per evento					//
//11. rapporto energia tracce associate al vertice/energia tracce nel jet 						//
//12. Per 5 eventi stampo info dettagliate su quali sono le particelle a MC associate a quali vertici ricostruiti	//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TLorentzVector.h"

//SELECTION FOR CC
int selection_cc(int mcid){

	if(abs(mcid)==411 || abs(mcid)== 421 || abs(mcid)==10411 || abs(mcid)==10421 || abs(mcid)==413 ||
		   abs(mcid)==423 || abs(mcid)==10413 || abs(mcid)==10423 || abs(mcid)==20413 || abs(mcid)==20423 ||
     		   abs(mcid)==415 || abs(mcid)==425 || abs(mcid)==431 || abs(mcid)==10431 || abs(mcid)==433 || 			   abs(mcid)==10433 ||abs(mcid)==20433 || abs(mcid)==435 || abs(mcid)==4122 || abs(mcid)==4222 || abs(mcid)==4212 || abs(mcid)==4112 || abs(mcid)==4224 || abs(mcid)==4214 || abs(mcid)==4114 || abs(mcid)==4232 || abs(mcid)==4132 || abs(mcid)==4322|| abs(mcid)==4312 || abs(mcid)==4324 || abs(mcid)==4314 || abs(mcid)==4332 || abs(mcid)==4334 || abs(mcid)==4412 || abs(mcid)==4422 || abs(mcid)==4414 || abs(mcid)==4424 || abs(mcid)==4432 || abs(mcid)==4434 || abs(mcid)==4444 || abs(mcid)==441 || abs(mcid)==10441 || abs(mcid)==100441 || abs(mcid)==443 || abs(mcid)==10443 || abs(mcid)==20443 || abs(mcid)==100443 || abs(mcid)==30443 || abs(mcid)==9000443 || abs(mcid)==9010443 || abs(mcid)==9020443 || abs(mcid)==445 || abs(mcid)==100445 ){ return 0;}

	else return -1;

}


//SELECTION FOR BB
int selection_bb(int mcid){

	if(abs(mcid)==511 || abs(mcid)==521 || abs(mcid)==10511 || abs(mcid)==10521 || abs(mcid)==513 || abs(mcid)==523 || abs(mcid)==10513 || abs(mcid)==10523 || abs(mcid)==20513 || abs(mcid)==20523 || abs(mcid)==515 || abs(mcid)==525 || abs(mcid)==531 || abs(mcid)==10531 || abs(mcid)==533 || abs(mcid)==10533 || abs(mcid)==20533 || abs(mcid)==535 || abs(mcid)==541 || abs(mcid)==10541 || abs(mcid)==543 || abs(mcid)==10543 || abs(mcid)==20543 || abs(mcid)==545 || abs(mcid)==551 || abs(mcid)==10551 || abs(mcid)==100551 || abs(mcid)==110551 || abs(mcid)==200551 || abs(mcid)==210551 || abs(mcid)==553 || abs(mcid)==10553 || abs(mcid)==20553 || abs(mcid)==30553 || abs(mcid)==100553 || abs(mcid)==110553 || abs(mcid)==120553 || abs(mcid)==130553 || abs(mcid)==200553 || abs(mcid)==210553 || abs(mcid)==220553 || abs(mcid)==300553 || abs(mcid)==9000553 || abs(mcid)==9010553 || abs(mcid)==555 || abs(mcid)==10555 || abs(mcid)==20555 || abs(mcid)==100555 || abs(mcid)==110555 || abs(mcid)==120555 || abs(mcid)==200555 || abs(mcid)==557 || abs(mcid)==100557 || abs(mcid)==5122 || abs(mcid)==5112 || abs(mcid)==5212 || abs(mcid)==5222 || abs(mcid)==5114 || abs(mcid)==5214 || abs(mcid)==5224 || abs(mcid)==5132 || abs(mcid)==5232 || abs(mcid)==5312 || abs(mcid)==5322 || abs(mcid)==5314 || abs(mcid)==5324 || abs(mcid)==5332 || abs(mcid)==5334 || abs(mcid)==5142 || abs(mcid)==5242 || abs(mcid)==5412 || abs(mcid)==5422 || abs(mcid)==5414 || abs(mcid)==5424 || abs(mcid)==5342 || abs(mcid)==5432 || abs(mcid)==5434 || abs(mcid)==5442 || abs(mcid)==5444 || abs(mcid)==5512 || abs(mcid)==5522 || abs(mcid)==5514 || abs(mcid)==5524 || abs(mcid)==5532 || abs(mcid)==5534 || abs(mcid)==5542 || abs(mcid)==5544 || abs(mcid)==5554){ return 0;}

	else return -1;
}


void vertex_analysis(){
	
	gStyle->SetOptStat("nemrou");
	
	Int_t nmcp,nrec,r2mnrel,nvt,nj,ntrk;
	
	//mc particles
	int num=2000000;
	Float_t *mcvtx,*mcvty,*mcvtz,*mcmox,*mcmoy,*mcmoz,*mcepx,*mcepy,*mcepz;
	Int_t *mcpdg,*mcpa0,*mcda0,*mcda1,*mcda2,*mcda3,*mcda4;
	mcpdg = (int*) malloc(sizeof(int)*num);
	mcpa0 = (int*) malloc(sizeof(int)*num);
	mcda0 = (int*) malloc(sizeof(int)*num);
	mcda1 = (int*) malloc(sizeof(int)*num);
	mcda2 = (int*) malloc(sizeof(int)*num);
	mcda3 = (int*) malloc(sizeof(int)*num);
	mcda4 = (int*) malloc(sizeof(int)*num);
	mcvtx = (float*) malloc(sizeof(float)*num);
	mcvty = (float*) malloc(sizeof(float)*num);
	mcvtz = (float*) malloc(sizeof(float)*num);
	mcepx = (float*) malloc(sizeof(float)*num);
	mcepy = (float*) malloc(sizeof(float)*num);
	mcepz = (float*) malloc(sizeof(float)*num);
	mcmox = (float*) malloc(sizeof(float)*num);
	mcmoy = (float*) malloc(sizeof(float)*num);
	mcmoz = (float*) malloc(sizeof(float)*num);
	
	//relations
	int num2=100000;
	Int_t *r2mt,*r2mf,*rcftr,*rcvts;
	Float_t *r2mw;
	r2mt = (int*) malloc(sizeof(int)*num2);
	r2mf = (int*) malloc(sizeof(int)*num2);
	r2mw = (float*) malloc(sizeof(float)*num2);
	rcftr = (int*) malloc(sizeof(int)*num2);
	rcvts = (int*) malloc(sizeof(int)*num2);
	
	//reco particles
	Float_t *rcmox,*rcmoy,*rcmoz,*rcene;
	Int_t *rccha;
	rcmox = (float*) malloc(sizeof(float)*num2);
	rcmoy = (float*) malloc(sizeof(float)*num2);
	rcmoz = (float*) malloc(sizeof(float)*num2);
	rcene = (float*) malloc(sizeof(float)*num2);
	rccha = (int*) malloc(sizeof(int)*num2);
	
	//tracks
	Float_t *trk_z0,*trk_d0;
	Int_t *trk_atIP;
	trk_z0 = (float*) malloc(sizeof(float)*num2);
	trk_d0 = (float*) malloc(sizeof(float)*num2);
	trk_atIP = (int*) malloc(sizeof(int)*num2);
	Float_t cov[100000][15];
	
	//vertices
	Float_t x[1000],y[1000],z[1000],chi[1000],vtcov[1000][6];
	Int_t id[100];
	
	//jets
	Float_t jmox[100],jmoy[100],jmoz[100],jene[100];
	
	//WHEN YOU CHANGE THE PATH CHANGE ALSO THIS FLAG
	char choice[]="o";

	TChain* fChain = new TChain("fChain");
   	//fChain->Add("/home/paola/Scrivania/xml_aggiornati/Irriducible_bkg/ntuple_cc_nohiggs_1000.root/MyLCTuple");
	//fChain->Add("/home/paola/Scrivania/xml_aggiornati/Fixed_version/10k/ntuple_cc_10000_fix.root/MyLCTuple");
	fChain->Add("/home/paola/Scrivania/xml_aggiornati/Bkg/ntuple_udsg_1000.root/MyLCTuple");
	//fChain->Add("/home/paola/Scrivania/xml_aggiornati/Fixed_version/10k/ntuple_bb_10000_fix.root/MyLCTuple");
	//fChain->Add("/home/paola/Scrivania/xml_aggiornati/0.7R/ntuple_bb_10000_fix_0.7DELTAR.root/MyLCTuple");
	
	fChain->SetBranchAddress("nmcp", &nmcp);
	fChain->SetBranchAddress("mcpdg", mcpdg);
	fChain->SetBranchAddress("mcpa0", mcpa0);
	fChain->SetBranchAddress("mcda1", mcda1);
	fChain->SetBranchAddress("mcda2", mcda2);
	fChain->SetBranchAddress("mcda3", mcda3);
	fChain->SetBranchAddress("mcda4", mcda4);
	fChain->SetBranchAddress("mcda0", mcda0);
	fChain->SetBranchAddress("mcvtx", mcvtx);
	fChain->SetBranchAddress("mcvty", mcvty);
	fChain->SetBranchAddress("mcvtz", mcvtz);
	fChain->SetBranchAddress("mcepx", mcepx);
	fChain->SetBranchAddress("mcepy", mcepy);
	fChain->SetBranchAddress("mcepz", mcepz);
	fChain->SetBranchAddress("mcmox", mcmox);
	fChain->SetBranchAddress("mcmoy", mcmoy);
	fChain->SetBranchAddress("mcmoz", mcmoz);
	
	
	fChain->SetBranchAddress("r2mnrel", &r2mnrel);
	fChain->SetBranchAddress("r2mt", r2mt);
	fChain->SetBranchAddress("r2mf", r2mf);	
        fChain->SetBranchAddress("r2mw", r2mw);
        
        fChain->SetBranchAddress("nrec", &nrec);   
	fChain->SetBranchAddress("rcftr", rcftr);
	fChain->SetBranchAddress("rcvts", rcvts);
	fChain->SetBranchAddress("rcmox", rcmox);
  	fChain->SetBranchAddress("rcmoy", rcmoy);
  	fChain->SetBranchAddress("rcmoz", rcmoz);
	fChain->SetBranchAddress("rcene", rcene);
	fChain->SetBranchAddress("rccha", rccha);
	
	fChain->SetBranchAddress("ntrk",&ntrk);
	fChain->SetBranchAddress("tszze",trk_z0);
	fChain->SetBranchAddress("trsip", trk_atIP);
        fChain->SetBranchAddress("tsdze", trk_d0);
        fChain->SetBranchAddress("tscov", cov);
        
        fChain->SetBranchAddress("nvt", &nvt);
        fChain->SetBranchAddress("vtori", id);
	fChain->SetBranchAddress("vtxxx", x);
	fChain->SetBranchAddress("vtyyy", y);
	fChain->SetBranchAddress("vtzzz", z);
	fChain->SetBranchAddress("vtchi", chi);
	fChain->SetBranchAddress("vtcov", vtcov);
	
	fChain->SetBranchAddress("nj", &nj);
   	fChain->SetBranchAddress("jmox", jmox);
	fChain->SetBranchAddress("jmoy", jmoy);
	fChain->SetBranchAddress("jmoz", jmoz);
	fChain->SetBranchAddress("jene", jene);
	
	/// PLOT 1 //////////////////////////////
	TH1F* SV_reco_perevent=new TH1F("SV_reco_perevent","SV_reco_perevent",12,-0.5,11.5);
	SV_reco_perevent->SetTitle("Number of SVs reconstructed per event");
   	SV_reco_perevent->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_perevent->SetLineWidth(2);
	/// PLOT 2 //////////////////////////////
	TH1F* SV_reco_perjet=new TH1F("SV_reco_perjet","SV_reco_perjet",12,-0.5,11.5);
	SV_reco_perjet->SetTitle("Number of SVs reconstructed per jet");
   	SV_reco_perjet->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_perjet->SetLineWidth(2);
	/// PLOT 2b //////////////////////////////
	TH1F* jet_fraction_vertex_cat=new TH1F("jet_fraction_vertex_cat","jet_fraction_vertex_cat",2,-0.5,1.5);
	jet_fraction_vertex_cat->SetTitle("Jet fraction per vertex category");
	jet_fraction_vertex_cat->SetLineWidth(2);
	/// PLOT 3 //////////////////////////////
	TH1F* SV_to_PV_distance=new TH1F("SV_to_PV_distance","SV_to_PV_distance",150,0.,1.5);
	SV_to_PV_distance->SetTitle("SV flight distance");
   	SV_to_PV_distance->GetXaxis()->SetTitle("distance[mm]");
	SV_to_PV_distance->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_sig=new TH1F("SV_to_PV_distance_sig","SV_to_PV_distance_sig",100,0.,100);
	SV_to_PV_distance_sig->SetTitle("SV flight distance significance");
   	SV_to_PV_distance_sig->GetXaxis()->SetTitle("significance");
	SV_to_PV_distance_sig->SetLineWidth(2);
	/// PLOT 3b //////////////////////////////
	TH1F* SV_to_PV_distance_j=new TH1F("SV_to_PV_distance_j","SV_to_PV_distance_j",200,0.,50);
	SV_to_PV_distance_j->SetTitle("SV flight distance (jet-matched only)");
   	SV_to_PV_distance_j->GetXaxis()->SetTitle("distance[mm]");
   	SV_to_PV_distance_j->SetLineColor(kRed);
	SV_to_PV_distance_j->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_sig_j=new TH1F("SV_to_PV_distance_sig_j","SV_to_PV_distance_sig_j",150,0.,150);
	SV_to_PV_distance_sig_j->SetTitle("SV flight distance significance (jet-matched only)");
   	SV_to_PV_distance_sig_j->GetXaxis()->SetTitle("significance");
   	SV_to_PV_distance_sig_j->SetLineColor(kRed);
	SV_to_PV_distance_sig_j->SetLineWidth(2);
	/// PLOT 3c //////////////////////////////
	TH1F* SV_to_PV_distance_j_B=new TH1F("SV_to_PV_distance_j_B","SV_to_PV_distance_j_B",150,0.,5);
	SV_to_PV_distance_j_B->SetTitle("SV flight distance (jet-matched only)");
   	SV_to_PV_distance_j_B->GetXaxis()->SetTitle("distance[mm]");
   	SV_to_PV_distance_j_B->SetLineColor(kBlue);
	SV_to_PV_distance_j_B->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_sig_j_B=new TH1F("SV_to_PV_distance_sig_j_B","SV_to_PV_distance_sig_j_B",100,0.,50);
	SV_to_PV_distance_sig_j_B->SetTitle("SV flight distance significance (jet-matched only)");
   	SV_to_PV_distance_sig_j_B->GetXaxis()->SetTitle("significance");
   	SV_to_PV_distance_sig_j_B->SetLineColor(kBlue);
	SV_to_PV_distance_sig_j_B->SetLineWidth(2);
	//////////////////////////////////////////
	TH1F* SV_to_PV_distance_j_C=new TH1F("SV_to_PV_distance_j_C","SV_to_PV_distance_j_C",150,0.,5);
	SV_to_PV_distance_j_C->SetTitle("SV flight distance (jet-matched only)");
   	SV_to_PV_distance_j_C->GetXaxis()->SetTitle("distance[mm]");
   	SV_to_PV_distance_j_C->SetLineColor(kRed);
	SV_to_PV_distance_j_C->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_sig_j_C=new TH1F("SV_to_PV_distance_sig_j_C","SV_to_PV_distance_sig_j_C",100,0.,50);
	SV_to_PV_distance_sig_j_C->SetTitle("SV flight distance significance (jet-matched only)");
   	SV_to_PV_distance_sig_j_C->GetXaxis()->SetTitle("significance");
   	SV_to_PV_distance_sig_j_C->SetLineColor(kRed);
	SV_to_PV_distance_sig_j_C->SetLineWidth(2);
	//////////////////////////////////////////
	TH1F* SV_to_PV_distance_j_BC=new TH1F("SV_to_PV_distance_j_BC","SV_to_PV_distance_j_BC",150,0.,5);
	SV_to_PV_distance_j_BC->SetTitle("SV flight distance (jet-matched only)");
   	SV_to_PV_distance_j_BC->GetXaxis()->SetTitle("distance[mm]");
   	SV_to_PV_distance_j_BC->SetLineColor(kGreen+3);
	SV_to_PV_distance_j_BC->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_sig_j_BC=new TH1F("SV_to_PV_distance_sig_j_BC","SV_to_PV_distance_sig_j_BC",100,0.,50);
	SV_to_PV_distance_sig_j_BC->SetTitle("SV flight distance significance (jet-matched only)");
   	SV_to_PV_distance_sig_j_BC->GetXaxis()->SetTitle("significance");
   	SV_to_PV_distance_sig_j_BC->SetLineColor(kGreen+3);
	SV_to_PV_distance_sig_j_BC->SetLineWidth(2);
	//////////////////////////////////////////
	TH1F* SV_to_PV_distance_j_noBC=new TH1F("SV_to_PV_distance_j_noBC","SV_to_PV_distance_j_noBC",150,0.,5);
	SV_to_PV_distance_j_noBC->SetTitle("SV flight distance (jet-matched only)");
   	SV_to_PV_distance_j_noBC->GetXaxis()->SetTitle("distance[mm]");
   	SV_to_PV_distance_j_noBC->SetLineColor(kMagenta+1);
	SV_to_PV_distance_j_noBC->SetLineWidth(2);
	
	TH1F* SV_to_PV_distance_sig_j_noBC=new TH1F("SV_to_PV_distance_sig_j_noBC","SV_to_PV_distance_sig_j_noBC",100,0.,50);
	SV_to_PV_distance_sig_j_noBC->SetTitle("SV flight distance significance (jet-matched only)");
   	SV_to_PV_distance_sig_j_noBC->GetXaxis()->SetTitle("significance");
   	SV_to_PV_distance_sig_j_noBC->SetLineColor(kMagenta+1);
	SV_to_PV_distance_sig_j_noBC->SetLineWidth(2);
	/// PLOT 4 //////////////////////////////
	TH1F* deltaR_SV_nearestJet=new TH1F("deltaR_SV_nearestJet","deltaR_SV_nearestJet",100,0.,10);
	deltaR_SV_nearestJet->SetTitle("#DeltaR between SV flight direction and its nearest jet axis");
   	deltaR_SV_nearestJet->GetXaxis()->SetTitle("#DeltaR");
	deltaR_SV_nearestJet->SetLineWidth(2);
	/// PLOT 5 //////////////////////////////
	TH1F* vtx_invariant_mass=new TH1F("vtx_invariant_mass","vtx_invariant_mass",100,0.,12);
	vtx_invariant_mass->SetTitle("Invariant mass of tracks associated to SV");
   	vtx_invariant_mass->GetXaxis()->SetTitle("M_{INV}[GeV]");
	vtx_invariant_mass->SetLineWidth(2);
	
	TH1F* vtx_corrected_mass=new TH1F("vtx_corrected_mass","vtx_corrected_mass",50,0.,10);
	vtx_corrected_mass->SetTitle("Corrected SV mass (jet-matched only)");
   	vtx_corrected_mass->GetXaxis()->SetTitle("M[GeV]");
	vtx_corrected_mass->SetLineWidth(2);
	vtx_corrected_mass->SetLineColor(kRed);
	
	TH1F* vtx_mom=new TH1F("vtx_mom","vtx_mom",200,0.,100);
	vtx_mom->SetTitle("SV momentum");
   	vtx_mom->GetXaxis()->SetTitle("P[GeV/c]");
	vtx_mom->SetLineWidth(2);
	
	TH1F* vtx_theta=new TH1F("vtx_theta","vtx_theta",50,0.,3.3);
	vtx_theta->SetTitle("Angle between SV momentum and SV flight direction");
   	vtx_theta->GetXaxis()->SetTitle("#theta[rad]");
	vtx_theta->SetLineWidth(2);
	
	/// PLOT 5b //////////////////////////////
	TH1F* vtx_invariant_mass_j=new TH1F("vtx_invariant_mass_j","vtx_invariant_mass_j",100,0.,12);
	vtx_invariant_mass_j->SetTitle("Invariant mass of tracks associated to SV (jet-matched only)");
   	vtx_invariant_mass_j->GetXaxis()->SetTitle("M_{INV}[GeV]");
   	vtx_invariant_mass_j->SetLineColor(kRed);
	vtx_invariant_mass_j->SetLineWidth(2);
	
	TH1F* vtx_mom_j=new TH1F("vtx_mom_j","vtx_mom_j",200,0.,100);
	vtx_mom_j->SetTitle("SV momentum (jet-matched only)");
	vtx_mom_j->SetLineColor(kRed);
   	vtx_mom_j->GetXaxis()->SetTitle("P[GeV/c]");
	vtx_mom_j->SetLineWidth(2);
	
	TH1F* vtx_theta_j=new TH1F("vtx_theta_j","vtx_theta_j",50,0.,3.3);
	vtx_theta_j->SetTitle("Angle between SV momentum and SV flight direction(jet-matched only)");
   	vtx_theta_j->GetXaxis()->SetTitle("#theta[rad]");
   	vtx_theta_j->SetLineColor(kRed);
	vtx_theta_j->SetLineWidth(2);
	
	//only for the corrected mass I differentiate between B/C/other SV
	double bins[26]={0};
	int in;
	for(in=0;in<16;in++){
		bins[in]=in*1;
		}
	for(in=16;in<26;in++){
		bins[in]=15+(in-15)*2;
		}
	TH1F* vtx_corrected_mass_j_B=new TH1F("vtx_corrected_mass_j_B","vtx_corrected_mass_j_B",50,0.,10);
	vtx_corrected_mass_j_B->SetTitle("Corrected SV mass (jet-matched only)");
   	vtx_corrected_mass_j_B->GetXaxis()->SetTitle("M[GeV]");
	vtx_corrected_mass_j_B->SetLineWidth(2);
	
	TH1F* vtx_corrected_mass_j_C=new TH1F("vtx_corrected_mass_j_C","vtx_corrected_mass_j_C",50,0.,10);
	vtx_corrected_mass_j_C->SetTitle("Corrected SV mass (jet-matched only)");
   	vtx_corrected_mass_j_C->GetXaxis()->SetTitle("M[GeV]");
   	vtx_corrected_mass_j_C->SetLineColor(kRed);
	vtx_corrected_mass_j_C->SetLineWidth(2);
	
	TH1F* vtx_corrected_mass_j_noBC=new TH1F("vtx_corrected_mass_j_noBC","vtx_corrected_mass_j_noBC",50,0.,10);
	vtx_corrected_mass_j_noBC->SetTitle("Corrected SV mass (jet-matched only)");
   	vtx_corrected_mass_j_noBC->GetXaxis()->SetTitle("M[GeV]");
   	vtx_corrected_mass_j_noBC->SetLineColor(kMagenta+1);
	vtx_corrected_mass_j_noBC->SetLineWidth(2);
	
	TH1F* vtx_corrected_mass_j_BC=new TH1F("vtx_corrected_mass_j_BC","vtx_corrected_mass_j_BC",50,0.,10);
	vtx_corrected_mass_j_BC->SetTitle("Corrected SV mass (jet-matched only)");
   	vtx_corrected_mass_j_BC->GetXaxis()->SetTitle("M[GeV]");
   	vtx_corrected_mass_j_BC->SetLineColor(kGreen+3);
	vtx_corrected_mass_j_BC->SetLineWidth(2);
	
	// PLOT 6 ///////////////////////////////
	TH1F* num_trk_forPV=new TH1F("num_trk_forPV","num_trk_forPV",30,-0.5,29.5);
	num_trk_forPV->SetTitle("Number of tracks used to fit PV");
   	num_trk_forPV->GetXaxis()->SetTitle("N_{trk}");
	num_trk_forPV->SetLineWidth(2);
	
	TH1F* num_trk_forSV=new TH1F("num_trk_forSV","num_trk_forSV",17,-0.5,16.5);
	num_trk_forSV->SetTitle("Number of tracks used to fit SV");
   	num_trk_forSV->GetXaxis()->SetTitle("N_{trk}");
	num_trk_forSV->SetLineWidth(2);
	// PLOT 7 ///////////////////////////////
	/*TH2F* chi2_VS_num_trk_PV=new TH2F("chi2_VS_num_trk_PV","chi2_VS_num_trk_PV",30,-0.5,29.5,55,0,55);
	chi2_VS_num_trk_PV->SetTitle("PV #chi^{2} VS number of tracks used to fit PV");
   	chi2_VS_num_trk_PV->GetXaxis()->SetTitle("N_{trk}");
   	chi2_VS_num_trk_PV->GetYaxis()->SetTitle("#chi^{2}");
	
	TH2F* chi2_VS_num_trk_SV=new TH2F("chi2_VS_num_trk_SV","chi2_VS_num_trk_SV",17,-0.5,16.5,40,0,20);
	chi2_VS_num_trk_SV->SetTitle("SV #chi^{2} VS number of tracks used to fit SV");
   	chi2_VS_num_trk_SV->GetXaxis()->SetTitle("N_{trk}");
   	chi2_VS_num_trk_SV->GetYaxis()->SetTitle("#chi^{2}");*/
   	TH1F* chi2_per_num_trk_PV=new TH1F("chi2_per_num_trk_PV","chi2_per_num_trk_PV",50,0,10);
   	chi2_per_num_trk_PV->SetTitle("PV #chi^{2} per number of tracks used to fit PV");
   	chi2_per_num_trk_PV->GetYaxis()->SetTitle("#chi^{2}/N_{trk}");
   	
   	TH1F* chi2_per_num_trk_SV=new TH1F("chi2_per_num_trk_SV","chi2_per_num_trk_SV",50,0,10);
   	chi2_per_num_trk_SV->SetTitle("SV #chi^{2} per number of tracks used to fit SV (jet-matched only)");
   	chi2_per_num_trk_SV->GetYaxis()->SetTitle("#chi^{2}/N_{trk}");
   	
	// PLOT 8 ///////////////////////////////
	TH1F* from_hh_wrong_to_PV=new TH1F("from_hh_wrong_to_PV","from_hh_wrong_to_PV",10,-0.5,9.5);
	from_hh_wrong_to_PV->SetTitle("Number of tracks from B/C decay wrongly associated to PV per event");
   	from_hh_wrong_to_PV->GetXaxis()->SetTitle("N_{trk}");
	from_hh_wrong_to_PV->SetLineWidth(2);
	// PLOT 9 ///////////////////////////////
	TH1F* from_IP_wrong_to_SV=new TH1F("from_IP_wrong_to_SV","from_IP_wrong_to_SV",10,-0.5,9.5);
	from_IP_wrong_to_SV->SetTitle("Number of tracks from hadronization wrongly associated to SV per event");
   	from_IP_wrong_to_SV->GetXaxis()->SetTitle("N_{trk}");
	from_IP_wrong_to_SV->SetLineWidth(2);
	// PLOT 10 ///////////////////////////////
	TH1F* SV_reco_B=new TH1F("SV_reco_B","SV_reco_B",10,-0.5,9.5);
	SV_reco_B->SetTitle("Number of SV (per event) recostructed with at least a track from B decay");
   	SV_reco_B->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_B->SetLineWidth(2);
	
	TH1F* SV_reco_C=new TH1F("SV_reco_C","SV_reco_C",10,-0.5,9.5);
	SV_reco_C->SetTitle("Number of SV (per event) recostructed with at least a track from C decay");
   	SV_reco_C->GetXaxis()->SetTitle("N_{SV}");
   	SV_reco_C->SetLineColor(kRed);
	SV_reco_C->SetLineWidth(2);
	
	TH1F* SV_reco_BC=new TH1F("SV_reco_BC","SV_reco_BC",10,-0.5,9.5);
	SV_reco_BC->SetTitle("Number of SV (per event) recostructed with at least a track both from B and C decay");
   	SV_reco_BC->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_BC->SetLineWidth(2);
	SV_reco_BC->SetLineColor(kGreen+3);
	
	TH1F* SV_reco_noBC=new TH1F("SV_reco_noBC","SV_reco_noBC",10,-0.5,9.5);
	SV_reco_noBC->SetTitle("Number of SV (per event) recostructed with no tracks from B/C decay");
   	SV_reco_noBC->GetXaxis()->SetTitle("N_{SV}");
	SV_reco_noBC->SetLineWidth(2);
	SV_reco_noBC->SetLineColor(kMagenta+1);
	// PLOT 11 //////////////////////////////
	TH1F* SV_energy_ratio=new TH1F("SV_energy_ratio","SV_energy_ratio",100,0.,3.);
	SV_energy_ratio->SetTitle("Energy ratio between tracks associated to SV and all the tracks in the jet");
   	SV_energy_ratio->GetXaxis()->SetTitle("Energy ratio");
	SV_energy_ratio->SetLineWidth(2);
	
	//addition
  	TH2F *add=new TH2F("add","add",50,0.,3.3,100,0.,100);
	add->GetXaxis()->SetTitle("theta");
	add->GetYaxis()->SetTitle("SV flight distance significance");
	
	TH1F *errx_PV=new TH1F("errx_PV","errx_PV",50,0.,0.01);
	TH1F *erry_PV=new TH1F("erry_PV","erry_PV",50,0.,0.01);
	TH1F *errz_PV=new TH1F("errz_PV","errz_PV",50,0.,0.1);
	
	TH1F *errx_SV=new TH1F("errx_SV","errx_SV",100,0.,1);
	TH1F *erry_SV=new TH1F("erry_SV","erry_SV",100,0.,1);
	TH1F *errz_SV=new TH1F("errz_SV","errz_SV",100,0.,1);
	
	TH1F *errxy_PV=new TH1F("errxy_PV","errxy_PV",100,-0.01,0.01);
	TH1F *errxz_PV=new TH1F("errxz_PV","errxz_PV",100,-0.01,0.01);
	TH1F *erryz_PV=new TH1F("erryz_PV","erryz_PV",100,-0.01,0.01);
	
	TH1F *errxy_SV=new TH1F("errxy_SV","errxy_SV",200,-1,1);
	TH1F *errxz_SV=new TH1F("errxz_SV","errxz_SV",200,-1,1);
	TH1F *erryz_SV=new TH1F("erryz_SV","erryz_SV",200,-1,1);
	
	//////TRUE FLAVOUR TAGGING////////////////////////////////////
	TH1F *b_pt=new TH1F("b_pt","b_pt",200,10,150);
	b_pt->SetLineColor(kBlue);
	TH1F *c_pt=new TH1F("c_pt","c_pt",200,10,150);
	c_pt->SetLineColor(kRed);
	TH1F *other_pt=new TH1F("other_pt","other_pt",200,10,150);
	other_pt->SetLineColor(kGreen+3);
	
	TH1F *num_jets=new TH1F("num_jets","num_jets",8,-0.5,7.5);
	TH1F *num_b_jets=new TH1F("num_b_jets","num_b_jets",8,-0.5,7.5);
	num_b_jets->SetLineColor(kRed);
	TH1F *num_c_jets=new TH1F("num_c_jets","num_c_jets",8,-0.5,7.5);
	num_c_jets->SetLineColor(kRed);
	
	long int ientry,k,i,n,m,l,j;
	int num_SV_perevent=0,num_SV_perjet=0,index_of_jet=-1,vtx_trk_vector[100]={0},vtx_trk_matrix[100][100]={0},wrong_PV=0,wrong_SV=0;
	int num_SV_B=0,num_SV_C=0,num_SV_BC=0,num_SV_noBC=0,num_trk_B=0,num_trk_C=0,num_trk_BC=0,v_i=0;
	int is_matched_with_vertex[100]={0},is_matched_with_jet[100]={-1},vtx_jet_matrix[100][100]={0};
	int jet_to_trk_matrix[100][500]={0},jet_to_trk_vector[100]={0};
	int j_b[10]={-1}, j_c[10]={-1}, j_other[10]={-1}, j_light[10]={-1}, j_index[10]={-1};
	int nj_b=0,nj_other=0,nj_c=0,nj_light=0,is_equal=0,nj_set=0;
	const double pT_cut=0.;
	double delta_R=0,delta_R2=0,delta_R_min=10000,corrected_mass=0,mom=0,theta=0,err_flight_dist,termX,termY,termZ,termXY,termXZ,termYZ;
	TLorentzVector vtx_momentum,reco,all_tracks_injet;
	TVector3 PV,SV,err_PV,err_SV,err_offd_PV,err_offd_SV,flight_dir,jet_mo,jet_axis,reco_v;
	TVector3  B_hadron, C_hadron;
	TVector3 START, END;
	
	int dumb=0,dumb1=0,dumb2=0;
	int num_of_trk_wrg_to_PV=0;
	int num_of_trk_right_to_SV=0;
	int trk_qg_to_fit_SV=0, trk_B_to_fit_SV=0, trk_C_to_fit_SV=0;
	
	//PERFORMANCE
	int prev_parent, same_parent;
	int tot_trk_b=0,tot_trk_c=0,tot_trk_other=0,tot_trk_ip=0;
	int trk_SV_b=0,trk_SV_c=0,trk_SV_other=0,trk_SV_ip=0;
	int trk_SV_same_parent_b=0,trk_SV_same_parent_c=0,trk_SV_same_parent_other=0;
	
	for(k=0;k<100;k++){ is_matched_with_jet[k]=-1;}
	for(k=0;k<10;k++){ j_b[k]=-1;j_c[k]=-1;j_other[k]=-1;j_light[k]=-1;}
	//for(k=0;k<100;k++){for(i=0;i<100;i++){cout<< jet_to_trk_matrix[k][i]<<endl;}}
	//fChain->GetEntries()
	
	//entry cycle
	for(ientry=0; ientry< fChain->GetEntries(); ++ientry){
  		fChain->GetEntry(ientry);
  		
  		
  	////////////////////////////////////////////////////////////////////////////////////////////
  	//////////////////////////// TRUE FLAVOUR TAG //////////////////////////////////////////////
  	////////////////////////////////////////////////////////////////////////////////////////////
  		// cycle on mc particle to find a B-hadron which has not a B-hadron daughter
  		for(k=0;k<nmcp;k++){
  			if(selection_bb(mcpdg[k])==0 && selection_bb(mcpdg[mcda0[k]])!=0 && selection_bb(mcpdg[mcda1[k]])!=0 && 			selection_bb(mcpdg[mcda2[k]])!=0 && selection_bb(mcpdg[mcda3[k]])!=0 && selection_bb(mcpdg[mcda4[k]])!=0){
  				cout<<"B"<<endl;
  				B_hadron.SetXYZ(mcmox[k],mcmoy[k],mcmoz[k]);
  				delta_R_min=10000;
  				//find nearest jet
  				for(i=0;i<nj;i++){
  					jet_mo.SetXYZ(jmox[i],jmoy[i],jmoz[i]);
					jet_axis=(jet_mo);
					
					delta_R=jet_axis.DeltaR(B_hadron);
					cout<<" index jet: "<<i<<" deltaR "<<delta_R<<endl;
					
					if(delta_R<delta_R_min && jet_mo.Pt()>pT_cut){
						delta_R_min=delta_R;
						index_of_jet=i;
						}
					}
				if(delta_R_min<1){
					is_equal=0;
					for(n=0;n<nj_b;n++){
						if(j_b[n]==index_of_jet){is_equal=1;}
						}
					if(is_equal==0){
						//jet is classified as a b-jet 
						j_b[nj_b]=index_of_jet;
						nj_b++;
						}
					}
				}
			}
		//cycle to identify jets not associated to B-hadrons and fill j_other
		for(i=0;i<nj;i++){
			is_equal=0;
			for(n=0;n<nj_b;n++){
				if(j_b[n]==i){is_equal=1;}
				}
			if(is_equal==0){
				j_other[nj_other]=i;
				nj_other++;
				}
			}
		// cycle on mc particle to find a C-hadron which has not a C-hadron daughter
		for(k=0;k<nmcp;k++){
			if(selection_cc(mcpdg[k])==0 && selection_cc(mcpdg[mcda0[k]])!=0 && selection_cc(mcpdg[mcda1[k]])!=0 && 			selection_cc(mcpdg[mcda2[k]])!=0 && selection_cc(mcpdg[mcda3[k]])!=0 && selection_cc(mcpdg[mcda4[k]])!=0){
				cout<<"C"<<endl;
				C_hadron.SetXYZ(mcmox[k],mcmoy[k],mcmoz[k]);
  				delta_R_min=10000;
  				//find nearest jet between not b-jets
  				for(i=0;i<nj_other;i++){
  					jet_mo.SetXYZ(jmox[j_other[i]],jmoy[j_other[i]],jmoz[j_other[i]]);
					jet_axis=(jet_mo);
					
					delta_R=jet_axis.DeltaR(C_hadron);
					cout<<" index jet: "<<j_other[i]<<" deltaR "<<delta_R<<endl;
					
					if(delta_R<delta_R_min && jet_mo.Pt()>pT_cut){
						delta_R_min=delta_R;
						index_of_jet=j_other[i];
						}
					}
				if(delta_R_min<1){
					is_equal=0;
					for(n=0;n<nj_c;n++){
						if(j_c[n]==index_of_jet){is_equal=1;}
						}
					if(is_equal==0){
						//jet is classified as a b-jet 
						j_c[nj_c]=index_of_jet;
						nj_c++;
						}
					}
				}
			}
			
		//cycle to identify jets not associated to B/C-hadrons and fill j_light
		for(i=0;i<nj_other;i++){
			jet_mo.SetXYZ(jmox[j_other[i]],jmoy[j_other[i]],jmoz[j_other[i]]);
			is_equal=0;
			for(n=0;n<nj_c;n++){
				if(j_c[n]==j_other[i]){is_equal=1;}
				}
			if(is_equal==0 && jet_mo.Pt()>pT_cut){
				j_light[nj_light]=j_other[i];
				nj_light++;
				}
			}
			
	cout<<"--------------------------------------------"<<endl;
	cout<<" In this event there are "<< nj<<" jets"<<endl;
	cout<<" In this event there are "<< nj_b<<" b-jets"<<endl;
	cout<< "Their indeces are: "<<endl;
	for(i=0;i<nj_b;i++){
		cout<<"index: "<<j_b[i]<<endl;
		}
	cout<<" In this event there are "<< nj_c<<" c-jets"<<endl;
	cout<< "Their indeces are: "<<endl;
	for(i=0;i<nj_c;i++){
		cout<<"index: "<<j_c[i]<<endl;
		}
	cout<<" In this event there are "<< nj_light<<" light-flavour-jets"<<endl;
	cout<< "Their indeces are: "<<endl;
	for(i=0;i<nj_light;i++){
		cout<<"index: "<<j_light[i]<<endl;
		}
		
	for(i=0;i<nj_light;i++){
		jet_mo.SetXYZ(jmox[j_light[i]],jmoy[j_light[i]],jmoz[j_light[i]]);
		other_pt->Fill(jet_mo.Pt());
		}
	for(i=0;i<nj_c;i++){
		jet_mo.SetXYZ(jmox[j_c[i]],jmoy[j_c[i]],jmoz[j_c[i]]);
		c_pt->Fill(jet_mo.Pt());
		}
	for(i=0;i<nj_b;i++){
		jet_mo.SetXYZ(jmox[j_b[i]],jmoy[j_b[i]],jmoz[j_b[i]]);
		b_pt->Fill(jet_mo.Pt());
		}	
	num_b_jets->Fill(nj_b);
	num_c_jets->Fill(nj_c);
	num_jets->Fill(nj);
	cout<<"END EVENT"<<endl;
	
	cout<<"--------------------------------------------"<<endl;
	
	//if the sample is bb take just b jets
	//if it is cc just the c jets
	
	if(choice[0]=='b'){
		nj_set=nj_b;
		for(i=0;i<nj_b;i++){
			j_index[i]=j_b[i];
			}
		}
		
	if(choice[0]=='c'){
		nj_set=nj_c;
		for(i=0;i<nj_c;i++){
			j_index[i]=j_c[i];
			}
		}
	if(choice[0]=='o'){
		nj_set=nj_light;
		for(i=0;i<nj_light;i++){
			j_index[i]=j_light[i];
			}
		}
	
		
  	////////////////////////////////////////////////////////////////////////////////////////////
  	////////////////////////////////////////////////////////////////////////////////////////////
  	
  		//cycle on vertices to find the PV (id=103) of the event and count SVs(id=104)
  		 for (k=0;k<nvt; k++){
  		 	if(id[k]==103){
  		 		PV.SetXYZ(x[k],y[k],z[k]);
  		 		err_PV.SetXYZ(sqrt(vtcov[k][0]),sqrt(vtcov[k][2]),sqrt(vtcov[k][5]));
  		 		err_offd_PV.SetXYZ(vtcov[k][1],vtcov[k][3],vtcov[k][4]);
  		 		//cout<<"PV =  X: "<<x[k]<<" Y: "<<y[k]<<" Z: "<<z[k]<<endl;
  		 		break;
  		 		}			 	
  		 	}
  		 
  		 //cycle on vertices to count SVs(id=104), evaluate flight distance, find nearest jet axis and plot deltaR	
  		 for (k=0;k<nvt; k++){
  		 	if(id[k]==104){
  		 		num_SV_perevent++;
  		 		//SV flight direction
	  	 		SV.SetXYZ(x[k],y[k],z[k]);
	  	 		err_SV.SetXYZ(sqrt(vtcov[k][0]),sqrt(vtcov[k][2]),sqrt(vtcov[k][5]));
	  	 		err_offd_SV.SetXYZ(vtcov[k][1],vtcov[k][3],vtcov[k][4]);
	  	 		flight_dir=SV-PV;
	  	 		//error on PV is negligible
	  	 		termX=pow(SV.X()-PV.X(),2)*pow(err_SV.X(),2);
	  	 		termY=pow(SV.Y()-PV.Y(),2)*pow(err_SV.Y(),2);
	  	 		termZ=pow(SV.Z()-PV.Z(),2)*pow(err_SV.Z(),2);
	  	 		termXY=2*(SV.X()-PV.X())*(SV.Y()-PV.Y())*err_offd_SV.X();
	  	 		termXZ=2*(SV.X()-PV.X())*(SV.Z()-PV.Z())*err_offd_SV.Y();
	  	 		termYZ=2*(SV.Y()-PV.Y())*(SV.Z()-PV.Z())*err_offd_SV.Z();
	  	 		err_flight_dist=sqrt(termX+termY+termZ+termXY+termXZ+termYZ)/flight_dir.Mag();
  		 		//fill histo with SV flight distance//////////////////////////////////////////////
  		 		SV_to_PV_distance->Fill(flight_dir.Mag());///////////// PLOT 3 ///////////////////
  		 		SV_to_PV_distance_sig->Fill(flight_dir.Mag()/err_flight_dist);////////////////////
  		 		//////////////////////////////////////////////////////////////////////////////////
  		 		
  		 		//find the nearest jet axis and plot deltaR
  		 		delta_R_min=10000;
  		 		for(i=0;i<nj_set;i++){
					//jet axis
					jet_mo.SetXYZ(jmox[j_index[i]],jmoy[j_index[i]],jmoz[j_index[i]]);
					jet_axis=(jet_mo);
					delta_R=jet_axis.DeltaR(flight_dir);
					if(delta_R<delta_R_min){
						delta_R_min=delta_R;
						index_of_jet=j_index[i];
						}
					}
				if(delta_R_min<1.){
					dumb++;
					is_matched_with_jet[k]=index_of_jet;
					//it contains the indeces of the verteces related to the jets
					vtx_jet_matrix[index_of_jet][is_matched_with_vertex[index_of_jet]]=k;
					is_matched_with_vertex[index_of_jet]++;
					if(index_of_jet>nj){cout<<"Number on jet in the event: "<<nj<<" index of jet:" <<index_of_jet<<endl;}
					}
				
				//fill histo with deltaR min//////////////////////////////////////////////////////
  		 		deltaR_SV_nearestJet->Fill(delta_R_min);///////////// PLOT 4//////////////////////
  		 		//////////////////////////////////////////////////////////////////////////////////
  		 		}
  		 	}//vertex cycle
  		 	
  		 	//TEST
				
				for(l=0;l<nj_set;l++){for(j=0;j<is_matched_with_vertex[j_index[l]];j++){dumb2++;}}
				
  		 	
  		for(i=0;i<nj_set;i++){	
	  		//fill histo with number of SV per jet//////////////////////////////////////////////////////////
	  	 	SV_reco_perjet->Fill(is_matched_with_vertex[j_index[i]]);//////////////////// PLOT 2 /////////////////////
	  		////////////////////////////////////////////////////////////////////////////////////////////////
  			if(is_matched_with_vertex[j_index[i]]==0){        ///////////////////////////////////////////////////////
  				jet_fraction_vertex_cat->Fill(1);///////////////////////////////////////////
  				}				 ///////////////// PLOT 2b /////////////////////////////
  			else if(is_matched_with_vertex[j_index[i]]>0){	 ///////////////////////////////////////////////////////
  				jet_fraction_vertex_cat->Fill(0);///////////////////////////////////////////
  				}				 ///////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////
  		 	}
  		//fill histo with number of SV per event//////////////////////////////////////////////////////////	
  		SV_reco_perevent->Fill(num_SV_perevent);//////////////////// PLOT 1 //////////////////////////////
  		//////////////////////////////////////////////////////////////////////////////////////////////////
  		
  		num_SV_perevent=0;
  		
  		for (k=0;k<nvt; k++){
  		 	if(id[k]==104){
  		 		//SV flight direction
	  	 		SV.SetXYZ(x[k],y[k],z[k]);
	  	 		err_SV.SetXYZ(sqrt(vtcov[k][0]),sqrt(vtcov[k][2]),sqrt(vtcov[k][5]));
	  	 		err_offd_SV.SetXYZ(vtcov[k][1],vtcov[k][3],vtcov[k][4]);
	  	 		flight_dir=SV-PV;
	  	 		//error on PV is negligible
	  	 		termX=pow(SV.X()-PV.X(),2)*pow(err_SV.X(),2);
	  	 		termY=pow(SV.Y()-PV.Y(),2)*pow(err_SV.Y(),2);
	  	 		termZ=pow(SV.Z()-PV.Z(),2)*pow(err_SV.Z(),2);
	  	 		termXY=2*(SV.X()-PV.X())*(SV.Y()-PV.Y())*err_offd_SV.X();
	  	 		termXZ=2*(SV.X()-PV.X())*(SV.Z()-PV.Z())*err_offd_SV.Y();
	  	 		termYZ=2*(SV.Y()-PV.Y())*(SV.Z()-PV.Z())*err_offd_SV.Z();
	  	 		err_flight_dist=sqrt(termX+termY+termZ+termXY+termXZ+termYZ)/flight_dir.Mag();
  		 		if(is_matched_with_jet[k]!=-1){
  		 			//fill histo with SV flight distance///////////////////////////////////////////////
	  		 		SV_to_PV_distance_j->Fill(flight_dir.Mag());///////////// PLOT 3b /////////////////
	  		 		SV_to_PV_distance_sig_j->Fill(flight_dir.Mag()/err_flight_dist);//////////////////
	  		 		//////////////////////////////////////////////////////////////////////////////////
	  		 		}
	  		 	}
	  		 }
  		
  		//cycle on reco-particles
  		for(i=0;i<nrec;i++){
  		
  			//cycle on relations to find the related MC
  			for(n=0;n<r2mnrel;n++){
  				if(r2mf[n]==i){
  					break;
  					}
  				}
  			
  			//if the particle has a track associated to it
  			
  			if(rcftr[i]!=-1 && rccha[i]!=0){
  					
  			///////////////////////////////////////////////////////////////
  			//////////// PERFORMANCE //////////////////////////////////////
  			///////////////////////////////////////////////////////////////
  			//tot tracks from IP (start point less than 1e -3)
  			START.SetXYZ(mcvtx[r2mt[n]],mcvty[r2mt[n]],mcvtz[r2mt[n]]);
  			END.SetXYZ(mcvtx[mcpa0[r2mt[n]]],mcvty[mcpa0[r2mt[n]]],mcvtz[mcpa0[r2mt[n]]]);
  			//if(abs(mcpdg[mcpa0[r2mt[n]]])==5||abs(mcpdg[mcpa0[r2mt[n]]])==4 ||abs(mcpdg[mcpa0[r2mt[n]]])==3 || abs(mcpdg[mcpa0[r2mt[n]]])==2 || abs(mcpdg[mcpa0[r2mt[n]]])==1 || abs(mcpdg[mcpa0[r2mt[n]]])==21){
  			if(START.Mag()<1e-3){
  				tot_trk_ip++;
  				}
  			//tot tracks from non-zero-lifetime bottom 
  			
  			else if(selection_bb(mcpdg[mcpa0[r2mt[n]]])==0 && (START-END).Mag()>0){
  				tot_trk_b++;
  				}
  			//tot tracks from non-zero-lifetime charm
  			else if(selection_cc(mcpdg[mcpa0[r2mt[n]]])==0 && (START-END).Mag()>0){
  				tot_trk_c++;
  				}
  			//all other tracks
  			else if((START-END).Mag()>0){
  			//else if(abs(mcpdg[mcpa0[r2mt[n]]])==22 || abs(mcpdg[mcpa0[r2mt[n]]])==310 || abs(mcpdg[mcpa0[r2mt[n]]])==3122){
  			//if(abs(mcpdg[mcpa0[r2mt[n]]])==310 || abs(mcpdg[mcpa0[r2mt[n]]])==130 ||abs(mcpdg[mcpa0[r2mt[n]]])==321 || abs(mcpdg[mcpa0[r2mt[n]]])==311 ||abs(mcpdg[mcpa0[r2mt[n]]])== 2212 ||abs(mcpdg[mcpa0[r2mt[n]]])== 2112||abs(mcpdg[mcpa0[r2mt[n]]])== 3122 ||abs(mcpdg[mcpa0[r2mt[n]]])==3222 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3212 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3112 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3322 ||abs(mcpdg[mcpa0[r2mt[n]]])==3312 ||abs(mcpdg[mcpa0[r2mt[n]]])==3334||abs(mcpdg[mcpa0[r2mt[n]]])==15||abs(mcpdg[mcpa0[r2mt[n]]])==22 || abs(mcpdg[mcpa0[r2mt[n]]])==2212 || abs(mcpdg[mcpa0[r2mt[n]]])==2112 ||abs(mcpdg[mcpa0[r2mt[n]]])==211 ||abs(mcpdg[mcpa0[r2mt[n]]])==111){
  				tot_trk_other++;
  				}
  			//else{ cout<<"ID: "<<mcpdg[mcpa0[r2mt[n]]]<<endl;}
  			
  			}

  			//if it has a vertex
  			if(rcvts[i]!=-1){
  				//row index is the id number of the vertex, the content is the id number of the tracks related to that vertex
  				vtx_trk_matrix[rcvts[i]][vtx_trk_vector[rcvts[i]]]=i;
  				//array index is the id number of the vertex, the content is the number of tracks related to that vertex 
  				vtx_trk_vector[rcvts[i]]++;
  				
  					
  				  	///////////////////////////////////////////////////////////////
  					//////////// PERFORMANCE //////////////////////////////////////
  					///////////////////////////////////////////////////////////////
  					
  					//trk in SV from IP
  					START.SetXYZ(mcvtx[r2mt[n]],mcvty[r2mt[n]],mcvtz[r2mt[n]]);
  					END.SetXYZ(mcvtx[mcpa0[r2mt[n]]],mcvty[mcpa0[r2mt[n]]],mcvtz[mcpa0[r2mt[n]]]);
  					//if((abs(mcpdg[mcpa0[r2mt[n]]])==5||abs(mcpdg[mcpa0[r2mt[n]]])==4 ||abs(mcpdg[mcpa0[r2mt[n]]])==3 || abs(mcpdg[mcpa0[r2mt[n]]])==2 || abs(mcpdg[mcpa0[r2mt[n]]])==1 || abs(mcpdg[mcpa0[r2mt[n]]])==21)&&id[rcvts[i]]==104){	
  					if(START.Mag()<1e-3 && id[rcvts[i]]==104){
  						trk_SV_ip++;
  						}
  					//trk in SV from non-zero-lifetime bottom
  					else if(selection_bb(mcpdg[mcpa0[r2mt[n]]])==0 && (START-END).Mag()>0 && id[rcvts[i]]==104){
  							
  						trk_SV_b++;
  						}
  					//trk in SV from non-zero-lifetime charm
  					else if(selection_cc(mcpdg[mcpa0[r2mt[n]]])==0&&(START-END).Mag()>0 &&id[rcvts[i]]==104){	
  						trk_SV_c++;
  						}
  					//all other tracks
  					else if(id[rcvts[i]]==104 && (START-END).Mag()>0 ){
  					//else if((abs(mcpdg[mcpa0[r2mt[n]]])==22 || abs(mcpdg[mcpa0[r2mt[n]]])==310 || abs(mcpdg[mcpa0[r2mt[n]]])==3122) && id[rcvts[i]]==104 ){
  					//selection_bb(mcpdg[mcpa0[r2mt[n]]])!=0){
  					//else if((abs(mcpdg[mcpa0[r2mt[n]]])==310 || abs(mcpdg[mcpa0[r2mt[n]]])==130 ||abs(mcpdg[mcpa0[r2mt[n]]])==321 || abs(mcpdg[mcpa0[r2mt[n]]])==311 ||abs(mcpdg[mcpa0[r2mt[n]]])== 2212 ||abs(mcpdg[mcpa0[r2mt[n]]])== 2112||abs(mcpdg[mcpa0[r2mt[n]]])== 3122 ||abs(mcpdg[mcpa0[r2mt[n]]])==3222 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3212 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3112 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3322 ||abs(mcpdg[mcpa0[r2mt[n]]])==3312 ||abs(mcpdg[mcpa0[r2mt[n]]])==3334||abs(mcpdg[mcpa0[r2mt[n]]])==15||abs(mcpdg[mcpa0[r2mt[n]]])==22||abs(mcpdg[mcpa0[r2mt[n]]])==2212 || abs(mcpdg[mcpa0[r2mt[n]]])==2112 || abs(mcpdg[mcpa0[r2mt[n]]])==211 ||abs(mcpdg[mcpa0[r2mt[n]]])==111)&&id[rcvts[i]]==104){
  						trk_SV_other++;
  						}
  					////////////////////////////////////////////////////////////////////////////////////////
  					////////////////////////////////////////////////////////////////////////////////////////
  				//if(r2mw[n]>0.9){

  					//tracks from B/C wrongly associated to PV
  					if((selection_cc(mcpdg[mcpa0[r2mt[n]]])==0||selection_bb(mcpdg[mcpa0[r2mt[n]]])==0)
  					&& (START-END).Mag()>0. && id[rcvts[i]]==103){
  						wrong_PV++;
  						}
  					//tracks from B/C right associated to SV
  					else if((selection_cc(mcpdg[mcpa0[r2mt[n]]])==0||selection_bb(mcpdg[mcpa0[r2mt[n]]])==0)
  					&& (START-END).Mag()>0. && id[rcvts[i]]==104){	
  						num_of_trk_right_to_SV++;
  						}
  					//tracks from hadronization wrongly associated to SV
  					//if((abs(mcpdg[mcpa0[r2mt[n]]])==5||abs(mcpdg[mcpa0[r2mt[n]]])==4 ||abs(mcpdg[mcpa0[r2mt[n]]])==3 || abs(mcpdg[mcpa0[r2mt[n]]])==2 || abs(mcpdg[mcpa0[r2mt[n]]])==1 || abs(mcpdg[mcpa0[r2mt[n]]])==21) &&id[rcvts[i]]==104){
  					else if(START.Mag()<1e-3 && id[rcvts[i]]==104){
  						wrong_SV++;
  						}
  					//}
  				}
  		 	}//reco-p cycle
  		////////////////////////////////////////////////////////////////////////////////////////
  		from_hh_wrong_to_PV->Fill(wrong_PV);////////////// PLOT 9///////////////////////////////
  		////////////////////////////////////////////////////////////////////////////////////////
  		from_IP_wrong_to_SV->Fill(wrong_SV);////////////// PLOT 8 //////////////////////////////
  		////////////////////////////////////////////////////////////////////////////////////////
  		trk_qg_to_fit_SV=wrong_SV+trk_qg_to_fit_SV;
  		num_of_trk_wrg_to_PV=num_of_trk_wrg_to_PV+wrong_PV;
  		wrong_PV=0;
  		wrong_SV=0;
  		
  		for(k=0;k<nvt;k++){
  			if(id[k]==103){
  				////////////////////////////////////////////////////////////////////////////////////////
  				num_trk_forPV->Fill(vtx_trk_vector[k]); ////////////// PLOT 6 //////////////////////////
  				////////////////////////////////////////////////////////////////////////////////////////
  				chi2_per_num_trk_PV->Fill(chi[k]/vtx_trk_vector[k]); ////////////// PLOT 7 //////////////
  				////////////////////////////////////////////////////////////////////////////////////////
  				//addition
  				errx_PV->Fill(sqrt(vtcov[k][0]));
  				erry_PV->Fill(sqrt(vtcov[k][2]));
  				errz_PV->Fill(sqrt(vtcov[k][5]));
  				errxy_PV->Fill(vtcov[k][1]);
  				errxz_PV->Fill(vtcov[k][3]);
  				erryz_PV->Fill(vtcov[k][4]);
  				}
  			else if(id[k]==104){
  			
  				//addition
  				errx_SV->Fill(sqrt(vtcov[k][0]));
  				erry_SV->Fill(sqrt(vtcov[k][2]));
  				errz_SV->Fill(sqrt(vtcov[k][5]));
  				errxy_SV->Fill(vtcov[k][1]);
  				errxz_SV->Fill(vtcov[k][3]);
  				erryz_SV->Fill(vtcov[k][4]);
  				
  				num_trk_forSV->Fill(vtx_trk_vector[k]);
  				
  				//cycle on number of tracks associated to the vertex
  				for(n=0;n<vtx_trk_vector[k];n++){
  					reco.SetPxPyPzE(rcmox[vtx_trk_matrix[k][n]],rcmoy[vtx_trk_matrix[k][n]],rcmoz[vtx_trk_matrix[k][n]],rcene[vtx_trk_matrix[k][n]]);
  					vtx_momentum=vtx_momentum+reco;
  					}
  				//SV flight direction
	  	 		SV.SetXYZ(x[k],y[k],z[k]);
	  	 		flight_dir=SV-PV;
	  	 		mom=sqrt(pow(vtx_momentum.Px(),2)+pow(vtx_momentum.Py(),2)+pow(vtx_momentum.Pz(),2));
	  	 		theta=flight_dir.Angle(vtx_momentum.Vect());
	  	 		corrected_mass=sqrt(pow(vtx_momentum.M(),2)+pow(mom*sin(theta),2))+mom*sin(theta);
	  	 		////////////////////////////////////////////////////////////////////////////////////////
  				vtx_invariant_mass->Fill(vtx_momentum.M());/////////////// PLOT 5 //////////////////////
  				vtx_mom->Fill(mom); ////////////////////////////////////////////////////////////////////
  				vtx_theta->Fill(theta);  ///////////////////////////////////////////////////////////////
  				///////////////////////////////////////////////
  				////////////////////////////////////////////////////////////////////////////////////////
  				if(is_matched_with_jet[k]!=-1){
	  				////////////////////////////////////////////////////////////////////////////////////////
	  				vtx_invariant_mass_j->Fill(vtx_momentum.M());/////////////// PLOT 5b ///////////////////
	  				vtx_mom_j->Fill(mom); //////////////////////////////////////////////////////////////////
	  				vtx_theta_j->Fill(theta); //////////////////////////////////////////////////////////////
	  				vtx_corrected_mass->Fill(corrected_mass);///////////////////////////////////////////////
	  				////////////////////////////////////////////////////////////////////////////////////////
	  				chi2_per_num_trk_SV->Fill(chi[k]/vtx_trk_vector[k]);
	  				//addition
	  				/*cout<<" "<<endl;
	  				cout<<"sigma xx: "<<sqrt(vtcov[k][0])<<endl;
	  				cout<<"sigma xy: "<<vtcov[k][1]<<endl;
	  				cout<<"sigma yy: "<<sqrt(vtcov[k][2])<<endl;
	  				cout<<"sigma zx: "<<vtcov[k][3]<<endl;
	  				cout<<"sigma zy: "<<vtcov[k][4]<<endl;
	  				cout<<"sigma zz: "<<sqrt(vtcov[k][5])<<endl;
	  				cout<<"FINE"<<endl;*/
	  				//error on PV is negligible
		  	 		termX=pow(SV.X()-PV.X(),2)*pow(err_SV.X(),2);
		  	 		termY=pow(SV.Y()-PV.Y(),2)*pow(err_SV.Y(),2);
		  	 		termZ=pow(SV.Z()-PV.Z(),2)*pow(err_SV.Z(),2);
		  	 		termXY=2*(SV.X()-PV.X())*(SV.Y()-PV.Y())*err_offd_SV.X();
		  	 		termXZ=2*(SV.X()-PV.X())*(SV.Z()-PV.Z())*err_offd_SV.Y();
		  	 		termYZ=2*(SV.Y()-PV.Y())*(SV.Z()-PV.Z())*err_offd_SV.Z();
		  	 		err_flight_dist=sqrt(termX+termY+termZ+termXY+termXZ+termYZ)/flight_dir.Mag();
	  				add->Fill(theta,flight_dir.Mag()/err_flight_dist);
	  				}
  				
  				vtx_momentum.SetPxPyPzE(0.,0.,0.,0.);  
  				prev_parent=-1;
  				same_parent=1;
  				//cycle on number of tracks associated to the vertex
  				for(n=0;n<vtx_trk_vector[k];n++){
  					//cycle on relations to find the related MC
  					for(m=0;m<r2mnrel;m++){
  						if(r2mf[m]==vtx_trk_matrix[k][n]){
  							break;
  							}
  						}
  						
  					///PERFORMANCE
  					if(prev_parent==-1){
  						prev_parent=mcpa0[r2mt[m]];
  						}
  					else if(prev_parent!=mcpa0[r2mt[m]] && same_parent==1){
  						same_parent=0;
  						}	
  					/////////////
  					
  					
  					//if the parent is a B/C hadron
  					if(selection_cc(mcpdg[mcpa0[r2mt[m]]])==0){
  						num_trk_C++;	
  						}
  					else if(selection_bb(mcpdg[mcpa0[r2mt[m]]])==0){
  						num_trk_B++;
  						}
  					}
  					
  				///PERFORMANCE
  				if(same_parent==1){
  					START.SetXYZ(mcepx[prev_parent],mcepy[prev_parent],mcepz[prev_parent]);
  					END.SetXYZ(mcvtx[prev_parent],mcvty[prev_parent],mcvtz[prev_parent]);
  					////trk in SV from bottom requiring that all SV tracks come from the same parent
  					if(selection_bb(mcpdg[prev_parent])==0 && (START-END).Mag()>1e-7){
  						trk_SV_same_parent_b=trk_SV_same_parent_b+vtx_trk_vector[k];
  						}
  					////trk in SV from charm requiring that all SV tracks come from the same parent
  					else if(selection_cc(mcpdg[prev_parent])==0 && (START-END).Mag()>1e-7){
  						trk_SV_same_parent_c=trk_SV_same_parent_c+vtx_trk_vector[k];
  						}
  					else if(abs(mcpdg[mcpa0[r2mt[n]]])==310 || abs(mcpdg[mcpa0[r2mt[n]]])==130 ||abs(mcpdg[mcpa0[r2mt[n]]])==321 || abs(mcpdg[mcpa0[r2mt[n]]])==311 ||abs(mcpdg[mcpa0[r2mt[n]]])== 2212 ||abs(mcpdg[mcpa0[r2mt[n]]])== 2112||abs(mcpdg[mcpa0[r2mt[n]]])== 3122 ||abs(mcpdg[mcpa0[r2mt[n]]])==3222 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3212 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3112 ||abs(mcpdg[mcpa0[r2mt[n]]])== 3322 ||abs(mcpdg[mcpa0[r2mt[n]]])==3312 ||abs(mcpdg[mcpa0[r2mt[n]]])==3334||abs(mcpdg[mcpa0[r2mt[n]]])==15||abs(mcpdg[mcpa0[r2mt[n]]])==22){
  						trk_SV_same_parent_other=trk_SV_same_parent_other+vtx_trk_vector[k];
  						}
  					}
  				/////////////////////////
  						
  				trk_B_to_fit_SV=trk_B_to_fit_SV+num_trk_B;
  				trk_C_to_fit_SV=trk_C_to_fit_SV+num_trk_C;
  				//only jet-matched
  				if(is_matched_with_jet[k]!=-1){
  				if(num_trk_B>0&&num_trk_C==0){
  					num_SV_B++;
  					if(is_matched_with_jet[k]!=-1){
  						vtx_corrected_mass_j_B->Fill(corrected_mass);
  						SV_to_PV_distance_j_B->Fill(flight_dir.Mag());
  						SV_to_PV_distance_sig_j_B->Fill(flight_dir.Mag()/err_flight_dist);
  						}
  					}
  				if(num_trk_C>0&&num_trk_B==0){
  					num_SV_C++;
  					if(is_matched_with_jet[k]!=-1){
  						vtx_corrected_mass_j_C->Fill(corrected_mass);
  						SV_to_PV_distance_j_C->Fill(flight_dir.Mag());
  						SV_to_PV_distance_sig_j_C->Fill(flight_dir.Mag()/err_flight_dist);
  						}
  					}
  				if(num_trk_C>0 && num_trk_B>0){
  					num_SV_BC++;
  					if(is_matched_with_jet[k]!=-1){
  						vtx_corrected_mass_j_BC->Fill(corrected_mass);
  						SV_to_PV_distance_j_BC->Fill(flight_dir.Mag());
  						SV_to_PV_distance_sig_j_BC->Fill(flight_dir.Mag()/err_flight_dist);
  						}
  					}
  				if(num_trk_C==0 && num_trk_B==0){
  					num_SV_noBC++;
  					if(is_matched_with_jet[k]!=-1){
  						vtx_corrected_mass_j_noBC->Fill(corrected_mass);
  						SV_to_PV_distance_j_noBC->Fill(flight_dir.Mag());
  						SV_to_PV_distance_sig_j_noBC->Fill(flight_dir.Mag()/err_flight_dist);
  						}
  					}
  				num_trk_C=0;
  				num_trk_B=0;
  				}
  				}
  			}//vertex cycle
  		////////////////////////////////////////////////////////////////////////////////////
  		SV_reco_B->Fill(num_SV_B); /////////////////////////////////////////////////////////
  		SV_reco_C->Fill(num_SV_C); /////////////////////// PLOT 10 /////////////////////////
  		SV_reco_BC->Fill(num_SV_BC); ///////////////////////////////////////////////////////
  		SV_reco_noBC->Fill(num_SV_noBC);  //////////////////////////////////////////////////
  		////////////////////////////////////////////////////////////////////////////////////
  		num_SV_C=0;num_SV_B=0;num_SV_BC=0;num_SV_noBC=0;
  		
  		//PLOT 11: fist cycle on recoParticle to find the tracks associated to the jet and evaluate the total momentum carried by these
  		//second cycle on vertex associated to jet and evaluate the total energy of the tracks used to fit that vertex
  		for(i=0;i<nrec;i++){
  			if(rcftr[i]!=-1){
  				reco.SetPxPyPzE(rcmox[i],rcmoy[i],rcmoz[i],rcene[i]);
  				reco_v=reco.Vect();
  				delta_R_min=10000;
  				for(k=0; k<nj_set; k++){
  					jet_mo.SetXYZ(jmox[j_index[k]],jmoy[j_index[k]],jmoz[j_index[k]]);
					jet_axis=(jet_mo);
  					delta_R=jet_axis.DeltaR(reco_v);
					if(delta_R<delta_R_min){
						delta_R_min=delta_R;
						index_of_jet=j_index[k];
						}
					
					}
				if(delta_R_min<1.){
					jet_to_trk_matrix[index_of_jet][jet_to_trk_vector[index_of_jet]]=i;
					jet_to_trk_vector[index_of_jet]++;
					}
				}
			}

  		//jet cycle
  		for(k=0; k<nj_set; k++){
  			jet_mo.SetXYZ(jmox[j_index[k]],jmoy[j_index[k]],jmoz[j_index[k]]);
			jet_axis=(jet_mo);
  			if(is_matched_with_vertex[j_index[k]]>0){
  				for(i=0;i<jet_to_trk_vector[j_index[k]];i++){
  					reco.SetPxPyPzE(rcmox[jet_to_trk_matrix[j_index[k]][i]],rcmoy[jet_to_trk_matrix[j_index[k]][i]],rcmoz[jet_to_trk_matrix[j_index[k]][i]],rcene[jet_to_trk_matrix[j_index[k]][i]]);
  					all_tracks_injet=all_tracks_injet+reco;
			
					}
				for(n=0;n<is_matched_with_vertex[j_index[k]];n++){
					v_i=vtx_jet_matrix[j_index[k]][n];	
					for(m=0;m<vtx_trk_vector[v_i];m++){
		  				reco.SetPxPyPzE(rcmox[vtx_trk_matrix[v_i][m]],rcmoy[vtx_trk_matrix[v_i][m]],rcmoz[vtx_trk_matrix[v_i][m]],rcene[vtx_trk_matrix[v_i][m]]);
		  				vtx_momentum=vtx_momentum+reco;
		  				//TEST
		  				reco_v=reco.Vect();
						delta_R=jet_axis.DeltaR(reco_v);
						
		  				//if(delta_R>1){cout<<"THIS TRACK IS NOT ASSOCIATED TO THE JET"<<endl;}
		  				}
		  			
		  			dumb1++;
		  			SV_energy_ratio->Fill(vtx_momentum.E()/all_tracks_injet.E());
		  			vtx_momentum.SetPxPyPzE(0.,0.,0.,0.);
		  	 		}
  	 			all_tracks_injet.SetPxPyPzE(0.,0.,0.,0.);
  	 			}
  	 		}//jet cycle
			
					
  		
  		if(ientry>990){
  			cout<<"VERTICI SECONDARI RICOSTRUITI NELL'EVENTO: "<<ientry<<endl;
  			for(k=0;k<nvt;k++){
  				if(id[k]==104){
  				cout<<""<<endl;
  				cout<<"Vertice secondario con id number "<<k<<" di coordinate: X= "<<x[k]<<" Y: "<<y[k]<<" Z: "<<z[k]<<endl;
  				cout <<"Viene ricostruito con "<<vtx_trk_vector[k]<<" tracce: "<<endl;
  					for(n=0;n<vtx_trk_vector[k];n++){
  						i=0;
						while( r2mf[i]!=vtx_trk_matrix[k][n] ){i++;}
  						cout<< "Particella: "<<mcpdg[r2mt[i]]<<" che ha come parent: "<<mcpdg[mcpa0[r2mt[i]]]<<endl;
  						}
  					}
  				}
  			//TEST//
  			/*for(k=0; k<nrec;k++){
				//se la particella ha un vertice secondario associato
  				if(rcvts[k]!=-1&&id[rcvts[k]]==104){
  					cout<< "La particella con id number "<<k<<"proviene dal vertice (id number "<<rcvts[k]<<")  di tipo "<<id[rcvts[k]]<<"e di coordinate X: "<<x[rcvts[k]]<<" Y: "<<y[rcvts[k]]<<" Z: "<< z[rcvts[k]] <<endl;
  					cout<<" "<<endl;
  					}
  				}*/
  					
  			cout<<"FINE EVENTO"<<endl;
  			}
  			
  		for(k=0;k<100;k++){ vtx_trk_vector[k]=0; is_matched_with_jet[k]=-1; is_matched_with_vertex[k]=0;jet_to_trk_vector[k]=0;for(n=0;n<100;n++){vtx_trk_matrix[k][n]=-1;vtx_jet_matrix[k][n]=-1;jet_to_trk_matrix[k][n]=-1;}}
  		nj_b=0;nj_c=0;nj_other=0;nj_light=0;
  		for(k=0;k<10;k++){ j_b[k]=-1;j_c[k]=-1;j_other[k]=-1;j_light[k]=-1;j_index[k]=-1;}
  			
  		}//entry cycle
  	
	cout<<"END"<<endl;
	cout<<""<<endl;
	
	cout<<"------------------------------------------"<<endl;
	cout<<"  PERFORMANCE "<<endl;
	cout<<"tot_trk_ip: "<<tot_trk_ip<<endl;
	cout<<"tot_trk_b: "<<tot_trk_b<<endl;
	cout<<"tot_trk_c: "<<tot_trk_c<<endl;
	cout<<"tot_trk_other: "<<tot_trk_other<<endl;
	cout<<" "<<endl;
	cout<<"TRACKS IN SV :"<<endl;
	cout<<" Primary: "<< double(trk_SV_ip)/double(tot_trk_ip)*100<<endl;
	cout<<" Bottom: "<< double(trk_SV_b)/double(tot_trk_b)*100<<endl;
	cout<<" Charm: "<< double(trk_SV_c)/double(tot_trk_c)*100<<endl;
	cout<<" Other: "<< double(trk_SV_other)/double(tot_trk_other)*100<<endl;
	cout<<" "<<endl;
	cout<<"TRACKS IN SV from the same parent particle:"<<endl;
	cout<<" Bottom: "<< double(trk_SV_same_parent_b)/double(tot_trk_b)*100<<endl;
	cout<<" Charm: "<< double(trk_SV_same_parent_c)/double(tot_trk_c)*100<<endl;
	cout<<" Other: "<< double(trk_SV_same_parent_other)/double(tot_trk_other)*100<<endl;
	cout<<" ------------------------------------------ "<<endl;
	
  	
  	TCanvas *cdum=new TCanvas();
  	auto legend0 = new TLegend(0.63,0.75,0.97,0.96);
	legend0->AddEntry(b_pt, "b-jets" , "l");
	legend0->AddEntry(c_pt, "c-jets" , "l");
	legend0->AddEntry(other_pt, "light-flavour-jets" , "l");
	legend0->SetTextSize(0.03);
	if(choice[0]=='b'){
		b_pt->SetTitle("bb sample");
		b_pt->Draw("s hist");
  		c_pt->Draw("hist sames");
  		other_pt->Draw("hist sames");
		}
	else if(choice[0]=='c'){
		c_pt->SetTitle("cc sample");
		c_pt->Draw("s hist");
  		b_pt->Draw("hist sames");
  		other_pt->Draw("hist sames");
		}
	else if(choice[0]=='o'){
		other_pt->SetTitle("light jet sample");
		other_pt->Draw("s hist");
  		b_pt->Draw("hist sames");
  		c_pt->Draw("hist sames");
		}
  	
  	legend0->Draw();
  	
  	TCanvas *cjets=new TCanvas();
  	cjets->Divide(2,1);
  	cjets->cd(1);
  	num_b_jets->Draw("");
  	num_jets->Draw("same");
  	cjets->cd(2);
  	num_c_jets->Draw("");
  	num_jets->Draw("same");
  	
  		
  	TCanvas *c1=new TCanvas();
  	SV_reco_perevent->Scale(1./SV_reco_perevent->Integral());
  	SV_reco_perevent->Draw("s hist");	
  	
  	TCanvas *c2=new TCanvas();
  	SV_reco_perjet->Scale(1./SV_reco_perjet->Integral());
  	SV_reco_perjet->Draw("s hist");	
  	
  	TCanvas *c2b=new TCanvas();
  	jet_fraction_vertex_cat->Scale(1./jet_fraction_vertex_cat->Integral());
  	jet_fraction_vertex_cat->SetStats(0);
  	jet_fraction_vertex_cat->GetXaxis()->SetBinLabel(1,"RecoVertex");
  	jet_fraction_vertex_cat->GetXaxis()->SetBinLabel(2,"NoVertex");
  	jet_fraction_vertex_cat->GetYaxis()->SetRangeUser(0,1);
  	jet_fraction_vertex_cat->Draw("hist");	
  	
  	cout<<"reco vertex: "<<jet_fraction_vertex_cat->GetBinContent(1)<<endl;
  	cout<<"no vertex: "<<jet_fraction_vertex_cat->GetBinContent(2)<<endl;
  	
  	TCanvas *c3=new TCanvas();
  	c3->Divide(2,1);
  	c3->cd(1);
  	SV_to_PV_distance->Scale(1./SV_to_PV_distance->Integral());
  	SV_to_PV_distance->Draw("s hist");
  	c3->cd(2);
  	SV_to_PV_distance_sig->Scale(1./SV_to_PV_distance_sig->Integral());
  	SV_to_PV_distance_sig->Draw("s hist");	
  	
  	TCanvas *c3b=new TCanvas();
  	c3b->Divide(2,1);
  	c3b->cd(1);
  	SV_to_PV_distance_j->Scale(1./SV_to_PV_distance_j->Integral());
  	SV_to_PV_distance_j->Draw("s hist");
  	c3b->cd(2);
  	SV_to_PV_distance_sig_j->Scale(1./SV_to_PV_distance_sig_j->Integral());
  	SV_to_PV_distance_sig_j->Draw("s hist");
  	
  	TCanvas *c4=new TCanvas();
  	deltaR_SV_nearestJet->Scale(1./deltaR_SV_nearestJet->Integral());
  	deltaR_SV_nearestJet->Draw("s hist");
  	
  	TCanvas *c5=new TCanvas();
  	c5->Divide(2,1);
  	c5->cd(1);
  	num_trk_forPV->Scale(1./num_trk_forPV->Integral());
  	num_trk_forPV->Draw("s hist");
	c5->cd(2);
  	num_trk_forSV->Scale(1./num_trk_forSV->Integral());
  	num_trk_forSV->Draw("s hist");
  	
  	TCanvas *c7=new TCanvas();
  	c7->Divide(2,1);
  	c7->cd(1);
  	//chi2_VS_num_trk_PV->Draw("colz");
  	chi2_per_num_trk_PV->Scale(1./chi2_per_num_trk_PV->Integral());
  	chi2_per_num_trk_PV->Draw("s hist");
  	c7->cd(2);
  	chi2_per_num_trk_SV->Scale(1./chi2_per_num_trk_SV->Integral());
  	chi2_per_num_trk_SV->Draw("s hist");
  	//chi2_VS_num_trk_SV->Draw("colz");
  	
  	TCanvas *c8=new TCanvas();
  	SV_energy_ratio->Scale(1./SV_energy_ratio->Integral());
  	SV_energy_ratio->Draw("s hist");
  	
  	TCanvas *c9=new TCanvas();
  	c9->Divide(2,1);
  	c9->cd(1);
  	vtx_invariant_mass->Scale(1./vtx_invariant_mass->Integral());
  	vtx_invariant_mass->Draw("s hist");
  	c9->cd(2);
  	vtx_invariant_mass_j->Scale(1./vtx_invariant_mass_j->Integral());
  	vtx_invariant_mass_j->Draw("s hist");
  	
  	TCanvas *c10_0=new TCanvas();
  	vtx_corrected_mass->Scale(1./vtx_corrected_mass->Integral());
  	vtx_corrected_mass->Draw("s hist");
  	
  	auto legend1 = new TLegend(0.63,0.75,0.97,0.96);
	legend1->AddEntry(SV_reco_C, "#splitline{SV with at least a track}{from C decay}" , "l");
	legend1->AddEntry(SV_reco_noBC, "SV with no tracks from C decay" , "l");
	legend1->SetTextSize(0.03);
  	
  	auto legend2 = new TLegend(0.63,0.75,0.97,0.96);
	legend2->AddEntry(SV_reco_B, "#splitline{SV with at least a track}{from B decay (no tracks from C)}", "l");
	legend2->AddEntry(SV_reco_C, "#splitline{SV with at least a track}{from C decay (no tracks from B)}" , "l");
	legend2->AddEntry(SV_reco_BC, "#splitline{SV with at least a track}{from both B and C}" , "l");
	legend2->AddEntry(SV_reco_noBC, "SV with no tracks from B/C decay" , "l");
	legend2->SetTextSize(0.03);
	
	auto legend3 = new TLegend(0.63,0.75,0.97,0.96);
	legend3->AddEntry(SV_reco_B, "#splitline{SV with at least a track}{from B decay (no tracks from C)}", "l");
	legend3->AddEntry(SV_reco_C, "#splitline{SV with at least a track}{from C decay (no tracks from B)}" , "l");
	legend3->SetTextSize(0.03);
	
	auto legend4 = new TLegend(0.63,0.75,0.97,0.96);
	legend4->AddEntry(SV_reco_BC, "#splitline{SV with at least a track}{from both B and C}" , "l");
	legend4->AddEntry(SV_reco_noBC, "#splitline{SV with no tracks}{ from B/C decay}" , "l");
	legend4->SetTextSize(0.03);
	
	vtx_corrected_mass_j_BC->Scale(1./vtx_corrected_mass_j_BC->Integral());
  	vtx_corrected_mass_j_B->Scale(1./vtx_corrected_mass_j_B->Integral());
  	vtx_corrected_mass_j_C->Scale(1./vtx_corrected_mass_j_C->Integral());
  	vtx_corrected_mass_j_noBC->Scale(1./vtx_corrected_mass_j_noBC->Integral());
  	
  	SV_to_PV_distance_sig_j_BC->Scale(1./SV_to_PV_distance_sig_j_BC->Integral());
  	SV_to_PV_distance_sig_j_B->Scale(1./SV_to_PV_distance_sig_j_B->Integral());
  	SV_to_PV_distance_sig_j_C->Scale(1./SV_to_PV_distance_sig_j_C->Integral());
  	SV_to_PV_distance_sig_j_noBC->Scale(1./SV_to_PV_distance_sig_j_noBC->Integral());
  	
  	SV_to_PV_distance_j_BC->Scale(1./SV_to_PV_distance_j_BC->Integral());
  	SV_to_PV_distance_j_B->Scale(1./SV_to_PV_distance_j_B->Integral());
  	SV_to_PV_distance_j_C->Scale(1./SV_to_PV_distance_j_C->Integral());
  	SV_to_PV_distance_j_noBC->Scale(1./SV_to_PV_distance_j_noBC->Integral());
	
  	
  	
  	if(choice[0]=='b'){
	  	//for bb
	  	
	  	TCanvas *c10=new TCanvas();
	  	c10->Divide(2,1);
	  	c10->cd(1);
		
	  	vtx_corrected_mass_j_noBC->Draw("s hist");
	  	vtx_corrected_mass_j_BC->Draw("sames hist");
	  	legend4->Draw("");
	  	
	  	c10->cd(2);
	  	vtx_corrected_mass_j_C->Draw("s hist");
	  	vtx_corrected_mass_j_B->Draw("sames hist");
	  	legend3->Draw("");	
	  	
	  	TCanvas *c10_a=new TCanvas();
	  	c10_a->Divide(2,1);
	  	c10_a->cd(1);
		SV_to_PV_distance_sig_j_noBC->Draw("s hist");
	  	SV_to_PV_distance_sig_j_BC->Draw("sames hist");
	  	legend4->Draw("");

	  	c10_a->cd(2);
	  	SV_to_PV_distance_sig_j_B->Draw("s hist");
	  	SV_to_PV_distance_sig_j_C->Draw("sames hist");
		legend3->Draw("");
		
		TCanvas *c10_b=new TCanvas();
	  	c10_b->Divide(2,1);
	  	c10_b->cd(1);
	  	SV_to_PV_distance_j_noBC->Draw("s hist");
	  	SV_to_PV_distance_j_BC->Draw("sames hist");
	  	legend4->Draw("");
	  	
	  	c10_b->cd(2);
	  	SV_to_PV_distance_j_B->Draw("s hist");
	  	SV_to_PV_distance_j_C->Draw("sames hist");
	  	legend3->Draw("");
	  	
	  	TCanvas *c15=new TCanvas();
	  	SV_reco_B->SetStats(0);
	  	SV_reco_BC->SetStats(0);
	  	SV_reco_C->SetStats(0);
	  	SV_reco_noBC->SetStats(0);
		
	  	SV_reco_B->Scale(1./SV_reco_B->Integral());
	  	SV_reco_B->SetTitle("SV per event divided in categories (jet-matched only)");
	  	SV_reco_B->Draw("s hist");
	  	SV_reco_BC->Scale(1./SV_reco_BC->Integral());
	  	SV_reco_BC->Draw("sames hist");
	  	SV_reco_C->Scale(1./SV_reco_C->Integral());
	  	SV_reco_C->Draw("sames hist");
	  	SV_reco_noBC->Scale(1./SV_reco_noBC->Integral());
	  	SV_reco_noBC->Draw("sames hist");
	  	legend2->Draw();
  				
		}
	else if(choice[0]=='c'){
		//for cc
		TCanvas *c10_c=new TCanvas();
		c10_c->Divide(2,1);
	  	c10_c->cd(1);
	  	SV_to_PV_distance_j_noBC->Draw("s hist");
	  	SV_to_PV_distance_j_C->Draw("sames hist");
	  	legend1->Draw();
	  	
	  	c10_c->cd(2);
	  	SV_to_PV_distance_sig_j_noBC->Draw("s hist");
	  	SV_to_PV_distance_sig_j_C->Draw("sames hist");
		legend1->Draw();
		
		TCanvas *c10_d=new TCanvas();
	  	vtx_corrected_mass_j_C->Draw("s hist");
	  	vtx_corrected_mass_j_noBC->Draw("sames hist");
	  	legend1->Draw();
	  	
	  	TCanvas *c10_e=new TCanvas();
	  	SV_reco_C->SetStats(0);
	  	SV_reco_noBC->SetStats(0);
	  	SV_reco_noBC->SetTitle("SV per event divided in categories (jet-matched only)");
	  	SV_reco_noBC->Scale(1./SV_reco_noBC->Integral());
	  	SV_reco_noBC->Draw("s hist");
	  	SV_reco_C->Scale(1./SV_reco_C->Integral());
	  	SV_reco_C->Draw("sames hist");
	  	legend1->Draw();
	  	}
  	
  	TCanvas *c11=new TCanvas();
  	c11->Divide(2,1);
  	c11->cd(1);
  	vtx_theta->Scale(1./vtx_theta->Integral());
  	vtx_theta->Draw("s hist");
 	c11->cd(2);
  	vtx_theta_j->Scale(1./vtx_theta_j->Integral());
  	vtx_theta_j->Draw("s hist");
  	
  	TCanvas *c12=new TCanvas();
  	c12->Divide(2,1);
  	c12->cd(1);
  	vtx_mom->Scale(1./vtx_mom->Integral());
  	vtx_mom->Draw("s hist");
  	c12->cd(2);
  	vtx_mom_j->Scale(1./vtx_mom_j->Integral());
  	vtx_mom_j->Draw("s hist");
  	
  	TCanvas *c13=new TCanvas();
  	c13->Divide(2,1);
  	c13->cd(1);
  	from_hh_wrong_to_PV->Scale(1./from_hh_wrong_to_PV->Integral());
  	from_hh_wrong_to_PV->Draw("s hist");
  	c13->cd(2);
  	from_IP_wrong_to_SV->Scale(1./from_IP_wrong_to_SV->Integral());
  	from_IP_wrong_to_SV->Draw("s hist");
  	
  	
  	TCanvas *c16=new TCanvas();
  	add->Draw("colz");
  	
  	TCanvas *c17=new TCanvas();  
  	c17->Divide(3,2);
  	c17->cd(1);
  	errx_PV->Draw();
  	gPad->SetLogy() ;
  	c17->cd(2);
  	erry_PV->Draw();
  	gPad->SetLogy() ;
  	c17->cd(3);
  	errz_PV->Draw();
  	gPad->SetLogy() ;
  	c17->cd(4);
  	errx_SV->Draw();
  	gPad->SetLogy() ;
  	c17->cd(5);
  	erry_SV->Draw();
  	gPad->SetLogy() ;
  	c17->cd(6);
  	errz_SV->Draw();
  	gPad->SetLogy() ;
  	
  	TCanvas *c18=new TCanvas();  
  	c18->Divide(3,2);
  	c18->cd(1);
  	errxy_PV->Draw();
  	gPad->SetLogy() ;
  	c18->cd(2);
  	errxz_PV->Draw();
  	gPad->SetLogy() ;
  	c18->cd(3);
  	erryz_PV->Draw();
  	gPad->SetLogy() ;
  	c18->cd(4);
  	errxy_SV->Draw();
  	gPad->SetLogy() ;
  	c18->cd(5);
  	errxz_SV->Draw();
  	gPad->SetLogy() ;
  	c18->cd(6);
  	erryz_SV->Draw();
  	gPad->SetLogy() ;
  		
	free(mcvtx);free(mcvty);free(mcvtz);free(mcpdg);free(mcpa0);free(r2mt);free(r2mf);free(rcftr);free(r2mw);free(rcvts);
	free(rcmox);free(rcmoy);free(rcmoz);free(rcene);free(trk_z0);free(trk_d0);free(trk_atIP);	
  	
  	//save histo in a rootfile
  	TFile* f;
  	if(choice[0]=='c'){
  		f=new TFile("/home/paola/Scrivania/vertex_cc_1000.root","RECREATE");
  		}
  	else if(choice[0]=='b'){
  		f=new TFile("/home/paola/Scrivania/vertex_bb_1000.root","RECREATE");
  		}
  	else if(choice[0]=='o'){
  		f=new TFile("/home/paola/Scrivania/vertex_light_1000.root","RECREATE");
  		}
  	else{cout<<"Non valid option"<<endl;}
  	
	if ( f->IsOpen() ) {
		cout << "File opened successfully" << endl;
	  	SV_reco_perevent->Write();
	  	SV_reco_perjet->Write();
	  	jet_fraction_vertex_cat->Write();
	  	SV_to_PV_distance_j->Write();
	  	SV_to_PV_distance_sig_j->Write();
	  	deltaR_SV_nearestJet->Write();
	  	vtx_corrected_mass->Write();
	  	chi2_per_num_trk_SV->Write();
	  	from_hh_wrong_to_PV->Write();
	  	from_IP_wrong_to_SV->Write();
	  	num_trk_forPV->Write();
	  	num_trk_forSV->Write();
	  	}
	else{cout<<"Error "<<endl;}
	}
