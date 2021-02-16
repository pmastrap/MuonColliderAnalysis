//struct for Jets, Vertices and Tracks

typedef struct {
	bool selected; //True if the jet passes the selection, False otherwise 
	char flavour_tag[2]; //true flavour tag 'b' or 'c'
	TLorentzVector j; //4-momentum of jet
	char vertex_category[20]; //'RecoVertex' or 'NoVertex'
	char lepton_category[20]; //'NoSoftLepton' or 'SoftMuon' or 'SoftElectron'
	int num_SV; //num of SV associated to the jet
	int SV_indices[50]; //indices of SVs associated to jet
	int num_trk; //num of trk associated to the jet
	int trk_indices[500]; //indices of trks associated to jet
	}Jet_struct;
	
typedef struct {
	bool selected; //True if the SV passes the selection, False otherwise 
	TVector3 v; //x-y-z coordinates of the vertex
	TVector3 errv;
	double fd_2; //flight distance
	double fd_3;
	double fd_sig_2; //flight dist sig 
	double fd_sig_3;
	int id; //103 for pv, 104 for sv, 105 for v0
	double chi2;
	int jet_index; //index of jet associated to SV
	int num_trk; //num of trk associated to the SV
	int trk_indices[500]; //indices of trks associated to SV
	}Vertex_struct;
	
typedef struct {
	bool selected; //True if the trk passes the selection, False otherwise 
	TLorentzVector t;
	double d0;
	double z0;
	double chi2;
	int ndf;
	double cov[15];
	int sub_det_hits[12];
	int SV_index; //index of SV associated to trk
	int jet_index; //index of jet associated to trk
	}Track_struct;
