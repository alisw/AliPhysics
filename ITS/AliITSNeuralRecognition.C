// ARGUMENTS:
// 1. number of azymuthal sectors (it's better not to go under 8 or over 40)
// 2. the ROOT file to read (WITHOUT exstension)
// 3. event number
// 4. if specified a string, a fstream named like the argument is opened and
//    the elapsed CPU time is stored (not useful)

// the macro will save a file named, for example "galice_<nsecs>.root"
// containing may AliITSneuralTrack objects

void AliITSNeuralRecognition
(Int_t nsecs = 20,
 const char* rfile = "its_recpoints_v1.root")
{
	TStopwatch timer;
	Double_t CONVERT = TMath::Pi() / 180.0;
	//cout << "Reading file " << rfile << " and saving in " << wfile << endl;
	
// ==================================
// ==== VERTEX READING ==============
// ==================================

	// If a "its_vertex.txt" file is provided, 
	// the vertex position is read from it.
	fstream f_vert("its_vertex.txt", ios::in);
	Double_t Vx, Vy, Vz, dummy;
	f_vert >> dummy >> Vx >> Vy >> Vz;
	if (fabs(Vx) < 0.05) Vx = 0.0;
	if (fabs(Vy) < 0.05) Vy = 0.0;
	cout << "Vertex position (x, y, z): " << Vx << ' ' << Vy << ' ' << Vz << endl;
	
	
// ==================================
// ==== CURVATURE CUT DEFINITION ====
// ==================================

	// These values define the curvature cuts for all steps
	// within a sector.
	// For a greater clarity, the cuts are indicated in units
	// of transverse momentum (GeV/c) but these value have no
	// exact physical meaning, but are useful to understand
	// well what means a choice in the value of a certain
	// curvature constraint
	// NOTE: be careful to make sure that the 'ncuts' variable
	//       have the same value of the dimension of the allocated arrays

	Int_t ncuts;
	Double_t *p, *cut;

	ncuts = 6;
	p = new Double_t[6];
	cut = new Double_t[6];
	cut[ 0] = 2.0 * 0.0003;
	cut[ 1] = 2.0 * 0.0006;
	cut[ 2] = 2.0 * 0.0009;
	cut[ 3] = 2.0 * 0.0010;
	cut[ 4] = 2.0 * 0.0012;
	cut[ 5] = 2.0 * 0.0015;


// ==========================
// ==== OTHER PARAMETERS ====
// ==========================

	Double_t helix_min[5]   = { 0.000, 0.000, 0.000, 0.00, 0.0 };
	Double_t helix_max[5]   = { 0.0215, 0.0215, 0.0206, 0.75, 0.1 };
	
	Double_t theta2D_min[5] = { 0.0, 0.0, 0.0, 0.0,  0.0 };
	Double_t theta2D_max[5] = { 1.0, 0.7, 0.8, 3.0, 30.0 };
	
	Double_t theta3D_min[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	Double_t theta3D_max[5] = { 1.2, 1.2, 2.0, 5.0, 5.0 };
	
	Double_t temp   = 1.0;     // temperature parameter
	Double_t var    = 0.00001; // stabilization threshold

	Double_t exp    = 20.0;    // straight-line excitator
	Double_t gtoc   = 6.0;     // gain/cost contribution ratio

	Double_t min    = 0.4;     // minimum in random activation initialization
	Double_t max    = 0.6;     // maximum in random activation initialization
	Double_t actmin = 0.55;    // activation threshold for binary map conversion


// =========================
// ==== NEURAL TRACKING ====
// =========================

	AliITSNeuralTracker *ANN = new AliITSNeuralTracker();
	
	ANN->SetVertex(Vx, Vy, Vz);

	ANN->SetCurvatureCuts(ncuts, cut);
	ANN->SetTemperature(temp);
	ANN->SetVariationLimit(var);
	ANN->SetGainToCostRatio(gtoc);
	ANN->SetWeightExponent(exp);
	ANN->SetInitInterval(min, max);
	ANN->SetActThreshold(actmin);

	ANN->SetPolarInterval(45.0);
	ANN->SetThetaCuts2D(theta2D_min, theta2D_max);
	ANN->SetThetaCuts3D(theta3D_min, theta3D_max);
	ANN->SetHelixMatchCuts(helix_min, helix_max);

	TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(rfile);
	if (!file) file = new TFile("its_recpoints_v1.root");
	TTree *points = 0;
	points = (TTree*)file->Get("TreeP");
	if (!points) {
		cerr << "No points found!!!" << endl;
		return;
	}

	cout << "Storing points..." << endl;
	ANN->CreateArrayStructure(nsecs);
	if (!ANN->ArrangePoints(points)) {
		cout << "Problems occurred while storing points. Aborted" << endl;
		return;
	}
	delete points;
	file.Close();

	cout << "Matching points..." << endl;
	ANN->StoreAbsoluteMatches();
	//ANN->PrintMatches(0);

	TCanvas *c = 0;//new TCanvas("c", "c", 0, 0, 500, 500);
	//c->Range(-50, -50, 50, 50);

	timer.Start();
	ANN->NeuralTracking("its_chains.root", c);
	timer.Stop();
	timer.Print();
	
	delete ANN;
}
