// ARGUMENTS:
// 1. number of azymuthal sectors (it's better not to go under 8 or over 40)
// 2. the ROOT file to read (WITHOUT exstension)
// 3. event number
// 4. if specified a string, a fstream named like the argument is opened and
//    the elapsed CPU time is stored (not useful)

// the macro will save a file named, for example "galice_<nsecs>.root"
// containing may AliITSneuralTrack objects

void AliITSNeuralTracking
(Int_t nsecs = 12, const char* rfile = "galice", Int_t event = 0, const char* save = 0)
{
	TStopwatch timer;
	Double_t CONVERT = TMath::Pi() / 180.0;
	const char* wfile = Form("%s_%d.root", rfile, nsecs);
	cout << "Reading file " << rfile << ".root and saving in " << wfile << endl;
	
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
	// NOTE: becareful to make sure that the 'ncuts' variable
	//       have the same value of the dimension of the allocated arrays
	
	Int_t ncuts;
	Double_t *p, *cut;
	
	ncuts = 5;
	p = new Double_t[5];
	cut = new Double_t[5];	
	p[0] = 2.0;
	p[1] = 1.0;
	p[2] = 0.7;
	p[3] = 0.5;
	p[4] = 0.3;
	
	for (Int_t i = 0; i < ncuts; i++) cut[i] = 0.003 * 0.2 / p[i];
	

// ==========================
// ==== OTHER PARAMETERS ====
// ==========================
	    
	Bool_t   flag   = kFALSE; // for now, don't change this line, please...
	
	Double_t diff   = 0.02;   // helicoidal cut 
	Double_t dtheta = 1.0;    // delta-theta cut
	Double_t temp   = 1.0;    // temperature parameter
	Double_t var    = 0.0001; // stabilization threshold
	
	Double_t exp    = 7.0;    // straight-line excitator
	Double_t gtoc   = 3.0;    // gain/cost contribution ratio
	
	Double_t min    = 0.4;    // minimum in random activation initialization
	Double_t max    = 0.6;    // maximum in random activation initialization
	
	
// =========================
// ==== NEURAL TRACKING ====
// =========================
	
	AliITSneuralTracker *ANN = new AliITSneuralTracker(nsecs, ncuts, cut, CONVERT*dtheta);
		
	TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("%s.root", rfile));
	if (!file) file = new TFile(Form("%s.root", rfile),"UPDATE");

	//Double_t Xv = -0.001097;
	//Double_t Yv = -0.00347647;
	//Double_t Zv =  0.000631345;
	//ANN->SetVertex(Xv, Yv, Zv);
	// You should find the vertex with VertexMacro.C
	// and then put by hand the found values with
	// the above method.
	
	Int_t points = ANN->ReadFile(Form("%s.root", rfile), event);
		
	ANN->SetTemperature(temp);
	ANN->SetVariationLimit(var);
	ANN->SetGainToCostRatio(gtoc);
	ANN->SetExponent(exp);
	ANN->SetInitLimits(min, max);
	ANN->SetDiff(diff);
	
	cout << points << " points found " << endl << endl;
	
	TStopwatch timer;
	timer.Start();
	
	ANN->Go(wfile, flag);
	
	timer.Stop();
	cout << endl;
	timer.Print();
	
	if (save) {
		fstream ksave(save, ios::app);
		ksave << nsecs << " " << timer->CpuTime() << endl;
		ksave.close();
	}
	
//	delete gAlice;
//	gAlice = 0;
}
