#if !defined(__CINT__) || defined(__MAKECINT__)
  #include "Riostream.h"
  #include "TFile.h"
  #include "TTree.h"
  #include "TStopwatch.h"

  #include "AliRun.h"
  #include "AliMagF.h"
  #include "AliRunLoader.h"
  #include "AliTPCLoader.h"
  #include "AliITSLoader.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerANN.h"
  #include "AliESD.h"
#endif

extern AliRun *gAlice;

Int_t AliITSFindTracksANN 
(Int_t nsectors = 10,      // number of azimutal sectors
 Int_t nev = 5)            // number of events to process
{

	// ==================================
	// ==== EVENT READING ===============
	// ==================================

	// Remove an already existing Run Loader
	if (gAlice) {
		delete AliRunLoader::GetRunLoader();
		delete gAlice; 
		gAlice=0;
	}
	
	// Instance the new Run Loader
	AliRunLoader* rl = AliRunLoader::Open("galice.root");
	if (rl == 0x0) {
		cerr<<"AliITSFindTracks.C : Can not open session RL=NULL"<< endl;
		return 3;
	}
	
	// Instance the ITS Loader
	AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
	if (itsl == 0x0) {
		cerr<<"AliITSFindTracksV2.C : Can not get ITS loader"<<endl;
		return 4;
	}
	
	// Load the gAlice
	if (rl->LoadgAlice()) {
		cerr<<"AliITSFindTracksV2.C : LoadgAlice returned error"<<endl;
		delete rl;
		return 3;
	}
	
	// Set NECESSARY conversion constant for magnetic field
	AliKalmanTrack::SetConvConst
	(
		1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField()
	);
	
	// Get the ITS data & geometry & recpoints
	AliITS *dITS = (AliITS*)rl->GetAliRun()->GetDetector("ITS");
	if (!dITS) {
		cerr<<"AliITSFindClusters.C : Can not find the ITS detector !"<<endl;
		return 6;
	}
	AliITSgeom *geom = dITS->GetITSgeom();
	itsl->LoadRecPoints("read");
	
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
	cut[ 0] = 0.0006;
	cut[ 1] = 0.0010;
	cut[ 2] = 0.0013;
	cut[ 3] = 0.0018;
	cut[ 4] = 0.0020;
	cut[ 5] = 0.0030;


	// ==========================
	// ==== OTHER PARAMETERS ====
	// ==========================

	
	Double_t helix[5]   = { 0.03, 0.02, 0.01, 0.03, 0.2 };
	Double_t theta2D[5] = { 1.5, 1.0, 1.0, 3.0, 10.0 };
	Double_t theta3D[5] = { 200., 200., 200., 200.0, 200.0 };
	
	
	/*
	Double_t helix[5]   = { 0.3, 0.3, 0.3, 0.5, 0.5 };
	Double_t theta2D[5] = { 5.0, 5.0, 5.0, 5.0, 5.0 };
	Double_t theta3D[5] = { 5.0, 5.0, 5.0, 5.0, 5.0 };
	*/
	
	Double_t temp   = 1.0;     // temperature parameter
	Double_t var    = 0.00001; // stabilization threshold

	Double_t exp    = 20.0;    // straight-line excitator
	Double_t gtoc   = 6.0;     // gain/cost contribution ratio

	Double_t min    = 0.4;     // minimum in random activation initialization
	Double_t max    = 0.6;     // maximum in random activation initialization
	Double_t actmin = 0.55;    // activation threshold for binary map conversion
	
	// Instance the Tracker
	AliITStrackerANN tracker(geom, 2);
	
	// Set cuts
	tracker.SetVertex(0.0, 0.0, 0.0);

	tracker.SetCuts(ncuts, cut, theta2D, theta3D, helix);
	tracker.SetTemperature(temp);
	tracker.SetVariationLimit(var);
	tracker.SetGainToCostRatio(gtoc);
	tracker.SetWeightExponent(exp);
	tracker.SetInitInterval(min, max);
	tracker.SetActThreshold(actmin);
	
	tracker.SetPolarInterval(45.0); 
		
	// Set number of events
	if (nev > rl->GetNumberOfEvents()) nev = rl->GetNumberOfEvents();
	Int_t rc = 0;
	
	// Get ESD files
	TFile *itsf=TFile::Open("AliESDits.root","RECREATE");
	if ((!itsf)||(!itsf->IsOpen())) {
		cerr << "Can't AliESDits.root !\n"; 
		return 1;
	}
	TFile *tpcf=TFile::Open("AliESDtpc.root");
	if ((!tpcf)||(!tpcf->IsOpen())) {
		cerr<<"Can't AliESDtpc.root !\n";
		return 1;
	}
	AliESD* event = new AliESD;
	TTree* esdTree = (TTree*) tpcf->Get("esdTree");
	if (!esdTree) {
	  cerr<<"no ESD tree found !\n";
	  return 1;
	}
	esdTree->SetBranchAddress("ESD", &event);
	
	// Loop on events
	TStopwatch timer; 
	for (Int_t i = 0; i < nev; i++) {
		cerr << "Processing event number: " << i << endl;
		esdTree->GetEvent(i);
		
		rl->GetEvent(i);
		
		// Get clusters from file
		TTree *cTree=itsl->TreeR();
		if (!cTree) {
			cerr<<"AliITSFindTracksANN.C : NO clusters tree!" << endl;
			return 4;
		}
		
		// Load clusters in tracker
		tracker.LoadClusters(cTree);
		
		// Create array structure and arrange points in it
		tracker.CreateArrayStructure(nsectors);
		tracker.ArrangePoints("its.recpoints.txt");
		tracker.StoreOverallMatches();
		//tracker.PrintMatches(kFALSE);
		
		// ***************************
		//  NEURAL TRACKING OPERATION
		// ***************************
		
		Bool_t isStable = kFALSE;
		Int_t i, nUnits = 0, nLinks = 0, nTracks = 0;
		Int_t sectTracks = 0, totTracks = 0;

		// tracking through sectors
		cout << endl;
		Int_t sector, curv;
//		nUnits = tracker.CreateNetwork(0, ncuts - 1);
//		cout << endl << nUnits << " NEURONS CREATED" << endl;
//		return;
		
		for (sector = 0; sector < nsectors; sector++) {
			sectTracks = 0;
			for(curv = 0; curv < ncuts; curv++) {
				// units creation
				nUnits = tracker.CreateNetwork(sector, curv);
				if (!nUnits) {
					cout << "no neurons --> skipped" << endl;
					continue;
				}
				// units connection
				nLinks = tracker.LinkNeurons();
				if (!nLinks) {
					cout << "no connections --> skipped" << endl;
					continue;
				}
				// neural updating
				for (i = 0;; i++) {
					isStable = tracker.Update();
					if (isStable) break;
				}
				// tracks saving
				tracker.CleanNetwork();
				nTracks = tracker.StoreTracks();
				cout << nUnits << " units, ";
				cout << nLinks << " connections, ";
				cout << i << " cycles --> ";
				cout << nTracks << " tracks" << endl;
				sectTracks += nTracks;
			}
			cout << "\n >>> Total tracks found in sector: " << sectTracks << endl;
			totTracks += sectTracks;
		}
		tracker.ExportTracks("ITS.Neural.root");
		cout << "\n\n--- Totally found " << totTracks << " tracks ---\n\n" << endl;
		
		// ***************************
		
		// End of operations: unload clusters
		tracker.UnloadClusters();
		
	}
	timer.Stop(); 
	timer.Print();
	
	// Close files & delete objects
	delete event;
	tpcf->Close();
	itsf->Close();
	delete rl;
	return rc;
}


