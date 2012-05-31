#include <fstream.h>

#include <TClassTable.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"

Int_t AliITSStoreFindableTracksCompiled
(Int_t nMinClusters = 5, const Text_t *evname = "galice", Int_t evnum = 0)
{
	// Make sure that ALICE objects are loaded
	if (gAlice) {
		delete gAlice;
		gAlice = 0;
	}
	
	// Define the names of all involved files
	TString strEventFile(evname);
	TString strOutputFile(evname);
	strEventFile.Append(".root");
	strOutputFile.Append("_tracks_");
	strOutputFile += evnum;
	strOutputFile.Append(".root");
		
	// Connect the Root Galice file containing Geometry, Kine and Hits
	TFile *fileEvent = (TFile*)gROOT->GetListOfFiles()->FindObject(strEventFile);
	if (!fileEvent) fileEvent = new TFile(strEventFile,"UPDATE");
	
	// Get AliRun object from file
	gAlice = (AliRun*)fileEvent->Get("gAlice");
	if (gAlice) cout << "OK, found an AliRun object in file" << endl;
	
	// Get ITS related objects and data
	AliITS* ITS =(AliITS *)gAlice->GetDetector("ITS");
	if (!ITS) {
		cerr << "ITS object not found!" << endl;
		return 1;
	}
	AliITSgeom *geometry = ITS->GetITSgeom(); 
	if (!geometry) {
		cerr << "ITS geometry object not found!" << endl;
		return 2;
	}
	
	// Count the number of modules per layer
	Int_t nLadders, nDetectors, mod_min[6], mod_max[6];
	for(Int_t i = 0; i < 6; i++) {
		nLadders = geometry->GetNladders(i + 1);
		nDetectors = geometry->GetNdetectors(i + 1);
		mod_min[i] = geometry->GetModuleIndex(i + 1, 1, 1);
		mod_max[i] = geometry->GetModuleIndex(i + 1, nLadders, nDetectors);
	}
	
	// Load event and ITS recpoints
	Int_t nParticles = gAlice->GetEvent(evnum);
	cout << "Event number: " << evnum << endl;
	cout << "# particles : " << nParticles <<endl;
	if (nParticles <= 0) {
		cerr << "Can't have <= 0 particles!" << endl;
		return 3; 
	}
	AliITSRecPoint *recp = 0;
	TClonesArray  *recPoints = ITS->RecPoints();
	TObjArray     *particles = gAlice->Particles();
	Int_t nTracks = gAlice->GetNtrack(); //FCA correction
	Bool_t *hitITSLayer[6];
	for (Int_t i = 0; i < 6; i++) {
		hitITSLayer[i] = new Bool_t[nTracks];
		for (Int_t j = 0; j < nTracks; j++) hitITSLayer[i][j] = kFALSE;
	}
	
	// Load recpoints in event
	TTree *TR = gAlice->TreeR();
	if (!TR) {
		cerr << "TreeR object not found!" << endl;
		return 4;
	}
	
	// Scan recpoints and define findable tracks
	Int_t nModules = (Int_t)TR->GetEntries(), nPoints = 0, nEmpty = 0;
	cout << "Found " << nModules;
	cout << " entries in the TreeR (must be one per module!)" << endl;
	for (Int_t layer = 1; layer <= 6; layer++) {
		for (Int_t mod = mod_min[layer - 1]; mod <= mod_max[layer - 1]; mod++) {
			ITS->ResetRecPoints();
			TR->GetEntry(mod);
			nPoints = recPoints->GetEntries();
			if(!nPoints) {
				nEmpty++;
				continue;
			}
			for (Int_t point = 0; point < nPoints; point++) {
				recp = (AliITSRecPoint*)recPoints->UncheckedAt(point);
				for (Int_t it = 0; it < 3; it++) {
					Int_t track = recp->GetLabel(it);
					if(track < 0) continue;
					if(track > nTracks) {
						cout << "Found track index " << track;
						cout << " whilw gAlice->GetNtrack() = " << nTracks << endl;
						continue;
					}
					hitITSLayer[layer - 1][track] = kTRUE;
				} // loop over recpoint labels
			} //loop over points
		} //loop over modules
	} //loop over layers
	cout << "Found " << nEmpty << " empty modules" << endl;
	
	// Scan the file of tracks in TPC to retrieve the findable TPC tracks
	TString strLabelsTPC;
	Int_t label, pdg_code, nFindablesTPC = 0;
	Double_t dummy;
	ifstream tpc("good_tracks_tpc");
	while (tpc >> label >> pdg_code) {
		for (Int_t i = 0; i < 6; i++) tpc >> dummy;
		nFindablesTPC++;
		strLabelsTPC.Append(Form("[%d]", label));		
	}
			
	// Define the TTree with tracks data by means of a set of variables
	Int_t nFindablesITS = 0, nFindablesITSTPC = 0;
	Int_t nhits, tpc_ok, mother, entry = 0;
	Double_t vx, vy, vz;
	Double_t px, py, pz, pt;

	TTree *tree = new TTree("Tracks", "Findable tracks in ITS");
	
	tree->Branch("vx", &vx, "vx/D");
	tree->Branch("vy", &vy, "vy/D");
	tree->Branch("vz", &vz, "vz/D");
	tree->Branch("px", &px, "px/D");
	tree->Branch("py", &py, "py/D");
	tree->Branch("pz", &pz, "pz/D");
	tree->Branch("pt", &pt, "pt/D");
	tree->Branch("label", &label, "label/I");
	tree->Branch("entry", &entry, "entry/I");
	tree->Branch("mother", &mother, "mother/I");
	tree->Branch("pdg_code", &pdg_code, "pdg_code/I");
	tree->Branch("nhits", &nhits, "nhits/I");
	tree->Branch("tpc_ok", &tpc_ok, "tpc_ok/I");
	
	// Fill the tree
	cout << endl;
	TParticle *p = 0;
	for (Int_t i = 0; i < nTracks; i++) {
		nhits = 0;
		for (Int_t j = 0; j < 6; j++) if (hitITSLayer[j][i]) nhits++;
		if (nhits < nMinClusters) continue;
		p = gAlice->Particle(i);
		px = p->Px();
		py = p->Py();
		pz = p->Pz();
		pt = p->Pt();
		vx = p->Vx();
		vy = p->Vy();
		vz = p->Vz();
		mother = p->GetFirstMother();
		cout << "Track " << i << " stored\r" << flush;
		tpc_ok = (strLabelsTPC.Contains(Form("[%d]", i)));
		pdg_code = p->GetPdgCode();
		label = i;
		nFindablesITS++;
		if (tpc_ok) nFindablesITSTPC++;
		tree->Fill();
		entry++;
	}
	
	// Save into a file
	TFile *fileOutput = new TFile(strOutputFile, "recreate");
	tree->Write();
	fileOutput->Close();
	
	cout << "# findable tracks in TPC     : " << nFindablesTPC << endl;
	cout << "# findable tracks in ITS     : " << nFindablesITS << endl;
	cout << "# findable tracks in ITS+TPC : " << nFindablesITSTPC << endl;
	
	return 0;
}
