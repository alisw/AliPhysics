#include "iostream.h"

void AliITSComparisonV1
(const char *foundfile = "itstracks.root", const char *truefile = "galice_tracks_0.root", Int_t evnum = 0) {

  	// Make sure that ALICE objects are loaded
	if (gClassTable->GetID("AliRun") < 0) {
		gROOT->LoadMacro("loadlibs.C");
		loadlibs();
	}  
	else {
		delete gAlice;
		gAlice=0;
	}
	
	// Load the tree with track data
	TFile *fileTrueTracks = new TFile(truefile);
	TTree *treeTrueTracks = (TTree*)fileTrueTracks->Get("Tracks");
	TFile *fileFoundTracks = new TFile(foundfile);
	TTree *treeFoundTracks = (TTree*)fileFoundTracks->Get(Form("TreeT%d", evnum));
	
	Int_t nTrue = (Int_t)treeTrueTracks->GetEntries();
	Int_t *cnt = new Int_t[nTrue];
	for (Int_t i = 0; i < nTrue; i++) cnt[i] = 0;
		
	Int_t label, tpc_ok, entry, pdg_code;
	Double_t px, py, pz, pt;
	treeTrueTracks->SetBranchAddress("label", &label);
	treeTrueTracks->SetBranchAddress("entry", &entry);
	treeTrueTracks->SetBranchAddress("tpc_ok", &tpc_ok);
	treeTrueTracks->SetBranchAddress("pdg_code", &pdg_code);
	treeTrueTracks->SetBranchAddress("px", &px);
	treeTrueTracks->SetBranchAddress("py", &py);
	treeTrueTracks->SetBranchAddress("pz", &pz);
	treeTrueTracks->SetBranchAddress("pt", &pt);
	
	AliITSIOTrack *iotrack = 0;
	treeFoundTracks->SetBranchAddress("ITStracks", &iotrack);
	
	const Int_t npos = 36, nneg = 31;
	Int_t pos[npos] = {2212, -11, -13, 211, 321, 3222, 213, 323, 10323, 3224, 
	                   2224, 2214, -1114, -3112, -3312, 3224, -3114, -3314, 411, 
	                   431, 413, 433, -15, 4232, 4222, 4322, 4422, 4412, 4432, 
	                   4224, 4214, 4324, 4424, 4414, 4434, 4444};
	Int_t neg[nneg] = {2212, 11, 13, -211, -321, 3112, -213, -323, -10323, 3114, 
	                   1114, -2224, -2214, 33112, -3222, 3114, 3314, 3334, -3224,
	                   -411, -431, -413, -433, 15, -4422, -4432, -4214, -4324, 
	                   -4424, -4434, -444};
			   
	// Evaluation tree definition 
	Int_t labITS, labTPC, signC;
	Double_t difpt, diflambda, difphi, Dz, Dr, Dtot, ptg;
	TTree *treeEvaluation = new TTree("Eval", "Parameters for tracking evaluation");
	treeEvaluation->Branch("labITS"   , &labITS   , "labITS/I"   );
	treeEvaluation->Branch("labTPC"   , &labTPC   , "labTPC/I"   );
	treeEvaluation->Branch("difpt"    , &difpt    , "difpt/D"    ); 
	treeEvaluation->Branch("diflambda", &diflambda, "diflambda/D");
	treeEvaluation->Branch("difphi"   , &difphi   , "difphi/D"   );
	treeEvaluation->Branch("Dz"       , &Dz       , "Dz/D"       );
	treeEvaluation->Branch("Dr"       , &Dr       , "Dr/D"       );
	treeEvaluation->Branch("Dtot"     , &Dtot     , "Dtot/D"     );
	treeEvaluation->Branch("ptg"      , &ptg      , "ptg/D"      );
	treeEvaluation->Branch("signC"    , &signC    , "signC/I"    );
	
	// Make comparison
	Double_t *var = 0;
	Bool_t isGood;
	Int_t nOK, trueEntry;
	Double_t found_px, found_py, found_pz, found_pt;
	Double_t found_tgl, found_lambda, found_phi, found_curv;
	Double_t true_lambda, true_phi, true_px, true_py, true_pz, true_pt;
	Double_t duepi = 2.*TMath::Pi();	 
	Bool_t ispos, isneg;
	Int_t countpos = 0, countneg = 0, found_charge;
	for (Int_t i = 0; i < treeFoundTracks->GetEntries(); i++) {
		treeFoundTracks->GetEntry(i);
		labITS = iotrack->GetLabel();
		labTPC = iotrack->GetTPCLabel();
		isGood = (labITS >= 0);
		nOK = treeTrueTracks->Draw("px:py:pz", Form("label==%d && tpc_ok==1", abs(labITS)), "goff");
		if (!nOK) {
		
//                        cerr << "ITS label not found among findable tracks:";
//			cerr << "   labITS = " << labITS;
//			cerr << "   labTPC = " << labTPC;
//			cerr << endl;
			
			isGood = kFALSE;
		}
		if (nOK > 1) {
			cerr << "More than 1 tracks with label " << labITS << ": taking the first" << endl;
		}
		true_px = *treeTrueTracks->GetV1();
		true_py = *treeTrueTracks->GetV2();
		true_pz = *treeTrueTracks->GetV3();
		true_pt = TMath::Sqrt(true_px * true_px + true_py * true_py);
		true_phi = TMath::ATan2(true_py, true_px);
		if(true_phi < 0) true_phi += duepi;
		true_phi *= 1000.0;
		true_lambda = TMath::ATan(true_pz / true_pt) * 1000.0;
		ptg = true_pt;
		
		// checks if two found good tracks have the same label in ITS
		treeTrueTracks->Draw("entry", Form("label==%d && tpc_ok==1", abs(labITS)), "goff");
		trueEntry = (Int_t)*treeTrueTracks->GetV1();
		if (isGood && cnt[trueEntry] == 0) 
			cnt[trueEntry] = 1;
		else if (isGood) {
			cout << "Entry " << trueEntry << " already analyzed!" << endl;
			continue;	
		}
				
		// charge matching
		found_charge = (Int_t)iotrack->GetCharge();
		ispos = isneg = kFALSE;
		for (Int_t j = 0; j < npos; j++) ispos = (pdg_code == pos[j]);
		for (Int_t j = 0; j < nneg; j++) isneg = (pdg_code == neg[j]);
		if (ispos) countpos++;
		if (isneg) countneg++;
		
		// pt resolution (%)
		found_px = iotrack->GetPx();
		found_py = iotrack->GetPy();
		found_pz = iotrack->GetPz();
		found_pt = TMath::Sqrt(found_px*found_px + found_py*found_py);
		difpt = ((found_pt - true_pt) / true_pt) * 100.;
		//cout << found_pt << " " << true_pt << " " << difpt << endl;
			
		// lambda (mrad)
		found_tgl = iotrack->GetStateTgl();
		found_lambda = TMath::ATan(found_tgl) * 1000.0;
		diflambda = found_lambda - true_lambda;
//		cout << "lambda " << found_lambda << " " << true_lambda << " " << diflambda << endl;
		
		// phi (mrad)
		found_phi = TMath::ACos(found_px / found_pt);
		if(found_phi > duepi) found_phi -= duepi;
		if(found_phi < 0.) found_phi += duepi;
		found_phi *= 1000.0;      
		difphi = found_phi - true_phi;
//		cout << "phi " << found_phi << " " << true_phi << " " << difphi << endl;
		
		// impact parameters (microns)
		Dr = iotrack->GetStateD() * 1.e4;
		Dz = iotrack->GetDz() * 1.e4;
		Dtot = sqrt(Dr*Dr + Dz*Dz);
		
		// curvature sign
		found_curv = iotrack->GetStateC();
		signC = (found_curv > 0.) ? 1 : -1;

		// fill tree
		treeEvaluation->Fill();
	}

	
	TFile *outFile = new TFile("AliITSComparisonV1.root", "recreate");
	treeEvaluation->Write("Eval", TObject::kOverwrite);
	outFile->Close();
}

