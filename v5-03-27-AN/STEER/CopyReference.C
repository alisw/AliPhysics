//
// This macro translate original Tracks Reference tree 
// into new tree sorted according to track label.
//
// TClonesArray connected to a branch, is filled
// with Track References from one track only.
//
// Algorithm:
// 
// 1) Open files and initialise trees
// 2) Create and initialise Look Up Table
// 3) Loop over entries in References tree
//    - store entry number in a LUT
//    - store first and last array index of references for one track 
// 4) Check LUT's consistancy.
// 5) Fill up new tree with references ordered acording to label.
//
//
// Sylwester Radomski (S.Radomski@gsi.de)
// GSI, Jan 20, 2003
//

void CopyReference () {

  // create output file
  const char *outFileName = "trackRef.root";
  TFile *outFile = new TFile(outFileName, "RECREATE");

  // connect to gAlice object in input file
  const char *inFileName = "galice.root";
  TFile *inFile = new TFile(inFileName, "READ");
  
  gAlice = (AliRun*)inFile->Get("gAlice");
  gAlice->GetEvent(0);

  // declare detectors to translate
  const Int_t nDet = 4;
  const char *detectors[nDet] = {"ITS","TPC", "TRD", "TOF"};

  // Create LUT's
  
  Int_t nPart = gAlice->GetNtrack();
  Int_t *lutTrack = new Int_t[nPart * nDet];
  Int_t *lutStart = new Int_t[nPart * nDet];
  Int_t *lutStop = new Int_t[nPart * nDet];

  // clean LUT
  
  for(Int_t i=0; i<nPart; i++)
	 for (Int_t j=0; j<nDet; j++)
		lutTrack[i*nDet+j] = lutStart[i*nDet+j] = lutStop[i*nDet+j] = -1;
  

  // Declare tree, branches and CloneArrays

  TTree *refTree = (TTree*)inFile->Get("TreeTR0");
  TBranch *refBranch[nDet];
  TClonesArray *refArray[nDet];

  // Connect branches with arrays

  for (Int_t j=0; j<nDet; j++) {
	 refBranch[j] = refTree->GetBranch(detectors[j]);
	 refArray[j] = new TClonesArray("AliTrackReference", 100);
	 refBranch[j]->SetAddress(&(refArray+j));
  }

  Int_t nTracks = refBranch[0]->GetEntries();

  cout << "N Particles\t" << nPart << endl
		 << "N Tracks\t" << nTracks << endl;

  cout << "Filling Look Up Tables ..." << endl;

  for (Int_t i=0; i<nTracks; i++) {

	 cout << i << "\r"; 

	 for (Int_t j=0; j<nDet; j++) {
		
		Int_t lab, start;
		AliTrackReference *ref = 0;
		refBranch[j]->GetEvent(i);
		Int_t nObj = refArray[j]->GetEntries();
		
		if (!nObj) continue;
		
		lab = ((AliTrackReference*)refArray[j]->At(0))->GetTrack();
		start = 0;
		
		for (Int_t e=0; e<nObj-1; e++) {
		  
		  ref = (AliTrackReference*)refArray[j]->At(e+1);
		  
		  if (ref->GetTrack() != lab) {
			 lutTrack[lab*nDet + j] = i;
			 lutStart[lab*nDet + j] = start;
			 lutStop[lab*nDet + j] = e+1;
			 
			 start = e+1;
			 lab = ref->GetTrack(); 		  
		  }
		}
		
		lutTrack[lab*nDet + j] = i;
		lutStart[lab*nDet + j] = start;
		lutStop[lab*nDet + j] = nObj;
	 }
  }
  
  cout << endl << "done" << endl;
  cout << "Checking consistancy of LUTs ..." << endl;

  // check consistancy
  
  for(Int_t i=0; i<nPart; i++) {
	 
	 cout << i << "\r";
	 Int_t ctrlTrack = -1;
	 
	 for(Int_t j=0; j<nDet; j++) { 

		if (ctrlTrack == -1 && lutTrack[i*nDet+j] != -1) 
		  ctrlTrack = lutTrack[i*nDet+j];

		if (lutTrack[i*nDet+j] != -1 && lutTrack[i*nDet+j] != ctrlTrack)
		  cout << "Error: " <<  i << " " << j << endl;
	 }
  }
  
  cout << endl <<  "done" << endl;
  cout << "Writing to a new Tree ..." << endl;

  // Create a Tree

  outFile->cd();
  TTree *outTree = new TTree("TreeTR0_Sorted", "Sorted Tracks References");
  TBranch *outBranch[nDet];
  TClonesArray *outArray[nDet];

  // Create Branches 

  for(Int_t j=0; j<nDet; j++) {
	 outArray[j] = new TClonesArray("AliTrackReference", 200);
	 outBranch[j] = outTree->Branch(detectors[j], &(outArray+j), 64000, 3);
  }

  // Loop over particles
  // Fill Tree if entries exists

  for(Int_t i=0; i<nPart; i++) {
	 
	 cout << i << "\r"; 

	 for(Int_t j=0; j<nDet; j++) {
		
		outArray[j]->Clear();

		if (lutTrack[i*nDet+j] != -1) {
		  
		  refBranch[j]->GetEvent(lutTrack[i*nDet+j]);

		  // rewrite objects in clonearrays
		  for(Int_t k = lutStart[i*nDet+j], en = 0; k<lutStop[i*nDet+j]; k++ ) {		  
			 AliTrackReference *ref = (AliTrackReference*)refArray[j]->At(k);
			 new( (*(outArray[j]))[en++] ) AliTrackReference(*ref);
		  }
		}
	 }

	 outTree->Fill();
  }

  cout << endl << "done" << endl;
  
  // close and clean
  
  outTree->Write();
  inFile->Close();
  outFile->Close();
  
  //if (refTree) delete refTree;
  //if (outTree) delete outTree;  
  //if (inFile) delete inFile;
  //if (outFile) delete outFile;
}
