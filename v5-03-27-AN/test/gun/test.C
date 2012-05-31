// Usage in compiled mode
// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include");  
// gROOT->LoadMacro("test.C+");
// test()

#if !defined(__CINT__) || defined(__MAKECINT__)

// Root include files
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TParticle.h>

// AliRoot include files
#include "AliESDEvent.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"

#endif

void test(const char * sdir =".") {

  TStopwatch timer;
  timer.Start();

  TString name;

  // Signal file and tree
  name = sdir;
  name += "/AliESDs.root";
  TFile * fSig = TFile::Open(name.Data());
  TTree * tSig = (TTree*)fSig->Get("esdTree");

  AliESDEvent * esdSig = new AliESDEvent();// The signal ESD object is put here
  esdSig->ReadFromTree(tSig);

  // Run loader (signal events)
  name = sdir;
  name += "/galice.root";
  AliRunLoader* rlSig = AliRunLoader::Open(name.Data());

  // gAlice
  rlSig->LoadgAlice();
  gAlice = rlSig->GetAliRun();

  // Now load kinematics and event header
  rlSig->LoadKinematics();
  rlSig->LoadHeader();

  // Loop on events: check that MC and data contain the same number of events
  Long64_t nevSig = rlSig->GetNumberOfEvents();

  cout << nevSig << " signal events" << endl;

  Int_t lab[3]; // Labels from TOF
  Double_t mom[3]; // Track momentum

  for (Int_t iev=0; iev<nevSig; iev++) {
    cout << "---------- Signal event ----------" << iev << endl;

    // Get signal ESD
    tSig->GetEntry(iev);

    // Particle stack
    rlSig->GetEvent(iev);
    AliStack * stackSig = rlSig->Stack();
    stackSig->DumpPStack();
    Int_t nPartSig = stackSig->GetNtrack();

    Int_t nrec = esdSig->GetNumberOfTracks();
    cout << nrec << " reconstructed tracks" << endl;
    for(Int_t irec=0; irec<nrec; irec++) {
      AliESDtrack * track = esdSig->GetTrack(irec);
      cout << "Labels:" << endl;
      cout << "Global: "<< track->GetLabel() << endl;
      cout << "ITS: "<< track->GetITSLabel() << endl;
      cout << "TPC: "<< track->GetTPCLabel() << endl;
      cout << "TRD: "<< track->GetTRDLabel() << endl;
      track->GetTOFLabel(lab);
      cout << "TOF: "<< lab[0] <<" "<< lab[1] <<" "<< lab[2] << endl;
      UInt_t label = TMath::Abs(track->GetLabel());
      if (label>=10000000) {
	// Underlying event. 10000000 is the
	// value of fkMASKSTEP in AliRunDigitizer
	cout <<"Strange, there should be no underlying event"<<endl;
      }
      else {
	if (label>=nPartSig) {
	  cout <<"Strange, label outside the range"<< endl;
	  continue;
	}
	TParticle * part = stackSig->Particle(label);
	if(part) part->Print();
	track->GetPxPyPz(mom);
	cout <<"Momentum: "<< mom[0] <<" "<< mom[1] <<" "<< mom[2] <<endl; 
      }

    }

  }

  fSig->Close();

  timer.Stop();
  timer.Print();
}
