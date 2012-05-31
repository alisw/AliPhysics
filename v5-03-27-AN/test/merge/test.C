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

void test(const char * sdir ="signal",
	  const char * bdir ="backgr") {

  TStopwatch timer;
  timer.Start();

  TString name;

  // Signal file, tree, and branch
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

  // Run loader (underlying events)
  name = bdir;
  name += "/galice.root";
  AliRunLoader* rlUnd = AliRunLoader::Open(name.Data(),"Underlying");

  // gAlice
  rlSig->LoadgAlice();
  rlUnd->LoadgAlice();
  gAlice = rlSig->GetAliRun();

  // Now load kinematics and event header
  rlSig->LoadKinematics();
  rlSig->LoadHeader();
  rlUnd->LoadKinematics();
  rlUnd->LoadHeader();

  // Loop on events: check that MC and data contain the same number of events
  Long64_t nevSig = rlSig->GetNumberOfEvents();
  Long64_t nevUnd = rlUnd->GetNumberOfEvents();
  Long64_t nSigPerUnd = nevSig/nevUnd;

  cout << nevSig << " signal events" << endl;
  cout << nevUnd << " underlying events" << endl;
  cout << nSigPerUnd << " signal events per one underlying" << endl;

  for (Int_t iev=0; iev<nevSig; iev++) {
    cout << "Signal event " << iev << endl;
    Int_t ievUnd = iev/nSigPerUnd;
    cout << "Underlying event " << ievUnd << endl;

    // Get signal ESD
    tSig->GetEntry(iev);
    // Get signal kinematics
    rlSig->GetEvent(iev);
    // Get underlying kinematics
    rlUnd->GetEvent(ievUnd);

    // Particle stack
    AliStack * stackSig = rlSig->Stack();
    Int_t nPartSig = stackSig->GetNtrack();
    AliStack * stackUnd = rlUnd->Stack();
    Int_t nPartUnd = stackUnd->GetNtrack();

    Int_t nrec = esdSig->GetNumberOfTracks();
    cout << nrec << " reconstructed tracks" << endl;
    for(Int_t irec=0; irec<nrec; irec++) {
      AliESDtrack * track = esdSig->GetTrack(irec);
      UInt_t label = TMath::Abs(track->GetTPCLabel());
      if (label>=10000000) {
	// Underlying event. 10000000 is the
	// value of fkMASKSTEP in AliRunDigitizer
// 	cout << " Track from the underlying event" << endl;
	label %=10000000;
	if (label>=nPartUnd) continue;
	TParticle * part = stackUnd->Particle(label);
 	if(part) part->Print();
      }
      else {
	cout << " Track " << label << " from the signal event" << endl;
	if (label>=nPartSig) {
	  cout <<"Strange, label outside the range "<< endl;
	  continue;
	}
	TParticle * part = stackSig->Particle(label);
	if(part) part->Print();
      }

    }

  }

  fSig->Close();

  timer.Stop();
  timer.Print();
}
