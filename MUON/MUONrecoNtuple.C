/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */



// Macro MUONrecoNtuple.C (TO BE COMPILED)
// for testing the C++ reconstruction code
// and producing an Ntuple with reconstructed tracks
// in output file "MUONtrackReco.root".
// An example for using the Ntuple is in the macro MUONmassPlot.C

// Arguments:
//   FirstEvent (default 0)
//   LastEvent (default 0)
//   RecGeantHits (1 to reconstruct GEANT hits) (default 0)
//   FileName (for signal) (default "galice.root")
//   BkgGeantFileName (for background),
//      needed only if RecGeantHits = 1 and background to be added

// IMPORTANT NOTICE FOR USERS:
// under "root" or "root.exe", execute the following commands:
// 1. "gSystem->SetIncludePath("-I$ALICE_ROOT/MUON -I$ALICE_ROOT/STEER -I$ROOTSYS/include")" to get the right path at compilation time
// 2. ".x loadlibs.C" to load the shared libraries
// 3. ".x MUONrecoNtuple.C+()" with the right arguments, without forgetting the "+" which implies the compilation of the macro before its execution

#include <iostream.h>

#include <TClassTable.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TTree.h>

#include "AliHeader.h"
#include "AliRun.h"

#include "AliMUON.h"
#include "AliMUONEventReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackHit.h"
#include "AliMUONTrackParam.h"

// Classes for Ntuple ////////////////////////////////////////////////////

class AliMUONTrackRecNtuple : public TObject {
 public:
  // for direct access to members
  Int_t fCharge; // charge of reconstructed track (+/- 1)
  Float_t fPxRec; // Px of reconstructed track at vertex (GeV/c)
  Float_t fPyRec; // Py of reconstructed track at vertex (GeV/c)
  Float_t fPzRec; // Pz of reconstructed track at vertex (GeV/c)
  Float_t fZRec; // Z of reconstructed track at vertex (cm)
  Float_t fZRec1; // Z of reconstructed track at first hit (cm)
  Int_t fNHits; // number of hits
  Float_t fChi2; // chi2 of fit
  Float_t fPxGen; // Px of best compatible generated track at vertex (GeV/c)
  Float_t fPyGen; // Py of best compatible generated track at vertex (GeV/c)
  Float_t fPzGen; // Pz of best compatible generated track at vertex (GeV/c)
  AliMUONTrackRecNtuple(){;} // Constructor
  virtual ~AliMUONTrackRecNtuple(){;} // Destructor
 protected:
 private:
  ClassDef(AliMUONTrackRecNtuple, 1) // AliMUONTrackRecNtuple
    };

class AliMUONHeaderRecNtuple : public TObject {
 public:
  // for direct access to members
  Int_t fEvent; // event number
  AliMUONHeaderRecNtuple(){;} // Constructor
  virtual ~AliMUONHeaderRecNtuple(){;} // Destructor
 protected:
 private:
  ClassDef(AliMUONHeaderRecNtuple, 1) // AliMUONHeaderRecNtuple
    };

ClassImp(AliMUONTrackRecNtuple) // Class implementation in ROOT context
ClassImp(AliMUONHeaderRecNtuple) // Class implementation in ROOT context

  //__________________________________________________________________________
void AliMUONEventRecNtupleFill(AliMUONEventReconstructor *Reco, Int_t FillWrite = 0)
{
  // Fill Ntuple for reconstructed event pointed to by "Reco"
  // if "FillWrite" different from -1.
  // Ntuple is created automatically before filling first event.
  // If "FillWrite" = -1, write and close the file

  static Bool_t firstTime = kTRUE;
  static TTree *ntuple;
  static AliMUONHeaderRecNtuple *header = new AliMUONHeaderRecNtuple();
  static TClonesArray *recTracks = new TClonesArray("AliMUONTrackRecNtuple",5);

  Int_t trackIndex;
  AliMUONTrack *track;
  AliMUONTrackParam *trackParam;
  AliMUONTrackRecNtuple *recTrackNt;
  Double_t bendingSlope, nonBendingSlope, pYZ;

  if (FillWrite == -1) {
    printf(">>> Writing Ntuple of reconstructed tracks\n");
    // better to create the file before the Ntuple ????
    TFile *file = new TFile("MUONtrackReco.root","recreate");
    ntuple->Write();
    file->Write();
    file->Close();
  }

  if (firstTime) {
    firstTime = kFALSE;
    // first call: create tree for Ntuple...
    printf(">>> Creating Ntuple of reconstructed tracks\n");
    ntuple = new TTree("MUONtrackReco", "MUONtrackReco");
    ntuple->Branch("Header","AliMUONHeaderRecNtuple", &header);
    ntuple->Branch("Tracks", &recTracks);
  }

  // header
  
  header->fEvent = gAlice->GetHeader()->GetEvent();

  TClonesArray *recoTracksPtr = Reco->GetRecTracksPtr();
  recoTracksPtr->Compress(); // for simple loop without "Next" since no hole
  recTracks->Clear(); // to reset the TClonesArray of tracks to be put in the ntuple
  // Loop over reconstructed tracks
  for (trackIndex = 0; trackIndex < Reco->GetNRecTracks(); trackIndex++) {
    track = (AliMUONTrack*) ((*recoTracksPtr)[trackIndex]);
    recTrackNt = (AliMUONTrackRecNtuple*)
      new ((*recTracks)[trackIndex]) AliMUONTrackRecNtuple();
    // track parameters at Vertex
    trackParam = track->GetTrackParamAtVertex();
    recTrackNt->fCharge =
      Int_t(TMath::Sign(1., trackParam->GetInverseBendingMomentum()));
    bendingSlope = trackParam->GetBendingSlope();
    nonBendingSlope = trackParam->GetNonBendingSlope();
    pYZ = 1/TMath::Abs(trackParam->GetInverseBendingMomentum());
    recTrackNt->fPzRec = pYZ / TMath::Sqrt(1.0 + bendingSlope * bendingSlope);
    recTrackNt->fPxRec = recTrackNt->fPzRec * nonBendingSlope;
    recTrackNt->fPyRec = recTrackNt->fPzRec * bendingSlope;
    recTrackNt->fZRec = trackParam->GetZ();
    // track parameters at first hit
    trackParam = ((AliMUONTrackHit*)
		  (track->GetTrackHitsPtr()->First()))->GetTrackParam();
    recTrackNt->fZRec1 = trackParam->GetZ();
    // chi2
    recTrackNt->fChi2 = track->GetFitFMin();
    // number of hits
    recTrackNt->fNHits = track->GetNTrackHits();
    printf("test> Px %f Py %f Pz %f \n", recTrackNt->fPxRec,  recTrackNt->fPyRec,  recTrackNt->fPzRec);
    // track parameters at vertex of best compatible generated track:
    // in fact muon with the right charge
    TTree* mtreeK=gAlice->TreeK();
    TBranch *brparticle = mtreeK->GetBranch("Particles");
    Int_t nPart = brparticle->GetEntries();
    TParticle *particle = new TParticle();
    mtreeK->SetBranchAddress("Particles",&particle);
    for (Int_t iPart = 0; iPart < nPart; iPart++) {
       brparticle->GetEntry(iPart);
       //cout << "Code Particle: " << particle->GetPdgCode() << "\n";
       if ((particle->GetPdgCode() * recTrackNt->fCharge) == -13) {
 	recTrackNt->fPxGen = particle->Px();
 	recTrackNt->fPyGen = particle->Py();
	recTrackNt->fPzGen = particle->Pz();
	printf("Gen: Px %f Py %f Pz %f \n", recTrackNt->fPxGen, recTrackNt->fPyGen, recTrackNt->fPzGen);
       }  
   }
  } // for (trackIndex = 0;...

  printf(">>> Filling Ntuple of reconstructed tracks\n");
  ntuple->Fill();

  return;
}

void MUONrecoNtuple (Int_t FirstEvent = 0, Int_t LastEvent = 0, Int_t RecGeantHits = 0, Text_t *FileName = "galice.root", Text_t *BkgGeantFileName = "")
{
  //
  cout << "MUON_recoNtuple" << endl;
  cout << "FirstEvent " << FirstEvent << endl;
  cout << "LastEvent " << LastEvent << endl;
  cout << "RecGeantHits " << RecGeantHits << endl;
  cout << "FileName ``" << FileName << "''" << endl;
  cout << "BkgGeantFileName ``" << BkgGeantFileName << "''" << endl;
//   // Dynamically link some shared libs                    
//   if (gClassTable->GetID("AliRun") < 0) {
//     gROOT->LoadMacro("loadlibs.C");
//     loadlibs();
//   }


  // Creating Run Loader and openning file containing Hits, Digits and RecPoints
  AliRunLoader * RunLoader = AliRunLoader::Open(FileName,"Event","UPDATE");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",FileName);
    return;
  }
  // Loading AliRun master
  RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();
  RunLoader->LoadKinematics("READ");
  
  // Loading MUON subsystem
  AliMUON * MUON = (AliMUON *) gAlice->GetDetector("MUON");
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadHits("READ");
  MUONLoader->LoadRecPoints("READ");

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();

  // Initializations
  // AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  MUON->SetTreeAddress(); 
  AliMUONEventReconstructor *Reco = new AliMUONEventReconstructor();

  Reco->SetRecGeantHits(RecGeantHits);

  // The right place for changing AliMUONEventReconstructor parameters
  // with respect to the default ones
//   Reco->SetMaxSigma2Distance(100.0);
//  Reco->SetPrintLevel(20);
   Reco->SetPrintLevel(1);
//   Reco->SetBendingResolution(0.0);
//   Reco->SetNonBendingResolution(0.0);
  cout << "AliMUONEventReconstructor: actual parameters" << endl;
  Reco->Dump();
//   gObjectTable->Print();
  // Loop over events
  if (LastEvent>nevents) LastEvent = nevents;
  for (Int_t event = FirstEvent; event < LastEvent; event++) {
    cout << "Event: " << event << endl;
    RunLoader->GetEvent(event);   
    //     Int_t nparticles = gAlice->GetEvent(event);
    //      cout << "nparticles: " << nparticles << endl;
    // prepare background file and/or event if necessary
    if (RecGeantHits == 1) {
      if (event == FirstEvent) Reco->SetBkgGeantFile(BkgGeantFileName);
      if (Reco->GetBkgGeantFile())Reco->NextBkgGeantEvent();
    }
    Reco->EventReconstruct();
    // Dump current event
    Reco->EventDump();
    // Fill Ntuple
    AliMUONEventRecNtupleFill(Reco, 0);
//     gObjectTable->Print();
    MUON->ResetRawClusters();
  } // Event loop
    // Write Ntuple
    AliMUONEventRecNtupleFill(Reco, -1);
}
