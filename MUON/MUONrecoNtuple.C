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

#include "AliRun.h"

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
    // better to create the file before the Ntuple ????
    TFile *file = new TFile("MUONtrackReco.root","recreate");
    ntuple->Write();
    file->Write();
    file->Close();
  }

  if (firstTime) {
    firstTime = kFALSE;
    // first call: create tree for Ntuple...
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
    // track parameters at vertex of best compatible generated track:
    // in fact muon with the right charge
    for (int iPart = 0; iPart < gAlice->Particles()->GetEntriesFast(); iPart++) {
      TParticle *particle = (TParticle*) gAlice->Particles()->UncheckedAt(iPart);
      if ((particle->GetPdgCode() * recTrackNt->fCharge) == -13) {
	recTrackNt->fPxGen = particle->Px();
	recTrackNt->fPyGen = particle->Py();
	recTrackNt->fPzGen = particle->Pz();
      }
    }
  } // for (trackIndex = 0;...

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

  // Connect the Root Galice file containing Geometry, Kine, Hits
  // and eventually RawClusters
  TFile *file = (TFile*) gROOT->GetListOfFiles()->FindObject(FileName);
  if (!file) {
    printf("\n Creating file %s\n", FileName);
    file = new TFile(FileName);
  }
  else printf("\n File %s found in file list\n", FileName);

  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*) file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) {
      printf("\n Create new gAlice object");
      gAlice = new AliRun("gAlice","Alice test program");
    }
  }

  // Initializations
  // AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  AliMUONEventReconstructor *Reco = new AliMUONEventReconstructor();

  Reco->SetRecGeantHits(RecGeantHits);

  // The right place for changing AliMUONEventReconstructor parameters
  // with respect to the default ones
//   Reco->SetMaxSigma2Distance(100.0);
  Reco->SetPrintLevel(10);
//   Reco->SetPrintLevel(1);
//   Reco->SetBendingResolution(0.0);
//   Reco->SetNonBendingResolution(0.0);
  cout << "AliMUONEventReconstructor: actual parameters" << endl;
  Reco->Dump();
//   gObjectTable->Print();

  // Loop over events
  for (Int_t event = FirstEvent; event <= LastEvent; event++) {
    cout << "Event: " << event << endl;
//     AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
    Int_t nparticles = gAlice->GetEvent(event);
    cout << "nparticles: " << nparticles << endl;
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

  } // Event loop
    // Write Ntuple
    AliMUONEventRecNtupleFill(Reco, -1);
}
