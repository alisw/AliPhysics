/* $Id$ */

#include "AliMultiplicityMCSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TParticle.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliMultiplicity.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "AliPWG0Helper.h"
#include "dNdEta/AliMultiplicityCorrection.h"

ClassImp(AliMultiplicityMCSelector)

AliMultiplicityMCSelector::AliMultiplicityMCSelector() :
  AliSelectorRL(),
  fMultiplicity(0),
  fEsdTrackCuts(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliMultiplicityMCSelector::~AliMultiplicityMCSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliMultiplicityMCSelector::Begin(TTree* tree)
{
  // Begin function

  AliSelectorRL::Begin(tree);

  ReadUserObjects(tree);
}

void AliMultiplicityMCSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");
}

void AliMultiplicityMCSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  ReadUserObjects(tree);

  fMultiplicity = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
}

void AliMultiplicityMCSelector::Init(TTree* tree)
{
  // read the user objects

  AliSelectorRL::Init(tree);

  // enable only the needed branches
  if (tree)
  {
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fTriggerMask", 1);
    tree->SetBranchStatus("fSPDVertex*", 1);
    tree->SetBranchStatus("fSPDMult*", 1);

    //AliESDtrackCuts::EnableNeededBranches(tree);
  }
}

Bool_t AliMultiplicityMCSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fTree is the pointer to the TChain being processed,
  //  use fTree->GetTree()->GetEntry(entry).

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  if (!fEsdTrackCuts)
  {
    AliDebug(AliLog::kError, "fESDTrackCuts not available");
    return kFALSE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  if (AliPWG0Helper::IsEventTriggered(fESD) == kFALSE)
    return kTRUE;

  if (AliPWG0Helper::IsVertexReconstructed(fESD) == kFALSE)
    return kTRUE;

  // get the ESD vertex
  const AliESDVertex* vtxESD = fESD->GetVertex();
  Double_t vtx[3];
  vtxESD->GetXYZ(vtx);

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // get tracklets
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  if (!mult)
  {
    AliDebug(AliLog::kError, "AliMultiplicity not available");
    return kFALSE;
  }

  //const Float_t kPtCut = 0.3;

  // get number of tracks from MC
  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();

  // tracks in different eta ranges
  Int_t nMCTracks05 = 0;
  Int_t nMCTracks10 = 0;
  Int_t nMCTracks15 = 0;
  Int_t nMCTracks20 = 0;
  Int_t nMCTracksAll = 0;

  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    AliDebug(AliLog::kDebug+1, Form("MC Loop: Processing particle %d.", iMc));

    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
      continue;
    }

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    //if (particle->Pt() < kPtCut)
    //  continue;

    if (TMath::Abs(particle->Eta()) < 0.5)
      nMCTracks05++;

    if (TMath::Abs(particle->Eta()) < 1.0)
      nMCTracks10++;

    if (TMath::Abs(particle->Eta()) < 1.5)
      nMCTracks15++;

    if (TMath::Abs(particle->Eta()) < 2.0)
      nMCTracks20++;

    nMCTracksAll++;
  }// end of mc particle

  // FAKE: put back vtxMC[2]
  fMultiplicity->FillGenerated(vtxMC[2], nMCTracks05, nMCTracks10, nMCTracks15, nMCTracks20, nMCTracksAll);

  Int_t nESDTracks05 = 0;
  Int_t nESDTracks10 = 0;
  Int_t nESDTracks15 = 0;
  Int_t nESDTracks20 = 0;

  // get multiplicity from ITS tracklets
  for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
  {
    //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));

    // this removes non-tracklets. Very bad solution. SPD guys are working on better solution...
    if (mult->GetDeltaPhi(i) < -1000)
      continue;

    Float_t theta = mult->GetTheta(i);
    Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));

    if (TMath::Abs(eta) < 0.5)
      nESDTracks05++;

    if (TMath::Abs(eta) < 1.0)
      nESDTracks10++;

    if (TMath::Abs(eta) < 1.5)
      nESDTracks15++;

    if (TMath::Abs(eta) < 2.0)
      nESDTracks20++;
  }

  // get multiplicity from ESD tracks
  /*TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);
  Int_t nGoodTracks = list->GetEntries();
  // loop over esd tracks
  for (Int_t i=0; i<nGoodTracks; i++)
  {
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (list->At(i));
    if (!esdTrack)
    {
      AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", i));
      continue;
    }

    Double_t p[3];
    esdTrack->GetConstrainedPxPyPz(p); // ### TODO should be okay because we have a vertex, however GetInnerPxPyPy / GetOuterPxPyPy also exist
    TVector3 vector(p);

    Float_t theta = vector.Theta();
    Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));
    Float_t pt = vector.Pt();

    if (pt < kPtCut)
      continue;

    if (TMath::Abs(eta) < 0.5)
      nESDTracks05++;

    if (TMath::Abs(eta) < 1.0)
      nESDTracks10++;

    if (TMath::Abs(eta) < 1.5)
      nESDTracks15++;

    if (TMath::Abs(eta) < 2.0)
      nESDTracks20++;
  }*/

  fMultiplicity->FillMeasured(vtxMC[2], nESDTracks05, nESDTracks10, nESDTracks15, nESDTracks20);

  // TODO should this be vtx[2] or vtxMC[2] ?
  fMultiplicity->FillCorrection(vtxMC[2], nMCTracks05, nMCTracks10, nMCTracks15, nMCTracks20, nMCTracksAll, nESDTracks05, nESDTracks10, nESDTracks15, nESDTracks20);

  return kTRUE;
}

void AliMultiplicityMCSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelector::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  fOutput->Add(fMultiplicity);
}

void AliMultiplicityMCSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fMultiplicity = dynamic_cast<AliMultiplicityCorrection*> (fOutput->FindObject("Multiplicity"));

  if (!fMultiplicity)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p", (void*) fMultiplicity));
    return;
  }

  TFile* file = TFile::Open("multiplicityMC.root", "RECREATE");

  fMultiplicity->SaveHistograms();

  file->Close();

  fMultiplicity->DrawHistograms();
}
