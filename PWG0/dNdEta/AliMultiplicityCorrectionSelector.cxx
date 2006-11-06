/* $Id$ */

#include "AliMultiplicityCorrectionSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>

#include <TChain.h>
#include <TSelector.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliTracker.h>
#include <AliESDVertex.h>
#include <AliESD.h>
#include <AliESDtrack.h>
#include <AliRunLoader.h>
#include <AliStack.h>

#include <AliMultiplicity.h>

#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <../PYTHIA6/AliGenPythiaEventHeader.h>
#include <../EVGEN/AliGenCocktailEventHeader.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEta/AliMultiplicityCorrection.h"
#include "AliPWG0Helper.h"
#include "AliPWG0depHelper.h"

ClassImp(AliMultiplicityCorrectionSelector)

AliMultiplicityCorrectionSelector::AliMultiplicityCorrectionSelector() :
  AliSelectorRL()
{
  //
  // Constructor. Initialization of pointers
  //
}

AliMultiplicityCorrectionSelector::~AliMultiplicityCorrectionSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}


void AliMultiplicityCorrectionSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");
}

void AliMultiplicityCorrectionSelector::Begin(TTree * tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::Begin(tree);

  ReadUserObjects(tree);

  TString option = GetOption();
  AliInfo(Form("Called with option %s.", option.Data()));

}

void AliMultiplicityCorrectionSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::SlaveBegin(tree);

  ReadUserObjects(tree);

  fMultiplicityCorrection = new AliMultiplicityCorrection("mult_correction", "mult_correction");

}

void AliMultiplicityCorrectionSelector::Init(TTree* tree)
{
  // read the user objects

  AliSelectorRL::Init(tree);

  // Enable only the needed branches
  if (tree)
  {
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fTriggerMask", 1);
    tree->SetBranchStatus("fSPDVertex*", 1);
    tree->SetBranchStatus("fTracks.fLabel", 1);
    tree->SetBranchStatus("fTracks.fITSncls", 1);
    tree->SetBranchStatus("fTracks.fTPCncls", 1);
    tree->SetBranchStatus("fSPDMult*", 1);

    AliESDtrackCuts::EnableNeededBranches(tree);
  }
}

Bool_t AliMultiplicityCorrectionSelector::Process(Long64_t entry)
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

  AliDebug(AliLog::kDebug+1,"Processing event ...\n");

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // check prerequesites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }

  // getting the stack
  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  // getting the its multiplicity data
  AliMultiplicity* itsMult = (AliMultiplicity*)fESD->GetMultiplicity();  
  if (!itsMult) {
    AliDebug(AliLog::kError, "Multiplicity object not found in ESD");
    return kFALSE;
  }

  //  if (!fEsdTrackCuts)
  // {
  //  AliDebug(AliLog::kError, "fESDTrackCuts not available");
  //  return kFALSE;
  // }

  Bool_t vertexReconstructed = AliPWG0Helper::IsVertexReconstructed(fESD);

  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD);

  // only look at triggered events with a vertex
  if (!(vertexReconstructed && eventTriggered))
    return kTRUE;

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // vertex cut - should not be hard coded!!!
  if (TMath::Abs(vtxMC[2]>10))
    return kTRUE;

  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();

  for (Int_t iMc = 0; iMc < nPrim; ++iMc) {
    
    AliDebug(AliLog::kDebug+1, Form("MC Loop: Processing particle %d.", iMc));

    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
      continue;
    }

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    if (particle->GetPDG()->Charge()==0)
      continue;

    Float_t eta = particle->Eta();

    fMultiplicityCorrection->FillMeasuredMultHit(eta);

  }// end of mc particle

  // ########################################################
  // loop over spd tracklets

  
  for(Int_t i=0; i<itsMult->GetNumberOfTracklets(); i++) {
    Float_t theta  = itsMult->GetTheta(i);
    Float_t eta    = -TMath::Log(TMath::Tan(0.5*theta));
    
    fMultiplicityCorrection->FillTrueMultHit(eta);
  } 

  // very important!
  fMultiplicityCorrection->NewEvent();

  return kTRUE;
}

void AliMultiplicityCorrectionSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, "ERROR: Output list not initialized");
    return;
  }

  fOutput->Add(fMultiplicityCorrection);

  AliDebug(AliLog::kDebug+1,"Slave terminate ...\n");

}

void AliMultiplicityCorrectionSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliDebug(AliLog::kDebug+1,"Terminate ...\n");


  AliSelectorRL::Terminate();

  fMultiplicityCorrection = dynamic_cast<AliMultiplicityCorrection*> (fOutput->FindObject("mult_correction"));
  if (!fMultiplicityCorrection)
  {
    AliDebug(AliLog::kError, "Could not read object from output list");
    return;
  }

  fMultiplicityCorrection->Finish();

  TFile* fout = new TFile(Form("correction_mult%s.root", GetOption()), "RECREATE");

  //  if (fEsdTrackCuts)
  //  fEsdTrackCuts->SaveHistograms("esd_track_cuts");
  fMultiplicityCorrection->SaveHistograms();

  fout->Write();
  fout->Close();

  fMultiplicityCorrection->DrawHistograms();

}
