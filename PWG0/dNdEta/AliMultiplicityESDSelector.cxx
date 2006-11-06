/* $Id$ */

#include "AliMultiplicityESDSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include <AliLog.h>
#include <AliESD.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "AliPWG0Helper.h"

#ifdef ALISELECTOR_USEMONALISA
  #include <TMonaLisaWriter.h>
#endif

ClassImp(AliMultiplicityESDSelector)

AliMultiplicityESDSelector::AliMultiplicityESDSelector() :
  AliSelector(),
  fMultiplicity(0),
  fEsdTrackCuts(0)
#ifdef ALISELECTOR_USEMONALISA
  ,fMonaLisaWriter(0)
#endif
{
  //
  // Constructor. Initialization of pointers
  //
}

AliMultiplicityESDSelector::~AliMultiplicityESDSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliMultiplicityESDSelector::Begin(TTree* tree)
{
  // Begin function

  AliSelector::Begin(tree);

  ReadUserObjects(tree);
}

void AliMultiplicityESDSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");
}

void AliMultiplicityESDSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  ReadUserObjects(tree);

  fMultiplicity = new TH1F("fMultiplicity", "multiplicity", 201, 0.5, 200.5);

  #ifdef ALISELECTOR_USEMONALISA
    TNamed *nm = 0;
    if (fInput)
      nm = dynamic_cast<TNamed*> (fInput->FindObject("PROOF_QueryTag"));
    if (!nm)
    {
      AliDebug(AliLog::kError, "Query tag not found. Cannot enable monitoring");
      return;
    }

    TString option = GetOption();
    option.ReplaceAll("#+", "");

    TString id;
    id.Form("%s_%s%d", gSystem->HostName(), nm->GetTitle(), gSystem->GetPid());
    fMonaLisaWriter = new TMonaLisaWriter(option, id, "CAF", "aliendb6.cern.ch");
  #endif
}

Bool_t AliMultiplicityESDSelector::Process(Long64_t entry)
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

  if (AliSelector::Process(entry) == kFALSE)
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

  if (AliPWG0Helper::IsEventTriggered(fESD) == kFALSE)
    return kTRUE;

  if (AliPWG0Helper::IsVertexReconstructed(fESD) == kFALSE)
    return kTRUE;

  // get number of "good" tracks
  Int_t nGoodTracks = fEsdTrackCuts->CountAcceptedTracks(fESD);

  fMultiplicity->Fill(nGoodTracks);

  return kTRUE;
}

void AliMultiplicityESDSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelector::SlaveTerminate();

  #ifdef ALISELECTOR_USEMONALISA
    if (fMonaLisaWriter)
    {
      delete fMonaLisaWriter;
      fMonaLisaWriter = 0;
    }
  #endif

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  fOutput->Add(fMultiplicity);
}

void AliMultiplicityESDSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fMultiplicity = dynamic_cast<TH1F*> (fOutput->FindObject("fMultiplicity"));

  if (!fMultiplicity)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histogram not available %p", (void*) fMultiplicity));
    return;
  }

  TFile* file = TFile::Open("multiplicityESD.root", "RECREATE");
  fMultiplicity->Write();
  file->Close();
}
