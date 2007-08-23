/* $Id$ */

#include "AliGenTriggerMapSelector.h"

#include <TVector3.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliRunLoader.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>

#include <AliITSgeom.h>
#include <AliITSLoader.h>
#include <AliITSdigitSPD.h>
#include <AliITSRecPoint.h>

#include "AliPWG0Helper.h"

//
//

ClassImp(AliGenTriggerMapSelector)

AliGenTriggerMapSelector::AliGenTriggerMapSelector() :
  AliSelectorRL(),
  fChipsFired(0),
  fTracklets(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliGenTriggerMapSelector::~AliGenTriggerMapSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliGenTriggerMapSelector::SlaveBegin(TTree *tree)
{
  AliSelectorRL::SlaveBegin(tree);

  fChipsFired = new TH2F("fChipsFired", ";Module;Chip;Count", 240, -0.5, 239.5, 5, -0.5, 4.5);
  fTracklets = new TNtuple("fTracklets", "", "vertex:chip1:chip2:dirty");
}

void AliGenTriggerMapSelector::Init(TTree* tree)
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
  }
}

Bool_t AliGenTriggerMapSelector::Process(Long64_t entry)
{
  //
  // processing
  //

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD, AliPWG0Helper::kMB1);

  if (!eventTriggered)
  {
    AliDebug(AliLog::kDebug, "Event not triggered. Skipping.");
    return kTRUE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);
  Float_t zVertex = vtxMC[2];

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }
  Int_t nPrim  = stack->GetNprimary();

  AliRunLoader* runLoader = GetRunLoader();
  if (!runLoader)
  {
    AliDebug(AliLog::kError, "runloader not available");
    return kFALSE;
  }

  // TDirectory::TContext restores the current directory is restored when the scope ends.
  // This helps around ROOT bug #26025 and is good behaviour anyway
  TDirectory::TContext context(0);
  AliITSLoader* loader = (AliITSLoader*) runLoader->GetLoader( "ITSLoader" );
  loader->LoadDigits("READ");
  TTree* treeD = loader->TreeD();
  if (!treeD)
  {
    AliDebug(AliLog::kError, "Could not retrieve TreeD of ITS");
    return kFALSE;
  }

  treeD->SetBranchStatus("*", 0);
  treeD->SetBranchStatus("ITSDigitsSPD.fTracks*", 1);
  treeD->SetBranchStatus("ITSDigitsSPD.fCoord1", 1);

  TClonesArray* digits = 0;
  treeD->SetBranchAddress("ITSDigitsSPD", &digits);
  if (digits);
    digits->Clear();

  const Int_t startSPD = 0; //geom->GetStartSPD();
  const Int_t lastSPD  = 239; //geom->GetLastSPD();

  const AliMultiplicity* mult = fESD->GetMultiplicity();
  if (!mult)
  {
    AliDebug(AliLog::kError, "AliMultiplicity not available");
    return kFALSE;
  }

  // loop over tracklets
  for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
  {
    Float_t theta = mult->GetTheta(i);
    Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));

    Int_t label = mult->GetLabel(i);
    if (label >= nPrim)
    {
      AliDebug(AliLog::kDebug, Form("Skipping particle %d because it is not a primary (label %d)", i, label));
      continue;
    }
    if (label < 0)
    {
      AliDebug(AliLog::kDebug, Form("Skipping particle %d because its label is negative (label %d)", i, label));
      continue;
    }

    Int_t chip[2];
    Int_t nClusters[2];

    for (Int_t j=0; j<2; ++j)
    {
      chip[j] = -1;
      nClusters[j] = 0;
    }

    // now find clusters that belong to this label
    // loop over modules (ladders)
    for (Int_t moduleIndex=startSPD; moduleIndex<lastSPD+1; moduleIndex++)
    {
      Int_t currentLayer = 0;
      if (moduleIndex >= 80)
        currentLayer = 1;

      treeD->GetEvent(moduleIndex);

      // get number of digits in this module
      Int_t ndigitsInModule = digits->GetEntriesFast();

      // loop over digits in this module
      for (Int_t iDig=0; iDig<ndigitsInModule; iDig++)
      {
        AliITSdigitSPD* dp = (AliITSdigitSPD*) digits->At(iDig);

        // check if digit belongs to this tracklet
        Bool_t belongs = kFALSE;
        for (Int_t iTrack = 0; iTrack < dp->GetNTracks(); ++iTrack)
        {
          //Printf("%d %d %d %d: %d =? %d", i, moduleIndex, iDig, iTrack, label, dp->GetTrack(iTrack));
          if (dp->GetTrack(iTrack) == label)
            belongs = kTRUE;
        }

        if (!belongs)
          continue;

        Int_t column = dp->GetCoord1();
        Int_t isChip = column / 32;

        //printf("Digit %d has column %d which translates to chip %d\n", iDig, column, isChip);

        fChipsFired->Fill(moduleIndex, isChip);

        if (chip[currentLayer] == -1)
        {
          chip[currentLayer] = moduleIndex * 5 + isChip;
          nClusters[currentLayer]++;
        }
        else
        {
          if (chip[currentLayer] != moduleIndex * 5 + isChip)
            nClusters[currentLayer]++;
        }
      }
    }

    AliDebug(AliLog::kDebug, Form("Tracklet %d fired by chip %d and %d", i, chip[0], chip[1]));

    Bool_t dirty = kFALSE;
    if (nClusters[0] == 0 || nClusters[1] == 0)
    {
      AliDebug(AliLog::kDebug, Form("Not enough clusters found for tracklet %d", i));
      continue;
    }
    if (nClusters[0] > 1 || nClusters[1] > 1)
    {
      AliDebug(AliLog::kDebug, Form("Too many clusters found for tracklet %d; L1: %d; L2: %d", i, nClusters[0], nClusters[1]));
      dirty = kTRUE;
    }
    fTracklets->Fill(zVertex, chip[0], chip[1], dirty);
  }

  return kTRUE;
}

Bool_t AliGenTriggerMapSelector::Notify()
{
  AliRunLoader* runLoader = GetRunLoader();

  if (runLoader)
  {
    AliITSLoader* loader = (AliITSLoader* )runLoader->GetLoader( "ITSLoader" );
    if (loader)
    {
      loader->UnloadDigits();
      loader->UnloadRecPoints();
    }
  }

  return AliSelectorRL::Notify();
}

void AliGenTriggerMapSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, "ERROR: Output list not initialized.");
    return;
  }

  fOutput->Add(fChipsFired);
  fOutput->Add(fTracklets);
}

void AliGenTriggerMapSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelectorRL::Terminate();

  fChipsFired = dynamic_cast<TH2F*> (fOutput->FindObject("fChipsFired"));
  fTracklets = dynamic_cast<TNtuple*> (fOutput->FindObject("fTracklets"));

  if (!fTracklets || !fChipsFired)
  {
    AliError("Histograms not available");
    return;
  }

  WriteHistograms();
}

void AliGenTriggerMapSelector::WriteHistograms(const char* filename)
{
  TFile* file = TFile::Open(filename, "RECREATE");

  fChipsFired->Write();
  fTracklets->Write();

  file->Close();
}

void AliGenTriggerMapSelector::ReadHistograms(const char* filename)
{
  TFile* file = TFile::Open(filename);

  if (!file)
    return;

  fTracklets  = dynamic_cast<TNtuple*> (file->Get("fTracklets"));
  fChipsFired  = dynamic_cast<TH2F*> (file->Get("fChipsFired"));
}
