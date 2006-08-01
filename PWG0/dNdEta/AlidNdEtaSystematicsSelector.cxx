/* $Id$ */

#include "AlidNdEtaSystematicsSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3F.h>
#include <TH1F.h>
#include <TParticle.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliGenEventHeader.h>
#include <AliStack.h>
#include <AliHeader.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "AliPWG0Helper.h"
#include "AlidNdEtaCorrection.h"

ClassImp(AlidNdEtaSystematicsSelector)

AlidNdEtaSystematicsSelector::AlidNdEtaSystematicsSelector() :
  AliSelectorRL(),
  fSecondaries(0),
  fSigmaVertex(0),
  fEsdTrackCuts(0),
  fOverallPrimaries(0),
  fOverallSecondaries(0)
{
  //
  // Constructor. Initialization of pointers
  //

  for (Int_t i=0; i<4; ++i)
    fdNdEtaCorrection[i] = 0;
}

AlidNdEtaSystematicsSelector::~AlidNdEtaSystematicsSelector()
{
  //
  // Destructor
  //
}

void AlidNdEtaSystematicsSelector::Begin(TTree* tree)
{
  // Begin function

  AliSelectorRL::Begin(tree);

  ReadUserObjects(tree);
}

void AlidNdEtaSystematicsSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");
}

void AlidNdEtaSystematicsSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  ReadUserObjects(tree);

  TString option(GetOption());

  printf("Running AlidNdEtaSystematicsSelector with options %s\n", option.Data());

  if (option.Contains("secondaries"))
  {
    fSecondaries = new TH3F("fSecondaries", "fSecondaries;NacceptedTracks;CratioSecondaries;p_{T}", 2000, -0.5, 205.5, 16, 0.45, 2.05, 10, 0, 10);
  }

  if (option.Contains("particle-composition"))
  {
    for (Int_t i=0; i<4; ++i)
    {
      TString name;
      name.Form("correction_%d", i);
      fdNdEtaCorrection[i] = new AlidNdEtaCorrection(name, name);
    }
  }

  if (option.Contains("sigma-vertex"))
  {
    fSigmaVertex = new TH1F("fSigmaVertex", "fSigmaVertex;Nsigma2vertex;NacceptedTracks", 10, 0.25, 5.25);
    printf("WARNING: sigma-vertex analysis enabled. This will produce weird results in the AliESDtrackCuts histograms\n");
  }
}

Bool_t AlidNdEtaSystematicsSelector::Process(Long64_t entry)
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

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);

  if (fdNdEtaCorrection[0])
    FillCorrectionMaps(list);

  if (fSecondaries)
    FillSecondaries(list);

  if (fSigmaVertex)
    FillSigmaVertex();

  delete list;
  list = 0;

  return kTRUE;
}

void AlidNdEtaSystematicsSelector::FillCorrectionMaps(TObjArray* listOfTracks)
{
  // fills the correction maps for different particle species

  AliStack* stack = GetStack();
  AliHeader* header = GetHeader();

  Bool_t vertexReconstructed = AliPWG0Helper::IsVertexReconstructed(fESD);
  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD);

  for (Int_t i=0; i<4; ++i)
  {
    fdNdEtaCorrection[i]->IncreaseEventCount();
    if (eventTriggered)
      fdNdEtaCorrection[i]->IncreaseTriggeredEventCount();
  }

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();

  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
      continue;
    }

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    Float_t eta = particle->Eta();
    Float_t pt = particle->Pt();

    Int_t id = -1;
    switch (TMath::Abs(particle->GetPdgCode()))
    {
      case 211: id = 0; break;
      case 321: id = 1; break;
      case 2212: id = 2; break;
      default: id = 3; break;
    }

    if (vertexReconstructed)
      fdNdEtaCorrection[id]->FillParticle(vtxMC[2], eta, pt);

    fdNdEtaCorrection[id]->FillParticleAllEvents(eta, pt);
    if (eventTriggered)
      fdNdEtaCorrection[id]->FillParticleWhenEventTriggered(eta, pt);
  }// end of mc particle

  // loop over esd tracks
  TIterator* iter = listOfTracks->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next()) != 0)
  {
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (obj);
    if (!esdTrack)
      continue;

    // using the properties of the mc particle
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label));
      continue;
    }

    Int_t id = -1;
    switch (TMath::Abs(particle->GetPdgCode()))
    {
      case 211:  id = 0; break;
      case 321:  id = 1; break;
      case 2212: id = 2; break;
      default:   id = 3; break;
    }

    if (vertexReconstructed)
      fdNdEtaCorrection[id]->FillParticleWhenMeasuredTrack(vtxMC[2], particle->Eta(), particle->Pt());
  } // end of track loop

  delete iter;
  iter = 0;
}

void AlidNdEtaSystematicsSelector::FillSecondaries(TObjArray* listOfTracks)
{
  // fills the secondary histograms

  AliStack* stack = GetStack();

  TH1* nPrimaries = new TH1F("secondaries_primaries", "secondaries_primaries", fSecondaries->GetZaxis()->GetNbins(), fSecondaries->GetZaxis()->GetXmin(), fSecondaries->GetZaxis()->GetXmax());
  TH1* nSecondaries = new TH1F("secondaries_secondaries", "secondaries_secondaries", fSecondaries->GetZaxis()->GetNbins(), fSecondaries->GetZaxis()->GetXmin(), fSecondaries->GetZaxis()->GetXmax());

  TIterator* iter = listOfTracks->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next()) != 0)
  {
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (obj);
    if (!esdTrack)
      continue;

    Int_t label = TMath::Abs(esdTrack->GetLabel());
    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label));
      continue;
    }

    Int_t nPrim  = stack->GetNprimary();
    if (label < nPrim)
      nPrimaries->Fill(particle->Pt());
    else
      nSecondaries->Fill(particle->Pt());
  }

  for (Int_t i=1; i<=nPrimaries->GetNbinsX(); ++i)
  {
    Int_t primaries = (Int_t) nPrimaries->GetBinContent(i);
    Int_t secondaries = (Int_t) nSecondaries->GetBinContent(i);

    if (primaries + secondaries > 0)
    {
      AliDebug(AliLog::kDebug, Form("The ratio between primaries and secondaries is %d:%d = %f", primaries, secondaries, ((secondaries > 0) ? (Double_t) primaries / secondaries : -1)));

      for (Double_t factor = 0.5; factor < 2.01; factor += 0.1)
      {
        Double_t nTracks = (Double_t) primaries + (Double_t) secondaries * factor;
        fSecondaries->Fill(nTracks, factor, nPrimaries->GetBinCenter(i));
        //if (secondaries > 0) printf("We fill: %f %f %f\n", nTracks, factor, nPrimaries->GetBinCenter(i));
      }
    }
  }

  fOverallPrimaries += (Int_t) nPrimaries->Integral();
  fOverallSecondaries += (Int_t) nSecondaries->Integral();

  delete nPrimaries;
  nPrimaries = 0;

  delete nSecondaries;
  nSecondaries = 0;

  delete iter;
  iter = 0;
}

void AlidNdEtaSystematicsSelector::FillSigmaVertex()
{
  // fills the fSigmaVertex histogram

  // save the old value
  Float_t oldSigmaVertex = fEsdTrackCuts->GetMinNsigmaToVertex();

  // set to maximum
  fEsdTrackCuts->SetMinNsigmaToVertex(5);

  TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next()) != 0)
  {
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (obj);
    if (!esdTrack)
      continue;

    Float_t sigma2Vertex = fEsdTrackCuts->GetSigmaToVertex(esdTrack);

    for (Double_t nSigma = 0.5; nSigma < 5.1; nSigma += 0.5)
    {
      if (sigma2Vertex < nSigma)
        fSigmaVertex->Fill(nSigma);
    }
  }

  delete iter;
  iter = 0;

  delete list;
  list = 0;

  // set back the old value
  fEsdTrackCuts->SetMinNsigmaToVertex(oldSigmaVertex);
}

void AlidNdEtaSystematicsSelector::SlaveTerminate()
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

  if (fSecondaries)
    fOutput->Add(fSecondaries);

  for (Int_t i=0; i<4; ++i)
    if (fdNdEtaCorrection[i])
      fOutput->Add(fdNdEtaCorrection[i]);

  if (fSigmaVertex)
    fOutput->Add(fSigmaVertex);
}

void AlidNdEtaSystematicsSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fSecondaries = dynamic_cast<TH3F*> (fOutput->FindObject("fSecondaries"));
  for (Int_t i=0; i<4; ++i)
    fdNdEtaCorrection[i] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject(Form("correction_%d", i)));
  fSigmaVertex = dynamic_cast<TH1F*> (fOutput->FindObject("fSigmaVertex"));

  TFile* fout = TFile::Open("systematics.root", "RECREATE");

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_track_cuts");

  if (fSecondaries)
  {
    fSecondaries->Write();
    printf("We had %d primaries and %d secondaries.\n", (Int_t) fOverallPrimaries, (Int_t) fOverallSecondaries);
  }

  if (fSigmaVertex)
    fSigmaVertex->Write();

  for (Int_t i=0; i<4; ++i)
    if (fdNdEtaCorrection[i])
      fdNdEtaCorrection[i]->SaveHistograms();

  fout->Write();
  fout->Close();

  if (fSecondaries)
  {
    new TCanvas;
    fSecondaries->Draw();
  }

  if (fSigmaVertex)
  {
    new TCanvas;
    fSigmaVertex->Draw();
  }
}
