/* $Id$ */

#include "AlidNdEtaAnalysisMCSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TTree.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>
#include <AliStack.h>
#include <AliESDtrack.h>

#include "dNdEta/dNdEtaAnalysis.h"
#include <AliPWG0Helper.h>
#include <AliCorrection.h>
#include <AliCorrectionMatrix2D.h>
#include "esdTrackCuts/AliESDtrackCuts.h"


ClassImp(AlidNdEtaAnalysisMCSelector)

AlidNdEtaAnalysisMCSelector::AlidNdEtaAnalysisMCSelector() :
  AliSelectorRL(),
  fEsdTrackCuts(0),
  fdNdEtaAnalysis(0),
  fdNdEtaAnalysisTr(0),
  fdNdEtaAnalysisTrVtx(0),
  fdNdEtaAnalysisTracks(0),
  fVertex(0),
  fPartPt(0),
  fEvents(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AlidNdEtaAnalysisMCSelector::~AlidNdEtaAnalysisMCSelector()
{
  //
  // Destructor
  //
}

void AlidNdEtaAnalysisMCSelector::Begin(TTree* tree)
{
  // Begin function

  ReadUserObjects(tree);
}

void AlidNdEtaAnalysisMCSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");
}

void AlidNdEtaAnalysisMCSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::SlaveBegin(tree);

  ReadUserObjects(tree);

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysisTr = new dNdEtaAnalysis("dndetaTr", "dndetaTr");
  fdNdEtaAnalysisTrVtx = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysisTracks = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks");
  fVertex = new TH3F("vertex_check", "vertex_check", 50, -50, 50, 50, -50, 50, 50, -50, 50);
  for (Int_t i=0; i<3; ++i)
  {
    fPartEta[i] = new TH1F(Form("dndeta_check_%d", i), Form("dndeta_check_%d", i), 60, -6, 6);
    fPartEta[i]->Sumw2();
  }
  fPartPt =  new TH1F("dndeta_check_pt", "dndeta_check_pt", 1000, 0, 10);
  fPartPt->Sumw2();
  fEvents = new TH1F("dndeta_check_vertex", "dndeta_check_vertex", 40, -20, 20);
}

void AlidNdEtaAnalysisMCSelector::Init(TTree *tree)
{
  AliSelectorRL::Init(tree);

  if (!tree)
  {
    AliDebug(AliLog::kError, "tree not available");
    return;
  }

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("fTriggerMask", 1);
  tree->SetBranchStatus("fSPDVertex*", 1);
  tree->SetBranchStatus("fTracks.fLabel", 1);

  AliESDtrackCuts::EnableNeededBranches(tree);
}

Bool_t AlidNdEtaAnalysisMCSelector::Process(Long64_t entry)
{
  // fill the dNdEtaAnalysis class from the monte carlo

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }

  Bool_t vertexReconstructed = AliPWG0Helper::IsVertexReconstructed(fESD);
  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD);

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();

  Int_t nAcceptedParticles = 0;

  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    TParticle* particle = stack->Particle(iMc);

    if (!particle)
      continue;

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    AliDebug(AliLog::kDebug+1, Form("Accepted primary %d, unique ID: %d", iMc, particle->GetUniqueID()));
    ++nAcceptedParticles;

    fdNdEtaAnalysis->FillTrack(vtxMC[2], particle->Eta(), particle->Pt());
    fVertex->Fill(particle->Vx(), particle->Vy(), particle->Vz());

    if (eventTriggered)
    {
      fdNdEtaAnalysisTr->FillTrack(vtxMC[2], particle->Eta(), particle->Pt());
      if (vertexReconstructed)
        fdNdEtaAnalysisTrVtx->FillTrack(vtxMC[2], particle->Eta(), particle->Pt());
    }

    if (TMath::Abs(vtxMC[2]) < 20)
    {
      fPartEta[0]->Fill(particle->Eta());

      if (vtxMC[2] < 0)
        fPartEta[1]->Fill(particle->Eta());
      else
        fPartEta[2]->Fill(particle->Eta());
    }

    if (TMath::Abs(particle->Eta()) < 0.8)
      fPartPt->Fill(particle->Pt());
  }

  fEvents->Fill(vtxMC[2]);

  fdNdEtaAnalysis->FillEvent(vtxMC[2], nAcceptedParticles);
  if (eventTriggered)
  {
    fdNdEtaAnalysisTr->FillEvent(vtxMC[2], nAcceptedParticles);
    if (vertexReconstructed)
      fdNdEtaAnalysisTrVtx->FillEvent(vtxMC[2], nAcceptedParticles);
  }

  if (!eventTriggered || !vertexReconstructed)
    return kTRUE;

  // from tracks is only done for triggered and vertex reconstructed events

  // get number of "good" tracks
  TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);
  Int_t nGoodTracks = list->GetEntries();

  // loop over esd tracks
  for (Int_t t=0; t<nGoodTracks; t++)
  {
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (list->At(t));
    if (!esdTrack)
    {
      AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", t));
      continue;
    }

    Int_t label = TMath::Abs(esdTrack->GetLabel());

    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("ERROR: Could not retrieve particle %d.", esdTrack->GetLabel()));
      continue;
    }

    fdNdEtaAnalysisTracks->FillTrack(vtxMC[2], particle->Eta(), particle->Pt());
  } // end of track loop

  delete list;
  list = 0;

  // for event count per vertex
  fdNdEtaAnalysisTracks->FillEvent(vtxMC[2], nGoodTracks);

  return kTRUE;
}

void AlidNdEtaAnalysisMCSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  fOutput->Add(fdNdEtaAnalysis);
  fOutput->Add(fdNdEtaAnalysisTr);
  fOutput->Add(fdNdEtaAnalysisTrVtx);
  fOutput->Add(fdNdEtaAnalysisTracks);

  fOutput->Add(fPartPt);
  fOutput->Add(fEvents);
  for (Int_t i=0; i<3; ++i)
    fOutput->Add(fPartEta[i]);
}

void AlidNdEtaAnalysisMCSelector::Terminate()
{
  //

  AliSelectorRL::Terminate();

  fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));
  fdNdEtaAnalysisTr = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTr"));
  fdNdEtaAnalysisTrVtx = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTrVtx"));
  fdNdEtaAnalysisTracks = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTracks"));
  fPartPt = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_pt"));
  fEvents = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_vertex"));
  for (Int_t i=0; i<3; ++i)
    fPartEta[i] = dynamic_cast<TH1F*> (fOutput->FindObject(Form("dndeta_check_%d", i)));

  if (!fdNdEtaAnalysis || !fdNdEtaAnalysisTr || !fdNdEtaAnalysisTrVtx || !fPartPt || !fEvents)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p", (void*) fdNdEtaAnalysis, (void*) fPartPt));
    return;
  }

  fdNdEtaAnalysis->Finish(0, -1, AlidNdEtaCorrection::kNone);
  fdNdEtaAnalysisTr->Finish(0, -1, AlidNdEtaCorrection::kNone);
  fdNdEtaAnalysisTrVtx->Finish(0, -1, AlidNdEtaCorrection::kNone);
  fdNdEtaAnalysisTracks->Finish(0, -1, AlidNdEtaCorrection::kNone);

  Int_t events = (Int_t) fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetMeasuredHistogram()->Integral();
  fPartPt->Scale(1.0/events);
  fPartPt->Scale(1.0/fPartPt->GetBinWidth(1));

  if (fPartEta[0])
  {
    Int_t events1 = (Int_t) fEvents->Integral(fEvents->GetXaxis()->FindBin(-19.9), fEvents->GetXaxis()->FindBin(-0.001));
    Int_t events2 = (Int_t) fEvents->Integral(fEvents->GetXaxis()->FindBin(0.001), fEvents->GetXaxis()->FindBin(19.9));

    fPartEta[0]->Scale(1.0 / (events1 + events2));
    fPartEta[1]->Scale(1.0 / events1);
    fPartEta[2]->Scale(1.0 / events2);

    for (Int_t i=0; i<3; ++i)
      fPartEta[i]->Scale(1.0 / fPartEta[i]->GetBinWidth(1));

    new TCanvas("control", "control", 500, 500);
    for (Int_t i=0; i<3; ++i)
    {
      fPartEta[i]->SetLineColor(i+1);
      fPartEta[i]->Draw((i==0) ? "" : "SAME");
    }
  }

  TFile* fout = new TFile("analysis_mc.root","RECREATE");

  fdNdEtaAnalysis->SaveHistograms();
  fdNdEtaAnalysisTr->SaveHistograms();
  fdNdEtaAnalysisTrVtx->SaveHistograms();
  fdNdEtaAnalysisTracks->SaveHistograms();
  fPartPt->Write();

  fout->Write();
  fout->Close();

  if (fPartPt)
  {
    new TCanvas("control2", "control2", 500, 500);
    fPartPt->DrawCopy();
  }

  if (fEvents)
  {
    new TCanvas("control3", "control3", 500, 500);
    fEvents->DrawCopy();
  }
}
