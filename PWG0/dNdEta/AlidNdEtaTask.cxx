/* $Id$ */

#include "AlidNdEtaTask.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TObjString.h>
#include <TF1.h>

#include <AliLog.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliMultiplicity.h>
#include <AliAnalysisManager.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliESDInputHandler.h>

#include "AliESDtrackCuts.h"
#include "AliPWG0Helper.h"
#include "AliCorrection.h"
#include "AliCorrectionMatrix3D.h"
#include "dNdEta/dNdEtaAnalysis.h"

ClassImp(AlidNdEtaTask)

AlidNdEtaTask::AlidNdEtaTask(const char* opt) :
  AliAnalysisTask("AlidNdEtaTask", ""),
  fESD(0),
  fOutput(0),
  fOption(opt),
  fAnalysisMode(AliPWG0Helper::kTPC),
  fTrigger(AliPWG0Helper::kMB1),
  fReadMC(kFALSE),
  fUseMCVertex(kFALSE),
  fUseMCKine(kFALSE),
  fEsdTrackCuts(0),
  fdNdEtaAnalysisESD(0),
  fMult(0),
  fMultVtx(0),
  fEvents(0),
  fVertexResolution(0),
  fdNdEtaAnalysis(0),
  fdNdEtaAnalysisNSD(0),
  fdNdEtaAnalysisTr(0),
  fdNdEtaAnalysisTrVtx(0),
  fdNdEtaAnalysisTracks(0),
  fVertex(0),
  fPartPt(0),
  fDeltaPhi(0)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  AliLog::SetClassDebugLevel("AlidNdEtaTask", AliLog::kWarning);
}

AlidNdEtaTask::~AlidNdEtaTask()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

Bool_t AlidNdEtaTask::Notify()
{
  static Int_t count = 0;
  count++;
  Printf("Processing %d. file", count);
  return kTRUE;
}

//________________________________________________________________________
void AlidNdEtaTask::ConnectInputData(Option_t *)
{
  // Connect ESD
  // Called once

  Printf("AlidNdEtaTask::ConnectInputData called");

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
  } else {
    fESD = esdH->GetEvent();

    TTree* tree = esdH->GetTree();
    if (!tree) {
      Printf("ERROR: Could not read tree");
    } else {
      // Disable all branches and enable only the needed ones
      tree->SetBranchStatus("*", 0);

      tree->SetBranchStatus("AliESDHeader*", 1);
      tree->SetBranchStatus("*Vertex*", 1);

      if (fAnalysisMode == AliPWG0Helper::kSPD) {
        tree->SetBranchStatus("AliMultiplicity*", 1);
      }

      if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS) {
        //AliESDtrackCuts::EnableNeededBranches(tree);
        tree->SetBranchStatus("Tracks*", 1);
      }

    }
  }

  // disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);
}

void AlidNdEtaTask::CreateOutputObjects()
{
  // create result objects and add to output list

  Printf("AlidNdEtaTask::CreateOutputObjects");

  fOutput = new TList;
  fOutput->SetOwner();

  fdNdEtaAnalysisESD = new dNdEtaAnalysis("fdNdEtaAnalysisESD", "fdNdEtaAnalysisESD", fAnalysisMode);
  fOutput->Add(fdNdEtaAnalysisESD);

  fMult = new TH1F("fMult", "fMult;Ntracks;Count", 201, -0.5, 200.5);
  fOutput->Add(fMult);

  fMultVtx = new TH1F("fMultVtx", "fMultVtx;Ntracks;Count", 201, -0.5, 200.5);
  fOutput->Add(fMultVtx);

  for (Int_t i=0; i<3; ++i)
  {
    fPartEta[i] = new TH1F(Form("dndeta_check_%d", i), Form("dndeta_check_%d", i), 60, -6, 6);
    fPartEta[i]->Sumw2();
    fOutput->Add(fPartEta[i]);
  }

  fEvents = new TH1F("dndeta_check_vertex", "dndeta_check_vertex", 160, -40, 40);
  fOutput->Add(fEvents);

  fVertexResolution = new TH1F("dndeta_vertex_resolution_z", "dndeta_vertex_resolution_z", 1000, 0, 10);
  fOutput->Add(fVertexResolution);

  if (fReadMC)
  {
    fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysis);

    fdNdEtaAnalysisNSD = new dNdEtaAnalysis("dndetaNSD", "dndetaNSD", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisNSD);

    fdNdEtaAnalysisTr = new dNdEtaAnalysis("dndetaTr", "dndetaTr", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisTr);

    fdNdEtaAnalysisTrVtx = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisTrVtx);

    fdNdEtaAnalysisTracks = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisTracks);

    fVertex = new TH3F("vertex_check", "vertex_check", 50, -50, 50, 50, -50, 50, 50, -50, 50);
    fOutput->Add(fVertex);

    fPartPt =  new TH1F("dndeta_check_pt", "dndeta_check_pt", 1000, 0, 10);
    fPartPt->Sumw2();
    fOutput->Add(fPartPt);
  }

  if (fAnalysisMode == AliPWG0Helper::kSPD)
    fDeltaPhi = new TH1F("fDeltaPhi", "fDeltaPhi;#Delta #phi;Entries", 1000, -3.14, 3.14);
}

void AlidNdEtaTask::Exec(Option_t*)
{
  // process the event

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return;
  }

  // trigger definition
  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD->GetTriggerMask(), fTrigger);

  // get the ESD vertex
  const AliESDVertex* vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);
  
  Double_t vtx[3];
  if (vtxESD) {
    vtxESD->GetXYZ(vtx);
  }
  else
    Printf("No vertex");
  
  // fill z vertex resolution
  if (vtxESD)
    fVertexResolution->Fill(vtxESD->GetZRes());

  // post the data already here
  PostData(0, fOutput);

  // needed for syst. studies
  AliStack* stack = 0;
  
  if (fUseMCVertex || fUseMCKine) {
    if (!fReadMC) {
      Printf("ERROR: fUseMCVertex or fUseMCKine set without fReadMC set!");
      return;
    }

    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    AliHeader* header = mcEvent->Header();
    if (!header)
    {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }

    if (fUseMCVertex)
    {
      Printf("WARNING: Replacing vertex by MC vertex. This is for systematical checks only.");
      // get the MC vertex
      AliGenEventHeader* genHeader = header->GenEventHeader();
      TArrayF vtxMC(3);
      genHeader->PrimaryVertex(vtxMC);

      vtx[2] = vtxMC[2];
    }

    if (fUseMCKine)
    {
      stack = mcEvent->Stack();
      if (!stack)
      {
        AliDebug(AliLog::kError, "Stack not available");
        return;
      }
    }
  }

  // create list of (label, eta, pt) tuples
  Int_t inputCount = 0;
  Int_t* labelArr = 0;
  Float_t* etaArr = 0;
  Float_t* ptArr = 0;
  if (fAnalysisMode == AliPWG0Helper::kSPD)
  {
    // get tracklets
    const AliMultiplicity* mult = fESD->GetMultiplicity();
    if (!mult)
    {
      AliDebug(AliLog::kError, "AliMultiplicity not available");
      return;
    }

    labelArr = new Int_t[mult->GetNumberOfTracklets()];
    etaArr = new Float_t[mult->GetNumberOfTracklets()];
    ptArr = new Float_t[mult->GetNumberOfTracklets()];

    if (fUseMCKine && stack)
      Printf("Processing only primaries (MC information used). This is for systematical checks only.");

    // get multiplicity from ITS tracklets
    for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
    {
      //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));

      if (fUseMCKine && stack)
        if (mult->GetLabel(i, 0) < 0 || mult->GetLabel(i, 0) != mult->GetLabel(i, 1) || !stack->IsPhysicalPrimary(mult->GetLabel(i, 0)))
          continue;
      
      Float_t deltaPhi = mult->GetDeltaPhi(i);
      // prevent values to be shifted by 2 Pi()
      if (deltaPhi < -TMath::Pi())
        deltaPhi += TMath::Pi() * 2;
      if (deltaPhi > TMath::Pi())
        deltaPhi -= TMath::Pi() * 2;

      if (TMath::Abs(deltaPhi) > 1)
        printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), deltaPhi);

      fDeltaPhi->Fill(deltaPhi);

      etaArr[inputCount] = mult->GetEta(i);
      labelArr[inputCount] = mult->GetLabel(i, 0);
      ptArr[inputCount] = 0; // no pt for tracklets
      ++inputCount;
    }

    //Printf("Accepted %d tracklets", inputCount);
  }
  else if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS)
  {
    if (!fEsdTrackCuts)
    {
      AliDebug(AliLog::kError, "fESDTrackCuts not available");
      return;
    }

    // get multiplicity from ESD tracks
    TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);
    Int_t nGoodTracks = list->GetEntries();

    labelArr = new Int_t[nGoodTracks];
    etaArr = new Float_t[nGoodTracks];
    ptArr = new Float_t[nGoodTracks];

    // loop over esd tracks
    for (Int_t i=0; i<nGoodTracks; i++)
    {
      AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (list->At(i));
      if (!esdTrack)
      {
        AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", i));
        continue;
      }

      etaArr[inputCount] = esdTrack->Eta();
      labelArr[inputCount] = TMath::Abs(esdTrack->GetLabel());
      ptArr[inputCount] = esdTrack->Pt();
      ++inputCount;
    }

    Printf("Accepted %d tracks", nGoodTracks);

    delete list;
  }
  else
    return;

  // Processing of ESD information (always)
  if (eventTriggered)
  {
    // control hist
    fMult->Fill(inputCount);
    fdNdEtaAnalysisESD->FillTriggeredEvent(inputCount);

    if (vtxESD)
    {
      // control hist
      fMultVtx->Fill(inputCount);

      for (Int_t i=0; i<inputCount; ++i)
      {
        Float_t eta = etaArr[i];
        Float_t pt  = ptArr[i];

        fdNdEtaAnalysisESD->FillTrack(vtx[2], eta, pt);

        if (TMath::Abs(vtx[2]) < 20)
        {
          fPartEta[0]->Fill(eta);

          if (vtx[2] < 0)
            fPartEta[1]->Fill(eta);
          else
            fPartEta[2]->Fill(eta);
        }
      }

      // for event count per vertex
      fdNdEtaAnalysisESD->FillEvent(vtx[2], inputCount);

      // control hist
      fEvents->Fill(vtx[2]);
    }
  }

  if (fReadMC)   // Processing of MC information (optional)
  {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    stack = mcEvent->Stack();
    if (!stack)
    {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    AliHeader* header = mcEvent->Header();
    if (!header)
    {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }

    // get the MC vertex
    AliGenEventHeader* genHeader = header->GenEventHeader();
    if (!genHeader)
    {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }

    TArrayF vtxMC(3);
    genHeader->PrimaryVertex(vtxMC);

    // get process type; NB: this only works for Pythia
    Int_t processType = AliPWG0Helper::GetPythiaEventProcessType(header);
    AliDebug(AliLog::kDebug+1, Form("Found pythia process type %d", processType));

    if (processType<0)
      AliDebug(AliLog::kError, Form("Unknown Pythia process type %d.", processType));

    // loop over mc particles
    Int_t nPrim  = stack->GetNprimary();

    Int_t nAcceptedParticles = 0;

    for (Int_t iMc = 0; iMc < nPrim; ++iMc)
    {
      //Printf("Getting particle %d", iMc);
      TParticle* particle = stack->Particle(iMc);

      if (!particle)
        continue;

      if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
        continue;

      AliDebug(AliLog::kDebug+1, Form("Accepted primary %d, unique ID: %d", iMc, particle->GetUniqueID()));
      Float_t eta = particle->Eta();
      Float_t pt = particle->Pt();

       // make a rough eta cut (so that nAcceptedParticles is not too far off the true measured value (NB: this histograms are only gathered for comparison))
      if (TMath::Abs(eta) < 1.5 && pt > 0.3)
        nAcceptedParticles++;

      fdNdEtaAnalysis->FillTrack(vtxMC[2], eta, pt);
      fVertex->Fill(particle->Vx(), particle->Vy(), particle->Vz());

      if (processType != 92 && processType != 93)
        fdNdEtaAnalysisNSD->FillTrack(vtxMC[2], eta, pt);

      if (eventTriggered)
      {
        fdNdEtaAnalysisTr->FillTrack(vtxMC[2], eta, pt);
        if (vtxESD)
          fdNdEtaAnalysisTrVtx->FillTrack(vtxMC[2], eta, pt);
      }

      if (TMath::Abs(eta) < 0.8)
        fPartPt->Fill(pt);
    }

    fdNdEtaAnalysis->FillEvent(vtxMC[2], nAcceptedParticles);
    if (processType != 92 && processType != 93)
      fdNdEtaAnalysisNSD->FillEvent(vtxMC[2], nAcceptedParticles);

    if (eventTriggered)
    {
      fdNdEtaAnalysisTr->FillEvent(vtxMC[2], nAcceptedParticles);
      if (vtxESD)
        fdNdEtaAnalysisTrVtx->FillEvent(vtxMC[2], nAcceptedParticles);
    }

    if (eventTriggered && vtxESD)
    {
      // from tracks is only done for triggered and vertex reconstructed events

      for (Int_t i=0; i<inputCount; ++i)
      {
        Int_t label = labelArr[i];

        if (label < 0)
          continue;

        //Printf("Getting particle of track %d", label);
        TParticle* particle = stack->Particle(label);

        if (!particle)
        {
          AliDebug(AliLog::kError, Form("ERROR: Could not retrieve particle %d.", label));
          continue;
        }

        fdNdEtaAnalysisTracks->FillTrack(vtxMC[2], particle->Eta(), particle->Pt());
      } // end of track loop

      // for event count per vertex
      fdNdEtaAnalysisTracks->FillEvent(vtxMC[2], inputCount);
    }
  }

  delete[] etaArr;
  delete[] labelArr;
  delete[] ptArr;
}

void AlidNdEtaTask::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: fOutput not available");
    return;
  }

  fdNdEtaAnalysisESD = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("fdNdEtaAnalysisESD"));
  fMult = dynamic_cast<TH1F*> (fOutput->FindObject("fMult"));
  fMultVtx = dynamic_cast<TH1F*> (fOutput->FindObject("fMultVtx"));
  fEvents = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_vertex"));
  fVertexResolution = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_vertex_resolution_z"));
  for (Int_t i=0; i<3; ++i)
    fPartEta[i] = dynamic_cast<TH1F*> (fOutput->FindObject(Form("dndeta_check_%d", i)));

  if (!fdNdEtaAnalysisESD)
  {
    AliDebug(AliLog::kError, "ERROR: fdNdEtaAnalysisESD not available");
    return;
  }

  if (fMult && fMultVtx)
  {
    new TCanvas;
    fMult->Draw();
    fMultVtx->SetLineColor(2);
    fMultVtx->Draw("SAME");
  }

  if (fPartEta[0])
  {
    Int_t events1 = (Int_t) fEvents->Integral(fEvents->GetXaxis()->FindBin(-19.9), fEvents->GetXaxis()->FindBin(-0.001));
    Int_t events2 = (Int_t) fEvents->Integral(fEvents->GetXaxis()->FindBin(0.001), fEvents->GetXaxis()->FindBin(19.9));

    Printf("%d events with vertex used", events1 + events2);

    if (events1 > 0 && events2 > 0)
    {
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
  }

  if (fEvents)
  {
    new TCanvas("control3", "control3", 500, 500);
    fEvents->Draw();
  }

  TFile* fout = new TFile("analysis_esd_raw.root", "RECREATE");

  if (fdNdEtaAnalysisESD)
    fdNdEtaAnalysisESD->SaveHistograms();

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_tracks_cuts");

  if (fMult)
    fMult->Write();

  if (fMultVtx)
    fMultVtx->Write();

  for (Int_t i=0; i<3; ++i)
    if (fPartEta[i])
      fPartEta[i]->Write();

  if (fEvents)
    fEvents->Write();

  if (fVertexResolution)
    fVertexResolution->Write();

  if (fDeltaPhi)
    fDeltaPhi->Write();

  fout->Write();
  fout->Close();

  Printf("Writting result to analysis_esd_raw.root");

  if (fReadMC)
  {
    fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));
    fdNdEtaAnalysisNSD = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaNSD"));
    fdNdEtaAnalysisTr = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTr"));
    fdNdEtaAnalysisTrVtx = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTrVtx"));
    fdNdEtaAnalysisTracks = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTracks"));
    fPartPt = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_pt"));

    if (!fdNdEtaAnalysis || !fdNdEtaAnalysisTr || !fdNdEtaAnalysisTrVtx || !fPartPt)
    {
      AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p", (void*) fdNdEtaAnalysis, (void*) fPartPt));
      return;
    }

    fdNdEtaAnalysis->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisNSD->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisTr->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisTrVtx->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisTracks->Finish(0, -1, AlidNdEtaCorrection::kNone);

    Int_t events = (Int_t) fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetMeasuredHistogram()->Integral();
    fPartPt->Scale(1.0/events);
    fPartPt->Scale(1.0/fPartPt->GetBinWidth(1));

    fout = new TFile("analysis_mc.root","RECREATE");

    fdNdEtaAnalysis->SaveHistograms();
    fdNdEtaAnalysisNSD->SaveHistograms();
    fdNdEtaAnalysisTr->SaveHistograms();
    fdNdEtaAnalysisTrVtx->SaveHistograms();
    fdNdEtaAnalysisTracks->SaveHistograms();
    fPartPt->Write();

    fout->Write();
    fout->Close();

    Printf("Writting result to analysis_mc.root");

    if (fPartPt)
    {
      new TCanvas("control2", "control2", 500, 500);
      fPartPt->Draw();
    }
  }
}
