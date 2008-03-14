/* $Id$ */

#include "AlidNdEtaCorrectionTask.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TParticle.h>

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

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "AliPWG0Helper.h"
//#include "AliCorrection.h"
//#include "AliCorrectionMatrix3D.h"
#include "dNdEta/dNdEtaAnalysis.h"
#include "dNdEta/AlidNdEtaCorrection.h"
#include "AliVertexerTracks.h"

ClassImp(AlidNdEtaCorrectionTask)

AlidNdEtaCorrectionTask::AlidNdEtaCorrectionTask(const char* opt) :
  AliAnalysisTask("AlidNdEtaCorrectionTask", ""),
  fESD(0),
  fOutput(0),
  fOption(opt),
  fAnalysisMode(AliPWG0Helper::kTPC),
  fSignMode(0),
  fEsdTrackCuts(0),
  fdNdEtaCorrection(0),
  fdNdEtaAnalysisMC(0),
  fdNdEtaAnalysisESD(0),
  fPIDParticles(0),
  fPIDTracks(0),
  fVertexCorrelation(0),
  fVertexProfile(0),
  fVertexShiftNorm(0),
  fSigmaVertexTracks(0),
  fSigmaVertexPrim(0),
  fMultAll(0),
  fMultTr(0),
  fMultVtx(0)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  for (Int_t i=0; i<2; i++)
    fdNdEtaCorrectionProcessType[i] = 0;
  
  for (Int_t i=0; i<3; i++)
    fDeltaPhi[i] = 0;
}

AlidNdEtaCorrectionTask::~AlidNdEtaCorrectionTask()
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

//________________________________________________________________________
void AlidNdEtaCorrectionTask::ConnectInputData(Option_t *)
{
  // Connect ESD
  // Called once

  Printf("AlidNdEtaCorrectionTask::ConnectInputData called");

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read tree from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    //tree->SetBranchStatus("*", 0);

    tree->SetBranchStatus("fTriggerMask", 1);
    tree->SetBranchStatus("fSPDVertex*", 1);
    // PrimaryVertex

    if (fAnalysisMode == AliPWG0Helper::kSPD)
      tree->SetBranchStatus("fSPDMult*", 1);

    if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS) {
      AliESDtrackCuts::EnableNeededBranches(tree);
      tree->SetBranchStatus("fTracks.fLabel", 1);
    }

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }

  // disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);
}

void AlidNdEtaCorrectionTask::CreateOutputObjects()
{
  // create result objects and add to output list

  if (fOption.Contains("only-positive"))
  {
    Printf("INFO: Processing only positive particles.");
    fSignMode = 1;
  }
  else if (fOption.Contains("only-negative"))
  {
    Printf("INFO: Processing only negative particles.");
    fSignMode = -1;
  }

  fOutput = new TList;
  fOutput->SetOwner();

  fdNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction", fAnalysisMode);
  fOutput->Add(fdNdEtaCorrection);

  fPIDParticles = new TH1F("pid_particles", "PID of generated primary particles", 10001, -5000.5, 5000.5);
  fOutput->Add(fPIDParticles);

  fPIDTracks = new TH1F("pid_tracks", "MC PID of reconstructed tracks", 10001, -5000.5, 5000.5);
  fOutput->Add(fPIDTracks);

  fdNdEtaAnalysisMC = new dNdEtaAnalysis("dndetaMC", "dndetaMC", fAnalysisMode);
  fOutput->Add(fdNdEtaAnalysisMC);

  fdNdEtaAnalysisESD = new dNdEtaAnalysis("dndetaESD", "dndetaESD", fAnalysisMode);
  fOutput->Add(fdNdEtaAnalysisESD);

  if (fOption.Contains("process-types")) {
    fdNdEtaCorrectionProcessType[0] = new AlidNdEtaCorrection("dndeta_correction_ND", "dndeta_correction_ND", fAnalysisMode);
    fdNdEtaCorrectionProcessType[1] = new AlidNdEtaCorrection("dndeta_correction_SD", "dndeta_correction_SD", fAnalysisMode);
    fdNdEtaCorrectionProcessType[2] = new AlidNdEtaCorrection("dndeta_correction_DD", "dndeta_correction_DD", fAnalysisMode);

    fOutput->Add(fdNdEtaCorrectionProcessType[0]);
    fOutput->Add(fdNdEtaCorrectionProcessType[1]);
    fOutput->Add(fdNdEtaCorrectionProcessType[2]);
  }

  if (fOption.Contains("sigma-vertex"))
  {
    fSigmaVertexTracks = new TH1F("fSigmaVertexTracks", "fSigmaVertexTracks;Nsigma2vertex;NacceptedTracks", 60, 0.05, 6.05);
    fSigmaVertexPrim = new TH1F("fSigmaVertexPrim", "fSigmaVertexPrim;Nsigma2vertex;NacceptedPrimaries", 60, 0.05, 6.05);
    fOutput->Add(fSigmaVertexTracks);
    fOutput->Add(fSigmaVertexPrim);
    Printf("WARNING: sigma-vertex analysis enabled. This will produce weird results in the AliESDtrackCuts histograms");
  }

  fVertexCorrelation = new TH2F("fVertexCorrelation", "fVertexCorrelation;MC z-vtx;ESD z-vtx", 80, -20, 20, 80, -20, 20);
  fVertexProfile = new TProfile("fVertexProfile", "fVertexProfile;MC z-vtx;MC z-vtx - ESD z-vtx", 40, -20, 20);
  fVertexShiftNorm = new TH1F("fVertexShiftNorm", "fVertexShiftNorm;(MC z-vtx - ESD z-vtx) / #sigma_{ESD z-vtx};Entries", 200, -100, 100);
  
  fMultAll = new TH1F("fMultAll", "fMultAll", 500, -0.5, 499.5);
  fMultVtx = new TH1F("fMultVtx", "fMultVtx", 500, -0.5, 499.5);
  fMultTr = new TH1F("fMultTr", "fMultTr", 500, -0.5, 499.5);

  for (Int_t i=0; i<3; i++)
    fDeltaPhi[i] = new TH1F(Form("fDeltaPhi_%d", i), ";#Delta phi;Entries", 200, -0.5, 0.5);
}

void AlidNdEtaCorrectionTask::Exec(Option_t*)
{
  // process the event

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return;
  }

  // trigger definition
  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD->GetTriggerMask(), AliPWG0Helper::kMB2);

  // post the data already here
  PostData(0, fOutput);

  // MC info
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

  AliStack* stack = mcEvent->Stack();
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

  // get process type; NB: this only works for Pythia
  Int_t processType = AliPWG0Helper::GetPythiaEventProcessType(header);
  AliDebug(AliLog::kDebug+1, Form("Found pythia process type %d", processType));

  if (processType<0)
    AliDebug(AliLog::kError, Form("Unknown Pythia process type %d.", processType));

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // get the ESD vertex
  const AliESDVertex* vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);

  Bool_t eventVertex = kFALSE;
  Double_t vtx[3];
  if (vtxESD) 
  {
    vtxESD->GetXYZ(vtx);
    eventVertex = kTRUE;
    
    Double_t diff = vtxMC[2] - vtx[2];
    if (vtxESD->GetZRes() > 0) 
      fVertexShiftNorm->Fill(diff / vtxESD->GetZRes());
  } 
  else
    Printf("No vertex found");


  // create list of (label, eta, pt) tuples
  Int_t inputCount = 0;
  Int_t* labelArr = 0;
  Float_t* etaArr = 0;
  Float_t* ptArr = 0;
  Float_t* deltaPhiArr = 0;
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
    deltaPhiArr = new Float_t[mult->GetNumberOfTracklets()];

    // get multiplicity from ITS tracklets
    for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
    {
      //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));

      Float_t deltaPhi = mult->GetDeltaPhi(i);
      // prevent values to be shifted by 2 Pi()
      if (deltaPhi < -TMath::Pi())
        deltaPhi += TMath::Pi() * 2;
      if (deltaPhi > TMath::Pi())
        deltaPhi -= TMath::Pi() * 2;

      if (TMath::Abs(deltaPhi) > 1)
        printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), deltaPhi);

      etaArr[inputCount] = mult->GetEta(i);
      labelArr[inputCount] = mult->GetLabel(i);
      ptArr[inputCount] = 0; // no pt for tracklets
      deltaPhiArr[inputCount] = deltaPhi;
      ++inputCount;
    }
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

    Printf("Accepted %d tracks", nGoodTracks);

    labelArr = new Int_t[nGoodTracks];
    etaArr = new Float_t[nGoodTracks];
    ptArr = new Float_t[nGoodTracks];
    deltaPhiArr = new Float_t[nGoodTracks];

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
      deltaPhiArr[inputCount] = 0; // no delta phi for tracks
      ++inputCount;
    }

    delete list;
  }
  else
    return;

  if (eventTriggered && eventVertex)
  {
    fVertexCorrelation->Fill(vtxMC[2], vtx[2]);
    fVertexProfile->Fill(vtxMC[2], vtxMC[2] - vtx[2]);
  }


  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();
  Int_t nAccepted = 0;

  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    //Printf("Getting particle %d", iMc);
    TParticle* particle = stack->Particle(iMc);

    if (!particle)
      continue;

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    if (SignOK(particle->GetPDG()) == kFALSE)
      continue;

    Float_t eta = particle->Eta();
    Float_t pt = particle->Pt();

    fdNdEtaCorrection->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, eventVertex, processType);

    if (fdNdEtaCorrectionProcessType[0])
    {
      // non diffractive
      if (processType!=92 && processType!=93 && processType!=94)
        fdNdEtaCorrectionProcessType[0]->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, eventVertex, processType);

      // single diffractive
      if (processType==92 || processType==93)
        fdNdEtaCorrectionProcessType[1]->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, eventVertex, processType);

      // double diffractive
      if (processType==94)
        fdNdEtaCorrectionProcessType[2]->FillMCParticle(vtxMC[2], eta, pt, eventTriggered, eventVertex, processType);
    }

    if (eventTriggered)
      if (eventVertex)
        fdNdEtaAnalysisMC->FillTrack(vtxMC[2], eta, pt);

    if (TMath::Abs(eta) < 1 && pt > 0.2)
      nAccepted++;
  }

  fMultAll->Fill(nAccepted);
  if (eventTriggered) {
    fMultTr->Fill(nAccepted);
    if (eventVertex)
      fMultVtx->Fill(nAccepted);
  }

  for (Int_t i=0; i<inputCount; ++i)
  {
    Int_t label = labelArr[i];

    if (label < 0)
    {
      Printf("WARNING: cannot find corresponding mc part for track(let) %d with label %d.", i, label);
      fDeltaPhi[2]->Fill(deltaPhiArr[i]);
      continue;
    }

    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label));
      continue;
    }

    if (SignOK(particle->GetPDG()) == kFALSE)
        continue;

    if (eventTriggered && eventVertex)
    {
      fdNdEtaCorrection->FillTrackedParticle(vtxMC[2], particle->Eta(), particle->Pt());
      fdNdEtaAnalysisESD->FillTrack(vtxMC[2], particle->Eta(), particle->Pt());
      if (particle->Pt() > 0.1 && particle->Pt() < 0.2)
      {
        fPIDTracks->Fill(particle->GetPdgCode());
      }

      if (fdNdEtaCorrectionProcessType[0])
      {
        // non diffractive
        if (processType!=92 && processType!=93 && processType!=94)
          fdNdEtaCorrectionProcessType[0]->FillTrackedParticle(vtxMC[2], particle->Eta(), particle->Pt());

        // single diffractive
        if (processType==92 || processType==93)
          fdNdEtaCorrectionProcessType[1]->FillTrackedParticle(vtxMC[2], particle->Eta(), particle->Pt());

        // double diffractive
        if (processType==94)
          fdNdEtaCorrectionProcessType[2]->FillTrackedParticle(vtxMC[2], particle->Eta(), particle->Pt());
      }

      if (stack->IsPhysicalPrimary(label))
      {
        fDeltaPhi[0]->Fill(deltaPhiArr[i]);
      }
      else
        fDeltaPhi[1]->Fill(deltaPhiArr[i]);
    }
  }

  if (eventTriggered && eventVertex)
  {
    fdNdEtaAnalysisMC->FillEvent(vtxMC[2], inputCount);
    fdNdEtaAnalysisESD->FillEvent(vtxMC[2], inputCount);
  }

   // stuff regarding the vertex reco correction and trigger bias correction
  fdNdEtaCorrection->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);

  if (fdNdEtaCorrectionProcessType[0])
  {
    // non diffractive
    if (processType!=92 && processType!=93 && processType!=94)
      fdNdEtaCorrectionProcessType[0]->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);

    // single diffractive
    if (processType==92 || processType==93)
      fdNdEtaCorrectionProcessType[1]->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);

    // double diffractive
    if (processType==94)
      fdNdEtaCorrectionProcessType[2]->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);
  }

  delete[] etaArr;
  delete[] labelArr;
  delete[] ptArr;
  delete[] deltaPhiArr;

  // fills the fSigmaVertex histogram (systematic study)
  if (fSigmaVertexTracks)
  {
    // save the old value
    Float_t oldSigmaVertex = fEsdTrackCuts->GetMinNsigmaToVertex();

    // set to maximum
    fEsdTrackCuts->SetMinNsigmaToVertex(6);

    TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);
    Int_t nGoodTracks = list->GetEntries();

    // loop over esd tracks
    for (Int_t i=0; i<nGoodTracks; i++)
    {
      AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (list->At(i));
      if (!esdTrack)
      {
        AliError(Form("ERROR: Could not retrieve track %d.", i));
        continue;
      }

      Float_t sigma2Vertex = fEsdTrackCuts->GetSigmaToVertex(esdTrack);

      for (Double_t nSigma = 0.1; nSigma < 6.05; nSigma += 0.1)
      {
        if (sigma2Vertex < nSigma)
        {
          fSigmaVertexTracks->Fill(nSigma);

          Int_t label = TMath::Abs(esdTrack->GetLabel());
          TParticle* particle = stack->Particle(label);
          if (!particle || label >= nPrim)
            continue;

          if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim))
            fSigmaVertexPrim->Fill(nSigma);
        }
      }
    }

    delete list;
    list = 0;

    // set back the old value
    fEsdTrackCuts->SetMinNsigmaToVertex(oldSigmaVertex);
  }
}

void AlidNdEtaCorrectionTask::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: fOutput not available");
    return;
  }

  fdNdEtaCorrection = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction"));
  fdNdEtaAnalysisMC = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaMC"));
  fdNdEtaAnalysisESD = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaESD"));
  if (!fdNdEtaCorrection || !fdNdEtaAnalysisMC || !fdNdEtaAnalysisESD)
  {
    AliDebug(AliLog::kError, "Could not read object from output list");
    return;
  }

  fdNdEtaCorrection->Finish();

  TString fileName;
  fileName.Form("correction_map%s.root", fOption.Data());
  TFile* fout = new TFile(fileName, "RECREATE");

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_track_cuts");
  fdNdEtaCorrection->SaveHistograms();
  fdNdEtaAnalysisMC->SaveHistograms();
  fdNdEtaAnalysisESD->SaveHistograms();

  fdNdEtaCorrectionProcessType[0] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_ND"));
  fdNdEtaCorrectionProcessType[1] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_SD"));
  fdNdEtaCorrectionProcessType[2] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_DD"));
  for (Int_t i=0; i<3; ++i)
    if (fdNdEtaCorrectionProcessType[i])
      fdNdEtaCorrectionProcessType[i]->SaveHistograms();

  fSigmaVertexTracks = dynamic_cast<TH1F*> (fOutput->FindObject("fSigmaVertexTracks"));
  if (fSigmaVertexTracks)
    fSigmaVertexTracks->Write();
  fSigmaVertexPrim = dynamic_cast<TH1F*> (fOutput->FindObject("fSigmaVertexPrim"));
  if (fSigmaVertexPrim)
    fSigmaVertexPrim->Write();

  if (fVertexCorrelation)
    fVertexCorrelation->Write();
  if (fVertexProfile)
    fVertexProfile->Write();
  if (fVertexShiftNorm)
    fVertexShiftNorm->Write();
  if (fMultAll)
    fMultAll->Write();

  if (fMultTr)
    fMultTr->Write();
  
  if (fMultVtx)
    fMultVtx->Write();

  for (Int_t i=0; i<3; ++i)
    if (fDeltaPhi[i])
      fDeltaPhi[i]->Write();

  fout->Write();
  fout->Close();

  //fdNdEtaCorrection->DrawHistograms();

  Printf("Writting result to %s", fileName.Data());

  if (fPIDParticles && fPIDTracks)
  {
    new TCanvas("pidcanvas", "pidcanvas", 500, 500);

    fPIDParticles->Draw();
    fPIDTracks->SetLineColor(2);
    fPIDTracks->Draw("SAME");

    TDatabasePDG* pdgDB = new TDatabasePDG;

    for (Int_t i=0; i <= fPIDParticles->GetNbinsX()+1; ++i)
      if (fPIDParticles->GetBinContent(i) > 0)
        printf("PDG = %d (%s): generated: %d, reconstructed: %d, ratio: %f\n", (Int_t) fPIDParticles->GetBinCenter(i), pdgDB->GetParticle((Int_t) fPIDParticles->GetBinCenter(i))->GetName(), (Int_t) fPIDParticles->GetBinContent(i), (Int_t) fPIDTracks->GetBinContent(i), ((fPIDTracks->GetBinContent(i) > 0) ? fPIDParticles->GetBinContent(i) / fPIDTracks->GetBinContent(i) : -1));

    delete pdgDB;
    pdgDB = 0;
  }
}

Bool_t AlidNdEtaCorrectionTask::SignOK(TParticlePDG* particle)
{
  // returns if a particle with this sign should be counted
  // this is determined by the value of fSignMode, which should have the same sign
  // as the charge
  // if fSignMode is 0 all particles are counted

  if (fSignMode == 0)
    return kTRUE;

  if (!particle)
  {
    printf("WARNING: not counting a particle that does not have a pdg particle\n");
    return kFALSE;
  }

  Double_t charge = particle->Charge();

  if (fSignMode > 0)
    if (charge < 0)
      return kFALSE;

  if (fSignMode < 0)
    if (charge > 0)
      return kFALSE;

  return kTRUE;
}

