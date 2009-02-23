/* $Id$ */

#include "AlidNdEtaCorrectionTask.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TParticlePDG.h>

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
#include "dNdEta/dNdEtaAnalysis.h"
#include "dNdEta/AlidNdEtaCorrection.h"

ClassImp(AlidNdEtaCorrectionTask)

AlidNdEtaCorrectionTask::AlidNdEtaCorrectionTask() :
  AliAnalysisTask(),
  fESD(0),
  fOutput(0),
  fOption(),
  fAnalysisMode(AliPWG0Helper::kTPC),
  fTrigger(AliPWG0Helper::kMB1),
  fFillPhi(kFALSE),
  fDeltaPhiCut(-1),
  fSignMode(0),
  fOnlyPrimaries(kFALSE),
  fStatError(0),
  fEsdTrackCuts(0),
  fdNdEtaCorrection(0),
  fdNdEtaAnalysisMC(0),
  fdNdEtaAnalysisESD(0),
  fPIDParticles(0),
  fPIDTracks(0),
  fVertexCorrelation(0),
  fVertexCorrelationShift(0),
  fVertexProfile(0),
  fVertexShift(0),
  fVertexShiftNorm(0),
  fEtaCorrelation(0),
  fEtaCorrelationShift(0),
  fEtaProfile(0),
  fEtaResolution(0),
  fDeltaPhiCorrelation(0),
  fpTResolution(0),
  fEsdTrackCutsPrim(0),
  fEsdTrackCutsSec(0),
  fTemp1(0),
  fTemp2(0),
  fMultAll(0),
  fMultTr(0),
  fMultVtx(0),
  fEventStats(0)
{
  //
  // Constructor. Initialization of pointers
  //

  for (Int_t i=0; i<4; i++)
    fdNdEtaCorrectionSpecial[i] = 0;
  
  for (Int_t i=0; i<8; i++)
    fDeltaPhi[i] = 0;
}

AlidNdEtaCorrectionTask::AlidNdEtaCorrectionTask(const char* opt) :
  AliAnalysisTask("AlidNdEtaCorrectionTask", ""),
  fESD(0),
  fOutput(0),
  fOption(opt),
  fAnalysisMode(AliPWG0Helper::kTPC),
  fTrigger(AliPWG0Helper::kMB1),
  fFillPhi(kFALSE),
  fDeltaPhiCut(0),
  fSignMode(0),
  fOnlyPrimaries(kFALSE),
  fStatError(0),
  fEsdTrackCuts(0),
  fdNdEtaCorrection(0),
  fdNdEtaAnalysisMC(0),
  fdNdEtaAnalysisESD(0),
  fPIDParticles(0),
  fPIDTracks(0),
  fVertexCorrelation(0),
  fVertexCorrelationShift(0),
  fVertexProfile(0),
  fVertexShift(0),
  fVertexShiftNorm(0),
  fEtaCorrelation(0),
  fEtaCorrelationShift(0),
  fEtaProfile(0),
  fEtaResolution(0),
  fDeltaPhiCorrelation(0),
  fpTResolution(0),
  fEsdTrackCutsPrim(0),
  fEsdTrackCutsSec(0),
  fTemp1(0),
  fTemp2(0),
  fMultAll(0),
  fMultTr(0),
  fMultVtx(0),
  fEventStats(0)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  for (Int_t i=0; i<4; i++)
    fdNdEtaCorrectionSpecial[i] = 0;
  
  for (Int_t i=0; i<8; i++)
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

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
  } else {
    fESD = esdH->GetEvent();

    // Enable only the needed branches
    esdH->SetActiveBranches("AliESDHeader Vertex");

    if (fAnalysisMode == AliPWG0Helper::kSPD)
      esdH->SetActiveBranches("AliESDHeader Vertex AliMultiplicity");

    if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS) {
      esdH->SetActiveBranches("AliESDHeader Vertex Tracks");
    }
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
  
  if (fOption.Contains("stat-error-1"))
  {
    Printf("INFO: Evaluation statistical errors. Mode: 1.");
    fStatError = 1;
  }
  else if (fOption.Contains("stat-error-2"))
  {
    Printf("INFO: Evaluation statistical errors. Mode: 2.");
    fStatError = 2;
  }

  fOutput = new TList;
  fOutput->SetOwner();

  fdNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction", fAnalysisMode);
  fOutput->Add(fdNdEtaCorrection);

  fPIDParticles = new TH1F("fPIDParticles", "PID of generated primary particles", 10001, -5000.5, 5000.5);
  fOutput->Add(fPIDParticles);

  fPIDTracks = new TH1F("fPIDTracks", "MC PID of reconstructed tracks", 10001, -5000.5, 5000.5);
  fOutput->Add(fPIDTracks);

  fdNdEtaAnalysisMC = new dNdEtaAnalysis("dndetaMC", "dndetaMC", fAnalysisMode);
  fOutput->Add(fdNdEtaAnalysisMC);

  fdNdEtaAnalysisESD = new dNdEtaAnalysis("dndetaESD", "dndetaESD", fAnalysisMode);
  fOutput->Add(fdNdEtaAnalysisESD);

  if (fEsdTrackCuts)
  {
    fEsdTrackCutsPrim = dynamic_cast<AliESDtrackCuts*> (fEsdTrackCuts->Clone("fEsdTrackCutsPrim"));
    fEsdTrackCutsSec  = dynamic_cast<AliESDtrackCuts*> (fEsdTrackCuts->Clone("fEsdTrackCutsSec"));
    fOutput->Add(fEsdTrackCutsPrim);
    fOutput->Add(fEsdTrackCutsSec);
  }

  if (fOption.Contains("process-types")) {
    fdNdEtaCorrectionSpecial[0] = new AlidNdEtaCorrection("dndeta_correction_ND", "dndeta_correction_ND", fAnalysisMode);
    fdNdEtaCorrectionSpecial[1] = new AlidNdEtaCorrection("dndeta_correction_SD", "dndeta_correction_SD", fAnalysisMode);
    fdNdEtaCorrectionSpecial[2] = new AlidNdEtaCorrection("dndeta_correction_DD", "dndeta_correction_DD", fAnalysisMode);

    fOutput->Add(fdNdEtaCorrectionSpecial[0]);
    fOutput->Add(fdNdEtaCorrectionSpecial[1]);
    fOutput->Add(fdNdEtaCorrectionSpecial[2]);
  }
  
  if (fOption.Contains("particle-species")) {
    fdNdEtaCorrectionSpecial[0] = new AlidNdEtaCorrection("dndeta_correction_pi", "dndeta_correction_pi", fAnalysisMode);
    fdNdEtaCorrectionSpecial[1] = new AlidNdEtaCorrection("dndeta_correction_K", "dndeta_correction_K", fAnalysisMode);
    fdNdEtaCorrectionSpecial[2] = new AlidNdEtaCorrection("dndeta_correction_p", "dndeta_correction_p", fAnalysisMode);
    fdNdEtaCorrectionSpecial[3] = new AlidNdEtaCorrection("dndeta_correction_other", "dndeta_correction_other", fAnalysisMode);

    for (Int_t i=0; i<4; i++)
      fOutput->Add(fdNdEtaCorrectionSpecial[i]);
  }

  
  /*
  fTemp1 = new TH2F("fTemp1", "fTemp1", 200, -0.08, 0.08, 200, -0.08, 0.08);
  fOutput->Add(fTemp1);
  fTemp2 = new TH1F("fTemp2", "fTemp2", 2000, -5, 5);
  fOutput->Add(fTemp2);
  */

  fVertexCorrelation = new TH2F("fVertexCorrelation", "fVertexCorrelation;MC z-vtx;ESD z-vtx", 120, -30, 30, 120, -30, 30);
  fOutput->Add(fVertexCorrelation);
  fVertexCorrelationShift = new TH2F("fVertexCorrelationShift", "fVertexCorrelationShift;MC z-vtx;MC z-vtx - ESD z-vtx", 120, -30, 30, 100, -1, 1);
  fOutput->Add(fVertexCorrelationShift);
  fVertexProfile = new TProfile("fVertexProfile", "fVertexProfile;MC z-vtx;MC z-vtx - ESD z-vtx", 40, -20, 20);
  fOutput->Add(fVertexProfile);
  fVertexShift = new TH1F("fVertexShift", "fVertexShift;(MC z-vtx - ESD z-vtx);Entries", 201, -2, 2);
  fOutput->Add(fVertexShift);
  fVertexShiftNorm = new TH1F("fVertexShiftNorm", "fVertexShiftNorm;(MC z-vtx - ESD z-vtx) / #sigma_{ESD z-vtx};Entries", 200, -100, 100);
  fOutput->Add(fVertexShiftNorm);

  fEtaCorrelation = new TH2F("fEtaCorrelation", "fEtaCorrelation;MC #eta;ESD #eta", 120, -3, 3, 120, -3, 3);
  fOutput->Add(fEtaCorrelation);
  fEtaCorrelationShift = new TH2F("fEtaCorrelationShift", "fEtaCorrelationShift;MC #eta;MC #eta - ESD #eta", 120, -3, 3, 100, -0.1, 0.1);
  fOutput->Add(fEtaCorrelationShift);
  fEtaProfile = new TProfile("fEtaProfile", "fEtaProfile;MC #eta;MC #eta - ESD #eta", 120, -3, 3);
  fOutput->Add(fEtaProfile);
  fEtaResolution = new TH1F("fEtaResolution", "fEtaResolution;MC #eta - ESD #eta", 201, -0.2, 0.2);
  fOutput->Add(fEtaResolution);

  fpTResolution = new TH2F("fpTResolution", ";MC p_{T} (GeV/c);(MC p_{T} - ESD p_{T}) / MC p_{T}", 160, 0, 20, 201, -0.2, 0.2);
  fOutput->Add(fpTResolution);

  fMultAll = new TH1F("fMultAll", "fMultAll", 500, -0.5, 499.5);
  fOutput->Add(fMultAll);
  fMultVtx = new TH1F("fMultVtx", "fMultVtx", 500, -0.5, 499.5);
  fOutput->Add(fMultVtx);
  fMultTr = new TH1F("fMultTr", "fMultTr", 500, -0.5, 499.5);
  fOutput->Add(fMultTr);

  for (Int_t i=0; i<8; i++)
  {
    fDeltaPhi[i] = new TH2F(Form("fDeltaPhi_%d", i), ";#Delta phi;#phi;Entries", 2000, -0.1, 0.1, 18 * 5, 0, TMath::TwoPi());
    fOutput->Add(fDeltaPhi[i]);
  }

  fEventStats = new TH2F("fEventStats", "fEventStats;event type;status;count", 106, -5.5, 100.5, 4, -0.5, 3.5);
  fOutput->Add(fEventStats);
  fEventStats->GetXaxis()->SetBinLabel(1, "INEL"); // x = -5
  fEventStats->GetXaxis()->SetBinLabel(2, "NSD");  // x = -4
  fEventStats->GetXaxis()->SetBinLabel(3, "ND");   // x = -3
  fEventStats->GetXaxis()->SetBinLabel(4, "SD");   // x = -2
  fEventStats->GetXaxis()->SetBinLabel(5, "DD");   // x = -1
  
  for (Int_t i=0; i<100; i++)
    fEventStats->GetXaxis()->SetBinLabel(7+i, Form("%d", i)); 

  fEventStats->GetYaxis()->SetBinLabel(1, "nothing");
  fEventStats->GetYaxis()->SetBinLabel(2, "trg");
  fEventStats->GetYaxis()->SetBinLabel(3, "vtx");
  fEventStats->GetYaxis()->SetBinLabel(4, "trgvtx");

  if (fEsdTrackCuts)
  {
    fEsdTrackCuts->SetName("fEsdTrackCuts");
    fOutput->Add(fEsdTrackCuts);
  }
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

  if (fOnlyPrimaries)
    Printf("WARNING: Processing only primaries. For systematical studies only!");
    
  if (fStatError > 0)
    Printf("WARNING: Statistical error evaluation active!");

  // trigger definition
  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD, fTrigger);

  if (!eventTriggered)
    Printf("No trigger");

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

  // get process type;
  Int_t processType = AliPWG0Helper::GetEventProcessType(header);
  AliDebug(AliLog::kDebug+1, Form("Found process type %d", processType));

  if (processType == AliPWG0Helper::kInvalidProcess)
    AliDebug(AliLog::kError, "Unknown process type.");

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // get the ESD vertex
  const AliESDVertex* vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);
  Bool_t eventVertex = kFALSE;
  if (vtxESD)
  {
    Double_t vtx[3];
    vtxESD->GetXYZ(vtx);

    Double_t diff = vtxMC[2] - vtx[2];
    fVertexShift->Fill(diff);
    if (vtxESD->GetZRes() > 0)
        fVertexShiftNorm->Fill(diff / vtxESD->GetZRes());

    if (!AliPWG0Helper::TestVertex(vtxESD, fAnalysisMode))
    {
        vtxESD = 0;
    }
    else
    {
      eventVertex = kTRUE;

      if (eventTriggered)
      {
        fVertexCorrelation->Fill(vtxMC[2], vtx[2]);
        fVertexCorrelationShift->Fill(vtxMC[2], vtxMC[2] - vtx[2]);
        fVertexProfile->Fill(vtxMC[2], vtxMC[2] - vtx[2]);
      }
    }
  }

  // fill process type
  Int_t biny = (Int_t) eventTriggered + 2 * (Int_t) eventVertex;
  // INEL
  fEventStats->Fill(-5, biny);
  // NSD
  if (processType != AliPWG0Helper::kSD)
    fEventStats->Fill(-4, biny);
  // SD, ND, DD
  if (processType == AliPWG0Helper::kND)
    fEventStats->Fill(-3, biny);
  if (processType == AliPWG0Helper::kSD)
    fEventStats->Fill(-2, biny);
  if (processType == AliPWG0Helper::kDD)
    fEventStats->Fill(-1, biny);
    
  // create list of (label, eta, pt) tuples
  Int_t inputCount = 0;
  Int_t* labelArr = 0;
  Int_t* labelArr2 = 0; // only for case of SPD
  Float_t* etaArr = 0;
  Float_t* thirdDimArr = 0;
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
    labelArr2 = new Int_t[mult->GetNumberOfTracklets()];
    etaArr = new Float_t[mult->GetNumberOfTracklets()];
    thirdDimArr = new Float_t[mult->GetNumberOfTracklets()];
    deltaPhiArr = new Float_t[mult->GetNumberOfTracklets()];

    // get multiplicity from SPD tracklets
    for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
    {
      //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));

      Float_t phi = mult->GetPhi(i);
      if (phi < 0)
        phi += TMath::Pi() * 2;
      Float_t deltaPhi = mult->GetDeltaPhi(i);

      if (TMath::Abs(deltaPhi) > 1)
        printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), deltaPhi);

      if (fOnlyPrimaries)
        if (mult->GetLabel(i, 0) < 0 || mult->GetLabel(i, 0) != mult->GetLabel(i, 1) || !stack->IsPhysicalPrimary(mult->GetLabel(i, 0)))
          continue;

      if (fDeltaPhiCut > 0 && TMath::Abs(deltaPhi) > fDeltaPhiCut)
        continue;
      
      etaArr[inputCount] = mult->GetEta(i);
      labelArr[inputCount] = mult->GetLabel(i, 0);
      labelArr2[inputCount] = mult->GetLabel(i, 1);
      thirdDimArr[inputCount] = phi;
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

    if (vtxESD)
    {
      // get multiplicity from ESD tracks
      TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD, (fAnalysisMode == AliPWG0Helper::kTPC));
      Int_t nGoodTracks = list->GetEntries();
  
      Printf("Accepted %d tracks", nGoodTracks);
  
      labelArr = new Int_t[nGoodTracks];
      labelArr2 = new Int_t[nGoodTracks];
      etaArr = new Float_t[nGoodTracks];
      thirdDimArr = new Float_t[nGoodTracks];
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
        
        // TODO fOnlyPrimaries not implemented for TPC
  
        etaArr[inputCount] = esdTrack->Eta();
        labelArr[inputCount] = TMath::Abs(esdTrack->GetLabel());
        labelArr2[inputCount] = labelArr[inputCount]; // no second label for tracks
        thirdDimArr[inputCount] = esdTrack->Pt();
        deltaPhiArr[inputCount] = 0; // no delta phi for tracks
        ++inputCount;
      }

      delete list;

      if (eventTriggered)
      {
        // collect values for primaries and secondaries
        for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++)
        {
          AliESDtrack* track = 0;
  
          if (fAnalysisMode == AliPWG0Helper::kTPC)
            track = AliESDtrackCuts::GetTPCOnlyTrack(fESD, iTrack);
          else if (fAnalysisMode == AliPWG0Helper::kTPCITS)
            track = fESD->GetTrack(iTrack);
  
          if (!track)
            continue;
  
          Int_t label = TMath::Abs(track->GetLabel());
          if (!stack->Particle(label) || !stack->Particle(label)->GetPDG())
          {
            Printf("WARNING: No particle for %d", label);
            if (stack->Particle(label))
              stack->Particle(label)->Print();
            continue;
          }
  
          if (stack->Particle(label)->GetPDG()->Charge() == 0)
            continue;
  
          if (TMath::Abs(track->Eta()) < 1)
          {
            if (stack->IsPhysicalPrimary(label))
            {
              // primary
              if (fEsdTrackCutsPrim->AcceptTrack(track)) 
              {
                if (AliESDtrackCuts::GetSigmaToVertex(track) > 900)
                {
                  Printf("Track %d has nsigma of %f. Printing track and vertex...", iTrack, AliESDtrackCuts::GetSigmaToVertex(track));
                  Float_t b[2];
                  Float_t r[3];
                  track->GetImpactParameters(b, r);
                  Printf("Impact parameter %f %f and resolution: %f %f %f", b[0], b[1], r[0], r[1], r[2]);
                  track->Print("");
                  if (vtxESD)
                    vtxESD->Print();
                }
              }
            }
            else
            {
              // secondary
              fEsdTrackCutsSec->AcceptTrack(track);
            }
          }
  
          // TODO mem leak in the continue statements above
          if (fAnalysisMode == AliPWG0Helper::kTPC)
            delete track;
        }
      }
    }
  }
  else
    return;

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
    {
      //if (TMath::Abs(particle->GetPdgCode()) > 3000 && TMath::Abs(particle->Eta()) < 1.0)
      //  fPIDParticles->Fill(particle->GetPdgCode());
      continue;
    }

    if (SignOK(particle->GetPDG()) == kFALSE)
      continue;

    if (fPIDParticles && TMath::Abs(particle->Eta()) < 1.0)
      fPIDParticles->Fill(particle->GetPdgCode());

    Float_t eta = particle->Eta();
    
    Float_t thirdDim = -1;
    if (fAnalysisMode == AliPWG0Helper::kSPD)
    {
      if (fFillPhi)
      {
        thirdDim = particle->Phi();
      }
      else
        thirdDim = inputCount;
    }
    else
      thirdDim = particle->Pt();

    // calculate y
    //Float_t y = 0.5 * TMath::Log((particle->Energy() + particle->Pz()) / (particle->Energy() - particle->Pz()));
    //fTemp1->Fill(eta);
    //fTemp2->Fill(y);

    fdNdEtaCorrection->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType);

    if (fOption.Contains("process-types"))
    {
      // non diffractive
      if (processType==AliPWG0Helper::kND)
        fdNdEtaCorrectionSpecial[0]->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType);

      // single diffractive
      if (processType==AliPWG0Helper::kSD)
        fdNdEtaCorrectionSpecial[1]->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType);

      // double diffractive
      if (processType==AliPWG0Helper::kDD)
        fdNdEtaCorrectionSpecial[2]->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType);
    }
    
    if (fOption.Contains("particle-species"))
    {
      Int_t id = -1;
      switch (TMath::Abs(particle->GetPdgCode()))
      {
        case 211:   id = 0; break;
        case 321:   id = 1; break;
        case 2212:  id = 2; break;
        default:    id = 3; break;
      }
      fdNdEtaCorrectionSpecial[id]->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType);
    }

    if (eventTriggered)
      if (eventVertex)
        fdNdEtaAnalysisMC->FillTrack(vtxMC[2], eta, thirdDim);

    // TODO this value might be needed lower for the SPD study (only used for control histograms anyway)
    if (TMath::Abs(eta) < 1.5) // && pt > 0.2)
      nAccepted++;
  }

  fEventStats->Fill(AliPWG0Helper::GetLastProcessType(), biny);

  fMultAll->Fill(nAccepted);
  if (eventTriggered) {
    fMultTr->Fill(nAccepted);
    if (eventVertex)
      fMultVtx->Fill(nAccepted);
  }

  Int_t processed = 0;
  
  Bool_t* primCount = 0;
  if (fStatError > 0)
  {
    primCount = new Bool_t[nPrim];
    for (Int_t i=0; i<nPrim; ++i)
      primCount[i] = kFALSE;
  }

  for (Int_t i=0; i<inputCount; ++i)
  {
    Int_t label = labelArr[i];
    Int_t label2 = labelArr2[i];

    if (label < 0)
    {
      Printf("WARNING: cannot find corresponding mc particle for track(let) %d with label %d.", i, label);
      continue;
    }

    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label));
      continue;
    }

    fPIDTracks->Fill(particle->GetPdgCode());
    
    // find particle that is filled in the correction map
    // this should be the particle that has been reconstructed
    // for TPC this is clear and always identified by the track's label
    // for SPD the labels might not be identical. In this case the particle closest to the beam line that is a primary is taken: 
    //  L1 L2 (P = primary, S = secondary)
    //   P P' : --> P
    //   P S  : --> P
    //   S P  : --> P
    //   S S' : --> S

    if (label != label2 && label2 >= 0)
    {
      if (stack->IsPhysicalPrimary(label) == kFALSE && stack->IsPhysicalPrimary(label2))
      {
        particle = stack->Particle(label2);
        if (!particle)
        {
          AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label2));
          continue;
        }
      }
    }

    if (eventTriggered && eventVertex)
    {
      if (SignOK(particle->GetPDG()) == kFALSE)
        continue;

      processed++;

      // resolutions
      fEtaResolution->Fill(particle->Eta() - etaArr[i]);

      if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS)
        if (TMath::Abs(particle->Eta() < 0.9) && particle->Pt() > 0)
          fpTResolution->Fill(particle->Pt(), (particle->Pt() - thirdDimArr[i]) / particle->Pt());

      Float_t eta = -999;
      Float_t thirdDim = -1;

      Bool_t firstIsPrim = stack->IsPhysicalPrimary(label);
      // in case of primary the MC values are filled, otherwise (background) the reconstructed values
      if (label == label2 && firstIsPrim)
      {
        eta = particle->Eta();
        
        if (fAnalysisMode == AliPWG0Helper::kSPD)
        {
          if (fFillPhi)
          {
            thirdDim = particle->Phi();
          }
          else
            thirdDim = inputCount;
        }
        else
          thirdDim = particle->Pt();
      }
      else
      {
        if (fAnalysisMode == AliPWG0Helper::kSPD && !fFillPhi)
        {
          thirdDim = inputCount;
        }
        else
          thirdDim = thirdDimArr[i];

        eta = etaArr[i];
      }

      Bool_t fillTrack = kTRUE;

      // statistical error evaluation active?
      if (fStatError > 0)
      {
        Bool_t statErrorDecision = kFALSE;
        
        // primary and not yet counted
        if (label == label2 && firstIsPrim && !primCount[label])
        {
          statErrorDecision = kTRUE;
          primCount[label] = kTRUE;
        }
  
        // in case of 1 we count only unique primaries, in case of 2 all the rest
        if (fStatError == 1)
        {
          fillTrack = statErrorDecision;
        }
        else if (fStatError == 2)
          fillTrack = !statErrorDecision;
      }

      if (fillTrack)
        fdNdEtaCorrection->FillTrackedParticle(vtxMC[2], eta, thirdDim);

      // eta comparison for tracklets with the same label (others are background)
      if (label == label2)
      {
        fEtaProfile->Fill(particle->Eta(), particle->Eta() - etaArr[i]);
        fEtaCorrelation->Fill(etaArr[i], particle->Eta());
        fEtaCorrelationShift->Fill(particle->Eta(), particle->Eta() - etaArr[i]);
      }

      fdNdEtaAnalysisESD->FillTrack(vtxMC[2], particle->Eta(), thirdDim);

      if (fOption.Contains("process-types"))
      {
        // non diffractive
        if (processType == AliPWG0Helper::kND)
          fdNdEtaCorrectionSpecial[0]->FillTrackedParticle(vtxMC[2], eta, thirdDim);

        // single diffractive
        if (processType == AliPWG0Helper::kSD)
          fdNdEtaCorrectionSpecial[1]->FillTrackedParticle(vtxMC[2], eta, thirdDim);

        // double diffractive
        if (processType == AliPWG0Helper::kDD)
          fdNdEtaCorrectionSpecial[2]->FillTrackedParticle(vtxMC[2], eta, thirdDim);
      }

      if (fOption.Contains("particle-species"))
      {
        // find mother first
        TParticle* mother = AliPWG0Helper::FindPrimaryMother(stack, label);
        
        Int_t id = -1;
        switch (TMath::Abs(mother->GetPdgCode()))
        {
          case 211:   id = 0; break;
          case 321:   id = 1; break;
          case 2212:  id = 2; break;
          default:    id = 3; break;
        }
        fdNdEtaCorrectionSpecial[id]->FillTrackedParticle(vtxMC[2], eta, thirdDim);
      }

      // control histograms
      Int_t hist = -1;
      if (label == label2)
      {
        if (firstIsPrim)
        {
          hist = 0;
        }
        else
          hist = 1; 
      }
      else if (label2 >= 0)
      {
        Bool_t secondIsPrim = stack->IsPhysicalPrimary(label2);
        if (firstIsPrim && secondIsPrim)
        {
          hist = 2;
        }
        else if (firstIsPrim && !secondIsPrim)
        {
          hist = 3;

          // check if secondary is caused by the primary or it is a fake combination
          //Printf("PS case --> Label 1 is %d, label 2 is %d", label, label2);
          Int_t mother = label2;
          while (stack->Particle(mother) && stack->Particle(mother)->GetMother(0) >= 0)
          {
            //Printf("  %d created by %d, %d", mother, stack->Particle(mother)->GetMother(0), stack->Particle(mother)->GetMother(1));
            if (stack->Particle(mother)->GetMother(0) == label)
            {
              /*Printf("The untraceable decay was:");
              Printf("   from:");
              particle->Print();
              Printf("   to (status code %d):", stack->Particle(mother)->GetStatusCode());
              stack->Particle(mother)->Print();*/
              hist = 4;
            }
            mother = stack->Particle(mother)->GetMother(0);
          }
        }
        else if (!firstIsPrim && secondIsPrim)
        {
          hist = 5;
        }
        else if (!firstIsPrim && !secondIsPrim)
        {
          hist = 6;
        }

      }
      else
        hist = 7;
      
      fDeltaPhi[hist]->Fill(deltaPhiArr[i], thirdDimArr[i]);
    }
  }
  
  if (primCount)
  {
    delete[] primCount;
    primCount = 0;
  }

  if (processed < inputCount)
    Printf("Only %d out of %d track(let)s were used", processed, inputCount); 

  if (eventTriggered && eventVertex)
  {
    fdNdEtaAnalysisMC->FillEvent(vtxMC[2], inputCount);
    fdNdEtaAnalysisESD->FillEvent(vtxMC[2], inputCount);
  }

   // stuff regarding the vertex reco correction and trigger bias correction
  fdNdEtaCorrection->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);

  if (fOption.Contains("process-types"))
  {
    // non diffractive
    if (processType == AliPWG0Helper::kND )
      fdNdEtaCorrectionSpecial[0]->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);

    // single diffractive
    if (processType == AliPWG0Helper::kSD)
      fdNdEtaCorrectionSpecial[1]->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);

    // double diffractive
    if (processType == AliPWG0Helper::kDD)
      fdNdEtaCorrectionSpecial[2]->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);
  }
  
  if (fOption.Contains("particle-species"))
    for (Int_t id=0; id<4; id++)
      fdNdEtaCorrectionSpecial[id]->FillEvent(vtxMC[2], inputCount, eventTriggered, eventVertex, processType);

  if (etaArr)
    delete[] etaArr;
  if (labelArr)
    delete[] labelArr;
  if (labelArr2)
    delete[] labelArr2;
  if (thirdDimArr)
    delete[] thirdDimArr;
  if (deltaPhiArr)
    delete[] deltaPhiArr;
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

  fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCuts"));
  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_track_cuts");

  fEsdTrackCutsPrim = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsPrim"));
  if (fEsdTrackCutsPrim)
    fEsdTrackCutsPrim->SaveHistograms("esd_track_cuts_primaries");

  fEsdTrackCutsSec = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsSec"));
  if (fEsdTrackCutsSec)
    fEsdTrackCutsSec->SaveHistograms("esd_track_cuts_secondaries");

  fdNdEtaCorrection->SaveHistograms();
  fdNdEtaAnalysisMC->SaveHistograms();
  fdNdEtaAnalysisESD->SaveHistograms();

  if (fOutput->FindObject("dndeta_correction_ND"))
  {
    fdNdEtaCorrectionSpecial[0] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_ND"));
    fdNdEtaCorrectionSpecial[1] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_SD"));
    fdNdEtaCorrectionSpecial[2] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_DD"));
  }
  else
  {
    fdNdEtaCorrectionSpecial[0] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_pi"));
    fdNdEtaCorrectionSpecial[1] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_K"));
    fdNdEtaCorrectionSpecial[2] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_p"));
    fdNdEtaCorrectionSpecial[3] = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction_other"));
  }
  for (Int_t i=0; i<4; ++i)
    if (fdNdEtaCorrectionSpecial[i])
      fdNdEtaCorrectionSpecial[i]->SaveHistograms();

  fTemp1 = dynamic_cast<TH1*> (fOutput->FindObject("fTemp1"));
  if (fTemp1)
    fTemp1->Write();
  fTemp2 = dynamic_cast<TH1*> (fOutput->FindObject("fTemp2"));
  if (fTemp2)
    fTemp2->Write();

  fVertexCorrelation = dynamic_cast<TH2F*> (fOutput->FindObject("fVertexCorrelation"));
  if (fVertexCorrelation)
    fVertexCorrelation->Write();
  fVertexCorrelationShift = dynamic_cast<TH2F*> (fOutput->FindObject("fVertexCorrelationShift"));
  if (fVertexCorrelationShift)
    fVertexCorrelationShift->Write();
  fVertexProfile = dynamic_cast<TProfile*> (fOutput->FindObject("fVertexProfile"));
  if (fVertexProfile)
    fVertexProfile->Write();
  fVertexShift = dynamic_cast<TH1F*> (fOutput->FindObject("fVertexShift"));
  if (fVertexShift)
    fVertexShift->Write();
  fVertexShiftNorm = dynamic_cast<TH1F*> (fOutput->FindObject("fVertexShiftNorm"));
  if (fVertexShiftNorm)
    fVertexShiftNorm->Write();

  fEtaCorrelation = dynamic_cast<TH2F*> (fOutput->FindObject("fEtaCorrelation"));
  if (fEtaCorrelation)
    fEtaCorrelation->Write();
  fEtaCorrelationShift = dynamic_cast<TH2F*> (fOutput->FindObject("fEtaCorrelationShift"));
  if (fEtaCorrelationShift)
    fEtaCorrelationShift->Write();
  fEtaProfile = dynamic_cast<TProfile*> (fOutput->FindObject("fEtaProfile"));
  if (fEtaProfile)
    fEtaProfile->Write();
  fEtaResolution = dynamic_cast<TH1F*> (fOutput->FindObject("fEtaResolution"));
  if (fEtaResolution)
    fEtaResolution->Write();
  fpTResolution = dynamic_cast<TH2F*> (fOutput->FindObject("fpTResolution"));
  if (fpTResolution)
  {
    fpTResolution->FitSlicesY();
    fpTResolution->Write();
  }

  fMultAll = dynamic_cast<TH1F*> (fOutput->FindObject("fMultAll"));
  if (fMultAll)
    fMultAll->Write();

  fMultTr = dynamic_cast<TH1F*> (fOutput->FindObject("fMultTr"));
  if (fMultTr)
    fMultTr->Write();

  fMultVtx = dynamic_cast<TH1F*> (fOutput->FindObject("fMultVtx"));
  if (fMultVtx)
    fMultVtx->Write();

  for (Int_t i=0; i<8; ++i)
  {
    fDeltaPhi[i] = dynamic_cast<TH2*> (fOutput->FindObject(Form("fDeltaPhi_%d", i)));
    if (fDeltaPhi[i])
      fDeltaPhi[i]->Write();
  }

  fEventStats = dynamic_cast<TH2F*> (fOutput->FindObject("fEventStats"));
  if (fEventStats)
    fEventStats->Write();

  fPIDParticles = dynamic_cast<TH1F*> (fOutput->FindObject("fPIDParticles"));
  if (fPIDParticles)
    fPIDParticles->Write();

  fPIDTracks = dynamic_cast<TH1F*> (fOutput->FindObject("fPIDTracks"));
  if (fPIDTracks)
    fPIDTracks->Write();

 //fdNdEtaCorrection->DrawHistograms();

  Printf("Writting result to %s", fileName.Data());

  if (fPIDParticles && fPIDTracks)
  {
    TDatabasePDG* pdgDB = new TDatabasePDG;

    for (Int_t i=0; i <= fPIDParticles->GetNbinsX()+1; ++i)
      if (fPIDParticles->GetBinContent(i) > 0)
      {
        TObject* pdgParticle = pdgDB->GetParticle((Int_t) fPIDParticles->GetBinCenter(i));
        printf("PDG = %d (%s): generated: %d, reconstructed: %d, ratio: %f\n", (Int_t) fPIDParticles->GetBinCenter(i), (pdgParticle) ? pdgParticle->GetName() : "not found", (Int_t) fPIDParticles->GetBinContent(i), (Int_t) fPIDTracks->GetBinContent(i), ((fPIDTracks->GetBinContent(i) > 0) ? fPIDParticles->GetBinContent(i) / fPIDTracks->GetBinContent(i) : -1));
      }

    delete pdgDB;
    pdgDB = 0;
  }

  fout->Write();
  fout->Close();
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

