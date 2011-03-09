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
#include <TDatabasePDG.h>
#include <TRandom.h>

#include <AliLog.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliESDRun.h>
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
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"

ClassImp(AlidNdEtaCorrectionTask)

AlidNdEtaCorrectionTask::AlidNdEtaCorrectionTask() :
  AliAnalysisTask(),
  fESD(0),
  fOutput(0),
  fOption(),
  fAnalysisMode((AliPWG0Helper::AnalysisMode) (AliPWG0Helper::kTPC | AliPWG0Helper::kFieldOn)),
  fTrigger(AliTriggerAnalysis::kMB1),
  fFillPhi(kFALSE),
  fDeltaPhiCut(-1),
  fSymmetrize(kFALSE),
  fMultAxisEta1(kFALSE),
  fDiffTreatment(AliPWG0Helper::kMCFlags),
  fSignMode(0),
  fOnlyPrimaries(kFALSE),
  fStatError(0),
  fSystSkipParticles(kFALSE),
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
  fEventStats(0),
  fEsdTrackCutsCheck(0),
  fEtaCorrelationAllESD(0),
  fpTCorrelation(0),
  fpTCorrelationShift(0),
  fpTCorrelationAllESD(0),
  fpTCorrelationShiftAllESD(0),
  fPtMin(0.15),
  fPtMC(0),
  fEtaMC(0),
  fPtESD(0),
  fEtaESD(0),
  fVtxMC(0),
  fNumberEventMC(0),
  fNumberEvent(0),
  fEventNumber(-1),
  fWeightSecondaries(kFALSE)
{
  //
  // Constructor. Initialization of pointers
  //

  for (Int_t i=0; i<4; i++)
    fdNdEtaCorrectionSpecial[i] = 0;
  
  for (Int_t i=0; i<8; i++)
    fDeltaPhi[i] = 0;

  AliLog::SetClassDebugLevel("AlidNdEtaCorrectionTask", AliLog::kWarning);
}

AlidNdEtaCorrectionTask::AlidNdEtaCorrectionTask(const char* opt) :
  AliAnalysisTask("AlidNdEtaCorrectionTask", ""),
  fESD(0),
  fOutput(0),
  fOption(opt),
  fAnalysisMode((AliPWG0Helper::AnalysisMode) (AliPWG0Helper::kTPC | AliPWG0Helper::kFieldOn)),
  fTrigger(AliTriggerAnalysis::kMB1),
  fFillPhi(kFALSE),
  fDeltaPhiCut(0),
  fSymmetrize(kFALSE),
  fMultAxisEta1(kFALSE),
  fDiffTreatment(AliPWG0Helper::kMCFlags),
  fSignMode(0),
  fOnlyPrimaries(kFALSE),
  fStatError(0),
  fSystSkipParticles(kFALSE),
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
  fEventStats(0),
  fEsdTrackCutsCheck(0),
  fEtaCorrelationAllESD(0),
  fpTCorrelation(0),
  fpTCorrelationShift(0),
  fpTCorrelationAllESD(0),
  fpTCorrelationShiftAllESD(0),
  fPtMin(0.15),
  fPtMC(0),
  fEtaMC(0),
  fPtESD(0),
  fEtaESD(0),
  fVtxMC(0),
  fNumberEventMC(0),
  fNumberEvent(0),
  fEventNumber(-1),
  fWeightSecondaries(kFALSE)
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

  AliLog::SetClassDebugLevel("AlidNdEtaCorrectionTask", AliLog::kWarning);
}

AlidNdEtaCorrectionTask::~AlidNdEtaCorrectionTask()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
	
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
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

    if (fAnalysisMode & AliPWG0Helper::kSPD)
      esdH->SetActiveBranches("AliESDHeader Vertex AliMultiplicity");

    if (fAnalysisMode & AliPWG0Helper::kTPC || fAnalysisMode & AliPWG0Helper::kTPCITS) {
      esdH->SetActiveBranches("AliESDHeader Vertex Tracks");
    }
  }

  // disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);
}

void AlidNdEtaCorrectionTask::CreateOutputObjects()
{
  // create result objects and add to output list

  AliDebug(2,Form("*********************************** fOption = %s", fOption.Data()));
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
    fEsdTrackCutsPrim = static_cast<AliESDtrackCuts*> (fEsdTrackCuts->Clone("fEsdTrackCutsPrim"));
    fEsdTrackCutsSec  = static_cast<AliESDtrackCuts*> (fEsdTrackCuts->Clone("fEsdTrackCutsSec"));
    fEsdTrackCutsCheck  = static_cast<AliESDtrackCuts*> (fEsdTrackCuts->Clone("fEsdTrackCutsCheck"));
    fEsdTrackCutsCheck->SetPtRange(0.15);
    fEsdTrackCutsCheck->SetEtaRange(-1.,1.);
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

  
  //fTemp1 = new TH2F("fTemp1", "fTemp1", 4, 0.5, 4.5, 101, -1.5, 99.5); // nsd study
  fTemp1 = new TH2F("fTemp1", "fTemp1", 300, -15, 15, 80, -2.0, 2.0); 
  fOutput->Add(fTemp1);
  
  fTemp2 = new TH2F("fTemp2", "fTemp2", 300, -15, 15, 80, -2.0, 2.0); 
  fOutput->Add(fTemp2);

  fVertexCorrelation = new TH2F("fVertexCorrelation", "fVertexCorrelation;MC z-vtx;ESD z-vtx", 120, -30, 30, 120, -30, 30);
  fOutput->Add(fVertexCorrelation);
  fVertexCorrelationShift = new TH3F("fVertexCorrelationShift", "fVertexCorrelationShift;MC z-vtx;MC z-vtx - ESD z-vtx;rec. tracks", 120, -30, 30, 100, -1, 1, 100, -0.5, 99.5);
  fOutput->Add(fVertexCorrelationShift);
  fVertexProfile = new TProfile("fVertexProfile", "fVertexProfile;MC z-vtx;MC z-vtx - ESD z-vtx", 40, -20, 20);
  fOutput->Add(fVertexProfile);
  fVertexShift = new TH1F("fVertexShift", "fVertexShift;(MC z-vtx - ESD z-vtx);Entries", 201, -2, 2);
  fOutput->Add(fVertexShift);
  fVertexShiftNorm = new TH2F("fVertexShiftNorm", "fVertexShiftNorm;(MC z-vtx - ESD z-vtx);rec. tracks;Entries", 200, -100, 100, 100, -0.5, 99.5);
  fOutput->Add(fVertexShiftNorm);

  fEtaCorrelation = new TH2F("fEtaCorrelation", "fEtaCorrelation;MC #eta;ESD #eta", 120, -3, 3, 120, -3, 3);
  fOutput->Add(fEtaCorrelation);
  fEtaCorrelationAllESD = new TH2F("fEtaCorrelationAllESD", "fEtaCorrelationAllESD;MC #eta;ESD #eta", 120, -3, 3, 120, -3, 3);
  fOutput->Add(fEtaCorrelationAllESD);
  fEtaCorrelationShift = new TH2F("fEtaCorrelationShift", "fEtaCorrelationShift;MC #eta;MC #eta - ESD #eta", 120, -3, 3, 100, -0.1, 0.1);
  fOutput->Add(fEtaCorrelationShift);
  fEtaProfile = new TProfile("fEtaProfile", "fEtaProfile;MC #eta;MC #eta - ESD #eta", 120, -3, 3);
  fOutput->Add(fEtaProfile);
  fEtaResolution = new TH1F("fEtaResolution", "fEtaResolution;MC #eta - ESD #eta", 201, -0.2, 0.2);
  fOutput->Add(fEtaResolution);

  fpTResolution = new TH2F("fpTResolution", ";MC p_{T} (GeV/c);(MC p_{T} - ESD p_{T}) / MC p_{T}", 160, 0, 20, 201, -0.2, 0.2);
  fOutput->Add(fpTResolution);

  fpTCorrelation = new TH2F("fpTCorrelation", "fpTCorrelation;MC p_{T} (GeV/c);ESD p_{T}", 160, 0, 20, 160, 0, 20);
  fOutput->Add(fpTCorrelation);
  fpTCorrelationShift = new TH2F("fpTCorrelationShift", "fpTCorrelationShift;MC p_{T} (GeV/c);MC p_{T} - ESD p_{T}", 160, 0, 20, 100, -1, 1);
  fOutput->Add(fpTCorrelationShift);
  fpTCorrelationAllESD = new TH2F("fpTCorrelationAllESD", "fpTCorrelationAllESD;MC p_{T} (GeV/c);ESD p_{T}", 160, 0, 20, 160, 0, 20);
  fOutput->Add(fpTCorrelationAllESD);
  fpTCorrelationShiftAllESD = new TH2F("fpTCorrelationShiftAllESD", "fpTCorrelationShiftAllESD;MC p_{T} (GeV/c);MC p_{T} - ESD p_{T}", 160, 0, 20, 100, -1, 1);
  fOutput->Add(fpTCorrelationShiftAllESD);

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

  fEventStats = new TH2F("fEventStats", "fEventStats;event type;status;count", 109, -6.5, 102.5, 4, -0.5, 3.5);
  fOutput->Add(fEventStats);
  fEventStats->GetXaxis()->SetBinLabel(1, "INEL"); // x = -5
  fEventStats->GetXaxis()->SetBinLabel(2, "NSD");  // x = -4
  fEventStats->GetXaxis()->SetBinLabel(3, "ND");   // x = -3
  fEventStats->GetXaxis()->SetBinLabel(4, "SD");   // x = -2
  fEventStats->GetXaxis()->SetBinLabel(5, "DD");   // x = -1

  fEventStats->GetXaxis()->SetBinLabel(108, "INEL=0");   // x = -101
  fEventStats->GetXaxis()->SetBinLabel(109, "INEL>0");   // x = -102
  
  for (Int_t i=-1; i<100; i++)
    fEventStats->GetXaxis()->SetBinLabel(7+i, Form("%d", i)); 

  fEventStats->GetYaxis()->SetBinLabel(1, "nothing");
  fEventStats->GetYaxis()->SetBinLabel(2, "trg");
  fEventStats->GetYaxis()->SetBinLabel(3, "vtx");
  fEventStats->GetYaxis()->SetBinLabel(4, "trgvtx");

  if (fEsdTrackCuts)
  {
    fEsdTrackCuts->SetName("fEsdTrackCuts");
    // TODO like this we send an empty object back...
    fOutput->Add(fEsdTrackCuts->Clone());
  }
  fPtMC = new TH1F("fPtMC", "Pt from MC for selected tracks;MC p_{T} (GeV/c)", 160, 0, 20);
  fOutput->Add(fPtMC); 
  fEtaMC = new TH1F("fEtaMC", "Eta from MC for selected tracks;MC #eta", 120, -3, 3);
  fOutput->Add(fEtaMC);
  fPtESD = new TH1F("fPtESD", "Pt from ESD for selected tracks;ESD p_{T} (GeV/c)", 160, 0, 20);
  fOutput->Add(fPtESD);
  fEtaESD = new TH1F("fEtaESD", "Eta from ESD for selected tracks;ESD #eta", 120, -3, 3);
  fOutput->Add(fEtaESD);

  fVtxMC = new TH1F("fVtxMC", "Vtx,z from MC for all events;MC vtx_z (cm)", 100, -30, 30);
  fOutput->Add(fVtxMC); 

  fNumberEventMC = new TH1F("fNumberEventMC","Number of event accepted at MC level",600000,-0.5,600000-0.5);
  fOutput->Add(fNumberEventMC);

  fNumberEvent = new TH1F("fNumberEvent","Number of event accepted at Reco level",600000,-0.5,600000-0.5);
  fOutput->Add(fNumberEvent);
}

void AlidNdEtaCorrectionTask::Exec(Option_t*)
{
  // process the event

  fEventNumber++;
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
    
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
    
  Bool_t eventTriggered = inputHandler->IsEventSelected();

  static AliTriggerAnalysis* triggerAnalysis = 0;
  if (!triggerAnalysis)
  {
    AliPhysicsSelection* physicsSelection = dynamic_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if (physicsSelection)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
  }
    
  if (eventTriggered)
    eventTriggered = triggerAnalysis->IsTriggerFired(fESD, fTrigger);
    
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

  // get process type
  Int_t processType = AliPWG0Helper::GetEventProcessType(fESD, header, stack, fDiffTreatment);
  
  AliDebug(AliLog::kDebug+1, Form("Found process type %d", processType));

  if (processType == AliPWG0Helper::kInvalidProcess)
  {
    AliDebug(AliLog::kWarning, "Unknown process type. Setting to ND");
    processType = AliPWG0Helper::kND;
  }

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();
  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);
  fVtxMC->Fill(vtxMC[2]);
  AliDebug(2,Form("MC vtx: x = %.6f, y = %.6f, z = %.6f",vtxMC[0],vtxMC[1],vtxMC[2]));

  // get the ESD vertex
  Bool_t eventVertex = kFALSE;
  Double_t vtx[3];
  const AliESDVertex* vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);
  if (vtxESD && AliPWG0Helper::TestVertex(vtxESD, fAnalysisMode))
  {
    vtxESD->GetXYZ(vtx);
    eventVertex = kTRUE;

    // remove vertices outside +- 15 cm
    if (TMath::Abs(vtx[2]) > 15)
    {
      eventVertex = kFALSE;
      vtxESD = 0;
    }
  }
  else
    vtxESD = 0;
    
  // create list of (label, eta, pt) tuples
  Int_t inputCount = 0;
  Int_t* labelArr = 0;
  Int_t* labelArr2 = 0; // only for case of SPD
  Float_t* etaArr = 0;
  Float_t* thirdDimArr = 0;
  Float_t* deltaPhiArr = 0;
  if (fAnalysisMode & AliPWG0Helper::kSPD)
  {
    if (vtxESD)
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
  
      Bool_t foundInEta10 = kFALSE;
      
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
  
        if (fDeltaPhiCut > 0 && (TMath::Abs(deltaPhi) > fDeltaPhiCut || TMath::Abs(mult->GetDeltaTheta(i)) > fDeltaPhiCut / 0.08 * 0.025))
          continue;
          
        if (fSystSkipParticles && gRandom->Uniform() < 0.0153)
        {
          Printf("Skipped tracklet!");
          continue;
        }
        
        // TEST exclude potentially inefficient phi region
        //if (phi > 5.70 || phi < 0.06)
        //  continue;
            
        // we have to repeat the trigger here, because the tracklet might have been kicked out fSystSkipParticles
        if (TMath::Abs(mult->GetEta(i)) < 1)
          foundInEta10 = kTRUE;
        
        etaArr[inputCount] = mult->GetEta(i);
        if (fSymmetrize)
          etaArr[inputCount] = TMath::Abs(etaArr[inputCount]);
        labelArr[inputCount] = mult->GetLabel(i, 0);
        labelArr2[inputCount] = mult->GetLabel(i, 1);
        thirdDimArr[inputCount] = phi;
        deltaPhiArr[inputCount] = deltaPhi;
        ++inputCount;
      }
      
      /*
      for (Int_t i=0; i<mult->GetNumberOfSingleClusters(); ++i)
      {
        if (TMath::Abs(TMath::Log(TMath::Tan(mult->GetThetaSingle(i)/2.))) < 1);
        {
          foundInEta10 = kTRUE;
          break;
        }
      }
      */
      
      if (fSystSkipParticles && (fTrigger & AliTriggerAnalysis::kOneParticle) && !foundInEta10)
        eventTriggered = kFALSE;
    }
  }
  else if (fAnalysisMode & AliPWG0Helper::kTPC || fAnalysisMode & AliPWG0Helper::kTPCITS)
  {
    if (!fEsdTrackCuts)
    {
      AliDebug(AliLog::kError, "fESDTrackCuts not available");
      return;
    }
    
    Bool_t foundInEta10 = kFALSE;
    
    if (vtxESD)
    {
      // control histograms on pT
      Int_t nfromstack = stack->GetNtrack();
      AliDebug(3,Form(" from stack we have %d tracks\n",nfromstack));
      for (Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
	AliESDtrack* esdTrackcheck = dynamic_cast<AliESDtrack*> (fESD->GetTrack(itrack));
	if (!esdTrackcheck){
	  AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", itrack));
	  continue;
	}
	if (fOnlyPrimaries){
	  Int_t label = TMath::Abs(esdTrackcheck->GetLabel());
	  AliDebug(4,Form("label = %d\n",label));
	  if (label == 0 || label > nfromstack) continue;
	  if (stack->IsPhysicalPrimary(label) == kFALSE) continue;
	}
        
	Int_t label = TMath::Abs(esdTrackcheck->GetLabel());
	if (label == 0 || label > nfromstack) continue;
	if (!stack->Particle(label)){
	  AliDebug(4,Form("WARNING: No particle for %d", label));
	}
	else{	
	  if (!fEsdTrackCuts->AcceptTrack(esdTrackcheck)){
	    TParticle* particle = stack->Particle(label);
	    if (fEsdTrackCutsCheck->AcceptTrack(esdTrackcheck)){
	      //if (TMath::Abs(particle->Eta() < 0.8) && particle->Pt() > fPtMin){
	      Float_t ptMC = particle->Pt();
	      Float_t etaMC = particle->Eta();
	      Float_t ptESD = esdTrackcheck->Pt();
	      Float_t etaESD = esdTrackcheck->Eta();
	      fEtaCorrelationAllESD->Fill(etaMC,etaESD);
	      fpTCorrelationAllESD->Fill(ptMC,ptESD);
	      fpTCorrelationShiftAllESD->Fill(ptMC,ptMC-ptESD);
	    }
	  }
	}
      } // end loop over all ESDs

      // get multiplicity from ESD tracks
      TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD, (fAnalysisMode & AliPWG0Helper::kTPC));
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
        
	AliDebug(3,Form("particle %d: pt = %.6f, eta = %.6f",i,esdTrack->Pt(), esdTrack->Eta())); 
	// 2 Options for INEL>0 trigger - choose one
 	// 1. HL
	//if (esdTrack->Pt() < 0.15)
 	// foundInEta10 = kTRUE;
 	// 2. MB Working Group definition
        if (esdTrack->Pt() < fPtMin)
          continue;
        
        if (fOnlyPrimaries)
        {
          Int_t label = TMath::Abs(esdTrack->GetLabel());
          if (label == 0)
            continue;
          
          if (stack->IsPhysicalPrimary(label) == kFALSE)
            continue;
        }
        
	Int_t label = TMath::Abs(esdTrack->GetLabel());
	if (!stack->Particle(label)){
	  AliDebug(3,Form("WARNING: No particle for %d", label));
	}
	else{
	  TParticle* particle = stack->Particle(label);
	  Float_t ptMC = particle->Pt();
	  Float_t etaMC = particle->Eta();
	  fPtMC->Fill(ptMC);
	  fEtaMC->Fill(etaMC);
	  fPtESD->Fill(esdTrack->Pt());
	  fEtaESD->Fill(esdTrack->Eta());
	}

        // 2 Options for INEL>0 trigger - choose one
	// 1. HL
        //if (TMath::Abs(esdTrack->Eta()) < 1 && esdTrack->Pt() > 0.15)
	// foundInEta10 = kTRUE;
	// 2. MB Working Group definition
	if (TMath::Abs(esdTrack->Eta()) < 0.8 && esdTrack->Pt() > fPtMin)
          foundInEta10 = kTRUE;

        etaArr[inputCount] = esdTrack->Eta();
        if (fSymmetrize)
          etaArr[inputCount] = TMath::Abs(etaArr[inputCount]);
        labelArr[inputCount] = TMath::Abs(esdTrack->GetLabel());
        labelArr2[inputCount] = labelArr[inputCount]; // no second label for tracks
        thirdDimArr[inputCount] = esdTrack->Pt();
        deltaPhiArr[inputCount] = 0; // no delta phi for tracks
        ++inputCount;
      }

      delete list;

      // TODO this code crashes for TPCITS because particles are requested from the stack for some labels that are out of bound
      if (0 && eventTriggered)
      {
        // collect values for primaries and secondaries
        for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++)
        {
          AliESDtrack* track = 0;
  
          if (fAnalysisMode & AliPWG0Helper::kTPC)
            track = AliESDtrackCuts::GetTPCOnlyTrack(fESD, iTrack);
          else if (fAnalysisMode & AliPWG0Helper::kTPCITS)
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
  
          if (TMath::Abs(track->Eta()) < 0.8 && track->Pt() > 0.15)
          {
            if (stack->IsPhysicalPrimary(label))
            {
              // primary
              if (fEsdTrackCutsPrim->AcceptTrack(track)) 
              {
//                 if (AliESDtrackCuts::GetSigmaToVertex(track) > 900)
//                 {
//                   Printf("Track %d has nsigma of %f. Printing track and vertex...", iTrack, AliESDtrackCuts::GetSigmaToVertex(track));
//                   Float_t b[2];
//                   Float_t r[3];
//                   track->GetImpactParameters(b, r);
//                   Printf("Impact parameter %f %f and resolution: %f %f %f", b[0], b[1], r[0], r[1], r[2]);
//                   track->Print("");
//                   if (vtxESD)
//                     vtxESD->Print();
//                 }
              }
            }
            else
            {
              // secondary
              fEsdTrackCutsSec->AcceptTrack(track);
            }
          }
  
          // TODO mem leak in the continue statements above
          if (fAnalysisMode & AliPWG0Helper::kTPC)
            delete track;
        }
      }
    }
    
    if (!foundInEta10)
      eventTriggered = kFALSE;
    else{
     //Printf("The event triggered. Its number in file is %d",fESD->GetEventNumberInFile());
      fNumberEvent->Fill(fESD->GetEventNumberInFile());
    }
  }
  else
    return;

  // fill process type
  Int_t biny = (Int_t) eventTriggered + 2 * (Int_t) eventVertex;
  // INEL
  fEventStats->Fill(-6, biny);
  // NSD
  if (processType != AliPWG0Helper::kSD)
    fEventStats->Fill(-5, biny);
  // SD, ND, DD
  if (processType == AliPWG0Helper::kND)
    fEventStats->Fill(-4, biny);
  if (processType == AliPWG0Helper::kSD)
    fEventStats->Fill(-3, biny);
  if (processType == AliPWG0Helper::kDD)
    fEventStats->Fill(-2, biny);
  
  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();
  Int_t nAccepted = 0;

  Bool_t oneParticleEvent = kFALSE;
  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    //Printf("Getting particle %d", iMc);
    TParticle* particle = stack->Particle(iMc);

    if (!particle)
      continue;

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;
    
    // for INEL > 0, MB Working Group definition use the second option
    // 1. standard
    //if (TMath::Abs(particle->Eta()) < 1.0)
    // 2. MB Working Group definition
    if (TMath::Abs(particle->Eta()) < 0.8 && particle->Pt() > fPtMin)
    {
      oneParticleEvent = kTRUE;
      fNumberEventMC->Fill(fESD->GetEventNumberInFile());
      break;
    }
  }
  
  if (TMath::Abs(vtxMC[2]) < 5.5)
  {
    if (oneParticleEvent)
      fEventStats->Fill(102, biny);
    else
      fEventStats->Fill(101, biny);
  }
  
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
    
    // for INEL > 0, MB Working Group definition use the second option
    // 1. standard
    //if (fPIDParticles && TMath::Abs(particle->Eta()) < 1.0)
    // 2. MB Working Group definition
    if (fPIDParticles && TMath::Abs(particle->Eta()) < 0.8 && particle->Pt() > fPtMin)
      fPIDParticles->Fill(particle->GetPdgCode());

    Float_t eta = particle->Eta();
    if (fSymmetrize)
      eta = TMath::Abs(eta);
    
    Float_t thirdDim = -1;
    if (fAnalysisMode & AliPWG0Helper::kSPD)
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

    Int_t processType2 = processType;
    if (oneParticleEvent)
      processType2 |= AliPWG0Helper::kOnePart;

    fdNdEtaCorrection->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType2);

    if (fOption.Contains("process-types"))
    {
      // non diffractive
      if (processType==AliPWG0Helper::kND)
        fdNdEtaCorrectionSpecial[0]->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType2);

      // single diffractive
      if (processType==AliPWG0Helper::kSD)
        fdNdEtaCorrectionSpecial[1]->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType2);

      // double diffractive
      if (processType==AliPWG0Helper::kDD)
        fdNdEtaCorrectionSpecial[2]->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType2);
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
      fdNdEtaCorrectionSpecial[id]->FillMCParticle(vtxMC[2], eta, thirdDim, eventTriggered, eventVertex, processType2);
    }

    if (eventTriggered)
      if (eventVertex)
        fdNdEtaAnalysisMC->FillTrack(vtxMC[2], eta, thirdDim);

    // TODO this value might be needed lower for the SPD study (only used for control histograms anyway)
    if (TMath::Abs(eta) < 1.5) // && pt > 0.2)
      nAccepted++;
  }

  if (AliPWG0Helper::GetLastProcessType() >= -1)
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

  Int_t nEta05 = 0;
  Int_t nEta10 = 0;
  for (Int_t i=0; i<inputCount; ++i)
  {
    if (TMath::Abs(etaArr[i]) < 0.5)
      nEta05++;
    if (TMath::Abs(etaArr[i]) < 1.0)
      nEta10++;
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

    // for INEL > 0, MB Working Group definition use the second option
    // 1. standard
    //if (TMath::Abs(particle->Eta()) < 1.0)
    // 2. INEL > 0 MB Working Group definition
    if (TMath::Abs(particle->Eta()) < 0.8 && particle->Pt()>fPtMin)
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
      if (fSymmetrize)
        fEtaResolution->Fill(TMath::Abs(particle->Eta()) - etaArr[i]);
      else
        fEtaResolution->Fill(particle->Eta() - etaArr[i]);

      if (fAnalysisMode & AliPWG0Helper::kTPC || fAnalysisMode & AliPWG0Helper::kTPCITS){
	// for INEL > 0, MB Working Group definition use the second option
	// 1. standard
        //if (TMath::Abs(particle->Eta() < 0.9) && particle->Pt() > 0)
	// 2. INEL > 0 MB WOrking Group definition
	if (TMath::Abs(particle->Eta() < 0.8) && particle->Pt() > 0){
	  fpTResolution->Fill(particle->Pt(), (particle->Pt() - thirdDimArr[i]) / particle->Pt());
	  fpTCorrelation->Fill(particle->Pt(),thirdDimArr[i]);
	  fpTCorrelationShift->Fill(particle->Pt(),particle->Pt()-thirdDimArr[i]);
	}
      }

      Float_t eta = -999;
      Float_t thirdDim = -1;

      Bool_t firstIsPrim = stack->IsPhysicalPrimary(label);
      // in case of same label the MC values are filled, otherwise (background) the reconstructed values
      if (label == label2)
      {
        eta = particle->Eta();
        if (fSymmetrize)
          eta = TMath::Abs(eta);
        
        if (fAnalysisMode & AliPWG0Helper::kSPD)
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
        if (fAnalysisMode & AliPWG0Helper::kSPD && !fFillPhi)
        {
          thirdDim = (fMultAxisEta1) ? nEta10 : inputCount;
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
      {
        Double_t weight = 1.;
	if (fWeightSecondaries){
	  if (!firstIsPrim){
	    weight = GetSecondaryCorrection(thirdDim);
	  }
	}        
	fdNdEtaCorrection->FillTrackedParticle(vtxMC[2], eta, thirdDim, weight);
        fTemp2->Fill(vtxMC[2], eta);
      }
      
      // eta comparison for tracklets with the same label (others are background)
      if (label == label2)
      {
        Float_t eta2 = particle->Eta();
        if (fSymmetrize)
          eta2 = TMath::Abs(eta2);
        
        fEtaProfile->Fill(eta2, eta2 - etaArr[i]);
        fEtaCorrelation->Fill(eta2, etaArr[i]);
        fEtaCorrelationShift->Fill(eta2, eta2 - etaArr[i]);
      }

      if (fSymmetrize)
        fdNdEtaAnalysisESD->FillTrack(vtxMC[2], TMath::Abs(particle->Eta()), thirdDim);
      else
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
    Double_t diff = vtxMC[2] - vtx[2];
    fVertexShift->Fill(diff);
    
    fVertexCorrelation->Fill(vtxMC[2], vtx[2]);
    fVertexCorrelationShift->Fill(vtxMC[2], vtxMC[2] - vtx[2], inputCount);
    fVertexProfile->Fill(vtxMC[2], vtxMC[2] - vtx[2]);
  
    if (vtxESD->IsFromVertexerZ() && inputCount > 0)
      fVertexShiftNorm->Fill(diff, vtxESD->GetNContributors());
  }

  Int_t multAxis = inputCount;
  if (fMultAxisEta1)
    multAxis = nEta10;

  if (eventTriggered && eventVertex)
  {
    fdNdEtaAnalysisMC->FillEvent(vtxMC[2], multAxis);
    fdNdEtaAnalysisESD->FillEvent(vtxMC[2], multAxis);
  }

  Int_t processType2 = processType;
  if (oneParticleEvent)
    processType2 |= AliPWG0Helper::kOnePart;

  // stuff regarding the vertex reco correction and trigger bias correction
  fdNdEtaCorrection->FillEvent(vtxMC[2], multAxis, eventTriggered, eventVertex, processType2);

  if (fOption.Contains("process-types"))
  {
    // non diffractive
    if (processType == AliPWG0Helper::kND)
      fdNdEtaCorrectionSpecial[0]->FillEvent(vtxMC[2], multAxis, eventTriggered, eventVertex, processType2);

    // single diffractive
    if (processType == AliPWG0Helper::kSD)
      fdNdEtaCorrectionSpecial[1]->FillEvent(vtxMC[2], multAxis, eventTriggered, eventVertex, processType2);

    // double diffractive
    if (processType == AliPWG0Helper::kDD)
      fdNdEtaCorrectionSpecial[2]->FillEvent(vtxMC[2], multAxis, eventTriggered, eventVertex, processType2);
  }
  
  if (fOption.Contains("particle-species"))
    for (Int_t id=0; id<4; id++)
      fdNdEtaCorrectionSpecial[id]->FillEvent(vtxMC[2], multAxis, eventTriggered, eventVertex, processType2);

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
  fVertexCorrelationShift = dynamic_cast<TH3F*> (fOutput->FindObject("fVertexCorrelationShift"));
  if (fVertexCorrelationShift)
  {
    ((TH2*) fVertexCorrelationShift->Project3D("yx"))->FitSlicesY();
    fVertexCorrelationShift->Write();
  }
  fVertexProfile = dynamic_cast<TProfile*> (fOutput->FindObject("fVertexProfile"));
  if (fVertexProfile)
    fVertexProfile->Write();
  fVertexShift = dynamic_cast<TH1F*> (fOutput->FindObject("fVertexShift"));
  if (fVertexShift)
    fVertexShift->Write();
  fVertexShiftNorm = dynamic_cast<TH2F*> (fOutput->FindObject("fVertexShiftNorm"));
  if (fVertexShiftNorm)
  {
    fVertexShiftNorm->ProjectionX();
    fVertexShiftNorm->Write();
  }  

  fEtaCorrelation = dynamic_cast<TH2F*> (fOutput->FindObject("fEtaCorrelation"));
  if (fEtaCorrelation)
    fEtaCorrelation->Write();
  fEtaCorrelationAllESD = dynamic_cast<TH2F*> (fOutput->FindObject("fEtaCorrelationAllESD"));
  if (fEtaCorrelationAllESD)
    fEtaCorrelationAllESD->Write();
  fEtaCorrelationShift = dynamic_cast<TH2F*> (fOutput->FindObject("fEtaCorrelationShift"));
  if (fEtaCorrelationShift)
  {
    fEtaCorrelationShift->FitSlicesY();
    fEtaCorrelationShift->Write();
  }
  fEtaProfile = dynamic_cast<TProfile*> (fOutput->FindObject("fEtaProfile"));
  if (fEtaProfile)
    fEtaProfile->Write();
  fEtaResolution = dynamic_cast<TH1F*> (fOutput->FindObject("fEtaResolution"));
  if (fEtaResolution)
    fEtaResolution->Write();
  fpTCorrelation = dynamic_cast<TH2F*> (fOutput->FindObject("fpTCorrelation"));
  if (fpTCorrelation)
    fpTCorrelation->Write();
  fpTCorrelationShift = dynamic_cast<TH2F*> (fOutput->FindObject("fpTCorrelationShift"));
  if (fpTCorrelationShift)
  {
    fpTCorrelationShift->FitSlicesY();
    fpTCorrelationShift->Write();
  }
  fpTResolution = dynamic_cast<TH2F*> (fOutput->FindObject("fpTResolution"));
  if (fpTResolution)
  {
    fpTResolution->FitSlicesY();
    fpTResolution->Write();
  }

  fpTCorrelationAllESD = dynamic_cast<TH2F*> (fOutput->FindObject("fpTCorrelationAllESD"));
  if (fpTCorrelationAllESD)
    fpTCorrelationAllESD->Write();
  fpTCorrelationShiftAllESD = dynamic_cast<TH2F*> (fOutput->FindObject("fpTCorrelationShiftAllESD"));
  if (fpTCorrelationShiftAllESD)
  {
    fpTCorrelationShiftAllESD->FitSlicesY();
    fpTCorrelationShiftAllESD->Write();
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

  fPtMC = dynamic_cast<TH1F*> (fOutput->FindObject("fPtMC"));
  if (fPtMC)
    fPtMC->Write();
  fEtaMC = dynamic_cast<TH1F*> (fOutput->FindObject("fEtaMC"));
  if (fEtaMC)
    fEtaMC->Write();
  fPtESD = dynamic_cast<TH1F*> (fOutput->FindObject("fPtESD"));
  if (fPtESD)
    fPtESD->Write();
  fEtaESD = dynamic_cast<TH1F*> (fOutput->FindObject("fEtaESD"));
  if (fEtaESD)
    fEtaESD->Write();
  fVtxMC = dynamic_cast<TH1F*> (fOutput->FindObject("fVtxMC"));
  if (fVtxMC)
    fVtxMC->Write();

  fNumberEventMC = dynamic_cast<TH1F*> (fOutput->FindObject("fNumberEventMC"));
  if (fNumberEventMC)
    fNumberEventMC->Write();
  fNumberEvent = dynamic_cast<TH1F*> (fOutput->FindObject("fNumberEvent"));
  if (fNumberEvent)
    fNumberEvent->Write();

 //fdNdEtaCorrection->DrawHistograms();

  Printf("Writing result to %s", fileName.Data());

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

Double_t AlidNdEtaCorrectionTask::GetSecondaryCorrection(Double_t pt){

	// getting the data driven correction factor to correct for 
	// the underestimate of secondaries in Pythia
	// (A. Dainese + J. Otwinowski

	if (pt <= 0.17) return 1.0;
	if (pt <= 0.4) return GetLinearInterpolationValue(0.17,1.0,0.4,1.07, pt);
	if (pt <= 0.6) return GetLinearInterpolationValue(0.4,1.07,0.6,1.25, pt);
	if (pt <= 1.2) return GetLinearInterpolationValue(0.6,1.25,1.2,1.5,  pt);
	return 1.5;

}

Double_t AlidNdEtaCorrectionTask::GetLinearInterpolationValue(Double_t const x1, Double_t const y1, Double_t const x2, Double_t const y2, const Double_t pt)
{

	//
	// linear interpolation to be used to get the secondary correction (see AlidNdEtaCorrectionTask::GetSecondaryCorrection)
	//

	return ((y2-y1)/(x2-x1))*pt+(y2-(((y2-y1)/(x2-x1))*x2));
}

