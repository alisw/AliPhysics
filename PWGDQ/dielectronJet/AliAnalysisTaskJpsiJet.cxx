/**************************************************************************
 * Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TChain.h"
#include "THnSparse.h"

#include "AliAnalysisManager.h"

#include "AliAODEvent.h"
#include "AliAODCaloCluster.h"

#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"

#include "AliDielectronVarCuts.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronPairLegCuts.h"

#include "AliAnalysisTaskJpsiJet.h"

class AliGenEventHeader;

class AliDielectronVarCuts;
class AliDielectronTrackCuts;
class AliDielectronPairLegCuts;
class AliJetContainer;

class AliAnalysisTaskJpsiJet;

using namespace std;

ClassImp(AliAnalysisTaskJpsiJet)

AliAnalysisTaskJpsiJet::AliAnalysisTaskJpsiJet():
  AliAnalysisTaskSE(),
  fAOD(NULL),
  fDielectron(NULL),
  fPairs(NULL),
  fDaughters(NULL),
  fTracksWithPair(NULL),
  fJetTasks(NULL),
  fJets(NULL),
  fSelectedTrigger(0),
  fSelectedTriggerClasses(""),
  fFiredTriggerTag(""),
  fRejectPileup(kFALSE),
  fIsPileup(kFALSE),
  fIsTriggerQA(kFALSE),
  fIsCellQA(kFALSE),
  fIsJetFinder(kTRUE),
  fIsMC(kFALSE),
  fEMCalEth(5.0),
  fMCParticles(NULL),
  fMCHeader(NULL),
  fMCGenType(""),
  fEventFilter(NULL),
  fJpsiPair(NULL),
  fTaggedJet(NULL),
  fJpsiMC(NULL),
  fTaggedJetMC(NULL),
  fHistos(NULL),
  fHistosMC(NULL),
  fRunNo(""),
  fVars(NULL)
{
  // Constructor
}

AliAnalysisTaskJpsiJet::AliAnalysisTaskJpsiJet(const char* taskName):
  AliAnalysisTaskSE(taskName),
  fAOD(NULL),
  fDielectron(NULL),
  fPairs(NULL),
  fDaughters(NULL),
  fTracksWithPair(NULL),
  fJetTasks(NULL),
  fJets(NULL),
  fSelectedTrigger(0),
  fSelectedTriggerClasses(""),
  fFiredTriggerTag(""),
  fRejectPileup(kFALSE),
  fIsPileup(kFALSE),
  fIsTriggerQA(kFALSE),
  fIsCellQA(kFALSE),
  fIsJetFinder(kTRUE),
  fIsMC(kFALSE),
  fEMCalEth(5.0),
  fMCParticles(NULL),
  fMCHeader(NULL),
  fMCGenType(""),
  fEventFilter(NULL),
  fJpsiPair(NULL),
  fTaggedJet(NULL),
  fJpsiMC(NULL),
  fTaggedJetMC(NULL),
  fHistos(NULL),
  fHistosMC(NULL),
  fRunNo(""),
  fVars(NULL)
{
  // IO
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  AliInfo(Form("Init task : %s", taskName));
}

AliAnalysisTaskJpsiJet::~AliAnalysisTaskJpsiJet(){
  // Destructor
  if(fEventFilter) delete fEventFilter;
  // Histogram list from AliAnalysisTaskSE
  if(fHistosQA) delete fHistosQA;
  if(fHistos) delete fHistos;
  if(fHistosMC) delete fHistosMC;
}

void AliAnalysisTaskJpsiJet::UserCreateOutputObjects(){
  // PID response
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handler = (AliInputEventHandler*)(mgr->GetInputEventHandler());
  if(handler->GetPIDResponse()){
    AliDielectronVarManager::SetPIDResponse(handler->GetPIDResponse());
  }
  else{
    AliFatal("This task needs the PID response attached to the input event handler!");
    return;
  }
  // Histograms
  fHistos = new THistManager("JpsiJetQA");
  fHistosQA = fHistos->GetListOfHistograms();

  TH1* fHistEventStat = fHistos->CreateTH1("EventStats","Event statistics;Status;N_{event}",int(kEventStatusN),-0.5,float(kEventStatusN)-0.5);
  fHistEventStat->GetXaxis()->SetBinLabel(kAllInAOD + 1, "Before PS");
  fHistEventStat->GetXaxis()->SetBinLabel(kPhysSelected + 1, "After PS");
  fHistEventStat->GetXaxis()->SetBinLabel(kV0ANDtrigger + 1, "V0AND Trig. ");
  fHistEventStat->GetXaxis()->SetBinLabel(kTRDtrigger + 1, "TRD Trig.");
  fHistEventStat->GetXaxis()->SetBinLabel(kFiltered + 1, "After cuts");
  fHistEventStat->GetXaxis()->SetBinLabel(kAfterPileUp + 1, "Pileup Rejected");
  fHistEventStat->GetXaxis()->SetBinLabel(kWithSinglePair + 1, "N_{pair}==1");
  fHistEventStat->GetXaxis()->SetBinLabel(kWithMultiPair + 1, "N_{pair}>1");
  fHistEventStat->GetXaxis()->SetBinLabel(kWithPairInJet + 1, "e^{+}e^{-} in jet");

  InitHistogramsForEventQA("Event_ALL");
  InitHistogramsForEventQA("Event_beforeCuts");
  InitHistogramsForEventQA("Event_afterCuts");

  fVars = new Double_t[AliDielectronVarManager::kNMaxValues];
  for(int i = 0; i < AliDielectronVarManager::kNMaxValues; i++)
    fVars[i] = 0.0;
  InitHistogramsForRunwiseQA("Runwise");

  if(!fIsTriggerQA)
    InitHistogramsForClusterQA("Cluster");
  else{
    InitHistogramsForClusterQA("Cluster_MB");
    InitHistogramsForClusterQA("Cluster_EG1");
    InitHistogramsForClusterQA("Cluster_EG2");
    InitHistogramsForClusterQA("Cluster_DG1");
    InitHistogramsForClusterQA("Cluster_DG2");
  }

  // Init dielectron
  InitHistogramsForDielectron("Dielectron");
  fDielectron->SetDontClearArrays();
  fDielectron->Init();
  fHistosQA->Add(const_cast<THashList*>(fDielectron->GetHistogramList()));

  // Init jet finder tasks
  AliEmcalJetTask* jetFinder = NULL;
  TIter next(fJetTasks);
  while((jetFinder=(AliEmcalJetTask*)next()))
    jetFinder->CreateOutputObjects();
  InitHistogramsForJetQA();

  // Init pair in jet analysis
  InitHistogramsForTaggedJet("PairInJet");

  PostData(1, fHistosQA);

  if(fIsMC){
    fHistosMC = new THistManager("JpsiJetMC");
    InitHistogramsForMC();
    TList* mcHistos = fHistosMC->GetListOfHistograms();
    PostData(2, mcHistos);
  }
}

void AliAnalysisTaskJpsiJet::UserExec(Option_t*){
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fFiredTriggerTag = "";
  fRunNo = Form("%d", fAOD->GetRunNumber());
  fHistos->FillTH1("EventStats", kAllInAOD);
  fHistos->FillTH2("Runwise/NEvent", fRunNo.Data(), "ALL");
  FillHistogramsForEventQA("Event_ALL");

  //
  // Event Selection
  //
  UInt_t isSelected = AliVEvent::kAny;
  // Select trigger
  AliAODHeader* header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  UInt_t offlineTrigger = header->GetOfflineTrigger();
  isSelected &= (fSelectedTrigger & offlineTrigger);
  if(!isSelected) return;
  // Select trigger classes
  TString triggerClass = fAOD->GetFiredTriggerClasses();
  Bool_t isFired = kFALSE;
  if(fSelectedTriggerClasses != ""){
    TObjArray* selectedTrigClasses = fSelectedTriggerClasses.Tokenize(";");
    selectedTrigClasses->SetOwner(kTRUE);
    for(int i = 0; i < selectedTrigClasses->GetEntries(); i++){
      TString tag = selectedTrigClasses->At(i)->GetName();
      if(triggerClass.Contains(tag)){
        isFired = kTRUE;
        fFiredTriggerTag += tag + "_";
      }
    }
    if(fFiredTriggerTag != "")
      fFiredTriggerTag.Remove(fFiredTriggerTag.Length()-1);
  }
  if(!isFired) return;

  fHistos->FillTH1("EventStats", kPhysSelected);
  fHistos->FillTH2("Runwise/NEvent", fRunNo.Data(), "PS");
  FillHistogramsForEventQA("Event_beforeCuts");

  // Event cuts
  if(fEventFilter && !fEventFilter->IsSelected(fAOD)) return;
  fHistos->FillTH1("EventStats", kFiltered);
  fHistos->FillTH2("Runwise/NEvent", fRunNo.Data(), "After cuts");
  // Pileup
  fIsPileup = fAOD->IsPileupFromSPDInMultBins();
  if(fRejectPileup && fIsPileup) return;
  fHistos->FillTH1("EventStats", kAfterPileUp);
  fHistos->FillTH2("Runwise/NEvent", fRunNo.Data(), "Analysis");
  FillHistogramsForEventQA("Event_afterCuts");
  
  // QA
  if(!fIsTriggerQA)
    FillHistogramsForClusterQA("Cluster");
  else{
    if(offlineTrigger & AliVEvent::kINT7)
      FillHistogramsForClusterQA("Cluster_MB");
    if(offlineTrigger & AliVEvent::kEMCEGA){
      if(fFiredTriggerTag.Contains("EG1"))
        FillHistogramsForClusterQA("Cluster_EG1");
      if(fFiredTriggerTag.Contains("EG2"))
        FillHistogramsForClusterQA("Cluster_EG2");
      if(fFiredTriggerTag.Contains("DG1"))
        FillHistogramsForClusterQA("Cluster_DG1");
      if(fFiredTriggerTag.Contains("DG2"))
        FillHistogramsForClusterQA("Cluster_DG2");
    }
  }

  // Reset data members
  fJpsiPair = NULL;
  fTaggedJet = NULL;
  fJpsiMC = NULL;
  fTaggedJetMC = NULL;

  // Process Data
  RunDetectoreLevelAnalysis();

  // Process MC truth
  if(fIsMC && fHistosMC) RunParticleLevelAnalysis();

  // Fill event level variables
  AliDielectronVarManager::Fill(fAOD, fVars);
  FillTH2("Runwise/NPVtxZ", fRunNo.Data(), fVars[AliDielectronVarManager::kZvPrim]);
  FillTH2("Runwise/NVContrib", fRunNo.Data(), fVars[AliDielectronVarManager::kNVtxContrib]);
  FillTH2("Runwise/NSPDtracklets", fRunNo.Data(), fVars[AliDielectronVarManager::kNaccTrcklts10Corr]);
  FillTH2("Runwise/NTracksAll", fRunNo.Data(), fVars[AliDielectronVarManager::kNTrk]);
  FillTH2("Runwise/NElectrons", fRunNo.Data(), fVars[AliDielectronVarManager::kTracks]);
  FillTH2("Runwise/NPairs", fRunNo.Data(), fPairs->GetEntries());
  // Fill jet level variables
  auto jets = static_cast<AliJetContainer*>(fJets->At(0));
  FillTH2("Runwise/NJets", fRunNo.Data(), jets->GetNAcceptedJets());
  auto jetFinder = static_cast<AliEmcalJetTask*>(fJetTasks->At(0));
  auto trkCont = jetFinder->GetTrackContainer();
  FillTH2("Runwise/Jet_NTracks", fRunNo.Data(), trkCont->GetNAcceptedTracks());

  // Fill electron level variables
  auto trkArr = fDielectron->GetTrackArray(0);
  for(int iTrk = 0; iTrk < trkArr->GetEntries(); iTrk++){
    auto trk = (AliAODTrack*)(trkArr->At(iTrk));
    FillHistogramsForRunwiseElectronQA("Runwise", trk);
  }
  trkArr = fDielectron->GetTrackArray(1);
  for(int iTrk = 0; iTrk < trkArr->GetEntries(); iTrk++){
    auto trk = (AliAODTrack*)(trkArr->At(iTrk));
    FillHistogramsForRunwiseElectronQA("Runwise", trk);
  }
}

Bool_t AliAnalysisTaskJpsiJet::RunDetectoreLevelAnalysis(){
  // Run jet finder tasks
  if(RunJetFinder("Jet"))
    FillHistogramsForJetQA("Jet");

  // Run dielectron task
  AliKFParticle::SetField(fAOD->GetMagneticField());
  AliDielectronPID::SetCorrVal(fAOD->GetRunNumber());
  fDielectron->Process(fAOD);

  // Build tracks with dielectron pair
  Int_t nCandidates = fDielectron->GetPairArray(1)->GetEntriesFast();
  if(nCandidates > 0){
    fHistos->FillTH2("Runwise/NEvent", fRunNo.Data(), "N_{pair}>0");
    if(nCandidates == 1){
      fHistos->FillTH1("EventStats", kWithSinglePair);
      fHistos->FillTH2("Runwise/NEvent", fRunNo.Data(), "N_{pair}=1");
    }else if(nCandidates > 1){
      fHistos->FillTH1("EventStats", kWithMultiPair);
      fHistos->FillTH2("Runwise/NEvent", fRunNo.Data(), "N_{pair}>1");
    }
  }else{return kFALSE;}
  // Pairs and daughters
  fPairs->Clear("C");
  fDaughters->Clear("C");
  Int_t nD = 0;
  const TObjArray *candidates = fDielectron->GetPairArray(1);
  for (Int_t i = 0; i < nCandidates; i++)
  {
    AliDielectronPair *pair = (AliDielectronPair *)(candidates->UncheckedAt(i));
    new ((*fPairs)[i]) AliDielectronPair(*pair);
    const AliKFParticle &d1 = pair->GetKFFirstDaughter();
    new ((*fDaughters)[nD++]) TLorentzVector(d1.GetPx(), d1.GetPy(), d1.GetPz(), d1.GetE());
    const AliKFParticle &d2 = pair->GetKFSecondDaughter();
    new ((*fDaughters)[nD++]) TLorentzVector(d2.GetPx(), d2.GetPy(), d2.GetPz(), d2.GetE());
  }
  // Fill tracks
  fTracksWithPair->Clear("C");
  TIter nextTrk(fAOD->GetTracks());
  AliAODTrack* trk = NULL;
  AliAODTrack* trkTemplate = NULL;
  Int_t nDaughters = 0;
  Int_t nTracks = 0;
  while (( trk = static_cast<AliAODTrack*>(nextTrk()) )) {
    if(FindDaughters(trk)){
      nDaughters ++;
      if(!trkTemplate || trkTemplate->Pt() < trk->Pt())
        trkTemplate = trk;
      continue;
    }
    new ((*fTracksWithPair)[nTracks++]) AliAODTrack(*trk);
  }
  AliDebug(1, Form("Find dielectron pairs : %d, daughters : %d", nCandidates, nDaughters));
  // Insert pair as AOD track
  AddTrackFromPair(trkTemplate);
  AliDebug(1, Form("AOD tracks : %d, tracks with pair : %d", fAOD->GetNumberOfTracks(), fTracksWithPair->GetEntriesFast()));
  // Register in AOD event
  fAOD->AddObject(fTracksWithPair);
  
  // Run jet finder tasks on J/psi embedded track array
  if(RunJetFinder("JpsiJet"))
    FillHistogramsForJetQA("JpsiJet");

  if(FillHistogramsForTaggedJet("PairInJet")){
    fHistos->FillTH1("EventStats",kWithPairInJet);
    fHistos->FillTH2("Runwise/NEvent", fRunNo.Data(), "JetTagged");
  }
  
  return kTRUE;
}

void AliAnalysisTaskJpsiJet::InitHistogramsForRunwiseQA(const char* histClass){
  fHistos->CreateTH2(
    Form("%s/NEvent", histClass),
    "Runwise QA - Event Number",
    201, -0.5, 200.5,
    int(kEventStatusN)+1, -0.5, int(kEventStatusN)+0.5);
  fHistos->CreateTH2(
    Form("%s/NPVtxZ", histClass),
    "Runwise QA - Primary vertex Z",
    201, -0.5, 200.5,
    300, -15., 15.);
  fHistos->CreateTH2(
    Form("%s/NVContrib", histClass),
    "Runwise QA - Number of Vertex contributors",
    201, -0.5, 200.5,
    1001, -0.5, 1000.5);
  fHistos->CreateTH2(
    Form("%s/NSPDtracklets", histClass),
    "Runwise QA - Number of accepted SPD tracklets in |eta|<1.0",
    201, -0.5, 200.5,
    501, -0.5, 500.5);
  fHistos->CreateTH2(
    Form("%s/NTracksAll", histClass),
    "Runwise QA - Number of all tracks",
    201, -0.5, 200.5,
    4001, -0.5, 4000.5);
  fHistos->CreateTH2(
    Form("%s/NElectrons", histClass),
    "Runwise QA - Number of selected tracks/electron",
    201, -0.5, 200.5,
    51, -0.5, 50.5);
  fHistos->CreateTH2(
    Form("%s/NPairs", histClass),
    "Runwise QA - Number of selected tracks/electron",
    201, -0.5, 200.5,
    21, -0.5, 20.5);

  // Jet level
  fHistos->CreateTH2(
    Form("%s/NJets", histClass),
    "Runwise Jet QA - Number of charged jets (R=0.4, |#eta|<0.5)",
    201, -0.5, 200.5,
    101, -0.5, 100.5);
  fHistos->CreateTH2(
    Form("%s/Jet_NTracks", histClass),
    "Runwise Jet QA - Number of accepted tracks (Hybrid)",
    201, -0.5, 200.5,
    1001, -0.5, 1000.5);

  // Cluster level
  fHistos->CreateTH2(
    Form("%s/EMCal_NClusters", histClass),
    "Runwise EMCal QA - Number of clusters",
    201, -0.5, 200.5,
    1001, -0.5, 1000.5);
  fHistos->CreateTH2(
    Form("%s/EMCal_NClsMatched", histClass),
    "Runwise EMCal QA - Number of clusters (track matched)",
    201, -0.5, 200.5,
    1001, -0.5, 1000.5);
  fHistos->CreateTH2(
    Form("%s/EMCal_E", histClass),
    "Runwise EMCal QA - Cluster energy (All)",
    201, -0.5, 200.5,
    200, 0., 100.);
  fHistos->CreateTH2(
    Form("%s/EMCal_Etrk", histClass),
    "Runwise EMCal QA - Cluster energy (track matched)",
    201, -0.5, 200.5,
    200, 0., 100.);

  // Electron level
  fHistos->CreateTH2(
    Form("%s/Ele_Pt", histClass),
    "Runwise Electron QA - p_{T}",
    201, -0.5, 200.5,
    200, 0., 100.);
  fHistos->CreateTH2(
    Form("%s/Ele_Eta", histClass),
    "Runwise Electron QA - #eta",
    201, -0.5, 200.5,
    100, -1., 1.);
  fHistos->CreateTH2(
    Form("%s/Ele_Phi", histClass),
    "Runwise Electron QA - #phi",
    201, -0.5, 200.5,
    160, -1., 7.);
  fHistos->CreateTH2(
    Form("%s/Ele_DCAxy", histClass),
    "Runwise Electron QA - DCA_{xy}",
    201, -0.5, 200.5,
    1000, -5., 5.);
  fHistos->CreateTH2(
    Form("%s/Ele_DCAz", histClass),
    "Runwise Electron QA - DCA_{z}",
    201, -0.5, 200.5,
    1000, -5., 5.);
  fHistos->CreateTH2(
    Form("%s/Ele_TPCNcls", histClass),
    "Runwise Electron QA - Number of TPC clusters",
    201, -0.5, 200.5,
    161, -0.5, 160.5);
  fHistos->CreateTH2(
    Form("%s/Ele_TPCChi2", histClass),
    "Runwise Electron QA - TPC #chi^{2}/N_{cluster}",
    201, -0.5, 200.5,
    100, 0., 10.);
  fHistos->CreateTH2(
    Form("%s/Ele_TPCNsigmaEle", histClass),
    "Runwise Electron QA - PID n_{#sigma e}^{TPC}",
    201, -0.5, 200.5,
    100, -5., 5.);
  fHistos->CreateTH2(
    Form("%s/Ele_EMCalE", histClass),
    "Runwise Electron QA - Cluster energy",
    201, -0.5, 200.5,
    200, 0., 100.);
  fHistos->CreateTH2(
    Form("%s/Ele_EMCalEP", histClass),
    "Runwise Electron QA - EMCal E/p",
    201, -0.5, 200.5,
    100, 0., 2.);
  fHistos->CreateTH2(
    Form("%s/Ele_EMCalNsigmaEle", histClass),
    "Runwise Electron QA - PID n_{#sigma e}^{EMC}",
    201, -0.5, 200.5,
    100, -5., 5.);
}

void AliAnalysisTaskJpsiJet::FillTH2(const char* histName, const char* labelX, Double_t value, Double_t weight){
  TH2* h2 = (TH2*)(fHistos->FindObject(histName));
  if(!h2){
    AliFatal(Form("Fail to find %s in Hitogram Manager", histName));
    return;
  }
  h2->Fill(labelX, value, weight);
}

void AliAnalysisTaskJpsiJet::FillHistogramsForRunwiseElectronQA(const char* histClass, AliAODTrack* ele){
  Double_t trkValues[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(ele, trkValues);
  FillTH2(Form("%s/Ele_Pt", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kPt]);
  FillTH2(Form("%s/Ele_Eta", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kEta]);
  FillTH2(Form("%s/Ele_Phi", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kPhi]);
  FillTH2(Form("%s/Ele_DCAxy", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kImpactParXY]);
  FillTH2(Form("%s/Ele_DCAz", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kImpactParZ]);
  FillTH2(Form("%s/Ele_TPCNcls", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kNclsTPC]);
  FillTH2(Form("%s/Ele_TPCChi2", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kTPCchi2Cl]);
  FillTH2(Form("%s/Ele_TPCNsigmaEle", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kTPCnSigmaEle]);
  FillTH2(Form("%s/Ele_EMCalE", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kEMCALE]);
  FillTH2(Form("%s/Ele_EMCalEP", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kEMCALEoverP]);
  FillTH2(Form("%s/Ele_EMCalNsigmaEle", histClass), fRunNo.Data(), trkValues[AliDielectronVarManager::kEMCALnSigmaEle]);
}

// Create event QA histograms in output list
void AliAnalysisTaskJpsiJet::InitHistogramsForEventQA(const char* histClass){

  TH1* hTrigger = fHistos->CreateTH1(Form("%s/Trigger", histClass),"Number of event by offline triggers;Nbits in AliVEvent;N_{events}",32,-0.5,31.5);
  const char* labelOfflineTrigger[32] = {
    "MB/INT1", "INT7", "MUON", "HM", "EMC1", "INT5", "MUS", "MUSH7",
    "MUL7", "MUU7", "EMC7/8", "MUS7", "PHI1", "PHI7/8", "EMCEJE", "EMCEGA",
    "HMV0/CENT", "SemiCENT", "DG/5", "ZED", "SPI/7", "INT8", "MUSL8", "MUSH8",
    "MULL8", "MUUL8", "MUUL0", "UserDef", "TRD", "MUCalo", "FastOnly", ""};
  for(int i = 0; i < 32; i++)
    hTrigger->GetXaxis()->SetBinLabel(i+1, labelOfflineTrigger[i]);

  TH1* hTriggerClass = fHistos->CreateTH1(Form("%s/TriggerClass", histClass),"Number of event by fired trigger class;Trig. Descriptor;N_{events}",10,-0.5,9.5);

  fHistos->CreateTH1(Form("%s/FiredTag", histClass),"Number of event by fired trigger tag;Trig. Tag;N_{events}",32,-0.5,31.5);
}

void AliAnalysisTaskJpsiJet::FillHistogramsForEventQA(const char* histClass){
  // Offline Trigger
  AliAODHeader* header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  UInt_t offlineTrigger = header->GetOfflineTrigger();
  for(Short_t i = 0; i < 32; i++){
    if(offlineTrigger & BIT(i)) fHistos->FillTH1(Form("%s/Trigger", histClass), i);
  }

  // Trigger Classes
  TString triggerClass = fAOD->GetFiredTriggerClasses();
  TObjArray* tcArray = triggerClass.Tokenize(" ");
  for(Short_t i = 0; i < tcArray->GetEntries(); i++){
    TString strClass = tcArray->At(i)->GetName();
    TObjArray* tmp = strClass.Tokenize("-");
    strClass = tmp->At(0)->GetName();
    fHistos->FillTH1(Form("%s/TriggerClass", histClass), strClass.Data());
    tmp->SetOwner(kTRUE);
    delete tmp;
  }
  tcArray->SetOwner(kTRUE);
  delete tcArray;

  // Trigger Tag
  fHistos->FillTH1(Form("%s/FiredTag", histClass), fFiredTriggerTag.Data());
}


// Cluster QA - EMCal ONLY
void AliAnalysisTaskJpsiJet::InitHistogramsForClusterQA(const char* histClass){
  // Cluster energy, and with Ncell, Nmatched tracks, dispersion
  fHistos->CreateTH1(Form("%s/E", histClass), "Cluster energy distribution;Energy (GeV);N_{cls};", 500, 0., 100.);
  fHistos->CreateTH2(Form("%s/E_Ncell", histClass), "Cluster cell number vs energy;Energy (GeV);N_{cell};", 500, 0., 100., 101, -0.5, 100.5);
  fHistos->CreateTH2(Form("%s/E_Ntrk", histClass), "Cluster matched tracks vs energy;Energy (GeV);N_{trk};", 500, 0., 100., 21, -0.5, 20.5);
  // Shower shape parameters - Dispersion, M02, M20
  fHistos->CreateTH2(Form("%s/E_Dispersion", histClass), "Cluster dispersion vs energy;Energy (GeV);Dispersion;", 500, 0., 100., 1000, 0., 10.);
  fHistos->CreateTH2(Form("%s/M02_M20", histClass), "Cluster shape parameter - M02 vs M20;M02;M20;", 200, 0., 20., 500, 0., 5.);

  // Cell QA
    // NOTICE: Large memory would be allocated here
  if(fIsCellQA){
    // E vs Cell ID (0-18000 for EMCal)
    fHistos->CreateTH2(Form("%s/CellEnergy", histClass), "Cluster Cell Energy;CellID;E_{cell}", 18001, -0.5, 18000.5, 100, 0., 20.);
    // Time vs Cell ID
    fHistos->CreateTH2(Form("%s/CellTime", histClass), "Cluster Cell Time;CellID;T_{cell} (#mus);", 18001, -0.5, 18000.5, 300, -1.5, 1.5);
  }
}

void AliAnalysisTaskJpsiJet::FillHistogramsForClusterQA(const char* histClass){
  Int_t nClsMatched = 0;
  FillTH2("Runwise/EMCal_NClusters", fRunNo.Data(), fAOD->GetNumberOfCaloClusters());
  // Loop EMCal clusters
  for(int iCls = 0; iCls < fAOD->GetNumberOfCaloClusters(); iCls++){
    auto cls = (AliAODCaloCluster*)(fAOD->GetCaloCluster(iCls));
    if(!cls->IsEMCAL()) continue;
    Double_t energy      = cls->E();
    Float_t  pos[3]      = {0.}; cls->GetPosition(pos);
    Double_t nCell       = cls->GetNCells();
    Double_t nTracks     = cls->GetNTracksMatched();
    Double_t m02         = cls->GetM02();
    Double_t m20         = cls->GetM20();
    Double_t dispersion  = cls->GetDispersion();

    fHistos->FillTH1(Form("%s/E", histClass), energy);
    fHistos->FillTH2(Form("%s/E_Ncell", histClass), energy, nCell);
    fHistos->FillTH2(Form("%s/E_Ntrk", histClass), energy, nTracks);
    fHistos->FillTH2(Form("%s/E_Dispersion", histClass), energy, dispersion);
    fHistos->FillTH2(Form("%s/M02_M20", histClass), m02, m20);

    // Cell QA - Loop cells
    if(fIsCellQA){
      auto cells = fAOD->GetEMCALCells();
      for(int iCell = 0; iCell < nCell; iCell++){
        Double_t cellID    = cls->GetCellAbsId(iCell);
        Double_t cellE  = cells->GetCellAmplitude(cellID);
        Double_t cellT  = cells->GetCellTime(cellID);

        fHistos->FillTH1(Form("%s/CellEnergy", histClass), cellID, cellE);
        fHistos->FillTH1(Form("%s/CellTime", histClass), cellID, cellT * 1e6);
      }
    }// End - Loop cells

    // Runwise QA
    FillTH2("Runwise/EMCal_E", fRunNo.Data(), energy);
    if(nTracks > 0){
      FillTH2("Runwise/EMCal_Etrk", fRunNo.Data(), energy);
      nClsMatched ++;
    }
  }// End - Loop clusters
  FillTH2("Runwise/EMCal_NClsMatched", fRunNo.Data(), nClsMatched);
}


// Copy from AliEmcalJetTask::AddTaskEmcalJet
void AliAnalysisTaskJpsiJet::AddTaskEmcalJet(
  const TString nTracks, const TString nClusters,
  const AliJetContainer::EJetAlgo_t jetAlgo, const Double_t radius, const AliJetContainer::EJetType_t jetType,
  const Double_t minTrPt, const Double_t minClPt,
  const Double_t ghostArea, const AliJetContainer::ERecoScheme_t reco,
  const TString tag, const Double_t minJetPt,
  const Bool_t lockTask, const Bool_t bFillGhosts
){
  // Setup containers
  TString trackName(nTracks);
  TString clusName(nClusters);
  if (trackName == "usedefault") {
    trackName = "tracks";
  }
  if (clusName == "usedefault") {
      clusName = "caloClusters";
  }

  AliParticleContainer* partCont = 0;
  if (trackName.Contains("mcparticles")) {  // must be contains in order to allow for non-standard particle containers
    AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
    partCont = mcpartCont;
  }
  else if (trackName.Contains("tracks")) {
    AliTrackContainer* trackCont = new AliTrackContainer(trackName);
    partCont = trackCont;
  }
  else if (!trackName.IsNull()) {
    partCont = new AliParticleContainer(trackName);
  }
  if (partCont) partCont->SetParticlePtCut(minTrPt);

  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
    clusCont->SetClusHadCorrEnergyCut(minClPt);
    clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  switch (jetType) {
  case AliJetContainer::kChargedJet:
    if (partCont) partCont->SetCharge(AliParticleContainer::kCharged);
    break;
  case AliJetContainer::kNeutralJet:
    if (partCont) partCont->SetCharge(AliParticleContainer::kNeutral);
    break;
  default:
    break;
  }
  
  TString name = AliJetContainer::GenerateJetName(jetType, jetAlgo, reco, radius, partCont, clusCont, tag);
  ;
  if(fJetTasks->FindObject(name.Data())){
    AliWarning(Form("Jet task existed: %s", name.Data()));
    return;
  }
  AliInfo(Form("Jet task name : %s", name.Data()));

  AliEmcalJetTask* jetTask = new AliEmcalJetTask(name);
  fJetTasks->Add(jetTask);

  jetTask->SetJetType(jetType);
  jetTask->SetJetAlgo(jetAlgo);
  jetTask->SetRecombScheme(reco);
  jetTask->SetRadius(radius);
  if (partCont) jetTask->AdoptParticleContainer(partCont);
  if (clusCont) jetTask->AdoptClusterContainer(clusCont);
  jetTask->SetJetsName(tag);
  jetTask->SetMinJetPt(minJetPt);
  jetTask->SetGhostArea(ghostArea);

  if (bFillGhosts) jetTask->SetFillGhost();
  if (lockTask) jetTask->SetLocked();

  jetTask->SetForceBeamType(AliAnalysisTaskEmcal::kpp);
  jetTask->SelectCollisionCandidates(AliVEvent::kAny);
  jetTask->SetUseAliAnaUtils(kTRUE);
  jetTask->SetZvertexDiffValue(0.5);
  jetTask->SetNeedEmcalGeom(kFALSE);

  // Connect input
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  jetTask->ConnectInput(0, mgr->GetCommonInputContainer());

  // Add jet container
  AliJetContainer *cont = new AliJetContainer(jetType, jetAlgo, reco, radius, partCont, clusCont, tag);
  cont->SetPercAreaCut(0.6);
  if(jetType == AliJetContainer::kChargedJet)
    cont->SetJetAcceptanceType(AliJetContainer::kTPCfid);
  // else for neutral jet, EMCal/DCal/PHOS should be considered
  fJets->Add(cont);
}

void AliAnalysisTaskJpsiJet::InitJetFinders(){
  if(!fJetTasks) fJetTasks = new TObjArray;
  fJetTasks->SetOwner(kTRUE);
  if(!fJets) fJets = new TObjArray;
  fJets->SetOwner(kTRUE);

  AddTaskEmcalJet("usedefault", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.3, 0.01, AliJetContainer::pt_scheme, "Jet", 1., kFALSE, kFALSE);
  AddTaskEmcalJet("tracksWithPair", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.3, 0.01, AliJetContainer::pt_scheme, "JpsiJet", 1., kFALSE, kFALSE);

  if(fIsMC){
    AddTaskEmcalJet("mcparticles", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.3, 0.01, AliJetContainer::pt_scheme, "JetMC", 1., kFALSE, kFALSE);
    AddTaskEmcalJet("mcparticles", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.3, 0.01, AliJetContainer::pt_scheme, "JpsiJetMC", 1., kFALSE, kFALSE);
  }
}

Bool_t AliAnalysisTaskJpsiJet::RunJetFinder(const char* jetTag){
  if(!fIsJetFinder) return kTRUE;
  AliEmcalJetTask* jetFinder = NULL;
  TIter next(fJetTasks);
  while((jetFinder=(AliEmcalJetTask*)next())){
    TString jetName = jetFinder->GetName();
    // Find jet task by name
    if(jetName.BeginsWith(Form("%s_", jetTag))){
      // Prosess as AliAnalysisManager::StartAnalysis
      jetFinder->Reset();
      jetFinder->SetActive(kTRUE);
      jetFinder->Exec("");
      return kTRUE;
    }
  }
  return kFALSE;
}

void AliAnalysisTaskJpsiJet::InitHistogramsForJetQA(){

  AliJetContainer* jets = NULL;
  TIter next(fJets);
  while((jets = static_cast<AliJetContainer*>(next()))){
    TString jetName = jets->GetName();
    // Skip MC jets
    if(jetName.Contains("mcparticles")) continue;
    // THnSparse - pT, eta， phi
    Int_t nBins[3]   = {2000, 200, 100};
    Double_t xmin[3] = {0.,   -1., -2.};
    Double_t xmax[3] = {100.,  1.,  8.};
    fHistos->CreateTHnSparse(Form("%s/jetVars", jets->GetName()), "Jet kinetic variables (p_{T}-#eta-#phi);p_{T} (GeV/c);#eta;#phi;", 3, nBins, xmin, xmax);
      // Constituents
    fHistos->CreateTH2(Form("%s/Ntracks_pT",jets->GetName()), "Jet constituents - number vs p_{T}^{jet};p_{T}^{jet} (GeV/c);N_{tracks};",
      200, 0., 100., 100, -0.5, 99.5);
    fHistos->CreateTH2(Form("%s/Area_pT",jets->GetName()), "Jet clustering area;p_{T}^{jet} (GeV/c);A_{jet};",
      200, 0., 100., 100, 0., 2.);
  }
}

void AliAnalysisTaskJpsiJet::FillHistogramsForJetQA(const char* jetTag){
  
  AliJetContainer* jets = NULL;
  TIter next(fJets);
  while((jets = static_cast<AliJetContainer*>(next()))){
    TString jetName = jets->GetName();
    if (!jetName.BeginsWith(Form("%s_", jetTag)))
      continue;

    // Found jet with tag
    jets->NextEvent(fAOD);
    jets->SetArray(fAOD);

    AliDebug(1, Form("%d jets found in %s", jets->GetNJets(), jets->GetName()));
    for (auto jet : jets->all())
    {
      // Jet cuts
      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason))
      {
        AliDebug(2, Form("Jet was rejected for reason : %d, details : %s", rejectionReason, (jet->toString()).Data()));
        continue;
      }
      Double_t x[3] = {0.};
      x[0] = jet->Pt();
      x[1] = jet->Eta();
      x[2] = jet->Phi();
      fHistos->FillTHnSparse(Form("%s/jetVars", jets->GetName()), x, 1.0);
      fHistos->FillTH2(Form("%s/Ntracks_pT",jets->GetName()), jet->Pt(), jet->GetNumberOfTracks());
      fHistos->FillTH2(Form("%s/Area_pT",jets->GetName()), jet->Pt(), jet->Area());
    } // End - Loop jets
  }// End - Loop jet containers
}

void AliAnalysisTaskJpsiJet::LocalInit(){
  // Dielectron task
  InitDielectron();
  // Jet task
  InitJetFinders();

  // Init call as AliAnalysisManager
  AliEmcalJetTask* jetFinder = NULL;
  TIter next(fJetTasks);
  while((jetFinder=(AliEmcalJetTask*)next()))
    jetFinder->LocalInit();
}

void AliAnalysisTaskJpsiJet::Terminate(Option_t*){
  // Terminate call as AliAnalysisManager
  AliEmcalJetTask* jetFinder = NULL;
  TIter next(fJetTasks);
  while((jetFinder=(AliEmcalJetTask*)next()))
    jetFinder->Terminate();
}

void AliAnalysisTaskJpsiJet::InitDielectron(){

  fDielectron = new AliDielectron("Diele","Dielectron with EMCal triggered");
  fPairs = new TClonesArray("AliDielectronPair",10);
  fPairs->SetName("dielectrons");
  fDaughters = new TClonesArray("TLorentzVector",20);
  fDaughters->SetName("daughters");
  fTracksWithPair = new TClonesArray("AliAODTrack",500);
  fTracksWithPair->SetName("tracksWithPair");

/**
 *  Track cuts
**/
  // Good track
  AliDielectronTrackCuts *trackCuts = new AliDielectronTrackCuts("trackCuts", "trackCuts");
  trackCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetRequireITSRefit(kTRUE);
  fDielectron->GetTrackFilter().AddCuts(trackCuts);
  // Track cuts for electron PID
  AliDielectronVarCuts *ePID = new AliDielectronVarCuts("ePidCut", "Track cuts for electron PID");
  ePID->AddCut(AliDielectronVarManager::kKinkIndex0, 0.);
  // Track quality cuts
  ePID->AddCut(AliDielectronVarManager::kNclsTPC, 85., 160.);
  ePID->AddCut(AliDielectronVarManager::kEta, -0.9, 0.9);
  ePID->AddCut(AliDielectronVarManager::kImpactParXY, -1., 1.);
  ePID->AddCut(AliDielectronVarManager::kImpactParZ, -3., 3.);
  ePID->AddCut(AliDielectronVarManager::kITSLayerFirstCls, 0., 4.);
  ePID->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0., 4.);
  ePID->AddCut(AliDielectronVarManager::kPt, 1.0, 1e30);
  // Electron PID with TPC
  ePID->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -2.0, 3.0);
  // Exclude hadrons
  ePID->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -100.0, 1.0, kTRUE);
  ePID->AddCut(AliDielectronVarManager::kTPCnSigmaKao, -100.0, 3.0, kTRUE);
  ePID->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -100.0, 3.0, kTRUE);
  fDielectron->GetTrackFilter().AddCuts(ePID);

/**
 *  Pair cuts
**/
  // EMCal gamma/electron energy threshold
    // TODO: now the inteface ONLY used for ALL and MC
  Double_t EMCal_E_Threshold = fEMCalEth; // GeV, ADC setting
  if(fSelectedTriggerClasses == "EG1" || fSelectedTriggerClasses == "DG1")
    EMCal_E_Threshold = 10.0;
  if(fSelectedTriggerClasses == "EG2" || fSelectedTriggerClasses == "DG2")
    EMCal_E_Threshold = 5.0;

  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut = new AliDielectronVarCuts("jpsiCuts", "1<M<5 + |Y|<.9");
  pairCut->AddCut(AliDielectronVarManager::kM, 1.0, 5.0);
  pairCut->AddCut(AliDielectronVarManager::kY, -0.9, 0.9);
  pairCut->AddCut(AliDielectronVarManager::kPt, 1, 1e30);
  pairCut->AddCut(AliDielectronVarManager::kPt, EMCal_E_Threshold, 1e30);
  fDielectron->GetPairFilter().AddCuts(pairCut);
  // Leg cuts
  AliDielectronVarCuts *emcCut = new AliDielectronVarCuts("CutEMCAL", "Jpsi leg cuts for EMCal");
  emcCut->AddCut(AliDielectronVarManager::kEMCALEoverP, 0.8, 1.3);
  emcCut->AddCut(AliDielectronVarManager::kEMCALE, EMCal_E_Threshold, 1e30);
  //emcCut->AddCut(AliDielectronVarManager::kPhi, 4.377, 5.7071, kTRUE); // Exclude DCal
  //emcCut->AddCut(AliDielectronVarManager::kPhi, 1.396, 3.2637, kTRUE); // Exclude EMCal
  AliDielectronPairLegCuts *legVar = new AliDielectronPairLegCuts();
  legVar->GetLeg1Filter().AddCuts(emcCut);
  legVar->GetLeg2Filter().AddCuts(emcCut);
  legVar->SetCutType(AliDielectronPairLegCuts::kAnyLeg);
  fDielectron->GetPairFilter().AddCuts(legVar);
}

void AliAnalysisTaskJpsiJet::InitHistogramsForDielectron(const char* histMgrName)
{
  //Setup histogram Manager
  AliDielectronHistos *histos = new AliDielectronHistos(histMgrName, "Histograms for dielectron");
  fDielectron->SetHistogramManager(histos);

  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  //Track classes
  for (Int_t i = 0; i < 2; ++i)
  {
    histos->AddClass(Form("Track_%s", AliDielectron::TrackClassName(i)));
  }
  //Pair classes
  for (Int_t i = 0; i < 3; ++i)
    histos->AddClass(Form("Pair_%s", AliDielectron::PairClassName(i)));
  //Legs from pair
  for (Int_t i = 0; i < 3; ++i)
    histos->AddClass(Form("Track_Legs_%s", AliDielectron::PairClassName(i)));

/**
 *  Event
**/
  TString histClass = "Event";
  histos->AddClass(histClass);
  // Event Primary vertex and diamond (IP) stats.
  histos->UserHistogram(histClass, "VtxZ", "Vertex Z;Z[cm];#events", 1000, -50., 50., AliDielectronVarManager::kZvPrim);
  histos->UserHistogram(histClass, "VtxX", "Vertex X;X[cm];#events", 1000, -1.0, 1.0, AliDielectronVarManager::kXvPrim);
  histos->UserHistogram(histClass, "VtxY", "Vertex Y;Y[cm];#events", 2000, -1.0, 1.0, AliDielectronVarManager::kYvPrim);
  histos->UserHistogram(histClass, "NVContrib", "Number of Vertex contributors", 1001, -0.5, 1000.5, AliDielectronVarManager::kNVtxContrib);
  // Event track and SPD (tracklets) stats.
  histos->UserHistogram(histClass, "kNTrk", "Number of tracks;kNTrk;Entries", 4001, -0.5, 4000.5, AliDielectronVarManager::kNTrk);
  histos->UserHistogram(histClass, "kNaccTrcklts", "Number of accepted SPD tracklets in |eta|<1.6;kNaccTrcklts;Entries", 1001, -0.5, 1000.5, AliDielectronVarManager::kNaccTrcklts);
  histos->UserHistogram(histClass, "kNaccTrcklts10Corr", "kNaccTrcklts10Corr;kNaccTrcklts10Corr;Entries", 501, -0.5, 500.5, AliDielectronVarManager::kNaccTrcklts10Corr);
  histos->UserHistogram(histClass, "VtxZ_kNaccTrcklts10Corr", "VtxZ vs. kNaccTrcklts10Corr;VtxZ;kNaccTrcklts10Corr", 800, -40., 40., 501, -0.5, 500.5, AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kNaccTrcklts10Corr);
  //new multiplicity estimator: V0
  histos->UserHistogram(histClass, "kMultV0", "kMultV0;kMultV0;Entries", 1000, 0., 1000., AliDielectronVarManager::kMultV0);
  histos->UserHistogram(histClass, "kMultV0A", "kMultV0;kMultV0;Entries", 1000, 0., 1000., AliDielectronVarManager::kMultV0A);
  histos->UserHistogram(histClass, "kMultV0C", "kMultV0;kMultV0;Entries", 1000, 0., 1000., AliDielectronVarManager::kMultV0C);
  // 2D
  histos->UserHistogram(histClass, "kMultV0A_kMultV0C", "kMultV0A vs. kMultV0C;kMultV0A;kMultV0C", 1000, 0., 1000., 1000, 0., 1000., AliDielectronVarManager::kMultV0A, AliDielectronVarManager::kMultV0C);
  histos->UserHistogram(histClass, "kMultV0_kMultV0A", "kMultV0 vs. kMultV0A;kMultV0;kMultV0A", 1000, 0., 1000., 1000, 0., 1000., AliDielectronVarManager::kMultV0, AliDielectronVarManager::kMultV0A);
  histos->UserHistogram(histClass, "kMultV0_kMultV0C", "kMultV0 vs. kMultV0C;kMultV0;kMultV0C", 1000, 0., 1000., 1000, 0., 1000., AliDielectronVarManager::kMultV0, AliDielectronVarManager::kMultV0C);
  // vs Vertex Z
  histos->UserHistogram(histClass, "VtxZ_kMultV0A", "VtxZ vs. kMultV0A;VtxZ;kMultV0A", 300, -15., 15., 1000, 0., 1000., AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kMultV0A);
  histos->UserHistogram(histClass, "VtxZ_kMultV0C", "VtxZ vs. kMultV0C;VtxZ;kMultV0C", 300, -15., 15., 1000, 0., 1000., AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kMultV0C);
  histos->UserHistogram(histClass, "VtxZ_kMultV0", "VtxZ vs. kMultV0;VtxZ;kMultV0", 300, -15., 15., 1000, 0., 1000., AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kMultV0);

  // Dielectron info.
  histos->UserHistogram(histClass, "Nelectrons", "Number of tracks/electron selected by AliDielectron after cuts;N_{e};#events", 50, -0.5, 49.5, AliDielectronVarManager::kTracks);
  histos->UserHistogram(histClass, "Npairs", "Number of Ev1PM pair candidates after all cuts;J/#psi candidates;#events", 20, -0.5, 19.5, AliDielectronVarManager::kPairs);
  /*
	  Histogram for Track
	*/
  // Track kinetics parameter
  histos->UserHistogram("Track", "Pt", "Pt;Pt [GeV/c];#tracks", 2000, 0, 100, AliDielectronVarManager::kPt, kTRUE);
  histos->UserHistogram("Track", "Eta_Phi", "Eta Phi Map; Eta; Phi;#tracks",
                        100, -1, 1, 144, 0, TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi, kTRUE);
  histos->UserHistogram("Track", "dXY", "dXY;dXY [cm];#tracks", 1000, -50, 50, AliDielectronVarManager::kImpactParXY, kTRUE);
  histos->UserHistogram("Track", "dXY_Pt", "DCA_{xy} vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm);#tracks;",
                        200, 0, 100., 1000, -5., 5., AliDielectronVarManager::kPt, AliDielectronVarManager::kImpactParXY, kTRUE);
  histos->UserHistogram("Track", "dZ", "dZ;dZ [cm];#tracks", 1000, -50., 50., AliDielectronVarManager::kImpactParZ, kTRUE);
  histos->UserHistogram("Track", "dZ_Pt", "DCA_{z} vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm);#tracks;",
                        200, 0, 100., 1000, -5., 5., AliDielectronVarManager::kPt, AliDielectronVarManager::kImpactParZ, kTRUE);
  // Tracking quality
  histos->UserHistogram("Track", "ITS_FirstCls", "ITS First Cluster;Layer No. of ITS 1st cluster;#Entries", 6, 0.5, 6.5, AliDielectronVarManager::kITSLayerFirstCls, kTRUE);
  histos->UserHistogram("Track", "TPCnCls", "Number of Clusters TPC;TPC number clusteres;#tracks", 161, -0.5, 160.5, AliDielectronVarManager::kNclsTPC, kTRUE);
  histos->UserHistogram("Track", "TPCchi2Cl", "Chi-2/Clusters TPC;Chi2/ncls number clusteres;#tracks", 100, 0, 10, AliDielectronVarManager::kTPCchi2Cl, kTRUE);
  // PID - TPC
  histos->UserHistogram("Track", "dEdx_P", "dEdx vs PinTPC;P [GeV];TPC signal (a.u.);#tracks",
                        800, 0., 40., 800, 20., 200., AliDielectronVarManager::kPIn, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "dEdx_Pt", "dEdx vs Pt;Pt [GeV];TPC signal (a.u.);#tracks",
                        800, 0., 40., 800, 20., 200., AliDielectronVarManager::kPt, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "dEdx_Phi", "dEdx vs #phi;#phi [rad];TPC signal (a.u.);#tracks",
                        200, 0., 2 * TMath::Pi(), 800, 20., 200., AliDielectronVarManager::kPhi, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "dEdx_Eta", "dEdx vs #eta;#eta;TPC signal (a.u.);#tracks",
                        200, -1., 1., 800, 20., 200., AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_P", "n#sigma_{e}(TPC) vs P_{in} TPC;P_{in} [GeV];n#sigma_{e}(TPC);#tracks",
                        800, 0., 40., 800, -12., 12., AliDielectronVarManager::kPIn, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_Pt", "n#sigma_{e}(TPC) vs Pt;Pt [GeV];n#sigma_{e}(TPC);#tracks",
                        800, 0., 40., 800, -12., 12., AliDielectronVarManager::kPt, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_Phi", "n#sigma_{e}(TPC) vs #phi;#phi [rad];n#sigma_{e}(TPC);#tracks",
                        200, 0., 2 * TMath::Pi(), 800, -12., 12., AliDielectronVarManager::kPhi, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_Eta", "n#sigma_{e}(TPC) vs #eta;#eta;n#sigma_{e}(TPC);#tracks",
                        200, -1., 1., 800, -12., 12., AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "dEdx_nSigmaEMCal", "dEdx vs n#sigma_{e}(EMCAL);n#sigma_{e}(EMCAL);TPC signal (a.u.);#tracks",
                        200, -5., 5., 800, 20., 200., AliDielectronVarManager::kEMCALnSigmaEle, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "dEdx_TPCnSigmaEle", "dEdx vs n#sigma_{e}(TPC);n#sigma_{e}(TPC);TPC signal (a.u.);#tracks",
                        100, -10., 10., 800, 20., 200., AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kTPCsignal, kTRUE);
  // Track - EMCal
  histos->UserHistogram("Track", "EMCalE", "EmcalE;Cluster Energy [GeV];#Clusters",
                        200, 0., 100., AliDielectronVarManager::kEMCALE, kTRUE);
  histos->UserHistogram("Track", "EMCalE_P", "Cluster energy vs. p_{IN}; EMCal_E (GeV);p_{IN} (GeV/c);#tracks",
                        800, 0., 40, 200, 0., 40.,  AliDielectronVarManager::kPIn, AliDielectronVarManager::kEMCALE, kTRUE);
  histos->UserHistogram("Track", "EMCalE_Pt", "Cluster energy vs. pT; EMCal_E;pT;#tracks",
                        800, 0., 40, 200, 0., 40.,  AliDielectronVarManager::kPt, AliDielectronVarManager::kEMCALE, kTRUE);
  //Ecluster versus Phi to separate EMCal and DCal
  histos->UserHistogram("Track", "EMCalE_Phi", "Cluster energy vs. #phi; EMCal_E (GeV);Phi;#tracks",
                        200, 0., TMath::TwoPi(), 200, 0., 40., AliDielectronVarManager::kPhi, AliDielectronVarManager::kEMCALE, kTRUE);
  histos->UserHistogram("Track", "EMCalE_Eta", "Cluster energy vs. #eta; EMCal_E (GeV);Eta;#tracks",
                        200, -1.0, 1.0, 200, 0., 40., AliDielectronVarManager::kEta, AliDielectronVarManager::kEMCALE, kTRUE);
  // PID - EMCal
  // E/p ratio
  histos->UserHistogram("Track", "EoverP", "EMCal E/p ratio;E/p;#Clusters",
                        200, 0., 2., AliDielectronVarManager::kEMCALEoverP, kTRUE);
  histos->UserHistogram("Track", "EoverP_P", "E/p ratio vs P;P_{in} (GeV/c);E/p;#tracks",
                        200, 0., 40., 200, 0., 2., AliDielectronVarManager::kPIn, AliDielectronVarManager::kEMCALEoverP, kTRUE);  
  histos->UserHistogram("Track", "EoverP_Pt", "E/p ratio vs Pt;Pt (GeV/c);E/p;#tracks",
                        200, 0., 40., 200, 0., 2., AliDielectronVarManager::kPt, AliDielectronVarManager::kEMCALEoverP, kTRUE);
  histos->UserHistogram("Track", "EoverP_Phi", "E/p ratio vs #phi;Phi;E/p;#tracks",
                        200, 0., TMath::TwoPi(), 200, 0., 2., AliDielectronVarManager::kPhi, AliDielectronVarManager::kEMCALEoverP, kTRUE);
  histos->UserHistogram("Track", "EoverP_Eta", "E/p ratio vs #eta;Eta;E/p;#tracks",
                        200, -1.0, 1.0, 200, 0., 2., AliDielectronVarManager::kEta, AliDielectronVarManager::kEMCALEoverP, kTRUE);
  // EMCal nSigma electron
  histos->UserHistogram("Track", "EMCALnSigmaE_P", "n#sigma_{e} vs P_{IN};P_{IN} (GeV/c);n#sigma_{e};#tracks",
                        200, 0., 40., 200, -12, 12, AliDielectronVarManager::kPIn, AliDielectronVarManager::kEMCALnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "EMCALnSigmaE_Phi", "n#sigma_{e} vs #phi;Phi;n#sigma_{e};#tracks",
                        200, 0., TMath::TwoPi(), 200, -12, 12, AliDielectronVarManager::kPhi, AliDielectronVarManager::kEMCALnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "EMCALnSigmaE_Eta", "n#sigma_{e} vs #eta;Eta;n#sigma_{e};#tracks",
                        200, -1.0, 1.0, 200, 0., 2., AliDielectronVarManager::kEta, AliDielectronVarManager::kEMCALnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "EMCALnSigmaEle_EoverP", "n#sigma_{e}(EMCal) vs E/p;E/p;n#sigma_{e}(EMCal);#tracks",
                        200, 0., 2., 200, -12., 12., AliDielectronVarManager::kEMCALEoverP, AliDielectronVarManager::kEMCALnSigmaEle, kTRUE);
  // PID - TPC + EMCal
  histos->UserHistogram("Track", "dEdx_EoverP", "dEdx vs E/p;E/P;TPC signal (a.u.);#tracks",
                        200, 0., 2., 800, 20., 200., AliDielectronVarManager::kEMCALEoverP, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_EoverP", "n#sigma_{e}(TPC) vs E/p;E/p;n#sigma_{e}(TPC);#tracks",
                        200, 0., 2., 200, -12., 12., AliDielectronVarManager::kEMCALEoverP, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "dEdx_EMCALnSigmaE", "dEdx vs n#sigma_{e}(EMCAL);n#sigma_{e}(EMCAL);TPC signal (a.u.);#tracks",
                        200, -12., 12., 800, 20., 200., AliDielectronVarManager::kEMCALnSigmaEle, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "nSigmaTPC_EMCal", "n#sigma_{e}(TPC vs EMCAL);n#sigma_{e}(EMCAL);n#sigma_{e}(TPC);#tracks",
                        200, -5., 5., 200, -12., 12., AliDielectronVarManager::kEMCALnSigmaEle, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);

/**
 *  Histograms for Pair
**/
  histos->UserHistogram("Pair", "InvMass", "Inv.Mass;Inv. Mass (GeV/c^{2});#pairs/(40 MeV/c^{2})",
                        100, 1.0, 5.0, AliDielectronVarManager::kM);
  histos->UserHistogram("Pair", "pT", "Pt;Pt (GeV/c);#pairs",
                        2000, 0., 100.0, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair", "Eta_Phi", "#eta-#phi map of dielectron pairs;#eta_{ee};#phi_{ee};#pairs",
                        200, -1, 1, 200, 0., 10, AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
  histos->UserHistogram("Pair", "Rapidity", "Rapidity;Rapidity;#pairs",
                        200, -1., 1., AliDielectronVarManager::kY);
  histos->UserHistogram("Pair", "OpeningAngle", "Opening angle / rad;#pairs",
                        50, 0., 3.15, AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("Pair", "PseudoProperTime", "Pseudoproper decay length; pseudoproper-decay-length[cm];#pairs / 40#mum",
                        600, -0.3, 0.3, AliDielectronVarManager::kPseudoProperTime);
  histos->UserHistogram("Pair", "InvMass_Pt", "Inv. Mass vs Pt;Pt (GeV/c); Inv. Mass (GeV/c^{2})",
                        200, 0., 100., 100, 1.0, 5.0, AliDielectronVarManager::kPt, AliDielectronVarManager::kM);
  histos->UserHistogram("Pair", "OpeningAngle_Pt", "Opening angle vs p_{T} ;p_{T} (GeV/c); angle",
                        200, 0., 40., 200, 0, TMath::Pi(), AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
  //InvMass versus Proper time
  histos->UserHistogram("Pair", "InvMass_ProperTime", "InvMass vs. ProperTime;pseudoproper-decay-length[cm]; Inv. Mass [GeV]",
                        600, -0.3, 0.3, 100, 1.0, 5.0, AliDielectronVarManager::kPseudoProperTime, AliDielectronVarManager::kM);
}

Bool_t AliAnalysisTaskJpsiJet::FindDaughters(AliVTrack* trk){
  static const Double_t ERR_LIMIT = 1e-6;
  TIter nextEle(fDaughters);
  TLorentzVector* vec = NULL;
  while((vec = static_cast<TLorentzVector*>(nextEle()))){
    if(TMath::Abs(vec->Pt() - trk->Pt()) < ERR_LIMIT &&
       TMath::Abs(vec->Eta() - trk->Eta()) < ERR_LIMIT &&
       TMath::Abs(TVector2::Phi_0_2pi(vec->Phi()) - trk->Phi()) < ERR_LIMIT )
      return kTRUE;
  }
  return kFALSE;
}

Double_t AliAnalysisTaskJpsiJet::GetPseudoProperDecayTime(AliDielectronPair* pair){
  auto priv = fAOD->GetPrimaryVertex();
  Double_t errPseudoProperTime2 = 0.;
  AliKFParticle kfPair = pair->GetKFParticle();
  Double_t lxy = kfPair.GetPseudoProperDecayTime(*priv, TDatabasePDG::Instance()->GetParticle(443)->Mass(), &errPseudoProperTime2 );
  return lxy;
}

void AliAnalysisTaskJpsiJet::AddTrackFromPair(AliAODTrack* trkTemplate){
  AliAODTrack* trk = new ((*fTracksWithPair)[fTracksWithPair->GetEntriesFast()]) AliAODTrack(*trkTemplate);
  AliDielectronPair *pair = static_cast<AliDielectronPair*>(fPairs->UncheckedAt(0));

  trk->SetProdVertex(fAOD->GetPrimaryVertex());
  //trk->SetStatus(AliVTrack::kEmbedded);

  trk->SetPt(pair->Pt());
  trk->SetPhi(TVector2::Phi_0_2pi(pair->Phi()));
  trk->SetTheta(TMath::ACos(pair->Pz() / pair->P()));

  // Remove EMCal
  trk->ResetStatus(AliVTrack::kEMCALmatch);
  trk->SetEMCALcluster(AliVTrack::kEMCALNoMatch);

  // Reset reference
  trk->ResetBit(kIsReferenced);
  trk->SetUniqueID(0);

  // DEBUG - Pseudo-proper decay length
  trk->SetTrackPhiEtaPtOnEMCal(pair->M(), GetPseudoProperDecayTime(pair), 0.);
}

Bool_t AliAnalysisTaskJpsiJet::FindTrackInJet(AliEmcalJet* jet, AliVParticle* p, Int_t trackID){
  // Find track in jet with track ID
  // AliEmcalJetTask::fgkConstIndexShift = 100000
  // track_long_ID_for_jet_finder = track_index_in_array + fgkConstIndexShift * index_of_track_container
  // - More details in AliEmcalJetTask.cxx
  Int_t fgkConstIndexShift = 100000;

  // Get index of track container
  Int_t iCont = jet->TrackAt(0) / fgkConstIndexShift;
  // Get track index for jet finder
  Int_t trackLongID = trackID + iCont * fgkConstIndexShift;
  // Get track index in jet
  Int_t trackJetID = jet->ContainsTrack(trackLongID);
  if(trackJetID == -1) return kFALSE;
  // Cross-check
  const Double_t TRACK_TOLERANCE = 1e-3;
  AliVParticle* pJ = jet->Track(trackJetID);
  Bool_t passCrossCheck = kTRUE;
  if(TMath::Abs(pJ->Pt() - p->Pt()) > TRACK_TOLERANCE)
    passCrossCheck = kFALSE;
  if(TMath::Abs(pJ->Phi() - p->Phi()) > TRACK_TOLERANCE)
    passCrossCheck = kFALSE;
  if(TMath::Abs(pJ->Theta() - p->Theta()) > TRACK_TOLERANCE)
    passCrossCheck = kFALSE;
  if(!passCrossCheck){
    AliWarning(Form("Difference between tracks in more than %.1e", TRACK_TOLERANCE));
    pJ->Print();
    AliWarning(Form("Track pT:%.3f/%.3f, phi:%.3f/%.3f, theta:%.3f/%.3f", pJ->Pt(), p->Pt(), pJ->Phi(), p->Phi(), pJ->Theta(), p->Theta()));
  }
  return passCrossCheck;
}

void AliAnalysisTaskJpsiJet::InitHistogramsForTaggedJet(const char *histClass){
  // THnSparse - pT_pair, M, Lxy, z, \DeltaR, pT_jet
  TString histName = Form("%s/PairVars",histClass);
  Int_t nBins[6]   = { 200, 100,  600,  11,  10, 200};
  Double_t xmin[6] = {  0.,  1., -0.3,  0.,  0.,  0.};
  Double_t xmax[6] = {100.,  5.,  0.3, 1.1,  1.,100.};
  THnSparse *hs = fHistos->CreateTHnSparse(histName.Data(), "Dielectron pair in jet variables (p_{T}^{pair}-M_{e^{+}e^{-}}-L_{xy}-z-#DeltaR-p_{T}^{jet});p^{pair}_{T} (GeV/c);M_{e^{+}e^{-}} (GeV/c^{2});L_{xy} (cm);z(p_{T}^{pair}/p_{T}^{jet});#DeltaR;p_{T}^{jet} (GeV/c)", 6, nBins, xmin, xmax);

  // Constituents
  histName = Form("%s/Ntracks_pT",histClass);
  fHistos->CreateTH2(histName.Data(), "Jet constituents - number vs p_{T}^{jet};p_{T}^{jet} (GeV/c);N_{tracks};",
      200, 0., 100., 100, -0.5, 99.5);
  histName = Form("%s/Area_pT",histClass);
  fHistos->CreateTH2(histName.Data(), "Jet clustering area;p_{T}^{jet} (GeV/c);A_{jet};",
      200, 0., 100., 100, 0., 2.);
  
  // Fragmentation Function - Prompt and Non-prompt
  // Pre-defined cuts: pT_jet > 15, pT_pair > 5, Mpair\in[2.92, 3.16]
  // Prompt: |Lxy| < 0.01; Non-prompt: |Lxy| > 0.01
  histName = Form("%s/PromptFF",histClass);
  fHistos->CreateTH1(histName.Data(), "z #equiv p_{T} of track in jet - Prompt J/#psi;z;N_{pairs}",11,0.,1.1);
  histName = Form("%s/NonPromptFF",histClass);
  fHistos->CreateTH1(histName.Data(), "z #equiv p_{T} of track in jet - Non-Prompt J/#psi;z;N_{pairs}",11,0.,1.1);
}

Bool_t AliAnalysisTaskJpsiJet::FillHistogramsForTaggedJet(const char* histClass){
  // Tagged jet container
  AliJetContainer* jets = NULL;
  TIter next(fJets);
  while((jets = static_cast<AliJetContainer*>(next()))){
    TString jetName = jets->GetName();
    if(jetName.BeginsWith("JpsiJet_")) break;
  }
  if(!jets->GetNJets()) return kFALSE;

  // Get dielectron pair and track
  AliDielectronPair *pair = (AliDielectronPair*)(fPairs->At(0));
  Int_t pairTrackID = fTracksWithPair->GetEntriesFast() - 1;
  AliAODTrack *pairTrack = (AliAODTrack*)(fTracksWithPair->UncheckedAt(pairTrackID));
  
  // Find pair in jet
  AliEmcalJet *taggedJet = NULL;
  for(auto jet : jets->all()){
    if(FindTrackInJet(jet, pairTrack, pairTrackID)){
      taggedJet = jet;
      break;
    }
  }
  UInt_t rejectionReason = 0;
  if(!taggedJet || !jets->AcceptJet(taggedJet, rejectionReason)) return kFALSE;
  AliDebug(1, Form("Found pair (%.2f, %.2f, %.2f) in jet (%s)", pair->Pt(), pair->Eta(), TVector2::Phi_0_2pi(pair->Phi()), (taggedJet->toString()).Data()));

  // DEBUG: Register the tagged jet on detector level
  fTaggedJet = taggedJet;

  // THnSparse - pT_pair, M, Lxy, z, \DeltaR, pT_jet
  TString histName = Form("%s/PairVars",histClass);
  Double_t x[6] = {0.};
  x[0] = pair->Pt();
  x[1] = pair->M();
  x[2] = GetPseudoProperDecayTime(pair);
  x[3] = (taggedJet->GetNumberOfTracks() == 1 ? 1.0 : pair->Pt() / taggedJet->Pt());
  x[4] = 0.0;
  x[5] = taggedJet->Pt();

  fHistos->FillTHnSparse(histName.Data(), x, 1.0);

  // Constituents Info
  fHistos->FillTH2(Form("%s/Ntracks_pT",histClass), taggedJet->Pt(), taggedJet->GetNumberOfTracks());
  fHistos->FillTH2(Form("%s/Area_pT",histClass), taggedJet->Pt(), taggedJet->Area());

  // FF of J/psi candidate
  if(taggedJet->Pt() > 15. && 
      pair->M() > 2.92 && pair->M() < 3.16){
    if(TMath::Abs(x[2]) < 0.01){
      histName = Form("%s/PromptFF",histClass);
      fHistos->FillTH1(histName.Data(), x[3]);
    }else if(x[2] > 0.01){
      histName = Form("%s/NonPromptFF",histClass);
      fHistos->FillTH1(histName.Data(), x[3]);
    }
  }

  return kTRUE;
}

/**
 *  MC truth
**/
void AliAnalysisTaskJpsiJet::InitHistogramsForMC(){
  if(!fHistosMC || !fIsMC){
    AliWarning("MC input not initialized");
    return;
  }

  // Event
  TString histGroup = "Event";
  fHistosMC->CreateTH1(Form("%s/NParticles", histGroup.Data()), "MC particles - ALL produced", 5001, -0.5, 5000.5);
  fHistosMC->CreateTH1(Form("%s/NPhysPrim", histGroup.Data()), "MC particles - Physical Primary with |#eta|<1.0", 501, -0.5, 500.5);
  fHistosMC->CreateTH2(
      Form("%s/NPhysPrim_Ntracks", histGroup.Data()),
      "MC particles vs Detector tracks;N_{trakcs};N_{MC particle};",
      501, -0.5, 500.5,
      501, -0.5, 500.5);

  // Electron
  histGroup = "Electron";
    // eleVars - pT, eta， phi, E
  Int_t nBins[4]   = {500, 400,  80,  500};
  Double_t xmin[4] = { 0., -2., -1.,   0.};
  Double_t xmax[4] = {100., 2.,  7., 100.};
  fHistosMC->CreateTHnSparse(Form("%s/eleVars", histGroup.Data()), "Electron kinetic variables (p_{T}-#eta-#phi-E);p_{T} (GeV/c);#eta;#phi;E (GeV)", 4, nBins, xmin, xmax);
    // PID check - Pure, Wrong, Miss
  fHistosMC->CreateTH1(
      Form("%s/PID_pure", histGroup.Data()),
      "Electron PID - Pure;p_{T,ele} (GeV/c);N_{pure};",
      500, 0., 100.);
  fHistosMC->CreateTH1(
      Form("%s/PID_wrong", histGroup.Data()),
      "Electron PID - Wrong;p_{T,ele} (GeV/c);N_{pure};",
      500, 0., 100.);
  fHistosMC->CreateTH1(
      Form("%s/EMCal_all", histGroup.Data()),
      "Electron on EMCal/DCal (MC);p_{T,ele} (GeV/c);N_{ele,EMC};",
      500, 0., 100.);
  fHistosMC->CreateTH1(
      Form("%s/EMCal_det", histGroup.Data()),
      "Electron on EMCal/DCal (Detector);p_{T,ele} (GeV/c);N_{ele,EMC};",
      500, 0., 100.);
  // Jpsi
  InitHistogramsForJpsiMC("JpsiPrompt");
  InitHistogramsForJpsiMC("JpsiBdecay");

  // Jet
  InitHistogramsForJetMC();
}

void AliAnalysisTaskJpsiJet::InitHistogramsForJpsiMC(const char* histClass){
  // jpsiVars - pT, rapidity, phi, E
  Int_t nBins[4]   = {500, 400,  80,  500};
  Double_t xmin[4] = { 0., -2., -1.,   0.};
  Double_t xmax[4] = {100., 2.,  7., 100.};
  fHistosMC->CreateTHnSparse(Form("%s/jpsiVars", histClass), "J/#psi kinetic variables (p_{T}-Y-#phi-E);p_{T} (GeV/c);Rapidity;#phi;E (GeV)", 4, nBins, xmin, xmax);
  
  // Signal check - pairVars : pT, M, Lxy
  Int_t nBinsDet[3]   = {500, 100,  600};
  Double_t xminDet[3] = { 0.,  1., -0.3};
  Double_t xmaxDet[3] = {100., 5.,  0.3};
  fHistosMC->CreateTHnSparse(
      Form("%s/Reco_sig", histClass),
      "J/#psi Reconstruction - Signal;p_{T,J/#psi} (GeV/c);M_{e^{+}e^{-}} (GeV/c^{2});L_{xy} (cm);N_{pair};",
      3, nBinsDet, xminDet, xmaxDet);
  fHistosMC->CreateTHnSparse(
      Form("%s/Reco_bkg", histClass),
      "J/#psi Reconstruction - Background;p_{T,pair} (GeV/c);M_{e^{+}e^{-}} (GeV/c^2);L_{xy} (cm);N_{pair};",
      3, nBinsDet, xminDet, xmaxDet);
  // Tagged jet
  fHistosMC->CreateTH2(Form("%s/Jet_PtNtracks",histClass),
  "Jet constituents - number vs p_{T}^{jet};p_{T,jet}^{gen} (GeV/c);N_{tracks};",
      200, 0., 100., 100, -0.5, 99.5);
  fHistosMC->CreateTH2(Form("%s/Jet_PtZ",histClass),
  "p_{T}^{jet} vs fragmentation function;z;p_{T,jet}^{gen} (GeV/c);",
      11, 0., 1.1, 200, 0., 100.);
  fHistosMC->CreateTH2(Form("%s/Jet_PtArea",histClass),
  "p_{T}^{jet} vs clustering area;p_{T,jet}^{gen} (GeV/c);A_{jet};",
      200, 0., 100., 100, 0., 2.);
  // Detector response with tagged jet
    // z_det, z_gen, pT-det, pT-gen, dz, dPtJet, dPtJpsi
  Int_t nBinsDetZ[7]   = {11,   11,  100,  100, 200, 200, 200};
  Double_t xminDetZ[7] = { 0.,  0.,   0.,   0., -1., -1., -1.};
  Double_t xmaxDetZ[7] = {1.1, 1.1, 100., 100.,  1.,   1.,  1.};
  fHistosMC->CreateTHnSparse(
      Form("%s/Jet_DetResponse", histClass),
      "Detector response matrix - Fragmentation Function with J/#psi tagged jet p_{T};z_{det};z_{gen};p_{T,jet}^{det} (GeV/c);p_{T,jet}^{gen} (GeV/c);#delta z;#delta p_{T,jet};#delta p_{T,J/#psi};",
      7, nBinsDetZ, xminDetZ, xmaxDetZ);
}

Bool_t AliAnalysisTaskJpsiJet::ApplyEmcalCut(AliVParticle* par, Bool_t isMCTruth = kTRUE){
  // Kinetic variables
  Double_t eta = TMath::Abs(par->Eta());
  Double_t phi = par->Phi();
  Double_t E = par->E();
  if(!isMCTruth){
    Int_t clsID = ((AliAODTrack*)par)->GetEMCALcluster();
    if(clsID > -1)
      E = fAOD->GetCaloCluster(clsID)->E();
  }
  // EMCal
  Bool_t isEMCal = kTRUE;
  if(eta > 0.7) isEMCal = kFALSE;
  if(phi < 1.396 || phi > 3.264) isEMCal = kFALSE;
  // DCal
  Bool_t isDCal = kTRUE;
  if(eta < 0.22 || eta > 0.7) isDCal = kFALSE;
  if(phi < 4.377 || phi > 5.707) isDCal = kFALSE;
  // Energy threshold
  Bool_t isTriggered = (E > fEMCalEth);
  
  return (isTriggered && (isEMCal || isDCal));
}

Bool_t AliAnalysisTaskJpsiJet::RunParticleLevelAnalysis(){
  // MC info.
  fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::StdBranchName()));
  fMCParticles = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));

  if(!fMCHeader || !fMCParticles){
    AliFatal("Fail to retrieve MC objects");
    return kFALSE;
  }

  SetJpsiGeneratorType();

  fHistosMC->FillTH1("Event/NParticles", fMCParticles->GetEntriesFast());

  // Loop on MC particles
  Int_t nPhysPrim = 0;
  AliAODMCParticle* mcp = NULL;
  TIter nextParticle(fMCParticles);
  while((mcp = static_cast<AliAODMCParticle*>(nextParticle()))){
    Int_t pdg = TMath::Abs(mcp->GetPdgCode());
    if(mcp->IsPhysicalPrimary() && TMath::Abs(mcp->Eta()) < 1.0){
      nPhysPrim++;
      if(pdg == PDG_ELECTRON){
        FillHistogramsForParticle("Electron/eleVars", mcp);
        if(ApplyEmcalCut(mcp))
          fHistosMC->FillTH1("Electron/EMCal_all", mcp->Pt());
      }// Electron - MC
    }
    if(pdg == PDG_JPSI){
        Double_t x[4] = {0.};
        x[0] = mcp->Pt();
        x[1] = mcp->Y();
        x[2] = mcp->Phi();
        x[3] = mcp->E();
        // Re-check J/psi generator type
        if(mcp->GetMother() > -1) fMCGenType = "JpsiBdecay";
        fHistosMC->FillTHnSparse(Form("%s/jpsiVars", fMCGenType.Data()), x, 1.0);
    }// Jpsi
  }// End - MC particles

  fHistosMC->FillTH1("Event/NPhysPrim", nPhysPrim);
  fHistosMC->FillTH2("Event/NPhysPrim_Ntracks", fAOD->GetNumberOfTracks(), nPhysPrim);

  // Jet
  FillHistogramsForJetMC("Jet");
  if(RunJetFinder("JetMC"))
    FillHistogramsForJetMC("JetMC");

  // Particle vs Detector
    // Electron PID
  FillHistogramsForElectronPID(fDielectron->GetTrackArray(0));
  FillHistogramsForElectronPID(fDielectron->GetTrackArray(1));

    // J/psi Acceptance X Efficiency
  if(fDielectron->HasCandidates()){
    FillHistogramsForJetMC("JpsiJet");
    FillHistogramsForJpsiMC();
  }
  return kTRUE;
}

// Set J/psi generator type from MC header
void AliAnalysisTaskJpsiJet::SetJpsiGeneratorType(){
  // MC from HFjets - JIRA/ALIROOT-7974 
  if(fMCHeader->GetCocktailHeaders()->GetSize() == 1){
    fMCGenType = "JpsiPrompt";
    return;
  }
  // Generator type - Prompt/Jpsi2ee, Bdecay/B2Jpsi2ee
  TString genType = fMCHeader->GetCocktailHeader(1)->GetName();
  if(genType.Contains("Jpsi2ee")){
    fMCGenType = "JpsiPrompt";
  }
  else if(genType.Contains("B2Jpsi2ee")){
    fMCGenType = "JpsiBdecay";
  }else
  {
    AliWarning("Unknown generator type for J/psi MC");
    return;
  }
}

void AliAnalysisTaskJpsiJet::FillHistogramsForParticle(const char* histName, AliVParticle* par){
  Double_t x[4] = {0.};
  x[0] = par->Pt();
  x[1] = par->Eta();
  x[2] = par->Phi();
  x[3] = par->E();

  fHistosMC->FillTHnSparse(histName, x, 1.0);
}

void AliAnalysisTaskJpsiJet::FillHistogramsForElectronPID(const TObjArray* eleArray){

  AliAODTrack* eleTrk = NULL;
  TIter nextEle(eleArray);
  while((eleTrk = static_cast<AliAODTrack*>(nextEle()))){
    if(ApplyEmcalCut(eleTrk, kFALSE))
      fHistosMC->FillTH1("Electron/EMCal_det", eleTrk->Pt());
    Int_t mcID = TMath::Abs(eleTrk->GetLabel());
    auto mcp = static_cast<AliAODMCParticle*>(fMCParticles->At(mcID));
    if(!mcp){
      AliWarning("Can not found MC particle");
      eleTrk->Print();
      continue;
    }
    Int_t pdg = TMath::Abs(mcp->GetPdgCode());
    if(pdg == PDG_ELECTRON)
      fHistosMC->FillTH1("Electron/PID_pure", mcp->Pt());
    else
      fHistosMC->FillTH1("Electron/PID_wrong", mcp->Pt());
  }
}

// Add J/psi particle for jet finder
// Charge() is retrieved by PDG code, and it can not be modified from AliAODMCParticle. To avoid RejectionReason::kChargeCut, build new AliAODMCParticle as electron from AliMCParticle.
AliAODMCParticle* AliAnalysisTaskJpsiJet::AddParticleFromJpsi(AliAODMCParticle* jpsi){
  
  auto part = new TParticle(
    PDG_ELECTRON, 0, -1, -1, -1, -1,
    jpsi->Px(), jpsi->Py(), jpsi->Pz(), jpsi->E(),
    jpsi->Xv(), jpsi->Yv(), jpsi->Zv(), jpsi->T());

  auto mcPart = new AliMCParticle(part);
  
  Int_t jpsiAsEleID = fMCParticles->GetEntries();
  auto jpsiAsEle = new ((*fMCParticles)[jpsiAsEleID]) AliAODMCParticle(mcPart, jpsiAsEleID, jpsi->GetFlag());

  jpsiAsEle->SetPhysicalPrimary(kTRUE);

  delete part;
  delete mcPart;
  return jpsiAsEle;
}

// J/psi reconstruction vs MC truth
void AliAnalysisTaskJpsiJet::FillHistogramsForJpsiMC(){

  AliDielectronPair* pair = NULL;
  TIter nextPair(fDielectron->GetPairArray(1));
  while( (pair = static_cast<AliDielectronPair*>(nextPair()) ) ){

    AliVParticle* d1 = pair->GetFirstDaughterP();
    AliVParticle* d2 = pair->GetSecondDaughterP();

    if(!d1 || !d2){
      AliWarning("Can not found daughters of AliDielectronPair");
      continue;
    }

    // For real signal, both daughters should come from J/psi decay.
    Double_t x[3] = {0.};
    x[0] = pair->Pt();
    x[1] = pair->M();
    x[2] = GetPseudoProperDecayTime(pair);
    Int_t mother1 = CheckDielectronDaughter(d1);
    Int_t mother2 = CheckDielectronDaughter(d2);
    if( mother1 == mother2 && mother1 > -1)
      fHistosMC->FillTHnSparse(Form("%s/Reco_sig", fMCGenType.Data()), x, 1.0);
    else{
      fHistosMC->FillTHnSparse(Form("%s/Reco_bkg", fMCGenType.Data()), x, 1.0);
      continue;
    }

    // J/psi in Jet - Particle Level
    auto jpsi = static_cast<AliAODMCParticle *>(fMCParticles->At(mother1));
    auto mcD1 = static_cast<AliAODMCParticle *>(fMCParticles->At(TMath::Abs(d1->GetLabel())));
    auto mcD2 = static_cast<AliAODMCParticle *>(fMCParticles->At(TMath::Abs(d2->GetLabel())));
      // Remove daughters and add J/psi for jet finder
      // See: AliMCParticleContainer::AcceptMCParticle
    auto* jpsiAsEle = AddParticleFromJpsi(jpsi);
    Bool_t mcD1_status = mcD1->IsPhysicalPrimary();
    Bool_t mcD2_status = mcD2->IsPhysicalPrimary();
    mcD1->SetPhysicalPrimary(kFALSE);
    mcD2->SetPhysicalPrimary(kFALSE);

    if(RunJetFinder("JpsiJetMC")){
      FillHistogramsForJetMC("JpsiJetMC");
      // DEBUG
      fJpsiPair = pair;
      fJpsiMC = jpsi;
      FillHistogramsForTaggedJetMC(jpsiAsEle);
    }

    fMCParticles->Remove(jpsiAsEle);
    mcD1->SetPhysicalPrimary(mcD1_status);
    mcD2->SetPhysicalPrimary(mcD2_status);
  }// End - Loop dielectron pairs
}

// Found mother of dielectron daughter
// Return MC label, if mother is J/psi, else return -1
Int_t AliAnalysisTaskJpsiJet::CheckDielectronDaughter(AliVParticle *par)
{
  Int_t mcID = TMath::Abs(par->GetLabel());
  auto mcp = static_cast<AliAODMCParticle *>(fMCParticles->At(mcID));
  Int_t motherID = mcp->GetMother();
  if (motherID == -1) return -1;

  auto mcMother = static_cast<AliAODMCParticle *>(fMCParticles->At(motherID));
  Int_t pdg = TMath::Abs(mcMother->GetPdgCode());
  if (pdg != PDG_JPSI) return -1;
  
  return motherID;
}

void AliAnalysisTaskJpsiJet::InitHistogramsForJetMC(){
  AliJetContainer* jets = NULL;
  TIter next(fJets);
  while((jets = static_cast<AliJetContainer*>(next()))){
    // THnSparse - pT, eta， phi, MCPt
    Int_t nBins[4]   = {2000, 200, 100, 2000};
    Double_t xmin[4] = {0.,   -1., -2.,   0.};
    Double_t xmax[4] = {100.,  1.,  8., 100.};
    THnSparse* hs = fHistosMC->CreateTHnSparse(Form("%s/jetVars", jets->GetName()), "Jet kinetic variables (p_{T}-#eta-#phi-p_{T,MC});p_{T,reco} (GeV/c);#eta;#phi;p_{T,true} (GeV/c)", 4, nBins, xmin, xmax);

    // Constituents
    fHistosMC->CreateTH2(Form("%s/Ntracks_pT",jets->GetName()), "Jet constituents - number vs p_{T}^{jet};p_{T}^{jet} (GeV/c);N_{tracks};",
      200, 0., 100., 100, -0.5, 99.5);
    fHistosMC->CreateTH2(Form("%s/Area_pT",jets->GetName()), "Jet clustering area;p_{T}^{jet} (GeV/c);A_{jet};",
      200, 0., 100., 100, 0., 2.);
    // Detector response matrix
    fHistosMC->CreateTH2(
      Form("%s/detResponse", jets->GetName()),
      "Detector response matrix for jets;p_{T,reco} (GeV/c);p_{T,true} *(GeV/c)",
      100, 0., 100.,
      100, 0., 100.);
  }
}

// Calculate jet MC pT by constituents
Double_t AliAnalysisTaskJpsiJet::GetJetMCPt(AliEmcalJet* jet){
  Double_t mcPt = 0.;
  for(int i = 0; i < jet->GetNumberOfTracks(); i++){
    Int_t mcID = jet->Track(i)->GetLabel();
    if(mcID < 0) continue;
    auto mcp = (AliVParticle*)(fMCParticles->UncheckedAt(mcID));
    mcPt += mcp->Pt();
  }
  jet->SetMCPt(mcPt);
  return mcPt;
}

void AliAnalysisTaskJpsiJet::FillHistogramsForJetMC(const char* jetTag){
  
  AliJetContainer* jets = NULL;
  TIter next(fJets);
  while((jets = static_cast<AliJetContainer*>(next()))){
    TString jetName = jets->GetName();
    if (!jetName.BeginsWith(Form("%s_", jetTag)))
      continue;

    // Found jet with tag
    jets->NextEvent(fAOD);
    jets->SetArray(fAOD);

    AliDebug(1, Form("%d jets found in %s", jets->GetNJets(), jets->GetName()));
    for(auto jet : jets->all()){
      // Jet cuts
      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason)) {
        AliDebug(2,Form("Jet was rejected for reason : %d, details : %s", rejectionReason, (jet->toString()).Data()));
        continue;
      }
      Double_t x[4] = {0.};
      x[0] = jet->Pt();
      x[1] = jet->Eta();
      x[2] = jet->Phi();
      x[3] = GetJetMCPt(jet);
      fHistosMC->FillTHnSparse(Form("%s/jetVars", jetName.Data()),x,1.0);
      fHistosMC->FillTH2(Form("%s/Ntracks_pT", jetName.Data()), jet->Pt(), jet->GetNumberOfTracks());
      fHistosMC->FillTH2(Form("%s/Area_pT", jetName.Data()), jet->Pt(), jet->Area());
      fHistosMC->FillTH2(Form("%s/detResponse", jetName.Data()), jet->Pt(), jet->MCPt());
    }
  }
}

Bool_t AliAnalysisTaskJpsiJet::FillHistogramsForTaggedJetMC(AliAODMCParticle* jpsi){
  // Retrive jet container on particle level
  AliJetContainer* jets = NULL;
  TIter next(fJets);
  while((jets = static_cast<AliJetContainer*>(next()))){
    TString jetName = jets->GetName();
    if (jetName.BeginsWith("JpsiJetMC_"))
      break;
  }
  // Find tagged jet
  AliEmcalJet* taggedJet = NULL;
  for(auto jet : jets->all()){
    if(FindTrackInJet(jet, jpsi, fMCParticles->IndexOf(jpsi))){
      taggedJet = jet;
      break;
    }
  }
  UInt_t rejectionReason = 0;
  if(!taggedJet) return kFALSE;
  if(!jets->AcceptJet(taggedJet, rejectionReason)){
    AliDebug(1, Form("J/psi tagged jet was reject by %d", rejectionReason));
    return kFALSE;
  }
  fTaggedJetMC = taggedJet;

  Int_t nTracks = taggedJet->GetNumberOfTracks();
  Double_t z = jpsi->Pt() / taggedJet->Pt();
  if(nTracks == 1) z = 1.0;
  fHistosMC->FillTH2(Form("%s/Jet_PtNtracks", fMCGenType.Data()), taggedJet->Pt(), nTracks);
  fHistosMC->FillTH2(Form("%s/Jet_PtZ", fMCGenType.Data()), z, taggedJet->Pt());
  fHistosMC->FillTH2(Form("%s/Jet_PtArea", fMCGenType.Data()), taggedJet->Pt(), taggedJet->Area());

  // Detector response
  if(!fTaggedJet) return kTRUE;
  Double_t x[7] = {0.};
  // z-det
  x[0] = fJpsiPair->Pt() / fTaggedJet->Pt();
  if(fTaggedJet->GetNumberOfTracks() == 1) x[0] = 1.0;
  x[1] = z;  // z-gen
  x[2] = fTaggedJet->Pt(); // jet pT-det
  x[3] = fTaggedJetMC->Pt(); // jet pT-gen
  x[4] = (x[0] - x[1]) / x[1]; // dz
  x[5] = (x[2] - x[3]) / x[3]; // dPtJet
  x[6] = (fJpsiPair->Pt() - fJpsiMC->Pt()) / fJpsiMC->Pt();
  fHistosMC->FillTHnSparse(Form("%s/Jet_DetResponse", fMCGenType.Data()), x, 1.0);

  return kTRUE;
}
