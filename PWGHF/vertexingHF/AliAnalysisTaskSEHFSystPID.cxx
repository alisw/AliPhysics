/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. */

/////////////////////////////////////////////////////////////////////////////////////////
// \class AliAnalysisTaskSEHFSystPID                                                   //
// \brief analysis task for the study of PID systematic uncertainties of HF particles  //
// \author: A. M. Barbano, anastasia.maria.barbano@cern.ch                             //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                           //
/////////////////////////////////////////////////////////////////////////////////////////

#include <cstddef>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TChain.h>
#include <TFile.h>

#include "AliAnalysisTaskSEHFSystPID.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliESDtrack.h"
#include "AliRDHFCuts.h"

ClassImp(AliAnalysisTaskSEHFSystPID)

//________________________________________________________________________
AliAnalysisTaskSEHFSystPID::AliAnalysisTaskSEHFSystPID() :
AliAnalysisTaskSE("TaskNsigmaPID"),
fOutputList(nullptr),
fHistNEvents(nullptr),
fHistQtVsMassKinks(nullptr),
fHistPDaughterVsMotherKink(nullptr),
fHistdEdxVsPMotherKink(nullptr),
fHistOpeningAngleVsPMotherKink(nullptr),
fHistNTPCclsVsRadius(nullptr),
fPIDtree(nullptr),
fPTPC(0),
fPTOF(0),
fPHMPID(0),
fdEdxTPC(0),
fdEdxITS(0),
fToF(0),
fPt(0),
fTPCNcls(0),
fTPCNclsPID(0),
fTrackLength(0),
fStartTimeRes(0),
fTPCNcrossed(0),
fTPCFindable(0),
fITSclsMap(0),
fHMPIDsignal(0),
fHMPIDoccupancy(0),
fTrackInfoMap(0),
fOOBPileupMap(0),
fEta(-9999),
fPhi(9999),
fPDGcode(-1),
fTag(0),
fNsigmaMaxForKaonTag(0.2),
fNsigmaMaxForNucleiTag(3.),
fQtMinKinks(0.15),
fRMinKinks(120),
fRMaxKinks(210),
fDeadZoneWidth(3.),
fCutGeoNcrNclLength(130.),
fCutGeoNcrNclGeom1Pt(1.5),
fCutGeoNcrNclFractionNcr(0.85),
fCutGeoNcrNclFractionNcl(0.7),
fCentMin(0.),
fCentMax(100.),
fCentEstimator(kCentOff),
fTriggerClass(""),
fTriggerMask(AliVEvent::kINT7),
fIsMC(false),
fSystem(0),
fConversionFactordEdx(100.),
fESDtrackCuts(nullptr),
fAOD(nullptr),
fPIDresp(nullptr),
fV0cuts(nullptr),
fFillTreeWithPIDInfo(true),
fFillTreeWithNsigmaPIDOnly(false),
fFillTreeWithRawPIDOnly(false),
fFillTreeWithTrackQualityInfo(false),
fEnabledDownSampling(false),
fFracToKeepDownSampling(0.1),
fPtMaxDownSampling(1.5),
fDownSamplingOpt(0),
fAODProtection(0),
fRunNumberPrevEvent(-1),
fEnableNsigmaTPCDataCorr(false),
fSystNsigmaTPCDataCorr(AliAODPidHF::kNone),
fMeanNsigmaTPCPionData{},
fMeanNsigmaTPCKaonData{},
fMeanNsigmaTPCProtonData{},
fSigmaNsigmaTPCPionData{},
fSigmaNsigmaTPCKaonData{},
fSigmaNsigmaTPCProtonData{},
fPlimitsNsigmaTPCDataCorr{},
fNPbinsNsigmaTPCDataCorr(0),
fEtalimitsNsigmaTPCDataCorr{},
fNEtabinsNsigmaTPCDataCorr(0),
fUseAliEventCuts(false),
fAliEventCuts(),
fApplyPbPbOutOfBunchPileupCuts(0),
fApplyPbPbOutOfBunchPileupCutsITSTPC(0),
fKeepOnlyPbPbOutOfBunchPileupCutsITSTPC(false),
fUseTimeRangeCutForPbPb2018(true),
fTimeRangeCut()
{
  //
  // default constructur
  //

  for(int iHisto=0; iHisto<5; iHisto++) fHistArmenteroPlot[iHisto] = nullptr;

  for(int iDet=0; iDet<kNMaxDet; iDet++) {
    for(int iHypo=0; iHypo<kNMaxHypo; iHypo++) {
      fPIDNsigma[iDet][iHypo]      = numeric_limits<short>::min();
      fHistNsigmaVsPt[iDet][iHypo] = nullptr;
    }
  }

  fEnabledSpecies[kPion]     = true;
  fEnabledSpecies[kKaon]     = true;
  fEnabledSpecies[kProton]   = true;
  fEnabledSpecies[kElectron] = false;
  fEnabledSpecies[kDeuteron] = false;
  fEnabledSpecies[kTriton]   = false;
  fEnabledSpecies[kHe3]      = false;

  fEnabledDet[kITS]          = false;
  fEnabledDet[kTPC]          = true;
  fEnabledDet[kTOF]          = true;
  fEnabledDet[kHMPID]        = false;

  for(int iP=0; iP<=AliAODPidHF::kMaxPBins; iP++) {
    fPlimitsNsigmaTPCDataCorr[iP] = 0.;
  }
  for(int iEta=0; iEta<=AliAODPidHF::kMaxEtaBins; iEta++) {
    fEtalimitsNsigmaTPCDataCorr[iEta] = 0.;
  }

  fAliEventCuts.SetManualMode();
}

//________________________________________________________________________
AliAnalysisTaskSEHFSystPID::AliAnalysisTaskSEHFSystPID(const char *name, int system) :
AliAnalysisTaskSE(name),
fOutputList(nullptr),
fHistNEvents(nullptr),
fHistQtVsMassKinks(nullptr),
fHistPDaughterVsMotherKink(nullptr),
fHistdEdxVsPMotherKink(nullptr),
fHistOpeningAngleVsPMotherKink(nullptr),
fHistNTPCclsVsRadius(nullptr),
fPIDtree(nullptr),
fPTPC(0),
fPTOF(0),
fPHMPID(0),
fdEdxTPC(0),
fdEdxITS(0),
fToF(0),
fPt(0),
fTPCNcls(0),
fTPCNclsPID(0),
fTrackLength(0),
fStartTimeRes(0),
fTPCNcrossed(0),
fTPCFindable(0),
fITSclsMap(0),
fHMPIDsignal(0),
fHMPIDoccupancy(0),
fTrackInfoMap(0),
fOOBPileupMap(0),
fEta(-9999),
fPhi(9999),
fPDGcode(-1),
fTag(0),
fNsigmaMaxForKaonTag(0.2),
fNsigmaMaxForNucleiTag(3.),
fQtMinKinks(0.15),
fRMinKinks(120),
fRMaxKinks(210),
fDeadZoneWidth(3.),
fCutGeoNcrNclLength(130.),
fCutGeoNcrNclGeom1Pt(1.5),
fCutGeoNcrNclFractionNcr(0.85),
fCutGeoNcrNclFractionNcl(0.7),
fCentMin(0.),
fCentMax(100.),
fCentEstimator(kCentOff),
fTriggerClass(""),
fTriggerMask(AliVEvent::kINT7),
fIsMC(false),
fSystem(system),
fConversionFactordEdx(100.),
fESDtrackCuts(nullptr),
fAOD(nullptr),
fPIDresp(nullptr),
fV0cuts(nullptr),
fFillTreeWithPIDInfo(true),
fFillTreeWithNsigmaPIDOnly(false),
fFillTreeWithRawPIDOnly(false),
fFillTreeWithTrackQualityInfo(false),
fEnabledDownSampling(false),
fFracToKeepDownSampling(0.1),
fPtMaxDownSampling(1.5),
fDownSamplingOpt(0),
fAODProtection(0),
fRunNumberPrevEvent(-1),
fEnableNsigmaTPCDataCorr(false),
fSystNsigmaTPCDataCorr(AliAODPidHF::kNone),
fMeanNsigmaTPCPionData{},
fMeanNsigmaTPCKaonData{},
fMeanNsigmaTPCProtonData{},
fSigmaNsigmaTPCPionData{},
fSigmaNsigmaTPCKaonData{},
fSigmaNsigmaTPCProtonData{},
fPlimitsNsigmaTPCDataCorr{},
fNPbinsNsigmaTPCDataCorr(0),
fEtalimitsNsigmaTPCDataCorr{},
fNEtabinsNsigmaTPCDataCorr(0),
fUseAliEventCuts(false),
fAliEventCuts(),
fApplyPbPbOutOfBunchPileupCuts(),
fUseTimeRangeCutForPbPb2018(true),
fTimeRangeCut()
{
  //
  // standard constructur
  //

  for(int iHisto=0; iHisto<5; iHisto++) fHistArmenteroPlot[iHisto] = nullptr;

  for(int iDet=0; iDet<kNMaxDet; iDet++) {
    for(int iHypo=0; iHypo<kNMaxHypo; iHypo++) {
      fPIDNsigma[iDet][iHypo]      = numeric_limits<short>::min();
      fHistNsigmaVsPt[iDet][iHypo] = nullptr;
    }
  }

  fEnabledSpecies[kPion]     = true;
  fEnabledSpecies[kKaon]     = true;
  fEnabledSpecies[kProton]   = true;
  fEnabledSpecies[kElectron] = false;
  fEnabledSpecies[kDeuteron] = false;
  fEnabledSpecies[kTriton]   = false;
  fEnabledSpecies[kHe3]      = false;

  fEnabledDet[kITS]          = false;
  fEnabledDet[kTPC]          = true;
  fEnabledDet[kTOF]          = true;
  fEnabledDet[kHMPID]        = false;

  for(int iP=0; iP<=AliAODPidHF::kMaxPBins; iP++) {
    fPlimitsNsigmaTPCDataCorr[iP] = 0.;
  }
  for(int iEta=0; iEta<=AliAODPidHF::kMaxEtaBins; iEta++) {
    fEtalimitsNsigmaTPCDataCorr[iEta] = 0.;
  }

  fAliEventCuts.SetManualMode();

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEHFSystPID::~AliAnalysisTaskSEHFSystPID()
{
  // Destructor
  if (fOutputList) {
    delete fHistNEvents;
    for(int iHisto=0; iHisto<5; iHisto++) delete fHistArmenteroPlot[iHisto];
    delete fHistQtVsMassKinks;
    delete fHistPDaughterVsMotherKink;
    delete fHistOpeningAngleVsPMotherKink;
    delete fHistdEdxVsPMotherKink;
    delete fHistNTPCclsVsRadius;
    delete fOutputList;
  }

  if(fPIDtree) delete fPIDtree;
  if(fESDtrackCuts) delete fESDtrackCuts;
  if(fV0cuts) delete fV0cuts;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFSystPID::UserCreateOutputObjects()
{
  // create track cuts
  if(!fESDtrackCuts) {
    fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,0);
    fESDtrackCuts->SetEtaRange(-0.8, 0.8);
    fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  }

  // V0 Kine cuts
  fV0cuts = new AliAODv0KineCuts;
  fV0cuts->SetNoKinks(false);

  fOutputList = new TList();
  fOutputList->SetOwner(true);

  fHistNEvents = new TH1F("fHistNEvents","Number of processed events;;Number of events",14,-1.5,12.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"Read from AOD");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Mismatch AOD");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Pass Phys. Sel. + Trig");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"No vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Vertex contributors < 1");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Without SPD vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Error on zVertex>0.5");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"|zVertex|>10");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"Good Z vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"excluded time range");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"Cent corr cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(12,"V0mult vs. nTPC cls");
  fHistNEvents->GetXaxis()->SetBinLabel(13,"OOB ITS vs. nTPC cls");
  fHistNEvents->GetXaxis()->SetBinLabel(14,"Selected events");
  fOutputList->Add(fHistNEvents);

  TString armenteronames[5] = {"All","K0s","Lambda","AntiLambda","Gamma"};
  TString armenterolabels[5] = {"all","K_{s}^{0}","#Lambda","#bar{#Lambda}","#gamma"};
  for(int iHisto=0; iHisto<5; iHisto++) {
    fHistArmenteroPlot[iHisto] = new TH2F(Form("fHistArmenteroPlot%s",armenteronames[iHisto].Data()),Form("Armenteros-Podolanski plot - %s;#it{#alpha};#it{q}_{T} (GeV/#it{c})",armenterolabels[iHisto].Data()),200,-1.,1.,200,0.,0.4);
    fOutputList->Add(fHistArmenteroPlot[iHisto]);
  }

  fHistQtVsMassKinks = new TH2F("fHistQtVsMassKinks","#it{q}_{T} vs. #it{M}(#mu#nu) kink mother;#it{M}(#mu#nu) (GeV/#it{c}^{2});#it{q}_{T} (GeV/#it{c})",250,0.1,0.6,150,0.,0.3);
  fOutputList->Add(fHistQtVsMassKinks);

  fHistPDaughterVsMotherKink = new TH2F("fHistPDaughterVsMotherKink","#it{p} (daughter) vs. #it{p} (mother) kinks;#it{p}^{mother} (GeV/#it{c});#it{p}^{dau} (GeV/#it{c})",100,0,20,100,0,20);
  fOutputList->Add(fHistPDaughterVsMotherKink);

  fHistOpeningAngleVsPMotherKink = new TH2F("fHistOpeningAngleVsPMotherKink","#it{#alpha} vs. #it{p} (mother) kinks;#it{p}^{mother} (GeV/#it{c});#it{#alpha} (deg)",100,0,20,90,0,90);
  fOutputList->Add(fHistOpeningAngleVsPMotherKink);

  fHistdEdxVsPMotherKink = new TH2F("fHistdEdxVsPMotherKink","dE/dx vs. #it{p} (mother) kinks;#it{p}^{mother} (GeV/#it{c});d#it{E}/d#it{x} (a.u.)",100,0,20,125,0,250);
  fOutputList->Add(fHistdEdxVsPMotherKink);

  fHistNTPCclsVsRadius = new TH2F("fHistNTPCclsVsRadius","N TPC clusters (mother) vs. #it{R} kinks;#it{R} (cm);N TPC clusters (mother)",50,0,250,160,-0.5,159.5);
  fOutputList->Add(fHistNTPCclsVsRadius);

  TString detnames[kNMaxDet]       = {"ITS","TPC","TOF","HMPID"};
  TString partnameshort[kNMaxHypo] = {"pi","K","p","e","d","t","He3"};
  TString hyponames[kNMaxHypo]     = {"Pion","Kaon","Proton","Electron","Deuteron","Triton","He3"};

  if(fIsMC) {
    for(int iDet=0; iDet<kNMaxDet; iDet++) {
      if(!fEnabledDet[iDet])
        continue;
      for(int iHypo=0; iHypo<kNMaxHypo; iHypo++) {
        if(!fEnabledSpecies[iHypo])
          continue;
        fHistNsigmaVsPt[iDet][iHypo] = new TH2F(Form("fHistNsigma%svsPt_%s",detnames[iDet].Data(),hyponames[iHypo].Data()),Form(";#it{p}_{T} (GeV/#it{c});N_{#sigma}^{%s}(%s)",detnames[iDet].Data(),hyponames[iHypo].Data()),500,0,50,1000,-50,50);
        fOutputList->Add(fHistNsigmaVsPt[iDet][iHypo]);
      }
    }
  }

  fPIDtree = new TTree("fPIDtree","fPIDtree");
  fPIDtree->Branch("pT",&fPt,"pT/s");
  fPIDtree->Branch("eta",&fEta,"eta/S");
  fPIDtree->Branch("phi",&fPhi,"phi/s");
  if(fFillTreeWithPIDInfo) {
    if(fEnabledDet[kITS])
      fPIDtree->Branch("p",&fP,"p/s");
    if(fEnabledDet[kTPC])
      fPIDtree->Branch("pTPC",&fPTPC,"pTPC/s");
    if(fEnabledDet[kTOF])
      fPIDtree->Branch("pTOF",&fPTOF,"pTOF/s");
    if(fEnabledDet[kHMPID])
      fPIDtree->Branch("pHMPID",&fPHMPID,"pHMPID/s");
    if(!fFillTreeWithRawPIDOnly) {
      for(int iDet=0; iDet<kNMaxDet; iDet++) {
        if(!fEnabledDet[iDet])
          continue;
        for(int iHypo=0; iHypo<kNMaxHypo; iHypo++) {
          if(!fEnabledSpecies[iHypo])
            continue;
          fPIDtree->Branch(Form("n_sigma_%s_%s",detnames[iDet].Data(),partnameshort[iHypo].Data()),&fPIDNsigma[iDet][iHypo],Form("n_sigma_%s_%s/S",detnames[iDet].Data(),partnameshort[iHypo].Data()));
        }
      }
    }
    if(!fFillTreeWithNsigmaPIDOnly) {
      if(fEnabledDet[kITS]) {
        fPIDtree->Branch("dEdxITS",&fdEdxITS,"dEdxITS/s");
        fPIDtree->Branch("NclusterPIDTPC",&fTPCNclsPID,"NclusterPIDTPC/b");
      }
      if(fEnabledDet[kTPC]) {
        fPIDtree->Branch("dEdxTPC",&fdEdxTPC,"dEdxTPC/s");
        fPIDtree->Branch("ITSclsMap",&fITSclsMap,"ITSclsMap/b");
      }
      if(fEnabledDet[kTOF]) {
        fPIDtree->Branch("ToF",&fToF,"ToF/s");
        fPIDtree->Branch("TrackLength",&fTrackLength,"TrackLength/s");
        fPIDtree->Branch("StartTimeRes",&fStartTimeRes,"StartTimeRes/s");
      }
      if(fEnabledDet[kHMPID]) {
        fPIDtree->Branch("HMPIDsig",&fHMPIDsignal,"HMPIDsig/s");
        fPIDtree->Branch("HMPIDocc",&fHMPIDoccupancy,"HMPIDocc/s");
      }
    }
  }
  if(fFillTreeWithTrackQualityInfo) {
      fPIDtree->Branch("NclusterTPC",&fTPCNcls,"NclusterTPC/b");
      fPIDtree->Branch("NcrossedRowsTPC",&fTPCNcrossed,"NcrossedRowsTPC/b");
      fPIDtree->Branch("NFindableTPC",&fTPCFindable,"NFindableClustersTPC/b");
  }
  fPIDtree->Branch("trackbits",&fTrackInfoMap,"trackbits/b"); // basic track info always filled
  if(fSystem == 1)
    fPIDtree->Branch("OOBpileupbits",&fOOBPileupMap,"OOBpileupbits/b");

  fPIDtree->Branch("tag",&fTag,"tag/s");
  if(fIsMC) fPIDtree->Branch("PDGcode",&fPDGcode,"PDGcode/I");

  if(fUseAliEventCuts) { //add QA plots if event cuts used
    fAliEventCuts.AddQAplotsToList(fOutputList,true);
  }

  // post data
  PostData(1, fOutputList);
  PostData(2, fPIDtree);
}

//________________________________________________________________________
void AliAnalysisTaskSEHFSystPID::UserExec(Option_t */*option*/)
{
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
  }
  if (!fAOD) {
    AliWarning("AliAnalysisTaskSEHFSystPID::Exec(): bad AOD");
    PostData(1, fOutputList);
    return;
  }
  fHistNEvents->Fill(-1);

  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistNEvents->Fill(0);
      PostData(1, fOutputList);
      return;
    }
  }

  if(TMath::Abs(fAOD->GetMagneticField())<0.001) return;

  AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(!aodHandler) {
    AliWarning("AliAnalysisTaskSEHFSystPID::Exec(): No AliInputEventHandler!");
    return;
  }

  unsigned int maskPhysSel = aodHandler->IsEventSelected();
  TString firedTriggerClasses = fAOD->GetFiredTriggerClasses();
  if(!fIsMC && (fAOD->GetRunNumber()<136851 || fAOD->GetRunNumber()>139517)) {
    if(!(firedTriggerClasses.Contains(fTriggerClass.Data()))) return;
  }
  if((maskPhysSel & fTriggerMask)==0) return;

  if(!IsCentralitySelected()) {
    PostData(1, fOutputList);
    return;
  }

  fHistNEvents->Fill(1);

  if (!IsVertexAccepted()) {
    fHistNEvents->Fill(2);
    PostData(1, fOutputList);
    return;
  }

  int selEvCuts = 0;
  if(fUseAliEventCuts) {
    selEvCuts = IsEventSelectedWithAliEventCuts();
    fHistNEvents->Fill(selEvCuts);
  }
  else {
    if(fUseTimeRangeCutForPbPb2018 && !fIsMC){
      if(fAOD->GetRunNumber() != fRunNumberPrevEvent){
        fTimeRangeCut.InitFromRunNumber(fAOD->GetRunNumber());
      }
      if(fTimeRangeCut.CutEvent(fAOD)){
        fHistNEvents->Fill(8);
        PostData(1, fOutputList);
        return;
      }
    }
  }

  if(fKeepOnlyPbPbOutOfBunchPileupCutsITSTPC) {
    if(selEvCuts > 0 && selEvCuts != 11) {
      PostData(1, fOutputList);
      return;
    }
  }
  else {
    if(selEvCuts > 0) {
      PostData(1, fOutputList);
      return;
    }
  }

  fHistNEvents->Fill(12);

  // tag OOB pileup and store info in the tree for each track
  if(fSystem == 1)
    TagOOBPileUpEvent();

  // load MC particles
  TClonesArray *arrayMC=0;
  if(fIsMC){
    arrayMC =  (TClonesArray*)fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      AliWarning("AliAnalysisTaskSEDplus::UserExec: MC particles branch not found!\n");
      return;
    }
  }

  //get pid response
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!fPIDresp)
    fPIDresp = inputHandler->GetPIDResponse();

  //load data correction for NsigmaTPC, if enabled
  if(fEnableNsigmaTPCDataCorr) {
    if(fAOD->GetRunNumber()!=fRunNumberPrevEvent) {
      Bool_t isPass1 = kFALSE;
      TTree *treeAOD = inputHandler->GetTree();
      TString currentFile = treeAOD->GetCurrentFile()->GetName();
      if((currentFile.Contains("LHC18q") || currentFile.Contains("LHC18r")) && currentFile.Contains("pass1"))
        isPass1 = kTRUE;

      AliAODPidHF::SetNsigmaTPCDataDrivenCorrection(fAOD->GetRunNumber(), fSystNsigmaTPCDataCorr, fNPbinsNsigmaTPCDataCorr, fPlimitsNsigmaTPCDataCorr, fNEtabinsNsigmaTPCDataCorr, fEtalimitsNsigmaTPCDataCorr, fMeanNsigmaTPCPionData, fMeanNsigmaTPCKaonData, fMeanNsigmaTPCProtonData, fSigmaNsigmaTPCPionData, fSigmaNsigmaTPCKaonData, fSigmaNsigmaTPCProtonData, isPass1);
    }
  }

  // V0 selection
  if(fSystem == 0)
    fV0cuts->SetMode(AliAODv0KineCuts::kPurity, AliAODv0KineCuts::kPP);
  else if(fSystem == 1)
    fV0cuts->SetMode(AliAODv0KineCuts::kPurity, AliAODv0KineCuts::kPbPb);
  fV0cuts->SetEvent(fAOD);

  vector<short> idPionFromK0s;
  vector<short> idPionFromL;
  vector<short> idProtonFromL;
  vector<short> idElectronFromGamma;
  vector<short> idKaonFromKinks;
  GetTaggedV0s(idPionFromK0s, idPionFromL, idProtonFromL, idElectronFromGamma);
  GetTaggedKaonsFromKinks(idKaonFromKinks);

  vector<short>::iterator it;

  const int nTracks = fAOD->GetNumberOfTracks();
  //loop on tracks

  for(int iTrack=0; iTrack<nTracks; iTrack++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) continue;
    //applying ESDtrackCut
    if(!fESDtrackCuts->IsSelected(track)) continue;

    if(fEnabledDownSampling && fFracToKeepDownSampling<1. && track->Pt()<fPtMaxDownSampling) {
      double pseudoRand = track->Pt()*1000.-(long)(track->Pt()*1000);
      if(fDownSamplingOpt==0 && pseudoRand>fFracToKeepDownSampling) continue; //keep tracks with pseudorand < fFracToKeepDownSampling
      else if(fDownSamplingOpt==1 && pseudoRand<(1-fFracToKeepDownSampling)) continue; //keep tracks with pseudorand > 1-fFracToKeepDownSampling
    }

    fPt = ConvertFloatToUnsignedShort(track->Pt()*1000);
    fP = ConvertFloatToUnsignedShort(track->P()*1000);
    fPTPC = ConvertFloatToUnsignedShort(track->GetTPCmomentum()*1000);
    fPTOF = ConvertFloatToUnsignedShort(GetTOFmomentum(track)*1000);
    double momHMPID[3];
    bool momHMPIDok = track->GetOuterHmpPxPyPz(momHMPID);
    if(momHMPIDok)
      fPHMPID = ConvertFloatToUnsignedShort(TMath::Sqrt(momHMPID[0]*momHMPID[0]+momHMPID[1]*momHMPID[1]+momHMPID[2]*momHMPID[2])*1000);
    else
      fPHMPID = numeric_limits<short>::min();

    fEta = ConvertFloatToShort(track->Eta()*1000);
    fPhi = ConvertFloatToUnsignedShort(track->Phi()*1000);

    //charge
    if(track->Charge()>0) {
      fTag |= kPositiveTrack;
      fTag &= ~kNegativeTrack;
    }
    else if(track->Charge()<0) {
      fTag |= kNegativeTrack;
      fTag &= ~kPositiveTrack;
    }

    if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))
      fTrackInfoMap |= kHasSPDAny;
    else
      fTrackInfoMap &= ~kHasSPDAny;
    if(track->HasPointOnITSLayer(1))
      fTrackInfoMap |= kHasSPDFirst;
    else
      fTrackInfoMap &= ~kHasSPDFirst;
    if(track->GetStatus() & AliESDtrack::kITSrefit)
      fTrackInfoMap |= kHasITSrefit;
    else
      fTrackInfoMap &= ~kHasITSrefit;
    if(track->GetStatus() & AliESDtrack::kTPCrefit)
      fTrackInfoMap |= kHasTPCrefit;
    else
      fTrackInfoMap &= ~kHasTPCrefit;
    if(IsSelectedByGeometricalCut(track))
      fTrackInfoMap |= kPassGeomCut;
    else
      fTrackInfoMap &= ~kPassGeomCut;

    if(fFillTreeWithTrackQualityInfo) {
      fTPCNcls = static_cast<unsigned char>(track->GetTPCNcls());
      fTPCNcrossed = static_cast<unsigned char>(track->GetTPCNCrossedRows());
      fTPCFindable = static_cast<unsigned char>(track->GetTPCNclsF());
    }

    //PID
    track->SetTOFsignalTunedOnData(100000); // force tune-on-data to have latest development of tail parametrisation in old AODs
    bool isDetOk[kNMaxDet] = {false,false,false};
    if (fPIDresp->CheckPIDStatus(AliPIDResponse::kITS,track) == AliPIDResponse::kDetPidOk) {
      fTrackInfoMap &= ~kHasNoITS;
      isDetOk[kITS] = true;
    }
    else {
      fTrackInfoMap |= kHasNoITS;
    }
    if (fPIDresp->CheckPIDStatus(AliPIDResponse::kTPC,track) == AliPIDResponse::kDetPidOk) {
      fTrackInfoMap &= ~kHasNoTPC;
      isDetOk[kTPC] = true;
    }
    else {
      fTrackInfoMap &= ~kHasNoTOF;
    }
    if (fPIDresp->CheckPIDStatus(AliPIDResponse::kTOF,track) == AliPIDResponse::kDetPidOk) {
      fTrackInfoMap &= ~kHasNoTOF;
      isDetOk[kTOF] = true;
    }
    else {
      fTrackInfoMap |= kHasNoTOF;
    }
    if (fPIDresp->CheckPIDStatus(AliPIDResponse::kHMPID,track) == AliPIDResponse::kDetPidOk) {
      isDetOk[kHMPID] = true;
    }

    if(fFillTreeWithPIDInfo) {
      if(!fFillTreeWithNsigmaPIDOnly) { //raw variables
        //ITS variables
        fdEdxITS = ConvertFloatToUnsignedShort(track->GetITSsignal()*fConversionFactordEdx);
        fITSclsMap = track->GetITSClusterMap();

        //TPC variables
        fdEdxTPC = ConvertFloatToUnsignedShort(track->GetTPCsignal()*fConversionFactordEdx);
        fTPCNclsPID = static_cast<unsigned char>(track->GetTPCsignalN());

        //TOF variables
        fTrackLength = ConvertFloatToUnsignedShort(track->GetIntegratedLength()*10);
        fStartTimeRes = ConvertFloatToUnsignedShort(fPIDresp->GetTOFResponse().GetStartTimeRes(track->P())*100);

        if (!(track->GetStatus() & AliESDtrack::kTOFout) || !(track->GetStatus() & AliESDtrack::kTIME)) {
          fToF = 0;
        }
        else {
          if (fTrackLength < 3500) fToF = 0;
          else {
            float tof = track->GetTOFsignal();
            float time0 = fPIDresp->GetTOFResponse().GetStartTime(track->P());
            fToF = ConvertFloatToUnsignedShort((tof-time0)/10);
          }
        }

        //HMPID variables
        fHMPIDsignal = ConvertFloatToUnsignedShort(track->GetHMPIDsignal()*100);
        fHMPIDoccupancy = ConvertFloatToUnsignedShort(track->GetHMPIDoccupancy()*100);
      }
      if(!fFillTreeWithRawPIDOnly) { // nsigma
        for(int iDet=0; iDet<kNMaxDet; iDet++) {
          if(isDetOk[iDet]) {
            FillNsigma(iDet, track);
          }
          else {
            for(int iHypo=0; iHypo<kNMaxHypo; iHypo++) {
              fPIDNsigma[iDet][iHypo] = numeric_limits<short>::min();
            }
          }
        }
      }
    }

    short trackid = track->GetID();

    bool filltree = false;
    it = find(idPionFromK0s.begin(),idPionFromK0s.end(),trackid);
    if(it!=idPionFromK0s.end()) {
      fTag |= kIsPionFromK0s;
      filltree = true;
    }
    else
      fTag &= ~kIsPionFromK0s;

    it = find(idPionFromL.begin(),idPionFromL.end(),trackid);
    if(it!=idPionFromL.end()) {
      fTag |= kIsPionFromL;
      filltree = true;
    }
    else
      fTag &= ~kIsPionFromL;

    it = find(idProtonFromL.begin(),idProtonFromL.end(),trackid);
    if(it!=idProtonFromL.end()) {
      filltree = true;
      fTag |= kIsProtonFromL;
    }
    else
      fTag &= ~kIsProtonFromL;

    it = find(idElectronFromGamma.begin(),idElectronFromGamma.end(),trackid);
    if(it!=idElectronFromGamma.end()) {
      filltree = true;
      fTag |= kIsElectronFromGamma;
    }
    else
      fTag &= ~kIsElectronFromGamma;

    it = find(idKaonFromKinks.begin(),idKaonFromKinks.end(),trackid);
    if(it!=idKaonFromKinks.end()) {
      filltree = true;
      fTag |= kIsKaonFromKinks;
    }
    else
      fTag &= ~kIsKaonFromKinks;

    if(TMath::Abs(fPIDresp->NumberOfSigmasTOF(track,AliPID::kKaon))<fNsigmaMaxForKaonTag && isDetOk[kTOF]) {
      filltree = true;
      fTag |= kIsKaonFromTOF;
    }
    else
      fTag &= ~kIsKaonFromTOF;

    if(TMath::Abs(fPIDresp->NumberOfSigmasTPC(track,AliPID::kKaon))<fNsigmaMaxForKaonTag && isDetOk[kTPC]) {
      filltree = true;
      fTag |= kIsKaonFromTPC;
    }
    else
      fTag &= ~kIsKaonFromTPC;

    if(TMath::Abs(fPIDresp->NumberOfSigmasHMPID(track,AliPID::kKaon))<fNsigmaMaxForKaonTag && isDetOk[kHMPID]) {
      filltree = true;
      fTag |= kIsKaonFromHMPID;
    }
    else
      fTag &= ~kIsKaonFromHMPID;

    if(isDetOk[kTPC] && isDetOk[kTOF] && TMath::Abs(fPIDresp->NumberOfSigmasTPC(track,AliPID::kDeuteron))<fNsigmaMaxForNucleiTag && TMath::Abs(fPIDresp->NumberOfSigmasTOF(track,AliPID::kDeuteron))<fNsigmaMaxForNucleiTag) {
      filltree = true;
      fTag |= kIsDeuteronFromTPCTOF;
    }
    else
      fTag &= ~kIsDeuteronFromTPCTOF;

    if(isDetOk[kTPC] && isDetOk[kTOF] && TMath::Abs(fPIDresp->NumberOfSigmasTPC(track,AliPID::kTriton))<fNsigmaMaxForNucleiTag && TMath::Abs(fPIDresp->NumberOfSigmasTOF(track,AliPID::kTriton))<fNsigmaMaxForNucleiTag) {
      filltree = true;
      fTag |= kIsTritonFromTPCTOF;
    }
    else
      fTag &= ~kIsTritonFromTPCTOF;

    if(isDetOk[kTPC] && isDetOk[kTOF] && TMath::Abs(fPIDresp->NumberOfSigmasTPC(track,AliPID::kHe3))<fNsigmaMaxForNucleiTag && TMath::Abs(fPIDresp->NumberOfSigmasTOF(track,AliPID::kHe3))<fNsigmaMaxForNucleiTag) {
      filltree = true;
      fTag |= kIsHe3FromTPCTOF;
    }
    else
      fTag &= ~kIsHe3FromTPCTOF;

    if(fIsMC) {
      fPDGcode = GetPDGcodeFromMC(track,arrayMC);
      AliPIDResponse::EDetector det[4] = {AliPIDResponse::kITS,AliPIDResponse::kTPC,AliPIDResponse::kTOF,AliPIDResponse::kHMPID};
      for(int iDet=0; iDet<kNMaxDet; iDet++) {
        if(!fEnabledDet[iDet])
          continue;

        switch(fPDGcode) {
          case 211:
            if(fEnabledSpecies[kPion])
              fHistNsigmaVsPt[iDet][kPion]->Fill(track->Pt(),fPIDresp->NumberOfSigmas(det[iDet],track,AliPID::kPion));
            break;
          case 321:
            if(fEnabledSpecies[kKaon])
              fHistNsigmaVsPt[iDet][kKaon]->Fill(track->Pt(),fPIDresp->NumberOfSigmas(det[iDet],track,AliPID::kKaon));
            break;
          case 2212:
            if(fEnabledSpecies[kProton])
              fHistNsigmaVsPt[iDet][kProton]->Fill(track->Pt(),fPIDresp->NumberOfSigmas(det[iDet],track,AliPID::kProton));
            break;
          case 11:
            if(fEnabledSpecies[kElectron])
              fHistNsigmaVsPt[iDet][kElectron]->Fill(track->Pt(),fPIDresp->NumberOfSigmas(det[iDet],track,AliPID::kElectron));
            break;
          case 1000010020:
            if(fEnabledSpecies[kDeuteron])
              fHistNsigmaVsPt[iDet][kDeuteron]->Fill(track->Pt(),fPIDresp->NumberOfSigmas(det[iDet],track,AliPID::kDeuteron));
            break;
          case 1000010030:
            if(fEnabledSpecies[kTriton])
              fHistNsigmaVsPt[iDet][kTriton]->Fill(track->Pt(),fPIDresp->NumberOfSigmas(det[iDet],track,AliPID::kTriton));
            break;
          case 1000020030:
            if(fEnabledSpecies[kHe3])
              fHistNsigmaVsPt[iDet][kHe3]->Fill(track->Pt(),fPIDresp->NumberOfSigmas(det[iDet],track,AliPID::kHe3));
            break;
        }
      }
    }
    else fPDGcode = 0;

    if(filltree && ((fIsMC && fPDGcode>=0) || !fIsMC)) fPIDtree->Fill();

    fTag = 0;
    fTrackInfoMap = 0;

    // Post output data
    PostData(2, fPIDtree);
  }

  //clear vectors of ids
  idPionFromK0s.clear();
  idPionFromL.clear();
  idProtonFromL.clear();
  idElectronFromGamma.clear();
  idKaonFromKinks.clear();

  fRunNumberPrevEvent = fAOD->GetRunNumber();

  // Post output data
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskSEHFSystPID::EnableParticleSpecies(bool pi, bool kao, bool pr, bool el, bool deu, bool tr, bool He3)
{
  // function to enable/disable different particle species
  for(int iHypo=0; iHypo<kNMaxHypo; iHypo++)
    fEnabledSpecies[iHypo] = false;

  if(pi) fEnabledSpecies[kPion]      = true;
  if(kao) fEnabledSpecies[kKaon]     = true;
  if(pr) fEnabledSpecies[kProton]    = true;
  if(el) fEnabledSpecies[kElectron]  = true;
  if(deu) fEnabledSpecies[kDeuteron] = true;
  if(tr) fEnabledSpecies[kTriton]    = true;
  if(He3) fEnabledSpecies[kHe3]      = true;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFSystPID::EnableDetectors(bool ITS, bool TPC, bool TOF, bool HMPID)
{
  // function to enable/disable different PID detectors
  for(int iDet=0; iDet<kNMaxDet; iDet++)
    fEnabledDet[iDet] = false;

  if(ITS)   fEnabledDet[kITS]  = true;
  if(TPC)   fEnabledDet[kTPC]  = true;
  if(TOF)   fEnabledDet[kTOF]  = true;
  if(HMPID) fEnabledDet[kHMPID]  = true;
}

//________________________________________________________________________
bool AliAnalysisTaskSEHFSystPID::IsVertexAccepted()
{
  // function to check if a proper vertex is reconstructed and write z-position in vertexZ
  const AliAODVertex *vertex = fAOD->GetPrimaryVertex();
  if(!vertex){
    fHistNEvents->Fill(2);
    return false;
  }
  else{
    TString title=vertex->GetTitle();
    if(title.Contains("Z") || title.Contains("3D")) return false;
    if(vertex->GetNContributors()<1) {
      fHistNEvents->Fill(3);
      return false;
    }
  }

  const AliVVertex *vSPD = fAOD->GetPrimaryVertexSPD();
  if(!vSPD || (vSPD && vSPD->GetNContributors()<1)){
    fHistNEvents->Fill(4);
    return false;
  }
  else{
    double dz = vSPD->GetZ()-vertex->GetZ();
    if(TMath::Abs(dz)>0.5) {
      fHistNEvents->Fill(5);
      return false;
    }
    double covTrc[6],covSPD[6];
    vertex->GetCovarianceMatrix(covTrc);
    vSPD->GetCovarianceMatrix(covSPD);
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
    if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) {
      fHistNEvents->Fill(5);
      return false;
    }
  }

  if(TMath::Abs(vertex->GetZ())>10.) {
    fHistNEvents->Fill(6);
    return false;
  }
  fHistNEvents->Fill(7);

  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskSEHFSystPID::IsCentralitySelected() {

  if(fCentEstimator==kCentOff) return true;

  AliMultSelection *multSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if(!multSelection){
    AliWarning("AliMultSelection could not be found in the aod event list of objects");
    return false;
  }

  float cent=-999;
  if(fCentEstimator==kCentV0M)      cent=multSelection->GetMultiplicityPercentile("V0M");
  else if(fCentEstimator==kCentV0A) cent=multSelection->GetMultiplicityPercentile("V0A");
  else if(fCentEstimator==kCentZNA) cent=multSelection->GetMultiplicityPercentile("ZNA");
  else if(fCentEstimator==kCentCL1) cent=multSelection->GetMultiplicityPercentile("CL1");
  else if(fCentEstimator==kCentCL0) cent=multSelection->GetMultiplicityPercentile("CL0");
  else {
    AliWarning(Form("CENTRALITY ESTIMATE WITH ESTIMATOR %d NOT YET IMPLEMENTED FOR NEW FRAMEWORK",fCentEstimator));
    return false;
  }

  if(cent>=fCentMin && cent<=fCentMax) return true;
  return false;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFSystPID::GetTaggedV0s(vector<short> &idPionFromK0s, vector<short> &idPionFromL, vector<short> &idProtonFromL, vector<short> &idElectronFromGamma) {
  // tag tracks from V0 decays
  const int nV0s = fAOD->GetNumberOfV0s();
  AliAODv0 *V0=nullptr;

  for (int iV0=0; iV0<nV0s; iV0++){
    V0 = (AliAODv0*)fAOD->GetV0(iV0);
    if(!V0) continue;
    //if(V0->GetOnFlyStatus()) continue;

    AliAODTrack* pTrack=dynamic_cast<AliAODTrack*>(V0->GetDaughter(0));
    AliAODTrack* nTrack=dynamic_cast<AliAODTrack*>(V0->GetDaughter(1));
    if(!fESDtrackCuts->IsSelected(pTrack) || !fESDtrackCuts->IsSelected(nTrack)) continue;

    // Get the particle selection
    bool foundV0 = false;
    int pdgV0=-1, pdgP=-1, pdgN=-1;
    foundV0 = fV0cuts->ProcessV0(V0, pdgV0, pdgP, pdgN);
    if(!foundV0) continue;

    // v0 Armenteros plot (QA)
    fHistArmenteroPlot[0]->Fill(V0->AlphaV0(),V0->PtArmV0());

    short iTrackP = V0->GetPosID();  // positive track
    short iTrackN = V0->GetNegID();  // negative track

    if(pdgP==211 && pdgN==-211) {
      idPionFromK0s.push_back(iTrackP);
      idPionFromK0s.push_back(iTrackN);
      fHistArmenteroPlot[1]->Fill(V0->AlphaV0(),V0->PtArmV0());
    }
    else if(pdgP==2212 && pdgN==-211) {
      idProtonFromL.push_back(iTrackP);
      idPionFromL.push_back(iTrackN);
      fHistArmenteroPlot[2]->Fill(V0->AlphaV0(),V0->PtArmV0());
    }
    else if(pdgP==211 && pdgN==-2212) {
      idPionFromL.push_back(iTrackP);
      idProtonFromL.push_back(iTrackN);
      fHistArmenteroPlot[3]->Fill(V0->AlphaV0(),V0->PtArmV0());
    }
    else if(pdgP==-11 && pdgN==11) {
      idElectronFromGamma.push_back(iTrackP);
      idElectronFromGamma.push_back(iTrackN);
      fHistArmenteroPlot[4]->Fill(V0->AlphaV0(),V0->PtArmV0());
    }
  }
}

//________________________________________________________________________
int AliAnalysisTaskSEHFSystPID::GetPDGcodeFromMC(AliAODTrack* track, TClonesArray* arrayMC)
{
  // Get pdg code
  int pdg = -1;
  if(!track) return pdg;
  int label = track->GetLabel();
  if(label<0) return pdg;
  AliAODMCParticle* partMC = dynamic_cast<AliAODMCParticle*>(arrayMC->At(label));
  if(!partMC) return pdg;
  pdg = TMath::Abs(partMC->GetPdgCode());

  return pdg;
}

//___________________________________________________________________
AliAODTrack* AliAnalysisTaskSEHFSystPID::IsKinkDaughter(AliAODTrack* track)
{
  // Check if track is a kink daughter --> in this case returns the mother track
  AliAODVertex *maybeKink=track->GetProdVertex();
  if(!maybeKink) return nullptr;
  if(maybeKink->GetType()==AliAODVertex::kKink) return ((AliAODTrack *)maybeKink->GetParent());

  return nullptr;
}

//___________________________________________________________________
void AliAnalysisTaskSEHFSystPID::GetTaggedKaonsFromKinks(vector<short> &idKaonFromKinks)
{
  // Tag kink mother tracks

  const int nTracks = fAOD->GetNumberOfTracks();
  for(int iTrack=0; iTrack<nTracks; iTrack++) {

    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(track->GetTPCNcls()<30 || track->IsPrimaryCandidate()) continue;
    AliAODTrack* kinkmothertrack = IsKinkDaughter(track);
    if(!kinkmothertrack || !fESDtrackCuts->IsSelected(kinkmothertrack) || !kinkmothertrack->IsPrimaryCandidate()) continue;
    AliAODVertex *prodvtx=track->GetProdVertex();
    float R = TMath::Sqrt(prodvtx->GetX()*prodvtx->GetX()+prodvtx->GetY()*prodvtx->GetY());
    if(R<fRMinKinks || R>fRMaxKinks) continue;
    if(fPIDresp->NumberOfSigmasTPC(track,AliPID::kMuon)>3.5 && fPIDresp->NumberOfSigmasTOF(track,AliPID::kMuon)>3.5) continue;

    //kinks at vertex
    float bfield = fAOD->GetMagneticField();
    AliESDtrack kinkMotherESD(kinkmothertrack);
    AliESDtrack kinkDaughterESD(track);
    double xKinkMother=-999;
    double xKinkDaughter=-999;
    double minDist = kinkMotherESD.GetDCA(&kinkDaughterESD, bfield, xKinkMother, xKinkDaughter);
		if (minDist > 2.0) continue;
		double mother3MomentumKinkVtxArray[3] = {0, 0, 0};
		double daughter3MomentumKinkVtxArray[3] = {0, 0, 0};
		if (!kinkMotherESD.GetPxPyPzAt(xKinkMother, bfield, mother3MomentumKinkVtxArray)) continue;
		if (!kinkDaughterESD.GetPxPyPzAt(xKinkDaughter, bfield, daughter3MomentumKinkVtxArray)) continue;
		TVector3 mother3MomentumDCA(mother3MomentumKinkVtxArray);
		TVector3 daughter3MomentumDCA(daughter3MomentumKinkVtxArray);
		float qt = daughter3MomentumDCA.Perp(mother3MomentumDCA);
    float openingAngle = daughter3MomentumDCA.Angle(mother3MomentumDCA);
    openingAngle *= 180/TMath::Pi();
    if(qt<fQtMinKinks || qt>0.300) continue;
    if(openingAngle<2 || openingAngle>MaxOpeningAngleKnu(kinkmothertrack->P())) continue;

    //inv-mass
    float massmu = TDatabasePDG::Instance()->GetParticle(13)->Mass(); //muon
    TVector3 transferedMom = mother3MomentumDCA-daughter3MomentumDCA;
		float energyDaughterMu = TMath::Sqrt(daughter3MomentumDCA.Mag()*daughter3MomentumDCA.Mag()+massmu*massmu);
    float invmassMuNu = (energyDaughterMu+transferedMom.Mag())*(energyDaughterMu+transferedMom.Mag())-mother3MomentumDCA.Mag()*mother3MomentumDCA.Mag();
		if(invmassMuNu>0) invmassMuNu = TMath::Sqrt(invmassMuNu);
    if(invmassMuNu<0 || invmassMuNu>0.8) continue;

    //N-clusters vs. R
    int TPCnclsMother = kinkmothertrack->GetTPCNcls();
    int TPCnclsMax = -31.67+(11./12.)*R;
    int TPCnclsMin = -85.5 +(65./95.)*R;
    if(TPCnclsMother<TPCnclsMin || TPCnclsMother>TPCnclsMax) continue;

    idKaonFromKinks.push_back(kinkmothertrack->GetID());
    fHistQtVsMassKinks->Fill(invmassMuNu,qt);
    fHistPDaughterVsMotherKink->Fill(kinkmothertrack->P(),track->P());
    fHistOpeningAngleVsPMotherKink->Fill(kinkmothertrack->P(),openingAngle);
    fHistdEdxVsPMotherKink->Fill(kinkmothertrack->P(),kinkmothertrack->GetTPCsignal());
    fHistNTPCclsVsRadius->Fill(R,TPCnclsMother);
  }
}

//___________________________________________________________________
float AliAnalysisTaskSEHFSystPID::MaxOpeningAngleKnu(float p) {

  float par0 = 0.493677;
  float par1 = 0.9127037;
  float par2 = TMath::Pi();

  return TMath::ATan(par0*par1*1./TMath::Sqrt(p*p*(1-par1*par1)-(par0*par0*par1*par1)))*180/par2;
}

//________________________________________________________________
float AliAnalysisTaskSEHFSystPID::GetTOFmomentum(AliAODTrack* track)
{
  float t_d = fPIDresp->GetTOFResponse().GetExpectedSignal(track, AliPID::kTriton); //largest mass possible with Z=1
  float len = track->GetIntegratedLength();
  float beta_d = len / (t_d * kCSPEED);
  float mass = AliPID::ParticleMassZ(AliPID::kTriton); //largest mass possible with Z=1

  if(TMath::Abs(beta_d-1.) < 1.e-12) return track->GetTPCmomentum();
  else return mass*beta_d/sqrt(1.-(beta_d*beta_d));
}

//___________________________________________________________________
short AliAnalysisTaskSEHFSystPID::ConvertFloatToShort(float num) {
  if(num>=static_cast<float>(numeric_limits<short>::max())) return numeric_limits<short>::max();
  else if(num<=static_cast<float>(numeric_limits<short>::min())) return numeric_limits<short>::min();

  num = round(num);
  return static_cast<short>(num);
}

//___________________________________________________________________
unsigned short AliAnalysisTaskSEHFSystPID::ConvertFloatToUnsignedShort(float num) {
  if(num>=static_cast<float>(numeric_limits<unsigned short>::max())) return numeric_limits<unsigned short>::max();
  else if(num<=static_cast<float>(numeric_limits<unsigned short>::min())) return numeric_limits<unsigned short>::min();

  num = round(num);
  return static_cast<unsigned short>(num);
}

//________________________________________________________________
void AliAnalysisTaskSEHFSystPID::GetNsigmaTPCMeanSigmaData(float &mean, float &sigma, AliPID::EParticleType species, float pTPC, float eta) {

  int bin = TMath::BinarySearch(fNPbinsNsigmaTPCDataCorr,fPlimitsNsigmaTPCDataCorr,pTPC);
  if(bin<0) bin=0; //underflow --> equal to min value
  else if(bin>fNPbinsNsigmaTPCDataCorr-1) bin=fNPbinsNsigmaTPCDataCorr-1; //overflow --> equal to max value

  int etabin = TMath::BinarySearch(fNEtabinsNsigmaTPCDataCorr,fEtalimitsNsigmaTPCDataCorr,TMath::Abs(eta));
  if(etabin<0) etabin=0; //underflow --> equal to min value
  else if(etabin>fNEtabinsNsigmaTPCDataCorr-1) etabin=fNEtabinsNsigmaTPCDataCorr-1; //overflow --> equal to max value

  switch(species) {
    case AliPID::kPion:
    {
      mean = fMeanNsigmaTPCPionData[etabin][bin];
      sigma = fSigmaNsigmaTPCPionData[etabin][bin];
      break;
    }
    case AliPID::kKaon:
    {
      mean = fMeanNsigmaTPCKaonData[etabin][bin];
      sigma = fSigmaNsigmaTPCKaonData[etabin][bin];
      break;
    }
    case AliPID::kProton:
    {
      mean = fMeanNsigmaTPCProtonData[etabin][bin];
      sigma = fSigmaNsigmaTPCProtonData[etabin][bin];
      break;
    }
    default:
    {
      mean = 0.;
      sigma = 1.;
      break;
    }
  }
}

//________________________________________________________________
int AliAnalysisTaskSEHFSystPID::IsEventSelectedWithAliEventCuts() {

  int run = fAOD->GetRunNumber();
  if(fAOD->GetRunNumber()!=fRunNumberPrevEvent) {
    if(run >= 244917 && run <= 246994)
      fAliEventCuts.SetupRun2PbPb();
    else if(run >= 295369 && run <= 297624)
      fAliEventCuts.SetupPbPb2018();
  }

  // setup cuts
  if(fUseTimeRangeCutForPbPb2018)
    fAliEventCuts.UseTimeRangeCut();
  if(fApplyPbPbOutOfBunchPileupCutsITSTPC)
    fAliEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(true, fApplyPbPbOutOfBunchPileupCutsITSTPC);
  fAliEventCuts.AcceptEvent(fAOD); //for QA plots

  if(fUseTimeRangeCutForPbPb2018){
    if(!fAliEventCuts.PassedCut(AliEventCuts::kTriggerClasses))
      return 8;
  }

  // cut on correlations for out of bunch pileup in PbPb run2
  if(fApplyPbPbOutOfBunchPileupCuts==1){
    if(!fAliEventCuts.PassedCut(AliEventCuts::kCorrelations))
      return 9;
  }
  else if(fApplyPbPbOutOfBunchPileupCuts==2 && run >= 295369 && run <= 297624){
    // Ionut cut on V0multiplicity vs. n TPC clusters (Pb-Pb 2018)
    AliAODVZERO* v0data = (AliAODVZERO*)((AliAODEvent*)fAOD)->GetVZEROData();
    float mTotV0 = v0data->GetMTotV0A()+v0data->GetMTotV0C();
    int nTPCcls = fAOD->GetNumberOfTPCClusters();
    float mV0TPCclsCut = -2000.+(0.013*nTPCcls)+(1.25e-9*nTPCcls*nTPCcls);
    if(mTotV0 < mV0TPCclsCut){
      return 10;
    }
  }

  // cut on ITS-TPC multiplicity correlation for OOB TPC pileup
  // IMPORTANT: it must be the last cut to have the possibility to select events that are good for all the other requirements but not for OOB pileup
  if(fApplyPbPbOutOfBunchPileupCutsITSTPC){
    if(!fAliEventCuts.PassedCut(AliEventCuts::kTPCPileUp))
      return 11;
  }

  return 0;
}

//________________________________________________________________
bool AliAnalysisTaskSEHFSystPID::IsSelectedByGeometricalCut(AliAODTrack* track) {

  // convert to ESD track here
  AliESDtrack esdTrack(track);
  // set the TPC cluster info
  esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(track->GetTPCNclsF());

  float nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
  float lengthInActiveZoneTPC=esdTrack.GetLengthInActiveZone(0,fDeadZoneWidth,220.,fAOD->GetMagneticField());
  double cutGeoNcrNclLength=fCutGeoNcrNclLength-TMath::Power(TMath::Abs(esdTrack.GetSigned1Pt()),fCutGeoNcrNclGeom1Pt);
  if (lengthInActiveZoneTPC<cutGeoNcrNclLength)
    return false;
  if (nCrossedRowsTPC<fCutGeoNcrNclFractionNcr*cutGeoNcrNclLength)
    return false;
  if (esdTrack.GetTPCncls()<fCutGeoNcrNclFractionNcl*cutGeoNcrNclLength)
    return false;

  return true;
}

//________________________________________________________________
bool AliAnalysisTaskSEHFSystPID::FillNsigma(int iDet, AliAODTrack* track) {

  if(iDet>=kNMaxDet)
    return false;

    AliPID::EParticleType hypopidresp[7] = {AliPID::kPion,AliPID::kKaon,AliPID::kProton,AliPID::kElectron,AliPID::kDeuteron,AliPID::kTriton,AliPID::kHe3};

    switch(iDet) {
      case kITS:
      {
        for(int iHypo=0; iHypo<kNMaxHypo; iHypo++)
          fPIDNsigma[iDet][iHypo] = ConvertFloatToShort(fPIDresp->NumberOfSigmasITS(track,hypopidresp[iHypo])*100);
        break;
      }
      case kTPC:
      {
        float nSigmaTPC[kNMaxHypo];

        for(int iHypo=0; iHypo<kNMaxHypo; iHypo++)
          nSigmaTPC[iHypo] = fPIDresp->NumberOfSigmasTPC(track,hypopidresp[iHypo]);

        if(fEnableNsigmaTPCDataCorr) { //only pion, kaon and protons
          float mean[kProton+1], sigma[kProton+1];
          for(int iHypo=0; iHypo<=kProton; iHypo++) {
            GetNsigmaTPCMeanSigmaData(mean[iHypo], sigma[iHypo], hypopidresp[iHypo], track->GetTPCmomentum(), track->Eta());

            if(nSigmaTPC[iHypo]>-990.)
              nSigmaTPC[iHypo] = (nSigmaTPC[iHypo]-mean[iHypo]) / sigma[iHypo];
          }
        }

        for(int iHypo=0; iHypo<kNMaxHypo; iHypo++)
          fPIDNsigma[iDet][iHypo] = ConvertFloatToShort(nSigmaTPC[iHypo]*100);

        break;
      }
      case kTOF:
      {
        for(int iHypo=0; iHypo<kNMaxHypo; iHypo++)
          fPIDNsigma[iDet][iHypo] = ConvertFloatToShort(fPIDresp->NumberOfSigmasTOF(track,hypopidresp[iHypo])*100);

        break;
      }
      case kHMPID:
      {
        for(int iHypo=0; iHypo<kNMaxHypo; iHypo++)
          fPIDNsigma[iDet][iHypo] = ConvertFloatToShort(fPIDresp->NumberOfSigmasHMPID(track,hypopidresp[iHypo])*100);

        break;
      }
    }

  return true;
}

//________________________________________________________________
void AliAnalysisTaskSEHFSystPID::TagOOBPileUpEvent() {

  fOOBPileupMap = 0;
  AliVMultiplicity* mult = fAOD->GetMultiplicity();
  double nTPCcls = static_cast<double>(fAOD->GetNumberOfTPCClusters());
  int nITScls = 0;
  for(int iLay = 2; iLay < 6; iLay++)
    nITScls += mult->GetNumberOfITSClusters(iLay);

  if(nITScls < -16000.+0.0099*nTPCcls+9.426e-10*nTPCcls*nTPCcls)
    fOOBPileupMap |= kVeryLooseITSTPC;
  else
    fOOBPileupMap &= ~kVeryLooseITSTPC;
  if(nITScls < -12000.+0.0099*nTPCcls+9.426e-10*nTPCcls*nTPCcls)
    fOOBPileupMap |= kLooseITSTPC;
  else
    fOOBPileupMap &= ~kLooseITSTPC;
  if(nITScls < -8000.+0.0099*nTPCcls+9.426e-10*nTPCcls*nTPCcls)
    fOOBPileupMap |= kMediumITSTPC;
  else
    fOOBPileupMap &= ~kMediumITSTPC;
  if(nITScls < -3000.+0.0099*nTPCcls+9.426e-10*nTPCcls*nTPCcls)
    fOOBPileupMap |= kTightITSTPC;
  else
    fOOBPileupMap &= ~kTightITSTPC;
}
