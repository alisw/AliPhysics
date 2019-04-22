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
fdEdxTPC(0),
fToF(0),
fPt(0),
fTPCNcls(0),
fTPCNclsPID(0),
fTrackLength(0),
fStartTimeRes(0),
fEta(-9999),
fPDGcode(-1),
fTag(0),
fNsigmaMaxForTag(0.02),
fQtMinKinks(0.15),
fRMinKinks(120),
fRMaxKinks(210),
fCentMin(0.),
fCentMax(100.),
fCentEstimator(kCentOff),
fTriggerClass(""),
fTriggerMask(AliVEvent::kINT7),
fIsMC(false),
fSystem(0),
fESDtrackCuts(nullptr),
fAOD(nullptr),
fPIDresp(nullptr),
fV0cuts(nullptr),
fFillTreeWithNsigmaPIDOnly(false),
fEnabledDownSampling(false),
fFracToKeepDownSampling(0.1),
fPtMaxDownSampling(1.5),
fDownSamplingOpt(0),
fAODProtection(1),
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
fNEtabinsNsigmaTPCDataCorr(0)
{
  //
  // default constructur
  //

  for(int iVar=0; iVar<6; iVar++) fPIDNsigma[iVar] = -999.;
  for(int iHisto=0; iHisto<5; iHisto++) fHistArmenteroPlot[iHisto] = nullptr;
  for(int iHisto=0; iHisto<kNHypo; iHisto++) {
    fHistNsigmaTPCvsPt[iHisto] = nullptr;
    fHistNsigmaTOFvsPt[iHisto] = nullptr;
  }

  for(int iP=0; iP<AliAODPidHF::kMaxPBins; iP++) {
    for(int iEta=0; iEta<AliAODPidHF::kMaxEtaBins; iEta++) {
      fMeanNsigmaTPCPionData[iEta][iP] = 0.;
      fMeanNsigmaTPCKaonData[iEta][iP] = 0.;
      fMeanNsigmaTPCProtonData[iEta][iP] = 0.;
      fSigmaNsigmaTPCPionData[iEta][iP] = 1.;
      fSigmaNsigmaTPCKaonData[iEta][iP] = 1.;
      fSigmaNsigmaTPCProtonData[iEta][iP] = 1.;
    }
    fPlimitsNsigmaTPCDataCorr[iP] = 0.;
  }
  fPlimitsNsigmaTPCDataCorr[AliAODPidHF::kMaxPBins] = 0.;

  for(int iEta=0; iEta<=AliAODPidHF::kMaxEtaBins; iEta++) {
    fEtalimitsNsigmaTPCDataCorr[iEta] = 0.;
  }

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
fdEdxTPC(0),
fToF(0),
fPt(0),
fTPCNcls(0),
fTPCNclsPID(0),
fTrackLength(0),
fStartTimeRes(0),
fEta(-9999),
fPDGcode(-1),
fTag(0),
fNsigmaMaxForTag(0.02),
fQtMinKinks(0.15),
fRMinKinks(120),
fRMaxKinks(210),
fCentMin(0.),
fCentMax(100.),
fCentEstimator(kCentOff),
fTriggerClass(""),
fTriggerMask(AliVEvent::kINT7),
fIsMC(false),
fSystem(system),
fESDtrackCuts(nullptr),
fAOD(nullptr),
fPIDresp(nullptr),
fV0cuts(nullptr),
fFillTreeWithNsigmaPIDOnly(false),
fEnabledDownSampling(false),
fFracToKeepDownSampling(0.1),
fPtMaxDownSampling(1.5),
fDownSamplingOpt(0),
fAODProtection(1),
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
fNEtabinsNsigmaTPCDataCorr(0)
{
  //
  // standard constructur
  //

  for(int iVar=0; iVar<6; iVar++) fPIDNsigma[iVar] = -999.;
  for(int iHisto=0; iHisto<5; iHisto++) fHistArmenteroPlot[iHisto] = nullptr;
  for(int iHisto=0; iHisto<kNHypo; iHisto++) {
    fHistNsigmaTPCvsPt[iHisto] = nullptr;
    fHistNsigmaTOFvsPt[iHisto] = nullptr;
  }

  for(int iP=0; iP<AliAODPidHF::kMaxPBins; iP++) {
    for(int iEta=0; iEta<AliAODPidHF::kMaxEtaBins; iEta++) {
      fMeanNsigmaTPCPionData[iEta][iP] = 0.;
      fMeanNsigmaTPCKaonData[iEta][iP] = 0.;
      fMeanNsigmaTPCProtonData[iEta][iP] = 0.;
      fSigmaNsigmaTPCPionData[iEta][iP] = 1.;
      fSigmaNsigmaTPCKaonData[iEta][iP] = 1.;
      fSigmaNsigmaTPCProtonData[iEta][iP] = 1.;
    }
    fPlimitsNsigmaTPCDataCorr[iP] = 0.;
  }
  fPlimitsNsigmaTPCDataCorr[AliAODPidHF::kMaxPBins] = 0.;

  for(int iEta=0; iEta<=AliAODPidHF::kMaxEtaBins; iEta++) {
    fEtalimitsNsigmaTPCDataCorr[iEta] = 0.;
  }

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

  fHistNEvents = new TH1F("fHistNEvents","Number of processed events;;Number of events",8,-1.5,6.5);
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

  if(fIsMC) {
    for(int iHisto=0; iHisto<kNHypo; iHisto++) {
      fHistNsigmaTPCvsPt[iHisto] = new TH2F(Form("fHistNsigmaTPCvsPt_%s",hyponames[iHisto].Data()),Form(";#it{p}_{T} (GeV/#it{c});N_{#sigma}^{TPC}(%s)",hyponames[iHisto].Data()),500,0,50,1000,-50,50);
      fHistNsigmaTOFvsPt[iHisto] = new TH2F(Form("fHistNsigmaTOFvsPt_%s",hyponames[iHisto].Data()),Form(";#it{p}_{T} (GeV/#it{c});N_{#sigma}^{TOF}(%s)",hyponames[iHisto].Data()),500,0,50,1000,-50,50);
      fOutputList->Add(fHistNsigmaTPCvsPt[iHisto]);
      fOutputList->Add(fHistNsigmaTOFvsPt[iHisto]);
    }
  }

  fPIDtree = new TTree("fPIDtree","fPIDtree");
  TString PIDbranchnames[6] = {"n_sigma_TPC_pi","n_sigma_TPC_K","n_sigma_TPC_p","n_sigma_TOF_pi","n_sigma_TOF_K","n_sigma_TOF_p"};
  for(int iVar=0; iVar<6; iVar++) {
    fPIDtree->Branch(PIDbranchnames[iVar].Data(),&fPIDNsigma[iVar],Form("%s/S",PIDbranchnames[iVar].Data()));
  }
  fPIDtree->Branch("pT",&fPt,"pT/s");
  fPIDtree->Branch("pTPC",&fPTPC,"pTPC/s");
  fPIDtree->Branch("pTOF",&fPTOF,"pTOF/s");
  fPIDtree->Branch("eta",&fEta,"eta/S");
  if(!fFillTreeWithNsigmaPIDOnly) {
    fPIDtree->Branch("dEdx",&fdEdxTPC,"dEdx/s");
    fPIDtree->Branch("ToF",&fToF,"ToF/s");
    fPIDtree->Branch("NclusterTPC",&fTPCNcls,"NclusterTPC/b");
    fPIDtree->Branch("NclusterPIDTPC",&fTPCNclsPID,"NclusterPIDTPC/b");
    fPIDtree->Branch("TrackLength",&fTrackLength,"TrackLength/s");
    fPIDtree->Branch("StartTimeRes",&fStartTimeRes,"StartTimeRes/s");
  }
  fPIDtree->Branch("tag",&fTag,"tag/b");
  if(fIsMC) fPIDtree->Branch("PDGcode",&fPDGcode,"PDGcode/S");

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
  if(!fPIDresp) fPIDresp = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();

  //load data correction for NsigmaTPC, if enabled
  if(fEnableNsigmaTPCDataCorr) {
    if(fAOD->GetRunNumber()!=fRunNumberPrevEvent) {
      AliAODPidHF::SetNsigmaTPCDataDrivenCorrection(fAOD->GetRunNumber(), fSystNsigmaTPCDataCorr, fNPbinsNsigmaTPCDataCorr, fPlimitsNsigmaTPCDataCorr, fNEtabinsNsigmaTPCDataCorr, fEtalimitsNsigmaTPCDataCorr, fMeanNsigmaTPCPionData, fMeanNsigmaTPCKaonData, fMeanNsigmaTPCProtonData, fSigmaNsigmaTPCPionData, fSigmaNsigmaTPCKaonData, fSigmaNsigmaTPCProtonData);
    }
  }

  // V0 selection
  if(fSystem==0) fV0cuts->SetMode(AliAODv0KineCuts::kPurity,AliAODv0KineCuts::kPP);
  else if(fSystem==1) fV0cuts->SetMode(AliAODv0KineCuts::kPurity,AliAODv0KineCuts::kPbPb);
  fV0cuts->SetEvent(fAOD);

  vector<short> idPionFromK0s;
  vector<short> idPionFromL;
  vector<short> idProtonFromL;
  vector<short> idElectronFromGamma;
  vector<short> idKaonFromKinks;
  GetTaggedV0s(idPionFromK0s,idPionFromL,idProtonFromL,idElectronFromGamma);
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
    fPTPC = ConvertFloatToUnsignedShort(track->GetTPCmomentum()*1000);
    fPTOF = ConvertFloatToUnsignedShort(GetTOFmomentum(track)*1000);
    fEta = ConvertFloatToShort(track->Eta()*1000);

    if(!fFillTreeWithNsigmaPIDOnly) {
      //TPC variables
      fTPCNcls = static_cast<unsigned char>(track->GetTPCNcls());
      fdEdxTPC = ConvertFloatToUnsignedShort(track->GetTPCsignal()*100);
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
    }

    bool isTPCok = false;
    bool isTOFok = false;
    if (fPIDresp->CheckPIDStatus(AliPIDResponse::kTPC,track) == AliPIDResponse::kDetPidOk) isTPCok = true;
    if (fPIDresp->CheckPIDStatus(AliPIDResponse::kTOF,track) == AliPIDResponse::kDetPidOk) isTOFok = true;

    if(isTPCok) {
      float nSigmaTPCPion = fPIDresp->NumberOfSigmasTPC(track,AliPID::kPion);
      float nSigmaTPCKaon = fPIDresp->NumberOfSigmasTPC(track,AliPID::kKaon);
      float nSigmaTPCProton = fPIDresp->NumberOfSigmasTPC(track,AliPID::kProton);
      if(fEnableNsigmaTPCDataCorr) {
        float meanPion = 0., meanKaon = 0., meanProton = 0., sigmaPion = 1., sigmaKaon = 1., sigmaProton = 1.;
        GetNsigmaTPCMeanSigmaData(meanPion, sigmaPion, AliPID::kPion, track->GetTPCmomentum(), track->Eta());
        GetNsigmaTPCMeanSigmaData(meanKaon, sigmaKaon, AliPID::kKaon, track->GetTPCmomentum(), track->Eta());
        GetNsigmaTPCMeanSigmaData(meanProton, sigmaProton, AliPID::kProton, track->GetTPCmomentum(), track->Eta());

        if(nSigmaTPCPion>-990.)
          nSigmaTPCPion = (nSigmaTPCPion-meanPion) / sigmaPion;
        if(nSigmaTPCKaon>-990.)
          nSigmaTPCKaon = (nSigmaTPCKaon-meanKaon) / sigmaKaon;
        if(nSigmaTPCProton>-990.)
          nSigmaTPCProton = (nSigmaTPCProton-meanProton) / sigmaProton;
      }

      fPIDNsigma[0] = ConvertFloatToShort(nSigmaTPCPion*100);
      fPIDNsigma[1] = ConvertFloatToShort(nSigmaTPCKaon*100);
      fPIDNsigma[2] = ConvertFloatToShort(nSigmaTPCProton*100);
    }
    else for(int iVar=0; iVar<3; iVar++) fPIDNsigma[iVar] = numeric_limits<short>::min();
    if(isTOFok) {
      fPIDNsigma[3] = ConvertFloatToShort(fPIDresp->NumberOfSigmasTOF(track,AliPID::kPion)*100);
      fPIDNsigma[4] = ConvertFloatToShort(fPIDresp->NumberOfSigmasTOF(track,AliPID::kKaon)*100);
      fPIDNsigma[5] = ConvertFloatToShort(fPIDresp->NumberOfSigmasTOF(track,AliPID::kProton)*100);
    }
    else for(int iVar=3; iVar<6; iVar++) fPIDNsigma[iVar] = numeric_limits<short>::min();

    short trackid = track->GetID();

    it = find(idPionFromK0s.begin(),idPionFromK0s.end(),trackid);
    if(it!=idPionFromK0s.end()) fTag |= kIsPionFromK0s;
    else fTag &= ~kIsPionFromK0s;
    it = find(idPionFromL.begin(),idPionFromL.end(),trackid);
    if(it!=idPionFromL.end()) fTag |= kIsPionFromL;
    else fTag &= ~kIsPionFromL;
    it = find(idProtonFromL.begin(),idProtonFromL.end(),trackid);
    if(it!=idProtonFromL.end()) fTag |= kIsProtonFromL;
    else fTag &= ~kIsProtonFromL;
    it = find(idElectronFromGamma.begin(),idElectronFromGamma.end(),trackid);
    if(it!=idElectronFromGamma.end()) fTag |= kIsElectronFromGamma;
    else fTag &= ~kIsElectronFromGamma;
    it = find(idKaonFromKinks.begin(),idKaonFromKinks.end(),trackid);
    if(it!=idKaonFromKinks.end()) fTag |= kIsKaonFromKinks;
    else fTag &= ~kIsKaonFromKinks;
    if(TMath::Abs(fPIDresp->NumberOfSigmasTOF(track,AliPID::kKaon))<fNsigmaMaxForTag && isTOFok && isTPCok) fTag |= kIsKaonFromTOF;
    else fTag &= ~kIsKaonFromTOF;
    if(TMath::Abs(fPIDresp->NumberOfSigmasTPC(track,AliPID::kKaon))<fNsigmaMaxForTag && isTOFok && isTPCok) fTag |= kIsKaonFromTPC;
    else fTag &= ~kIsKaonFromTPC;
    
    if(fIsMC) {
      fPDGcode = GetPDGcodeFromMC(track,arrayMC);
      if(fPDGcode==211) {
        fHistNsigmaTPCvsPt[kPion]->Fill(track->Pt(),fPIDresp->NumberOfSigmasTPC(track,AliPID::kPion));
        fHistNsigmaTOFvsPt[kPion]->Fill(track->Pt(),fPIDresp->NumberOfSigmasTOF(track,AliPID::kPion));
      }
      else if(fPDGcode==321) {
        fHistNsigmaTPCvsPt[kKaon]->Fill(track->Pt(),fPIDresp->NumberOfSigmasTPC(track,AliPID::kKaon));
        fHistNsigmaTOFvsPt[kKaon]->Fill(track->Pt(),fPIDresp->NumberOfSigmasTOF(track,AliPID::kKaon));
      }
      else if(fPDGcode==2212) {
        fHistNsigmaTPCvsPt[kProton]->Fill(track->Pt(),fPIDresp->NumberOfSigmasTPC(track,AliPID::kProton));
        fHistNsigmaTOFvsPt[kProton]->Fill(track->Pt(),fPIDresp->NumberOfSigmasTOF(track,AliPID::kProton));
      }
    }
    else fPDGcode = 0;
    
    if(fTag!=0 && ((fIsMC && fPDGcode>=0) || !fIsMC)) fPIDtree->Fill();

    fTag = 0;

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
short AliAnalysisTaskSEHFSystPID::GetPDGcodeFromMC(AliAODTrack* track, TClonesArray* arrayMC)
{
  // Get pdg code
  short pdg = -1;
  int label = track->GetLabel();
  if(label<0) return pdg;
  AliAODMCParticle* partMC = dynamic_cast<AliAODMCParticle*>(arrayMC->At(label));
  if(!partMC) return pdg;
  pdg = TMath::Abs(partMC->GetPdgCode());
  if(partMC->GetPdgCode()>numeric_limits<short>::max()) pdg = -1;

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
 
  return static_cast<short>(num);
}

//___________________________________________________________________
unsigned short AliAnalysisTaskSEHFSystPID::ConvertFloatToUnsignedShort(float num) {
  if(num>=static_cast<float>(numeric_limits<unsigned short>::max())) return numeric_limits<unsigned short>::max();
  else if(num<=static_cast<float>(numeric_limits<unsigned short>::min())) return numeric_limits<unsigned short>::min();
 
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
