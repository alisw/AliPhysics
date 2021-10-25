//Class for <pt>-v2 correlations
//Author: Vytautas Vislavicius


#include "AliAnalysisTaskMeanPtV2Corr.h"
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"
#include "TList.h"
#include "TProfile.h"
#include "AliEventCuts.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "AliStack.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "AliGFWWeights.h"
#include "AliGFWFlowContainer.h"
#include "AliGFW.h"
#include "TClonesArray.h"
#include "AliGFWCuts.h"
#include "AliAODMCParticle.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TNamed.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODMCHeader.h"
#include "AliEffFDContainer.h"

ClassImp(AliAnalysisTaskMeanPtV2Corr);

AliAnalysisTaskMeanPtV2Corr::AliAnalysisTaskMeanPtV2Corr():
  AliAnalysisTaskSE(),
  fStageSwitch(0),
  fSystFlag(0),
  fEventCutFlag(0),
  fContSubfix(0),
  fCentEst(0),
  fExtendV0MAcceptance(kTRUE),
  fIsMC(kFALSE),
  fBypassTriggerAndEvetCuts(kFALSE),
  fMCEvent(0),
  fUseRecoNchForMC(kTRUE),
  fRndm(0),
  fNBootstrapProfiles(10),
  fPtAxis(0),
  fMultiAxis(0),
  fV0MMultiAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fUseNch(kFALSE),
  fUseWeightsOne(kFALSE),
  fEta(0.8),
  fEtaLow(-9999),
  fEtaNch(0.8),
  fEtaV2Sep(0.4),
  fPIDResponse(0),
  fBayesPID(0),
  fQAList(0),
  fMultiDist(0),
  fMultiVsV0MCorr(0),
  fNchTrueVsReco(0),
  fESDvsFB128(0),
  fNchVsMulti(0),
  fNchInBins(0),
  fptVarList(0),
  fCkCont(0),
  fCovList(0),
  fV2dPtList(0),
  fCovariance(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fWeightList(0),
  fWeights(0),
  fRunNo(0),
  fGFWSelection(0),
  fGFWNtotSelection(0),
  fFC(0),
  fGFW(0),
  fEfficiencyList(0),
  fEfficiency(0),
  fEfficiencies(0),
  fPseudoEfficiency(2.),
  fDCAxyVsPt_noChi2(0),
  fWithinDCAvsPt_withChi2(0),
  fDCAxyVsPt_withChi2(0),
  fWithinDCAvsPt_noChi2(0),
  fV0MMulti(0),
  fITSvsTPCMulti(0),
  fV2dPtMulti(0),
  fSPDCutPU(0),
  fV0CutPU(0),
  fCenCutLowPU(0),
  fCenCutHighPU(0),
  fMultCutPU(0),
  fStdTPCITS2011(0),
  fDisablePID(kTRUE),
  fConsistencyFlag(3),
  fRequireReloadOnRunChange(kFALSE),
  fEfFd(0)
{
};
AliAnalysisTaskMeanPtV2Corr::AliAnalysisTaskMeanPtV2Corr(const char *name, Bool_t IsMC, TString stageSwitch, TString ContSubfix):
  AliAnalysisTaskSE(name),
  fStageSwitch(0),
  fSystFlag(0),
  fEventCutFlag(0),
  fContSubfix(0),
  fCentEst(0),
  fExtendV0MAcceptance(kTRUE),
  fIsMC(IsMC),
  fBypassTriggerAndEvetCuts(kFALSE),
  fMCEvent(0),
  fUseRecoNchForMC(kTRUE),
  fNBootstrapProfiles(10),
  fRndm(0),
  fPtAxis(0),
  fMultiAxis(0),
  fV0MMultiAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fUseNch(kFALSE),
  fUseWeightsOne(kFALSE),
  fEta(0.8),
  fEtaLow(-9999),
  fEtaNch(0.8),
  fEtaV2Sep(0.4),
  fPIDResponse(0),
  fBayesPID(0),
  fQAList(0),
  fMultiDist(0),
  fMultiVsV0MCorr(0),
  fNchTrueVsReco(0),
  fESDvsFB128(0),
  fNchVsMulti(0),
  fNchInBins(0),
  fptVarList(0),
  fCkCont(0),
  fCovList(0),
  fV2dPtList(0),
  fCovariance(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fWeightList(0),
  fWeights(0),
  fRunNo(0),
  fGFWSelection(0),
  fGFWNtotSelection(0),
  fFC(0),
  fGFW(0),
  fEfficiencyList(0),
  fEfficiency(0),
  fEfficiencies(0),
  fPseudoEfficiency(2.),
  fDCAxyVsPt_noChi2(0),
  fWithinDCAvsPt_withChi2(0),
  fDCAxyVsPt_withChi2(0),
  fWithinDCAvsPt_noChi2(0),
  fV0MMulti(0),
  fITSvsTPCMulti(0),
  fV2dPtMulti(0),
  fSPDCutPU(0),
  fV0CutPU(0),
  fCenCutLowPU(0),
  fCenCutHighPU(0),
  fMultCutPU(0),
  fStdTPCITS2011(0),
  fDisablePID(kTRUE),
  fConsistencyFlag(3),
  fRequireReloadOnRunChange(kFALSE),
  fEfFd(0)
{
  fStageSwitch = GetStageSwitch(stageSwitch);
  SetContSubfix(ContSubfix);
  fCentEst = new TString("V0M");
  if(!fStageSwitch) AliFatal("Stage switch is 0, not sure what should be done!\n");
  if(fStageSwitch==1)
    DefineOutput(1,TList::Class());
  if(fStageSwitch==7) {
    DefineOutput(1,TList::Class());
  }
  if(fStageSwitch==9) {
    if(!fIsMC) { //Efficiency and NUA only important for data
      DefineInput(1,TList::Class()); //NUE weights; ultimately, should be combined with NUA, but don't want to rerun now
      DefineInput(2,TList::Class()); //NUA weights from other analysis; quickfix
    };
    DefineOutput(1,TList::Class());
    DefineOutput(2,AliGFWFlowContainer::Class());
    DefineOutput(3,TList::Class());
    DefineOutput(4,TList::Class());
  };
  if(fStageSwitch==10) {
    DefineOutput(1,AliEffFDContainer::Class());//AliEffFDContainer::Class());
  }
  SetNchCorrelationCut(1,0,kFALSE);
};
AliAnalysisTaskMeanPtV2Corr::~AliAnalysisTaskMeanPtV2Corr() {
  SetNchCorrelationCut(1,0,kFALSE);
};
void AliAnalysisTaskMeanPtV2Corr::UserCreateOutputObjects(){
  printf("Stage switch is %i\n\n\n",fStageSwitch);
  if(!fGFWSelection) SetSystFlag(0);
  fGFWSelection->PrintSetup();
  fSystFlag = fGFWSelection->GetSystFlagIndex();
  if(fGFWSelection->GetSystFlagIndex() == 20) SetCentralityEstimator("CL0");
  else if(fGFWSelection->GetSystFlagIndex() == 21) SetCentralityEstimator("CL1");
  OpenFile(1);
  const Int_t temp_NV0MBinsDefault=fExtendV0MAcceptance?11:10;
  Double_t temp_V0MBinsDefault[12] = {0,5,10,20,30,40,50,60,70,80,90,101}; //Last bin to include V0M beyond anchor point
  if(!fV0MMultiAxis) SetV0MBins(temp_NV0MBinsDefault,temp_V0MBinsDefault);
  Double_t *l_V0MBinsDefault=GetBinsFromAxis(fV0MMultiAxis);
  Int_t l_NV0MBinsDefault=fV0MMultiAxis->GetNbins();
  if(l_V0MBinsDefault[l_NV0MBinsDefault]>90) fExtendV0MAcceptance = kTRUE; //If V0M is beyond 90, then we need to extend the V0M acceptance!
  if(!fMultiAxis) SetMultiBins(l_NV0MBinsDefault,l_V0MBinsDefault);
  fMultiBins = GetBinsFromAxis(fMultiAxis);
  fNMultiBins = fMultiAxis->GetNbins();
  if(!fV2dPtMulti) {
    Double_t temp_bn[] = {0,1e6};
    SetV2dPtMultiBins(1,temp_bn);
  };
  const Int_t l_NPtBinsDefault = 25;
  Double_t l_PtBinsDefault[l_NPtBinsDefault+1] = {0.20, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
                     1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90,
                     2.00, 2.20, 2.40, 2.60, 2.80, 3.00};
  if(!fPtAxis) SetPtBins(l_NPtBinsDefault,l_PtBinsDefault);
  fPtBins = GetBinsFromAxis(fPtAxis);
  fNPtBins = fPtAxis->GetNbins();
  TString spNames[] = {"ch","pi","ka","pr"};
  if(fStageSwitch==1) {
    fRequireReloadOnRunChange = kFALSE;
    fWeightList = new TList();
    fWeightList->SetOwner(kTRUE);
      TString wNames[] = {"ch","pi","ka","pr"};
      fWeights = new AliGFWWeights*[4];
      for(Int_t i=0; i<4;i++) {
        fWeights[i] = new AliGFWWeights();
        fWeights[i]->SetPtBins(fNPtBins,fPtBins);
        // fWeights[i]->SetPtBins(NbinsPtForV2,binsPtForV2);
        fWeights[i]->SetName(Form("weight_%s",wNames[i].Data()));
        fWeights[i]->Init(!fIsMC,fIsMC);
        fWeightList->Add(fWeights[i]);
      }
      PostData(1,fWeightList);
  };
  if(fStageSwitch==7) {
    fRequireReloadOnRunChange = kFALSE;
    fEfficiencyList = new TList();
    fEfficiencyList->SetOwner(kTRUE);
    if(fIsMC) {
      Int_t nSpecies=fDisablePID?1:7; //Disabling PID part to save some output for now
      fEfficiency = new TH2D*[nSpecies*3];
      TString lNames[] = {"ch","pi","ka","pr","sigma","xi","om"};
      for(Int_t i=0;i<nSpecies;i++) {
        lNames[i].Prepend("Spectra_");
        fEfficiency[i] = new TH2D(lNames[i].Data(),lNames[i].Data(),fNPtBins,fPtBins,l_NV0MBinsDefault,l_V0MBinsDefault);
        lNames[i].Append("_Gen");
        fEfficiency[nSpecies+i] = new TH2D(lNames[i].Data(),lNames[i].Data(),fNPtBins,fPtBins,l_NV0MBinsDefault,l_V0MBinsDefault);
        lNames[i].Append("_Sec");
        fEfficiency[2*nSpecies+i] = new TH2D(lNames[i].Data(),lNames[i].Data(),fNPtBins,fPtBins,l_NV0MBinsDefault,l_V0MBinsDefault);
        fEfficiencyList->Add(fEfficiency[i]);
        fEfficiencyList->Add(fEfficiency[i+nSpecies]);
        fEfficiencyList->Add(fEfficiency[i+(2*nSpecies)]);
      }
      fNchTrueVsReco = new TH2D("NchTrueVsReco",";Nch (MC-true); Nch (MC-reco)",fNMultiBins,fMultiBins,fNMultiBins,fMultiBins);
      fEfficiencyList->Add(fNchTrueVsReco);
    }
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",l_NV0MBinsDefault,l_V0MBinsDefault);
    fEfficiencyList->Add(fV0MMulti);
    fStdTPCITS2011 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    //DCA histograms:
    Double_t binsDCA[61] = {-3.00, -2.90, -2.80, -2.70, -2.60, -2.50, -2.40, -2.30, -2.20, -2.10, -2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10, -1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00};
    Int_t NbinsDCA = 60;
    Double_t binsType[4] = {-0.5,0.5,1.5,2.5};
    Int_t NbinsType = 3;
    fDCAxyVsPt_noChi2 = new TH3D("DCAxy_noChi2","DCAxy_noChi2; pt; type; dcaxy", fNPtBins, fPtBins, NbinsType, binsType, NbinsDCA, binsDCA);
    fDCAxyVsPt_withChi2 = new TH3D("DCAxy_withChi2","DCAxy_withChi2; pt; type; dcaxy", fNPtBins, fPtBins, NbinsType, binsType, NbinsDCA, binsDCA);
    fWithinDCAvsPt_noChi2 = new TH2D("WithinDCA_noChi2","WithinDCA_noChi2; pt; type", fNPtBins, fPtBins, NbinsType, binsType);
    fWithinDCAvsPt_withChi2 = new TH2D("WithinDCA_withChi2","WithinDCA_withChi2; pt; type", fNPtBins, fPtBins, NbinsType, binsType);
    fEfficiencyList->Add(fDCAxyVsPt_noChi2);
    fEfficiencyList->Add(fWithinDCAvsPt_noChi2);
    fEfficiencyList->Add(fDCAxyVsPt_withChi2);
    fEfficiencyList->Add(fWithinDCAvsPt_withChi2);
    //ITS vs TPC tracklets cut for PU
    fITSvsTPCMulti = new TH2D("TPCvsITSclusters",";TPC clusters; ITS clusters",1000,0,10000,5000,0,50000);
    fEfficiencyList->Add(fITSvsTPCMulti);
    PostData(1,fEfficiencyList);
  }
  if(fStageSwitch==9) {
    fRndm = new TRandom(0);
    fRequireReloadOnRunChange = kFALSE;
    if(!fIsMC) { //Efficiencies and NUA are only for the data
      fEfficiencyList = (TList*)GetInputData(1);
      fEfficiencies = new TH1D*[l_NV0MBinsDefault];
      for(Int_t i=0;i<l_NV0MBinsDefault;i++) {
        fEfficiencies[i] = (TH1D*)fEfficiencyList->FindObject(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
        if(fEfficiencies[i] && fPseudoEfficiency<1) fEfficiencies[i]->Scale(fPseudoEfficiency);
        if(!fEfficiencies[i]) {
          if(!i) AliFatal("Could not fetch efficiency!\n");
          printf("Could not find efficiency for V0M bin no. %i! Cloning the previous efficiency instead...\n",i);
          fEfficiencies[i] = (TH1D*)fEfficiencies[i-1]->Clone(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
        };
      }
      fWeightList = (TList*)GetInputData(2);
      fWeights = new AliGFWWeights*[1];
    };
    // if(!LoadMyWeights(0)) return; //Loading run-avg NUA weights
    fptVarList = new TList();
    fptVarList->SetOwner(kTRUE);
    for(Int_t i=0;i<1;i++) {
      fCkCont = new AliCkContainer(Form("ckcont_%s",spNames[i].Data()),Form("ckcont_%s",spNames[i].Data()),fNMultiBins,fMultiBins);
      fptVarList->Add(fCkCont);
      if(fNBootstrapProfiles) fCkCont->InitializeSubsamples(fNBootstrapProfiles);
    }
    fMultiDist = new TH1D("MultiDistribution","Multiplicity distribution; #it{N}_{ch}; N(events)",fNMultiBins,fMultiBins);
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",l_NV0MBinsDefault,l_V0MBinsDefault);
    fptVarList->Add(fMultiDist);
    fptVarList->Add(fV0MMulti);
    fMultiVsV0MCorr = new TH2D*[2];
    fMultiVsV0MCorr[0] = new TH2D("MultVsV0M_BeforeConsistency","MultVsV0M_BeforeConsistency",103,0,103,fNMultiBins,fMultiBins[0],fMultiBins[fNMultiBins]);
    fMultiVsV0MCorr[1] = new TH2D("MultVsV0M_AfterConsistency","MultVsV0M_AfterConsistency",103,0,103,fNMultiBins,fMultiBins[0],fMultiBins[fNMultiBins]);
    fESDvsFB128 = new TH2D("ESDvsFB128","; N(FB128); N(ESD)",500,-0.5,4999.5,1500,-0.5,14999.5);
    fptVarList->Add(fMultiVsV0MCorr[0]);
    fptVarList->Add(fMultiVsV0MCorr[1]);
    fptVarList->Add(fESDvsFB128);
    //ITS vs TPC tracklets cut for PU
    fITSvsTPCMulti = new TH2D("TPCvsITSclusters",";TPC clusters; ITS clusters",1000,0,10000,5000,0,50000);
    fptVarList->Add(fITSvsTPCMulti);
    if(fIsMC) {
      fNchTrueVsReco = new TH2D("NchTrueVsReco",";Nch (MC-true); Nch (MC-reco)",fNMultiBins,fMultiBins,fNMultiBins,fMultiBins);
      fptVarList->Add(fNchTrueVsReco);
    }
    PostData(1,fptVarList);
    //Setting up the FlowContainer
    TObjArray *oba = new TObjArray();
    oba->Add(new TNamed("ChGap22","ChGap22")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24","ChGap24")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull22","ChFull22")); //no-gap case
    oba->Add(new TNamed("ChFull24","ChFull24")); //no-gap case
//adding v3n
    oba->Add(new TNamed("ChGap32","ChGap32")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap34","ChGap34")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull32","ChFull32")); //no-gap case
    oba->Add(new TNamed("ChFull34","ChFull34")); //no-gap case

    oba->Add(new TNamed("ChGap42","ChGap42")); //gap case

    oba->Add(new TNamed("LM22","LM22")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("MR22","MR22")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("LR22","LR22")); //for gap (|eta|>0.4) case

    oba->Add(new TNamed("LLMR24","LLMR24")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("LMMR24","LMMR24")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("LMRR24","LMRR24")); //for gap (|eta|>0.4) case

    oba->Add(new TNamed("ChSC234","ChSC234")); //for SC{2,3}
    oba->Add(new TNamed("ChSC244","ChSC244")); //for SC{2,3}

    fFC = new AliGFWFlowContainer();
    TString fcname("FlowContainer");
    if(!fContSubfix->IsNull()) fcname.Append(fContSubfix->Data());
    fFC->SetName(fcname.Data());
    fFC->Initialize(oba,fNMultiBins,fMultiBins,fNBootstrapProfiles);
    delete oba;
    PostData(2,fFC);
    Int_t pows[] = {3,0,2,2,3,3,3}; //5th harm. sum = 3, b/c {-2 -3}
    Int_t powsFull[] = {5,0,4,4,3,3,3};
    fGFW = new AliGFW();
    fGFW->AddRegion("refN",7,pows,-0.8,-fEtaV2Sep,1,1);
    fGFW->AddRegion("refP",7,pows,fEtaV2Sep,0.8,1,1);
    if(fEtaV2Sep>=0)
      fGFW->AddRegion("subMid",7,pows,-fEtaV2Sep,fEtaV2Sep,1,1);
    fGFW->AddRegion("mid",7,powsFull,-0.8,0.8,1,2);
    CreateCorrConfigs();
    //Covariance
    fCovList = new TList();
    fCovList->SetOwner(kTRUE);
    fCovariance = new AliProfileBS*[6];
    for(Int_t i=0;i<1;i++) {
      fCovList->Add(new AliProfileBS(Form("covmpt_%s",spNames[i].Data()),Form("covmpt_%s",spNames[i].Data()),fNMultiBins,fMultiBins));
      fCovariance[i] = (AliProfileBS*)fCovList->At(i);
      fCovList->Add(new AliProfileBS(Form("covnopt_%s",spNames[i].Data()),Form("covnopt_%s",spNames[i].Data()),fNMultiBins,fMultiBins));
      fCovariance[i+1] = (AliProfileBS*)fCovList->At(i+1);
      fCovList->Add(new AliProfileBS(Form("covmpt_v3_%s",spNames[i].Data()),Form("covmpt_v3_%s",spNames[i].Data()),fNMultiBins,fMultiBins));
      fCovariance[i+2] = (AliProfileBS*)fCovList->At(i+2);
      fCovList->Add(new AliProfileBS(Form("covnopt_v3_%s",spNames[i].Data()),Form("covnopt_v3_%s",spNames[i].Data()),fNMultiBins,fMultiBins));
      fCovariance[i+3] = (AliProfileBS*)fCovList->At(i+3);
      fCovList->Add(new AliProfileBS(Form("covmpt_v23_%s",spNames[0].Data()),Form("covmpt_v23_%s",spNames[0].Data()),fNMultiBins,fMultiBins));
      fCovariance[i+4] = (AliProfileBS*)fCovList->At(4);
      fCovList->Add(new AliProfileBS(Form("covnopt_v23_%s",spNames[0].Data()),Form("covnopt_v23_%s",spNames[0].Data()),fNMultiBins,fMultiBins));
      fCovariance[i+5] = (AliProfileBS*)fCovList->At(5);

    };
    if(fNBootstrapProfiles) for(Int_t i=0;i<6;i++) fCovariance[i]->InitializeSubsamples(fNBootstrapProfiles);
    PostData(3,fCovList);
    fQAList = new TList();
    fQAList->SetOwner(kTRUE);
    fEventCuts.AddQAplotsToList(fQAList,kTRUE);
    PostData(4,fQAList);
  }
  if(fStageSwitch==10) {
    fEfFd = new AliEffFDContainer("EffAndFD","EffAndFD",fIsMC);
    fEfFd->SetCentralityBins(l_NV0MBinsDefault,l_V0MBinsDefault);
    fEfFd->SetPtBins(fNPtBins,fPtBins);
    fEfFd->SetEta(fEta);
    AliESDtrackCuts *tc = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    fEfFd->AddCut(tc);
    PostData(1,fEfFd);
  };
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerType,true);
  if(fExtendV0MAcceptance) {
    fEventCuts.OverrideCentralityFramework(1);
    fEventCuts.SetCentralityEstimators("V0M","CL0");
    fEventCuts.SetCentralityRange(0.f,101.f);
  }
  //Creating cuts for 15o_pass2 and 18qr_pass3. 18qr_pass3 not implemented yet.
  //Would like to do that in a more elegant way, but not at this point, unfortunatelly
  if(fEventCutFlag) { //Only initialize them if necessary
    fSPDCutPU = new TF1("fSPDCutPU", "450. + 3.9*x", 0, 50000);
    if(!fV0CutPU) fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000); //Only if not initialized externally. Set to 0 for ESD MC, as that seems to be problematic?
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    if(fEventCutFlag==1 || fEventCutFlag==101) {
       Double_t parV0[8] = {33.4237, 0.953516, 0.0712137, 227.923, 8.9239, -0.00319679, 0.000306314, -7.6627e-07};
       fV0CutPU->SetParameters(parV0);
       Double_t parV0CL0[6] = {0.0193587, 0.975914, 0.675714, 0.0292263, -0.000549509, 5.86421e-06};
       fCenCutLowPU->SetParameters(parV0CL0);
       fCenCutHighPU->SetParameters(parV0CL0);
       Double_t parFB32[9] = {-812.822, 6.41796, 5421.83, -0.382601, 0.0299686, -26.6249, 321.388, -0.82615, 0.0167828};
       fMultCutPU->SetParameters(parFB32);
    }
  };
  fGFWNtotSelection = new AliGFWCuts();
  fGFWNtotSelection->SetupCuts(0);
  fGFWNtotSelection->SetEta(fEtaNch);
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fBayesPID = new AliPIDCombined();
  fBayesPID->SetDefaultTPCPriors();
  fBayesPID->SetSelectedSpecies(AliPID::kSPECIES);
  fBayesPID->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  LoadWeightAndMPT();
};
void AliAnalysisTaskMeanPtV2Corr::UserExec(Option_t*) {
  if(fStageSwitch == 7) { //Efficiencies on ESDs
    AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fESD) return;
    //Checking multiplicity & trigger
    AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
    if(!lMultSel) AliFatal("Mult selection not found!\n");
    Double_t l_Cent = lMultSel->GetMultiplicityPercentile(fCentEst->Data());
    if(!fBypassTriggerAndEvetCuts)
      if(!CheckTrigger(l_Cent)) return;
    if(fEventCutFlag) if(!AcceptCustomEvent(fESD)) return;
    Bool_t dummy = fEventCuts.AcceptEvent(fESD); //for QA
    if(!fEventCutFlag && !dummy) return;
    if(!fGFWSelection->AcceptVertex(fESD)) return;
    if(fIsMC) {
      fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
      if (!fMCEvent) return;
    }
    Double_t vtxXYZ[] = {0,0,0};
    Double_t vzt=0;
    FillTPCITSClusters(fESD);
    ProduceEfficiencies(fESD,vzt,l_Cent,vtxXYZ); //Disabling vz vertex
    return;
  };
  if(fStageSwitch == 10) { //Efficiencies on ESDs
    AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fESD) return;
    //Checking multiplicity & trigger
    AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
    if(!lMultSel) AliFatal("Mult selection not found!\n");
    Double_t l_Cent = lMultSel->GetMultiplicityPercentile(fCentEst->Data());
    if(!fBypassTriggerAndEvetCuts)
      if(!CheckTrigger(l_Cent)) return;
    if(fEventCutFlag) if(!AcceptCustomEvent(fESD)) return;
    Bool_t dummy = fEventCuts.AcceptEvent(fESD); //for QA
    if(!fEventCutFlag && !dummy) return;
    if(!fGFWSelection->AcceptVertex(fESD)) return;
    if(fIsMC) {
      fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
      if (!fMCEvent) return;
    }
    if(fIsMC) fEfFd->Fill(*fESD, *fMCEvent);
    else fEfFd->Fill(*fESD);
    PostData(1,fEfFd);
    return;
  }

  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCEvent) return;
  }
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t l_Cent = lMultSel->GetMultiplicityPercentile(fCentEst->Data());
  if(!fBypassTriggerAndEvetCuts)
    if(!CheckTrigger(l_Cent)) return;
  Double_t vtxXYZ[] = {0.,0.,0.};
  if(!AcceptAOD(fAOD, vtxXYZ)) return;
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  FillTPCITSClusters(fAOD);
  if(!fGFWSelection->AcceptVertex(fAOD)) return;
  if(fStageSwitch==1)
    fIsMC?FillWeightsMC(fAOD, vz,l_Cent,vtxXYZ):FillWeights(fAOD, vz,l_Cent,vtxXYZ);
  if(fStageSwitch==9)
    CovSkipMpt(fAOD,vz,l_Cent,vtxXYZ);
};
void AliAnalysisTaskMeanPtV2Corr::NotifyRun() {
  if(!fEventCutFlag || fEventCutFlag>100) { //Only relevant if we're using the standard AliEventCuts
    //Reinitialize AliEventCuts (done automatically on check):
    Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);

    //Then override PU cut if required:
    if(fGFWSelection->GetSystFlagIndex()==22)
      fEventCuts.fESDvsTPConlyLinearCut[0] = 1500.;
  };
}
void AliAnalysisTaskMeanPtV2Corr::Terminate(Option_t*) {
};
Bool_t AliAnalysisTaskMeanPtV2Corr::CheckTrigger(Double_t lCent) {
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  //Apparently, MB trigger can also mark special triggers, leaving depleted regions in multi. To avoid this, pass true, if MB has been triggered.
  //This would fail if spec. triggers would also flag MB trigger, which seems to NOT be the case.
  if(!(fTriggerType&fSelMask)) { return kFALSE; }; //printf("Returning from the generic check\n");
  if(fSelMask&(fTriggerType&(AliVEvent::kINT7+AliVEvent::kMB))) {return kTRUE; }; //printf("Passed by MB trigger!\n");
  if((fSelMask&fTriggerType&AliVEvent::kCentral) && lCent>10) {return kFALSE; }; //printf("Returnning from kCent case\n");
  if((fSelMask&fTriggerType&AliVEvent::kSemiCentral) && (lCent<30 || lCent>50)) {return kFALSE; }; //printf("Returning from kSC case\n");
  return kTRUE;
};
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptAOD(AliAODEvent *inEv, Double_t *lvtxXYZ) {
  if(!fBypassTriggerAndEvetCuts) {
    if(!fEventCutFlag) { if(!fEventCuts.AcceptEvent(inEv)) return 0; } //Don't perform AcceptEvent if not relevant
    else if(!AcceptCustomEvent(inEv)) return 0;
    if(fEventCutFlag>100) Bool_t dummy = fEventCuts.AcceptEvent(inEv); //if flag > 100, then also store QA output from AcceptEvent
  };
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertexSPD());
  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;
  const Double_t aodVtxZ = vtx->GetZ();
  if(TMath::Abs(aodVtxZ) > 10)
    return kFALSE;
  vtx->GetXYZ(lvtxXYZ);
  return kTRUE;
};
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp) {
  if(mtr->Pt()<ptMin) return kFALSE;
  if(mtr->Pt()>ptMax) return kFALSE;
  if(ltrackXYZ && vtxp) {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2]-vtxp[2];
  } else return kFALSE; //DCA cut is a must for now
  return fGFWSelection->AcceptTrack(mtr,fSystFlag==1?0:ltrackXYZ,0,kFALSE);//All complementary DCA track cuts for FB768 are disabled
};
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptESDTrack(AliESDtrack *mtr, UInt_t& primFlag, Double_t *ltrackXYZ, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp) {
  if(mtr->Pt()<ptMin) return kFALSE;
  if(mtr->Pt()>ptMax) return kFALSE;
  if(ltrackXYZ) {
    Float_t fD, fZ;
    mtr->GetImpactParameters(fD,fZ);
    ltrackXYZ[0] = fD;
    ltrackXYZ[1] = fZ;
  } else return kFALSE; //DCA cut is a must for now
  return fGFWSelection->AcceptTrack(mtr,fSystFlag==1?0:ltrackXYZ,0,primFlag);
};

Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp, Int_t &nTot) {
  if(mtr->Pt()<ptMin) return kFALSE;
  if(mtr->Pt()>ptMax) return kFALSE;
  if(ltrackXYZ && vtxp) {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2]-vtxp[2];
  } else return kFALSE; //DCA cut is a must for now
  if(fGFWNtotSelection->AcceptTrack(mtr,ltrackXYZ,0,kFALSE)) nTot++;
  return fGFWSelection->AcceptTrack(mtr,fSystFlag==1?0:ltrackXYZ,0,kFALSE); //All complementary DCA track cuts for FB768 are disabled
};
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptESDTrack(AliESDtrack *mtr, UInt_t& primFlag, Double_t *ltrackXYZ, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp, Int_t &nTot) {
  if(mtr->Pt()<ptMin) return kFALSE;
  if(mtr->Pt()>ptMax) return kFALSE;
  if(ltrackXYZ) {
    Float_t fD, fZ;
    mtr->GetImpactParameters(fD,fZ);
    ltrackXYZ[0] = fD;
    ltrackXYZ[1] = fZ;
  } else return kFALSE; //DCA cut is a must for now
  UInt_t dummy;
  if(fGFWNtotSelection->AcceptTrack(mtr,ltrackXYZ,0,dummy)) nTot++;
  return fGFWSelection->AcceptTrack(mtr,fSystFlag==1?0:ltrackXYZ,0,primFlag); //All complementary DCA track cuts for FB768 are disabled
};

Int_t AliAnalysisTaskMeanPtV2Corr::GetStageSwitch(TString instr) {
  if(instr.Contains("weights")) return 1;
  if(instr.Contains("Efficiency")) return 7;
  if(instr.Contains("CovSkipMpt")) return 9;
  if(instr.Contains("EfTest")) return 10;
  return 0;
}
void AliAnalysisTaskMeanPtV2Corr::FillWeightsMC(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp) {
  //MC generated
  AliAODTrack *lTrack;
  // AliVParticle *lPart;
  Double_t trackXYZ[3];
  Double_t dummyDouble[] = {0.,0.};
  Double_t ptMin = fPtBins[0];
  Double_t ptMax = fPtBins[fNPtBins];
  TClonesArray *tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
  Int_t nPrim = tca->GetEntries();
  AliVParticle *lPart;
  Int_t partNotFetched=0;
  for (Int_t ipart = 0; ipart < nPrim; ipart++) {
    lPart = (AliAODMCParticle*)tca->At(ipart);
    if (!lPart) { continue; };
    /* get particlePDG */
    Int_t pdgcode = TMath::Abs(lPart->PdgCode());
    if (!lPart->IsPhysicalPrimary()) continue;
    if (lPart->Charge()==0.) continue;
    if (TMath::Abs(lPart->Eta()) > fEta) continue;
    Double_t pt = lPart->Pt();
    if (pt<ptMin || pt>ptMax) continue;
    fWeights[0]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    if(fDisablePID) continue;
    if(pdgcode==211) fWeights[1]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    if(pdgcode==321) fWeights[2]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    if(pdgcode==2212) fWeights[3]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
  };

  //MC reconstructed
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!AcceptAODTrack(lTrack,trackXYZ,ptMin,ptMax,vtxp)) continue;
    lPart = (AliAODMCParticle*)tca->At(TMath::Abs(lTrack->GetLabel()));
    if(!lPart) continue;
    if(!lPart->IsPhysicalPrimary()) continue;
    if(TMath::Abs(lTrack->Eta())>fEta) continue;
    if(!fGFWSelection->AcceptTrack(lTrack,dummyDouble)) continue;
    fWeights[0]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
    if(fDisablePID) continue;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    if(PIDIndex) fWeights[PIDIndex]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
  };
  PostData(1,fWeightList);
}
void AliAnalysisTaskMeanPtV2Corr::FillWeights(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp) {
  AliAODTrack *lTrack;
  Double_t trackXYZ[3];
  Double_t ptMin = fPtBins[0];
  Double_t ptMax = fPtBins[fNPtBins];
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!AcceptAODTrack(lTrack,trackXYZ,ptMin,ptMax,vtxp)) continue;
    Double_t leta = lTrack->Eta();
    Double_t lpt = lTrack->Pt();
    ((AliGFWWeights*)fWeightList->At(0))->Fill(lTrack->Phi(),lTrack->Eta(),vz,lTrack->Pt(),l_Cent,0);
    if(fDisablePID) continue;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    if(PIDIndex) ((AliGFWWeights*)fWeightList->At(PIDIndex))->Fill(lTrack->Phi(),lTrack->Eta(),vz,lTrack->Pt(),l_Cent,0);

  };
  PostData(1,fWeightList);
}
Int_t AliAnalysisTaskMeanPtV2Corr::GetNtotTracks(AliAODEvent* lAOD, const Double_t &ptmin, const Double_t &ptmax, Double_t *vtxp) {
  Double_t ltrackXYZ[3];
  AliAODTrack *lTrack;
  Int_t nTotNoTracks=0;
  for(Int_t lTr=0;lTr<lAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)lAOD->GetTrack(lTr);
    if(!lTrack) continue;
    if(!AcceptAODTrack(lTrack,ltrackXYZ,ptmin,ptmax,vtxp,nTotNoTracks)) continue;
  };
  return nTotNoTracks;
}
void AliAnalysisTaskMeanPtV2Corr::FillWPCounter(Double_t inArr[5], Double_t w, Double_t p) {
  inArr[0] += w;       // = w1p0
  inArr[1] += w*p;     // = w1p1
  inArr[2] += w*w*p*p; // = w2p2
  inArr[3] += w*w*p;   // = w2p1
  inArr[4] += w*w;     // = w2p0
}

void AliAnalysisTaskMeanPtV2Corr::ProduceEfficiencies(AliESDEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp) {
  AliESDtrack *lTrack;
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  Double_t ptMin = fPtBins[0];
  Double_t ptMax = fPtBins[fNPtBins];
  Float_t dcaxy, dcaz;
  fV0MMulti->Fill(l_Cent);
  if(!fIsMC) { //If running on data, only the following
    for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
      lTrack = (AliESDtrack*)fAOD->GetTrack(lTr);
      if(!lTrack) continue;
      Double_t lpt = lTrack->Pt();
      if(lpt<0.1) continue; //if below 100 MeV, just ignore & proceed
      UInt_t primSelFlag;
      if(!AcceptESDTrack(lTrack,primSelFlag,trackXYZ,ptMin,ptMax,0)) continue;
      Bool_t passDCA =(primSelFlag&1);
      Bool_t passChi2=(primSelFlag&2);
      Bool_t passBoth=(passChi2&&passDCA);
      Double_t l_Eta = lTrack->Eta();
      if(fEtaLow>-999) { if((l_Eta<fEtaLow) || (l_Eta>fEta)) continue; }
      else if(TMath::Abs(l_Eta) > fEta) continue;
      // if(TMath::Abs(lTrack->Eta()) > fEta) continue;
      fDCAxyVsPt_noChi2->Fill(lTrack->Pt(),0.,trackXYZ[0]);
      if(passDCA) fWithinDCAvsPt_noChi2->Fill(lTrack->Pt(),0.);
      if(passChi2) {
        fDCAxyVsPt_withChi2->Fill(lTrack->Pt(),0.,trackXYZ[0]);
        if(passDCA) fWithinDCAvsPt_withChi2->Fill(lTrack->Pt(),0.);
      }
    }
    PostData(1,fEfficiencyList);
    return;
  };
  Int_t nPrim = fMCEvent->GetNumberOfTracks();
  AliMCParticle *lPart;
  Int_t partNotFetched=0;
  Int_t lNchGen=0;
  Int_t lNchRec=0;
  Int_t nSpecies=fDisablePID?1:7;
  AliMCSpectraWeightsHandler* mcWeightsHandler = static_cast<AliMCSpectraWeightsHandler*>(fAOD->FindListObject("fMCSpectraWeights"));
  if(!mcWeightsHandler) AliFatal("MC weight handler not found!\n");
  AliMCSpectraWeights *fMCSpectraWeights = (mcWeightsHandler) ? mcWeightsHandler->fMCSpectraWeight : nullptr;
  if(!fMCSpectraWeights) AliFatal("Spectra weights not found!\n");
  for (Int_t ipart = 0; ipart < nPrim; ipart++) {
    lPart = (AliMCParticle*)fMCEvent->GetTrack(ipart);//ca->At(ipart);
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(ipart, fMCEvent)) continue;
    if (!lPart) { partNotFetched++; continue; };
    /* get particlePDG */
    Int_t pdgcode = TMath::Abs(lPart->PdgCode());
    if (!lPart->IsPhysicalPrimary()) continue;
    if (lPart->Charge()==0.) continue;
    Double_t pt = lPart->Pt();
    Double_t lEta = lPart->Eta();
    if (pt<0.15 || pt>50.) continue;
    if(pt>0.2 && pt<3 && TMath::Abs(lEta)<fEtaNch) lNchGen++;
    if(fEtaLow>-999) { if((lEta<fEtaLow) || (lEta>fEta)) continue; }
    else if(TMath::Abs(lEta) > fEta) continue;
    Double_t CompWeight = 1;
    if(fMCSpectraWeights) {
      CompWeight = fMCSpectraWeights->GetMCSpectraWeightNominal(lPart->Particle());
    };
    fEfficiency[nSpecies]->Fill(pt,l_Cent,CompWeight);
    if(!fDisablePID) {
      Int_t pidind = GetPIDIndex(pdgcode);
      if(pidind) fEfficiency[nSpecies+pidind]->Fill(pt,l_Cent,CompWeight);
    };
  };
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliESDtrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t lpt = lTrack->Pt();
    if(lpt>0.2 && lpt<3 && TMath::Abs(lTrack->Eta())<fEtaNch) lNchRec++;
    UInt_t primSelFlag;
    if(!AcceptESDTrack(lTrack,primSelFlag,trackXYZ,ptMin,ptMax,0)) continue;
    Bool_t passDCA =(primSelFlag&1);
    Bool_t passChi2=(primSelFlag&2);
    Bool_t passBoth=(passChi2&&passDCA);
    // if((primSelFlag&3)!=3) continue; //1 for DCA cut passed, 2 for chi2TPCGlobalVsConstrained. Needed for DCAxy dist.
    Int_t fLabel = lTrack->GetLabel();
    Int_t index = TMath::Abs(fLabel);
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(index, fMCEvent)) continue;
    if (index < 0) continue;
    lPart = (AliMCParticle*)fMCEvent->GetTrack(index);//(AliVParticle*)tca->At(index);//fMCEvent->Particle(index);
    if(!lPart) continue;
    Double_t l_Eta = lPart->Eta();
    if(fEtaLow>-999) { if((l_Eta<fEtaLow) || (l_Eta>fEta)) continue; }
    else if(TMath::Abs(l_Eta) > fEta) continue;
    Double_t CompWeight = 1;
    if(fMCSpectraWeights) {
      CompWeight = fMCSpectraWeights->GetMCSpectraWeightNominal(lPart->Particle());
    };
    Int_t pdgcode = lPart->PdgCode();
    Int_t pidind = fDisablePID?0:GetPIDIndex(pdgcode);
    Int_t pidindBayes = fDisablePID?0:GetBayesPIDIndex(lTrack);
    if(pidind && (pidind!=pidindBayes+1)) pidind=0; //If the true PID does not correspond to bayesian, then drop this track
    if(lPart->IsPhysicalPrimary()) {
        fDCAxyVsPt_noChi2->Fill(lPart->Pt(),0.,trackXYZ[0],CompWeight);
        if(passDCA) fWithinDCAvsPt_noChi2->Fill(lPart->Pt(),0.,CompWeight);
        if(passChi2) {
          fDCAxyVsPt_withChi2->Fill(lPart->Pt(),0.,trackXYZ[0],CompWeight);
          if(passDCA) fWithinDCAvsPt_withChi2->Fill(lPart->Pt(),0.,CompWeight);
        }
        if(passBoth) {
          fEfficiency[0]->Fill(lPart->Pt(),l_Cent,CompWeight);
          if(!fDisablePID && pidind) fEfficiency[pidind]->Fill(lPart->Pt(),l_Cent,CompWeight);
        };
    } else if(lPart->IsSecondaryFromWeakDecay()) {
      fDCAxyVsPt_noChi2->Fill(lPart->Pt(),1.,trackXYZ[0],CompWeight);
      if(passDCA) fWithinDCAvsPt_noChi2->Fill(lPart->Pt(),1.,CompWeight);
      if(passChi2) {
        fDCAxyVsPt_withChi2->Fill(lPart->Pt(),1.,trackXYZ[0],CompWeight);
        if(passDCA) fWithinDCAvsPt_withChi2->Fill(lPart->Pt(),1.,CompWeight);
      }
      if(passBoth) {
        fEfficiency[2*nSpecies]->Fill(lPart->Pt(),l_Cent);
        if(!fDisablePID && pidind) fEfficiency[pidind+(2*nSpecies)]->Fill(lPart->Pt(),l_Cent,CompWeight);
      }
    } else if(lPart->IsSecondaryFromMaterial()) {
      fDCAxyVsPt_noChi2->Fill(lPart->Pt(),2.,trackXYZ[0],CompWeight);
      if(passDCA) fWithinDCAvsPt_noChi2->Fill(lPart->Pt(),2.,CompWeight);
      if(passChi2) {
        fDCAxyVsPt_withChi2->Fill(lPart->Pt(),2.,trackXYZ[0],CompWeight);
        if(passDCA) fWithinDCAvsPt_withChi2->Fill(lPart->Pt(),2.,CompWeight);
      }
      if(passBoth) {
        fEfficiency[2*nSpecies]->Fill(lPart->Pt(),l_Cent);
        if(!fDisablePID && pidind) fEfficiency[pidind+(2*nSpecies)]->Fill(lPart->Pt(),l_Cent,CompWeight);
      }
    };
  };
  if(fUseCorrCuts) if(!CheckNchCorrelation(lNchGen,lNchRec)) return;
  fNchTrueVsReco->Fill(lNchGen,lNchRec);
  PostData(1,fEfficiencyList);
}

void AliAnalysisTaskMeanPtV2Corr::CovSkipMpt(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp) {
  AliAODTrack *lTrack;
  Double_t wp[5] = {0,0,0,0,0}; //Initial values, [species][w*p]
  Double_t trackXYZ[3];
  fGFW->Clear();
  Int_t iCent = fV0MMulti->FindBin(l_Cent);
  if(!iCent || iCent>fV0MMulti->GetNbinsX()) return;
  iCent--;
  Int_t lPosCount=0, lNegCount=0, lMidCount=0;
  Double_t ptMin = fPtBins[0];
  Double_t ptMax = fPtBins[fNPtBins];
  Int_t nTotNoTracks=0;
  Int_t nTotTracksFB128=0;
  if(fIsMC) {
    Int_t nTotNoTracksMC=0;
    Int_t nTotNoTracksReco=0;
    if(fUseRecoNchForMC) nTotNoTracksReco = GetNtotTracks(fAOD,ptMin,ptMax,vtxp);
    TClonesArray *tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
    Int_t nPrim = tca->GetEntries();
    AliAODMCParticle *lPart;
    for(Int_t ipart = 0; ipart < nPrim; ipart++) {
      lPart = (AliAODMCParticle*)tca->At(ipart);
      if (!lPart->IsPhysicalPrimary()) continue;
      if (lPart->Charge()==0.) continue;
      //Hardcoded cuts to inhereted from AcceptAODTrack
      Double_t leta = lPart->Eta();
      if (TMath::Abs(leta) > 0.8) continue;
      Double_t pt = lPart->Pt();
      if (pt<0.2 || pt>3.) continue;
      if(leta<-fEtaV2Sep) lNegCount++;
      if(leta>fEtaV2Sep) lPosCount++;
      if(TMath::Abs(leta)<fEtaNch) nTotNoTracksMC++; //Nch calculated in EtaNch region
      Double_t lpt = lPart->Pt();
      if(TMath::Abs(leta)<fEta)  { //for mean pt, only consider -0.4-0.4 region
        FillWPCounter(wp,1,lpt); //weight = 1, naturally
      }  //Actually, no need for if() statememnt now since GFW knows about eta's, so I can fill it all the time
      fGFW->Fill(leta,1,lPart->Phi(),1,3); //filling both gap (bit mask 1) and full (bit mas 2). Since this is MC, weight is 1.
    };
    nTotNoTracks = fUseRecoNchForMC?nTotNoTracksReco:nTotNoTracksMC;
    if(fUseRecoNchForMC) fNchTrueVsReco->Fill(nTotNoTracksMC,nTotNoTracksReco);
  } else {
    if(!LoadMyWeights(fAOD->GetRunNumber())) return; //Only load wieghts for data
    Bool_t usingPseudoEff = (fPseudoEfficiency<1);
    nTotNoTracks=GetNtotTracks(fAOD,ptMin,ptMax,vtxp);
    for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
      if(usingPseudoEff) if(fRndm->Uniform()>fPseudoEfficiency) continue;
      lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
      if(!lTrack) continue;
      Double_t leta = lTrack->Eta();
      Double_t trackXYZ[] = {0.,0.,0.};
      //Counting FB128 for QA:
      if(lTrack->TestFilterBit(128)) nTotTracksFB128++;
      if(!AcceptAODTrack(lTrack,trackXYZ,ptMin,ptMax,vtxp)) continue;
      if(leta<-fEtaV2Sep) lNegCount++;
      if(leta>fEtaV2Sep) lPosCount++;
      if(fEtaV2Sep>0 && TMath::Abs(leta)<fEtaV2Sep) lMidCount++;
      Double_t p1 = lTrack->Pt();
      Double_t weff = fEfficiencies[iCent]->GetBinContent(fEfficiencies[iCent]->FindBin(p1));
      if(weff==0) continue;
      Double_t wacc = fWeights[0]->GetNUA(lTrack->Phi(),lTrack->Eta(),vz);
      weff = 1./weff;
      if(TMath::Abs(lTrack->Eta())<fEta)  { //for mean pt, only consider -0.4-0.4 region
        FillWPCounter(wp,weff,p1);
      }  //Actually, no need for if() statememnt now since GFW knows about eta's, so I can fill it all the time
      fGFW->Fill(lTrack->Eta(),1,lTrack->Phi(),wacc*weff,3); //filling both gap (bit mask 1) and full (bit mas 2)
    };
  };
  if(wp[0]==0) return; //if no single charged particles, then surely no PID either, no sense to continue
  fMultiVsV0MCorr[0]->Fill(l_Cent,nTotNoTracks);
  //here in principle one could use the GFW output to check if the values are calculated, but this is more efficient
  if(fConsistencyFlag&1) if(!lPosCount || !lNegCount) return; // only events where v2{2, gap} could be calculated
  if(fConsistencyFlag&2) if(nTotNoTracks<4) return; //only events where v2{4} can be calculated (assuming same region as nch)
  if(fConsistencyFlag&4) if(lPosCount<2 || lNegCount<2) return; //Only events where v2{4, gap} can be calculated
  if(fConsistencyFlag&8) if(lMidCount<2) return; //If less than 2 particles in mid, reject. Relevant, if calculating v24{3-sub}
  fMultiVsV0MCorr[1]->Fill(l_Cent,nTotNoTracks);
  //Filling pT variance
  Double_t l_Multi = fUseNch?(1.0*nTotNoTracks):l_Cent;
  //A check in case l_Multi is completely off the charts (in MC, sometimes it ends up being... -Xe-310???)
  if(fUseNch && l_Multi<1) return;
  //Fetching number of ESD tracks -> for QA. Only after all the events are/were rejected
  AliAODHeader *head = (AliAODHeader*)fAOD->GetHeader();
  Int_t nESD = head->GetNumberOfESDTracks();
  fESDvsFB128->Fill(nTotTracksFB128,nESD);
  Double_t l_Random = fRndm->Rndm();
  fCkCont->FillObs(wp,l_Multi,l_Random);
  fV0MMulti->Fill(l_Cent);
  fMultiDist->Fill(l_Multi);
  PostData(1,fptVarList);
  //Filling FCs
  for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
    Bool_t filled = FillFCs(corrconfigs.at(l_ind),l_Multi,l_Random);
  };
  PostData(2,fFC);
  Double_t mptev = wp[1]/wp[0];
  for(Int_t i=0;i<1;i++) {
    FillCovariance(fCovariance[0],corrconfigs.at(i*4),l_Multi,mptev,wp[0],l_Random);
    FillCovariance(fCovariance[1],corrconfigs.at(i*4),l_Multi,1,wp[0],l_Random);
    FillCovariance(fCovariance[2],corrconfigs.at((i+1)*4),l_Multi,mptev,wp[0],l_Random);
    FillCovariance(fCovariance[3],corrconfigs.at((i+1)*4),l_Multi,1,wp[0],l_Random);
    FillCovariance(fCovariance[4],corrconfigs.at(15),l_Multi,mptev,wp[0],l_Random);
    FillCovariance(fCovariance[5],corrconfigs.at(15),l_Multi,1,wp[0],l_Random);
    //following is not necessary since we don't have any POIs
  };
  PostData(3,fCovList);
}
Bool_t AliAnalysisTaskMeanPtV2Corr::FillFCs(const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &rndmn, const Bool_t debug) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(debug) printf("FillFCs: dnx = %f\n",dnx);
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(debug) printf("FillFCs: val = %f\n",val);
    if(TMath::Abs(val)<1)
      fFC->FillProfile(corconf.Head.Data(),cent,val,fUseWeightsOne?1:dnx,rndmn);
    return kTRUE;
  };
  return kTRUE;
};
Bool_t AliAnalysisTaskMeanPtV2Corr::Fillv2dPtFCs(const AliGFW::CorrConfig &corconf, const Double_t &dpt, const Double_t &rndmn, const Int_t index) {
  if(!index || index>fV2dPtList->GetEntries()) return kFALSE;
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      ((AliGFWFlowContainer*)fV2dPtList->At(index))->FillProfile(corconf.Head.Data(),dpt,val,fUseWeightsOne?1:dnx,rndmn);
    return kTRUE;
  };
  return kTRUE;
};

Bool_t AliAnalysisTaskMeanPtV2Corr::FillCovariance(AliProfileBS *target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt, const Double_t &l_rndm) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      target->FillProfile(cent,val*d_mpt,fUseWeightsOne?1:dnx*dw_mpt,l_rndm);
    return kTRUE;
  };
  return kTRUE;
};
void AliAnalysisTaskMeanPtV2Corr::CreateCorrConfigs() {

  corrconfigs.push_back(GetConf("ChGap22","refP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("ChGap24","refP {2 2} refN {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("ChFull22","mid {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("ChFull24","mid {2 2 -2 -2}", kFALSE));
//v3:
  corrconfigs.push_back(GetConf("ChGap32","refP {3} refN {-3}", kFALSE));
  corrconfigs.push_back(GetConf("ChGap34","refP {3 3} refN {-3 -3}", kFALSE));
  corrconfigs.push_back(GetConf("ChFull32","mid {3 -3}", kFALSE));
  corrconfigs.push_back(GetConf("ChFull34","mid {3 3 -3 -3}", kFALSE));

  corrconfigs.push_back(GetConf("ChGap42","refP {4} refN {-4}", kFALSE));
//v24 3-sub
  if(fEtaV2Sep<0) return; //if eta < 0, then pos & neg are w/o SE and thus doesn't make sense to calculate v24
  corrconfigs.push_back(GetConf("LM22","refP {2} subMid {-2}", kFALSE));
  corrconfigs.push_back(GetConf("MR22","subMid {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("LR22","refP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("LLMR24","refP {2 2} subMid {-2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("LMMR24","refP {2} subMid {-2 -2} refN {2}", kFALSE));
  corrconfigs.push_back(GetConf("LMRR24","refP {2} subMid {2} refN {-2 -2}", kFALSE));

  corrconfigs.push_back(GetConf("ChSC234","refP {2 3} refN {-2 -3}", kFALSE));
  corrconfigs.push_back(GetConf("ChSC244","refP {2 4} refN {-2 -4}", kFALSE));
  return;

  //ditch the last code for now, since we don't need PID
  corrconfigs.push_back(GetConf("ChPos22","chP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("ChNeg22","chN {2} refP {-2}", kFALSE));
  corrconfigs.push_back(GetConf("ChPos24","chP refP | OLchP {2 2} refN {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("ChNeg24","chN refN | OLchN {2 2} refP {-2 -2}", kFALSE));
//pi
  corrconfigs.push_back(GetConf("PiPos22","piP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("PiNeg22","piN {2} refP {-2}", kFALSE));
  corrconfigs.push_back(GetConf("PiPos24","piP refP | OLpiP {2 2} refN {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("PiNeg24","piN refN | OLpiN {2 2} refP {-2 -2}", kFALSE));
//ka
  corrconfigs.push_back(GetConf("KaPos22","kaP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("KaNeg22","kaN {2} refP {-2}", kFALSE));
  corrconfigs.push_back(GetConf("KaPos24","kaP refP | OLkaP {2 2} refN {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("KaNeg24","kaN refN | OLkaN {2 2} refP {-2 -2}", kFALSE));
//pr
  corrconfigs.push_back(GetConf("PrPos22","prP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("PrNeg22","prN {2} refP {-2}", kFALSE));
  corrconfigs.push_back(GetConf("PrPos24","prP refP | OLprP {2 2} refN {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("PrNeg24","prN refN | OLprN {2 2} refP {-2 -2}", kFALSE));


};
void AliAnalysisTaskMeanPtV2Corr::GetSingleWeightFromList(AliGFWWeights **inWeights, TString pf) {
  (*inWeights) = (AliGFWWeights*)fWeightList->FindObject(Form("weight_%s",pf.Data()));
  if(!(*inWeights)) AliFatal(Form("Could not find weight %s in weight list\n", pf.Data()));
  if(!(*inWeights)->CalculateIntegratedEff()) AliFatal("Could not calculate integrated efficiency!\n");
  (*inWeights)->CreateNUA();
};
void AliAnalysisTaskMeanPtV2Corr::LoadWeightAndMPT() {//AliAODEvent *inEv) {
  if(!fRequireReloadOnRunChange) return;
  if(!fWeightList) AliFatal("Weight list not set!\n");

  // Int_t l_RunNo = inEv->GetRunNumber();
  TString spNames[] = {"ch","pi","ka","pr"};
  fWeights = new AliGFWWeights*[4];
  for(Int_t i=0;i<4;i++) GetSingleWeightFromList(&fWeights[i],spNames[i]);
}
Int_t AliAnalysisTaskMeanPtV2Corr::GetBayesPIDIndex(AliVTrack *l_track) {
  Double_t l_Probs[AliPID::kSPECIES];
  Double_t l_MaxProb[] = {0.95,0.85,0.85};
  Bool_t l_TOFUsed = fBayesPID->ComputeProbabilities(l_track, fPIDResponse, l_Probs) & AliPIDResponse::kDetTOF;
  Int_t pidInd = 0;
  for(Int_t i=0;i<AliPID::kSPECIES; i++) pidInd=(l_Probs[i]>l_Probs[pidInd])?i:pidInd;
  Int_t retInd = pidInd-AliPID::kPion; //Not interested in e+mu, so realign to 0
  if(retInd<0 || retInd>2) return -1; //Shouldn't be larger than 2, but just to be safe
  if(l_Probs[pidInd] < l_MaxProb[retInd]) return -1;
  //check nsigma cuts
  if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(l_track,(AliPID::EParticleType)pidInd))>3) return -1;
  if(l_TOFUsed) if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(l_track,(AliPID::EParticleType)pidInd))>3) return -1;
  return retInd;
}
Bool_t AliAnalysisTaskMeanPtV2Corr::LoadMyWeights(const Int_t &lRunNo) {
  if(!fWeightList) AliFatal("NUA list not set or does not exist!\n");
  if(lRunNo && lRunNo == fRunNo) return kTRUE;
  TString lBase(""); //base
  TString lSubfix(""); //subfix
  if(fWeightSubfix.IsNull()) { //If none specified, then follow the usual procedure
    lBase = Form("w%i",lRunNo);
    lSubfix = fGFWSelection->NeedsExtraWeight()?fGFWSelection->GetSystPF():"";
  } else {
    Int_t delind = fWeightSubfix.Index(";");
    if(delind<0) {//Only base
      lBase = fWeightSubfix;
      lSubfix = fGFWSelection->NeedsExtraWeight()?fGFWSelection->GetSystPF():"";
    } else if(!delind) {//Standard base, override subfix
      lBase = Form("w%i",lRunNo);
      lSubfix = fWeightSubfix(1,fWeightSubfix.Length());
    } else {
      lBase = fWeightSubfix(0,delind);
      lSubfix = fWeightSubfix(delind+1,fWeightSubfix.Length());
    }
  }
  lBase+=lSubfix;
  fWeights[0] = (AliGFWWeights*)fWeightList->FindObject(lBase.Data());
  if(!fWeights[0]) AliFatal(Form("Weights %s not not found in the list provided!\n",lBase.Data()));
  fWeights[0]->CreateNUA();
  return kTRUE;
}
void AliAnalysisTaskMeanPtV2Corr::SetPtBins(Int_t nPtBins, Double_t *PtBins) {
  if(fPtAxis) delete fPtAxis;
  fPtAxis = new TAxis(nPtBins, PtBins);
}
void AliAnalysisTaskMeanPtV2Corr::SetMultiBins(Int_t nMultiBins, Double_t *multibins) {
  if(fMultiAxis) delete fMultiAxis;
  fMultiAxis = new TAxis(nMultiBins, multibins);
}
void AliAnalysisTaskMeanPtV2Corr::SetV0MBins(Int_t nMultiBins, Double_t *multibins) {
  if(fV0MMultiAxis) delete fV0MMultiAxis;
  fV0MMultiAxis = new TAxis(nMultiBins, multibins);
}
void AliAnalysisTaskMeanPtV2Corr::SetV2dPtMultiBins(Int_t nMultiBins, Double_t *multibins) {
  if(fV2dPtMulti) delete fV2dPtMulti;
  fV2dPtMulti = new TH1D("v2_vs_mpt_mbins","v2_vs_mpt_mbins",nMultiBins, multibins);
}
Double_t *AliAnalysisTaskMeanPtV2Corr::GetBinsFromAxis(TAxis *inax) {
  Int_t lBins = inax->GetNbins();
  Double_t *retBins = new Double_t[lBins+1];
  for(Int_t i=0;i<lBins;i++)
    retBins[i] = inax->GetBinLowEdge(i+1);
  retBins[lBins] = inax->GetBinUpEdge(lBins);
  return retBins;
}
Int_t AliAnalysisTaskMeanPtV2Corr::GetPIDIndex(const Int_t &pdgcode) {
  if(TMath::Abs(pdgcode)==211) return 1;
  if(TMath::Abs(pdgcode)==321) return 2;
  if(TMath::Abs(pdgcode)==2212) return 3;
  if(TMath::Abs(pdgcode)==3222 || TMath::Abs(pdgcode)==3112) return 4;
  if(TMath::Abs(pdgcode)==3312) return 5;
  if(TMath::Abs(pdgcode)==3334) return 6;
  return 0;
}
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptCustomEvent(AliAODEvent* fAOD) { //From Alex
  Float_t v0Centr    = -100.;
  Float_t cl1Centr   = -100.;
  Float_t cl0Centr   = -100.;
  AliMultSelection* MultSelection = 0x0;
  MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if(!MultSelection) {
    AliWarning("AliMultSelection object not found!");
    return kFALSE;
  } else {
    v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
    cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
    cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
  }
  if(v0Centr>=80.||v0Centr<0) return kFALSE; //This would have to be adjusted for vs. V0M
  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;
  AliAODTracklets *aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets(); //ESD: esd->GetMultiplicity()->GetNumberOfTracklets()
  const Int_t nTracks = fAOD->GetNumberOfTracks(); //ESD: est->GetNumberOfTracks()
  Int_t multTrk = 0;
  for (Int_t it = 0; it < nTracks; it++) {
    AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
    if(!aodTrk){
        delete aodTrk;
        continue;
    }
    if(aodTrk->TestFilterBit(32)) multTrk++; //GetStandardITSTPCTrackCuts2011()
  }
  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;
  //pile-up cuts
  if(cl0Centr<fCenCutLowPU->Eval(v0Centr)) return kFALSE;
  if (cl0Centr > fCenCutHighPU->Eval(v0Centr)) return kFALSE;
  if(Float_t(nITSCls)>fSPDCutPU->Eval(nITSTrkls)) return kFALSE;
  if(multV0On<fV0CutPU->Eval(multV0Tot)) return kFALSE;
  if(Float_t(multTrk)<fMultCutPU->Eval(v0Centr)) return kFALSE;
  if(((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08()<0) return kFALSE;
  if(fAOD->IsIncompleteDAQ()) return kFALSE;
  return kTRUE;
}
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptCustomEvent(AliESDEvent* fESD) { //From Alex
  Float_t v0Centr    = -100.;
  Float_t cl1Centr   = -100.;
  Float_t cl0Centr   = -100.;
  AliMultSelection* MultSelection = 0x0;
  MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
  if(!MultSelection) {
    AliWarning("AliMultSelection object not found!");
    return kFALSE;
  } else {
    v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
    cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
    cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
  }
  if(v0Centr>=80.||v0Centr<0) return kFALSE; //This would have to be adjusted for vs. V0M
  Int_t nITSClsLy0 = fESD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fESD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;
  Int_t nITSTrkls = fESD->GetMultiplicity()->GetNumberOfTracklets();
  const Int_t nTracks = fESD->GetNumberOfTracks();
  Int_t multTrk = 0;
  AliESDtrack *esdTrack;
  for (Int_t it = 0; it < nTracks; it++) {
    esdTrack = (AliESDtrack*)fESD->GetTrack(it);
    if(!esdTrack) continue;
    if(fStdTPCITS2011->AcceptTrack(esdTrack)) multTrk++;
  }
  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  Float_t multV0a = esdV0->GetMTotV0A();
  Float_t multV0c = esdV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = esdV0->GetTriggerChargeA();
  UShort_t multV0cOn = esdV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;
  //pile-up cuts
  if(cl0Centr<fCenCutLowPU->Eval(v0Centr)) return kFALSE;
  if (cl0Centr > fCenCutHighPU->Eval(v0Centr)) return kFALSE;
  if(Float_t(nITSCls)>fSPDCutPU->Eval(nITSTrkls)) return kFALSE;
  if(multV0On<fV0CutPU->Eval(multV0Tot)) return kFALSE; //Problematic for MC for whatever reason? On AODs work perfectly fine
  if(Float_t(multTrk)<fMultCutPU->Eval(v0Centr)) return kFALSE;
  AliESDtrackCuts::MultEstTrackType estType = fESD->GetPrimaryVertexTracks()->GetStatus() ? AliESDtrackCuts::kTrackletsITSTPC : AliESDtrackCuts::kTracklets;
  if(AliESDtrackCuts::GetReferenceMultiplicity(fESD,estType,0.8) < 0) return kFALSE;
  if(fESD->IsIncompleteDAQ()) return kFALSE;
  return kTRUE;
}
void AliAnalysisTaskMeanPtV2Corr::FillTPCITSClusters(AliVEvent* inev) {
  Double_t nTPCClusters = inev->GetNumberOfTPCClusters()/1000.;
  Int_t nITSClusters = 0;
  AliVMultiplicity *multiObj = inev->GetMultiplicity();
  for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
  fITSvsTPCMulti->Fill(nTPCClusters,nITSClusters);
}
