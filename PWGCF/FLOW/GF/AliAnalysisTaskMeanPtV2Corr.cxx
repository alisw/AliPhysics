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
  fEvNomFlag(0),
  fTrNomFlag(0),
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
  fEtaAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fEtaBins(0),
  fNEtaBins(0),
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
  fRequireReloadOnRunChange(kFALSE)
{
};
AliAnalysisTaskMeanPtV2Corr::AliAnalysisTaskMeanPtV2Corr(const char *name, Bool_t IsMC, TString stageSwitch, TString ContSubfix):
  AliAnalysisTaskSE(name),
  fStageSwitch(0),
  fSystFlag(0),
  fEventCutFlag(0),
  fEvNomFlag(0),
  fTrNomFlag(0),
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
  fEtaAxis(0),
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
  fRequireReloadOnRunChange(kFALSE)
{
  fStageSwitch = GetStageSwitch(stageSwitch);
  SetContSubfix(ContSubfix);
  fCentEst = new TString("V0M");
  if(!fStageSwitch) AliFatal("Stage switch is 0, not sure what should be done!\n");
  if(fStageSwitch==1)
    DefineOutput(1,TList::Class());
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
};
AliAnalysisTaskMeanPtV2Corr::~AliAnalysisTaskMeanPtV2Corr() {
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

  const Int_t l_NEtaBinsDefault = 1;
  Double_t l_EtaBinsDefault[l_NEtaBinsDefault+1] = {-0.8,0.8};
  if(!fEtaAxis) SetEtaEffBins(l_NEtaBinsDefault,l_EtaBinsDefault);
  fEtaBins = GetBinsFromAxis(fEtaAxis);
  fNEtaBins = fEtaAxis->GetNbins();

  TString spNames[] = {"ch","pi","ka","pr"};

  if(fStageSwitch==9) {
    fRndm = new TRandom(0);
    fRequireReloadOnRunChange = kFALSE;
    if(!fIsMC) { //Efficiencies and NUA are only for the data

      //For now, ditching multiplicity eff. in favour of eta efficiency + new prescription
      /*fEfficiencyList = (TList*)GetInputData(1);
      fEfficiencies = new TH1D*[l_NV0MBinsDefault];
      for(Int_t i=0;i<l_NV0MBinsDefault;i++) {
        fEfficiencies[i] = (TH1D*)fEfficiencyList->FindObject(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
        if(fEfficiencies[i] && fPseudoEfficiency<1) fEfficiencies[i]->Scale(fPseudoEfficiency);
        if(!fEfficiencies[i]) {
          if(!i) AliFatal("Could not fetch efficiency!\n");
          printf("Could not find efficiency for V0M bin no. %i! Cloning the previous efficiency instead...\n",i);
          fEfficiencies[i] = (TH1D*)fEfficiencies[i-1]->Clone(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
        };
      }*/
      TList *effFile = (TList*)GetInputData(1);
      TList *effList = (TList*)effFile->FindObject("EfficiencyMB");
      TList *fdList  = (TList*)effFile->FindObject("FeeddownMB");
      if(!effList || !fdList) AliFatal("Efficiency or feeddown were not picked up");
      TString tempName = effList->At(0)->GetName();
      TRegexp r_period("LHC[0-9,a-z]+_");
      TString period=tempName(r_period);
      TString sysPF = GetSystPF(BitIndex(fEvNomFlag), BitIndex(fTrNomFlag));
      // fEfficiencyList = new TList();
      // fEfficiencyList->SetOwner(kTRUE);
      fEfficiencies = new TH1D*[fNEtaBins];
      for(Int_t i=0;i<fNEtaBins;i++) {
        TString etaDig(Form("%2.1f%2.1f",TMath::Abs(fEtaBins[i]),TMath::Abs(fEtaBins[i+1])));
        etaDig.ReplaceAll(".","");
        TString l_name = period+"ch_Eta_"+etaDig+sysPF;
        fEfficiencies[i] = (TH1D*)effList->FindObject(l_name.Data());
        if(!fEfficiencies[i]) printf(Form("Could not find efficiency %s!\n",l_name.Data()));
        fEfficiencies[i] = (TH1D*)fEfficiencies[i]->Clone(Form("Efficiency_%s",etaDig.Data()));
        fEfficiencies[i]->SetDirectory(0);
        TH1D *hTemp = (TH1D*)fdList->FindObject(l_name.Data());
        if(!hTemp) printf(Form("Feeddown %s not found!\n",l_name.Data()));
        fEfficiencies[i]->Divide(hTemp);
      };
      //GetSystPF(BitIndex(fEvNomFlag), BitIndex(fTrNomFlag)).Data())

      //NUAs:
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
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerType,true);
  if(fExtendV0MAcceptance) {
    fEventCuts.OverrideCentralityFramework(1);
    fEventCuts.SetCentralityEstimators("V0M","CL0");
    fEventCuts.SetCentralityRange(0.f,101.f);
  }

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
};
void AliAnalysisTaskMeanPtV2Corr::UserExec(Option_t*) {
  //On AODs with GFWFilter
  AliGFWFlags *lFlags = (AliGFWFlags*)fInputEvent->FindListObject("GFWFlags");
  if(!lFlags) { printf("GFWFlags not found!\n"); return; };
  UInt_t gEventFlag = lFlags->GetEventFlags();
  if(!(gEventFlag&fEvNomFlag)) return;
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
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  /*
  if(!AcceptAOD(fAOD, vtxXYZ)) return;
  FillTPCITSClusters(fAOD);
  if(!fGFWSelection->AcceptVertex(fAOD)) return;*/
  if(fStageSwitch==9)
    CovSkipMpt(lFlags, fAOD,vz,l_Cent,vtxXYZ);
};
void AliAnalysisTaskMeanPtV2Corr::NotifyRun() {
  LoadWeights(fInputEvent->GetRunNumber());
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

Int_t AliAnalysisTaskMeanPtV2Corr::GetStageSwitch(TString instr) {
  // if(instr.Contains("weights")) return 1;
  // if(instr.Contains("Efficiency")) return 7;
  if(instr.Contains("CovSkipMpt")) return 9;
  // if(instr.Contains("EfTest")) return 10;
  return 0;
}

Int_t AliAnalysisTaskMeanPtV2Corr::GetNtotTracks(AliAODEvent* lAOD, const Double_t &ptmin, const Double_t &ptmax, Double_t *vtxp) {
  Double_t ltrackXYZ[3];
  AliAODTrack *lTrack;
  Int_t nTotNoTracks=0;
  for(Int_t lTr=0;lTr<lAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)lAOD->GetTrack(lTr);
    if(!lTrack) continue;
    // if(!AcceptAODTrack(lTrack,ltrackXYZ,ptmin,ptmax,vtxp,nTotNoTracks)) continue;
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


void AliAnalysisTaskMeanPtV2Corr::CovSkipMpt(AliGFWFlags *lFlags, AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp) {
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
    // if(!LoadMyWeights(fAOD->GetRunNumber())) return; //Only load wieghts for data
    Bool_t usingPseudoEff = (fPseudoEfficiency<1);

    //Event check performed already
    UInt_t gTrackFlags=0;
    Int_t iEta;
    for(Int_t lTr=0;lTr<lFlags->GetNFiltered();lTr++) {
      gTrackFlags = lFlags->GetTrackFlag(lTr);
      if(!(gTrackFlags&fTrNomFlag)) continue; //Check if we want to accept the track
      nTotNoTracks++;
      Int_t trInd = lFlags->GetTrackIndex(lTr);
      lTrack = (AliAODTrack*)fAOD->GetTrack(trInd);

      if(usingPseudoEff) if(fRndm->Uniform()>fPseudoEfficiency) continue;
      Double_t leta = lTrack->Eta();
      //Counting FB128 for QA:
      if(lTrack->TestFilterBit(128)) nTotTracksFB128++;
      if(leta<-fEtaV2Sep) lNegCount++;
      if(leta>fEtaV2Sep) lPosCount++;
      if(fEtaV2Sep>0 && TMath::Abs(leta)<fEtaV2Sep) lMidCount++;
      Double_t p1 = lTrack->Pt();
      iEta = fEtaAxis->FindBin(leta);
      if(!iEta || iEta>fEtaAxis->GetNbins()) continue; //We are outside of eta efficiency range
      Double_t weff = fEfficiencies[iEta-1]->GetBinContent(fEfficiencies[iEta-1]->FindBin(p1));
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
  dnx = fGFW->Calculate(corconf,0,kTRUE).real();
  if(debug) printf("FillFCs: dnx = %f\n",dnx);
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).real()/dnx;
    if(debug) printf("FillFCs: val = %f\n",val);
    if(TMath::Abs(val)<1)
      fFC->FillProfile(corconf.Head.c_str(),cent,val,fUseWeightsOne?1:dnx,rndmn);
    return kTRUE;
  };
  return kTRUE;
};
Bool_t AliAnalysisTaskMeanPtV2Corr::Fillv2dPtFCs(const AliGFW::CorrConfig &corconf, const Double_t &dpt, const Double_t &rndmn, const Int_t index) {
  if(!index || index>fV2dPtList->GetEntries()) return kFALSE;
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).real();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).real()/dnx;
    if(TMath::Abs(val)<1)
      ((AliGFWFlowContainer*)fV2dPtList->At(index))->FillProfile(corconf.Head.c_str(),dpt,val,fUseWeightsOne?1:dnx,rndmn);
    return kTRUE;
  };
  return kTRUE;
};

Bool_t AliAnalysisTaskMeanPtV2Corr::FillCovariance(AliProfileBS *target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt, const Double_t &l_rndm) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).real();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).real()/dnx;
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

Bool_t AliAnalysisTaskMeanPtV2Corr::LoadWeights(const Int_t &runno) { //Cannot be used when running on the trains
  if(fWeightList) {
    // fWeights[0] = (AliGFWWeights*)fWeightList->FindObject(Form("w%i",runno));
    fWeights[0] = (AliGFWWeights*)fWeightList->FindObject(Form("w%i%s",runno,GetSystPF(BitIndex(fEvNomFlag), BitIndex(fTrNomFlag)).Data()));
    if(!fWeights) {
      AliFatal("Weights could not be found in the list!\n");
      return kFALSE;
    };
    fWeights[0]->CreateNUA();
    fWeights[0]->CreateNUE();
    return kTRUE;
  } else {
    AliFatal("Weight list (for some reason) not set!\n");
    return kFALSE;
  };
};


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
void AliAnalysisTaskMeanPtV2Corr::SetEtaEffBins(Int_t nEtaBins, Double_t *etabins) {
  if(fEtaAxis) delete fEtaAxis;
  fEtaAxis = new TAxis(nEtaBins, etabins);
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
void AliAnalysisTaskMeanPtV2Corr::SetupFlagsByIndex(Int_t ind) {
  fEvNomFlag=1<<kNominal;
  fTrNomFlag=1<<kFB96;
  switch(ind) {
    default: // also 0
      break;
    //Event flags:
    case 1:
      fEvNomFlag = 1<<kVtx9;
      break;
    case 2:
      fEvNomFlag = 1<<kVtx7;
      break;
    case 3:
      fEvNomFlag = 1<<kVtx5;
      break;
    //Track flags:
    case 4:
      fTrNomFlag = 1<<kFB768;
      break;
    case 5:
      fTrNomFlag = 1<<kDCAz10;
      break;
    case 6:
      fTrNomFlag = 1<<kDCAz05;
      break;
    case 7:
      fTrNomFlag = 1<<kDCA4Sigma;
      break;
    case 8:
      fTrNomFlag = 1<<kDCA10Sigma;
      break;
    case 9:
      fTrNomFlag = 1<<kChiSq2;
      break;
    case 10:
      fTrNomFlag = 1<<kChiSq3;
      break;
    case 11:
      fTrNomFlag = 1<<kNTPC80;
      break;
    case 12:
      fTrNomFlag = 1<<kNTPC90;
      break;
    case 13:
      fTrNomFlag = 1<<kNTPC100;
      break;
    case 14:
      fTrNomFlag = 1<<kFB768Tuned;
      break;
    case 15:
      fTrNomFlag = 1<<kFB96Tuned;
      break;
    case 16:
      fTrNomFlag = 1<<kFB768DCAz;
      break;
    case 17:
      fTrNomFlag = 1<<kFB768DCAxyLow;
      break;
    case 18:
      fTrNomFlag = 1<<kFB768DCAxyHigh;
      break;
    case 19:
      fTrNomFlag = 1<<kFB768ChiSq2;
      break;
    case 20:
      fTrNomFlag = 1<<kFB768ChiSq3;
      break;
    case 21:
      fTrNomFlag = 1<<kFB768nTPC;
      break;
    case 22:
      fTrNomFlag = 1<<kFB96MergedDCA;
      break;
  }
}
