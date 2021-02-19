//Class for <pt>-v2 correlations
#include "AliAnalysisTaskDeform.h"
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliAODTrack.h"
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

ClassImp(AliAnalysisTaskDeform);

AliAnalysisTaskDeform::AliAnalysisTaskDeform():
  AliAnalysisTaskSE(),
  fStageSwitch(0),
  fSystSwitch(0),
  fCentEst(0),
  fExtendV0MAcceptance(kTRUE),
  fIsMC(kFALSE),
  fIsPP(kFALSE),
  fIsXeXe(kFALSE),
  fMCEvent(0),
  fPtAxis(0),
  fMultiAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fUseNch(kFALSE),
  fUseWeightsOne(kFALSE),
  fEta(0.8),
  fEtaNch(0.8),
  fEtaV2Sep(0.4),
  fPIDResponse(0),
  fBayesPID(0),
  fMPTList(0),
  fmPT(0),
  fMultiDist(0),
  fNchVsMulti(0),
  fNchInBins(0),
  fptVarList(0),
  fptvar(0),
  fCovList(0),
  fV2dPtList(0),
  fNmptBinModifier(1),
  fCovariance(0),
  fQAList(0),
  fhCentvsNch(0),
  fmptSet(kFALSE),
  fTriggerType(AliVEvent::kMB),
  fWeightList(0),
  fWeights(0),
  fNUAList(0),
  fNUAHist(0),
  fRunNo(0),
  fGFWSelection(0),
  fFC(0),
  fGFW(0),
  fSpectraList(0),
  fSpectra(0),
  fEfficiencyList(0),
  fEfficiency(0),
  fEfficiencies(0),
  fV0MMulti(0),
  fV2dPtMulti(0),
  fFilterBit(96),
  fDisablePID(kFALSE),
  fRequireReloadOnRunChange(kFALSE)
{
};
AliAnalysisTaskDeform::AliAnalysisTaskDeform(const char *name, Bool_t IsMC, TString stageSwitch):
  AliAnalysisTaskSE(name),
  fStageSwitch(0),
  fSystSwitch(0),
  fCentEst(0),
  fExtendV0MAcceptance(kTRUE),
  fIsMC(IsMC),
  fIsPP(kFALSE),
  fIsXeXe(kFALSE),
  fMCEvent(0),
  fPtAxis(0),
  fMultiAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fUseNch(kFALSE),
  fUseWeightsOne(kFALSE),
  fEta(0.8),
  fEtaNch(0.8),
  fEtaV2Sep(0.4),
  fPIDResponse(0),
  fBayesPID(0),
  fMPTList(0),
  fmPT(0),
  fmptSet(kFALSE),
  fMultiDist(0),
  fNchVsMulti(0),
  fNchInBins(0),
  fptVarList(0),
  fptvar(0),
  fCovList(0),
  fV2dPtList(0),
  fNmptBinModifier(1),
  fCovariance(0),
  fQAList(0),
  fhCentvsNch(0),
  fTriggerType(AliVEvent::kMB),
  fWeightList(0),
  fWeights(0),
  fNUAList(0),
  fNUAHist(0),
  fRunNo(0),
  fGFWSelection(0),
  fFC(0),
  fGFW(0),
  fSpectraList(0),
  fSpectra(0),
  fEfficiencyList(0),
  fEfficiency(0),
  fEfficiencies(0),
  fV0MMulti(0),
  fV2dPtMulti(0),
  fFilterBit(96),
  fDisablePID(kFALSE),
  fRequireReloadOnRunChange(kFALSE)
{
  fStageSwitch = GetStageSwitch(stageSwitch);
  if(!fStageSwitch) AliFatal("Stage switch is 0, not sure what should be done!\n");
  if(fStageSwitch==1)
    DefineOutput(1,TList::Class());
  if(fStageSwitch==2) {
    if(!fIsMC) DefineInput(1,TList::Class());
    DefineOutput(1,TList::Class());
    DefineOutput(2,TH1D::Class());
  };
  if(fStageSwitch==3) {
    DefineInput(1,TList::Class()); //Mean Pt, should be rerun with Bayes PID
    if(!fIsMC) { //Efficiency and NUA only important for data
      DefineInput(2,TList::Class()); //NUE weights; ultimately, should be combined with NUA, but don't want to rerun now
      DefineInput(3,TList::Class()); //NUA weights from other analysis; quickfix
    };
    DefineOutput(1,TList::Class());
    DefineOutput(2,AliGFWFlowContainer::Class());
    DefineOutput(3,TList::Class());
    DefineOutput(4,TList::Class());
    DefineOutput(5,TList::Class());
  }
  if(fStageSwitch==4) {
    DefineOutput(1,TList::Class());
  }
  if(fStageSwitch==5) {
    DefineInput(1,TList::Class());
    DefineOutput(1,TList::Class());
  }
  if(fStageSwitch==6) {
    DefineOutput(1,TList::Class());
  }
  if(fStageSwitch==7) {
    DefineOutput(1,TList::Class());
  }
};
AliAnalysisTaskDeform::~AliAnalysisTaskDeform() {
};
void AliAnalysisTaskDeform::UserCreateOutputObjects(){
  printf("Stage switch is %i\n\n\n",fStageSwitch);
  if(!fGFWSelection) SetSystFlag(0);
  fGFWSelection->PrintSetup();
  if(fGFWSelection->GetSystFlagIndex() == 13) fCentEst = new TString("CL0");
  else if(fGFWSelection->GetSystFlagIndex() == 14) fCentEst = new TString("CL1");
  else fCentEst = new TString("V0M");
  OpenFile(1);
  if(fIsPP) { fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true); }
  // const Int_t nMultiBins = 300;
  // Double_t lMultiBins[nMultiBins+1];
  // for(Int_t i=0;i<=nMultiBins;i++) lMultiBins[i] = i*10;
  const Int_t l_NV0MBinsDefault = fExtendV0MAcceptance?11:(fIsXeXe?8:10);
  Double_t l_V0MBinsDefault[12] = {0,5,10,20,30,40,50,60,70,80,90,101}; //Last bin to include V0M beyond anchor point
  if(!fMultiAxis) SetMultiBins(l_NV0MBinsDefault,l_V0MBinsDefault);
  fMultiBins = GetBinsFromAxis(fMultiAxis);
  fNMultiBins = fMultiAxis->GetNbins();
  if(!fV2dPtMulti) {
    Double_t temp_bn[] = {0,1e6};
    SetV2dPtMultiBins(1,temp_bn);
  };
  const Int_t l_NPtBinsDefault = 25;
  Double_t l_PtBinsDefault[l_NPtBinsDefault+1] = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
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
        fWeights[i]->Init(kFALSE,kTRUE);
        fWeightList->Add(fWeights[i]);
      }
      PostData(1,fWeightList);
  };
  if(fStageSwitch==2) {
    fRequireReloadOnRunChange=kFALSE;
    if(!fIsMC) {
      fEfficiencyList = (TList*)GetInputData(1);
      fEfficiencies = new TH1D*[l_NV0MBinsDefault];
      for(Int_t i=0;i<l_NV0MBinsDefault;i++) {
        fEfficiencies[i] = (TH1D*)fEfficiencyList->FindObject(Form("EffRescaled_Cent%i",i));
        if(!fEfficiencies[i]) AliFatal("Could not fetch efficiency!\n");
      }
    };
    fMPTList = new TList();
    fMPTList->SetOwner(kTRUE);
    fmPT = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fmPT[i] = new TProfile(Form("MeanPt_%s",spNames[i].Data()),Form("MeanPt_%s",spNames[i].Data()),fNMultiBins,fMultiBins);
      fMPTList->Add(fmPT[i]);
    }
    fMultiDist = new TH1D("MultiDistribution","Multiplicity distribution; #it{N}_{ch}; N(events)",fNMultiBins,fMultiBins);
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",l_NV0MBinsDefault,l_V0MBinsDefault);
    fMPTList->Add(fMultiDist);
    fMPTList->Add(fV0MMulti);
    PostData(1,fMPTList);
  };
  if(fStageSwitch==3) {
    fRequireReloadOnRunChange = kFALSE;
    fMPTList = (TList*)GetInputData(1);
    if(!fMPTList) AliFatal("Could not fetch input mean pT list!\n");
    fmPT = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fmPT[i] = (TProfile*)fMPTList->FindObject(Form("MeanPt_%s",spNames[i].Data()));
      if(!fmPT[i]) AliFatal("Could not fetch mean pt!\n");
    }
    if(!fIsMC) { //Efficiencies and NUA are only for the data
      fEfficiencyList = (TList*)GetInputData(2);
      fEfficiencies = new TH1D*[l_NV0MBinsDefault];
      for(Int_t i=0;i<l_NV0MBinsDefault;i++) {
        fEfficiencies[i] = (TH1D*)fEfficiencyList->FindObject(Form("EffRescaled_Cent%i",i));
        if(!fEfficiencies[i]) AliFatal("Could not fetch efficiency!\n");
      }
      fWeightList = (TList*)GetInputData(3);
      fWeights = new AliGFWWeights*[1];
    };
    // if(!LoadMyWeights(0)) return; //Loading run-avg NUA weights
    fptVarList = new TList();
    fptVarList->SetOwner(kTRUE);
    fptvar = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fptVarList->Add(new TProfile(Form("varpt_%s",spNames[i].Data()),Form("varpt_%s",spNames[i].Data()),fNMultiBins,fMultiBins));
      fptvar[i] = (TProfile*)fptVarList->At(i);
    }
    fMultiDist = new TH1D("MultiDistribution","Multiplicity distribution; #it{N}_{ch}; N(events)",fNMultiBins,fMultiBins);
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",l_NV0MBinsDefault,l_V0MBinsDefault);
    fptVarList->Add(fMultiDist);
    fptVarList->Add(fV0MMulti);
    PostData(1,fptVarList);
    //Setting up the FlowContainer
    TObjArray *oba = new TObjArray();
    oba->Add(new TNamed("ChGap22","ChGap22")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24","ChGap24")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull22","ChFull22")); //no-gap case
    oba->Add(new TNamed("ChFull24","ChFull24")); //no-gap case

    //Following is for PID. Let's remove it for now to save some memory
/*    oba->Add(new TNamed("ChPos22","ChPos22"));
    oba->Add(new TNamed("ChPos24","ChPos24"));
    oba->Add(new TNamed("PiPos22","PiPos22"));
    oba->Add(new TNamed("PiPos24","PiPos24"));
    oba->Add(new TNamed("KaPos22","KaPos22"));
    oba->Add(new TNamed("KaPos24","KaPos24"));
    oba->Add(new TNamed("PrPos22","PrPos22"));
    oba->Add(new TNamed("PrPos24","PrPos24"));
    oba->Add(new TNamed("ChNeg22","ChNeg22"));
    oba->Add(new TNamed("ChNeg24","ChNeg24"));
    oba->Add(new TNamed("PiNeg22","PiNeg22"));
    oba->Add(new TNamed("PiNeg24","PiNeg24"));
    oba->Add(new TNamed("KaNeg22","KaNeg22"));
    oba->Add(new TNamed("KaNeg24","KaNeg24"));
    oba->Add(new TNamed("PrNeg22","PrNeg22"));
    oba->Add(new TNamed("PrNeg24","PrNeg24"));*/
    fFC = new AliGFWFlowContainer();
    TString fcname("FlowContainer");
    if(fSystSwitch) fcname.Append(Form("_%i",fSystSwitch));
    fcname.Append(fGFWSelection->GetSystPF());
    fFC->SetName(fcname.Data());
    fFC->Initialize(oba,fNMultiBins,fMultiBins);
    delete oba;
    PostData(2,fFC);
    //Initializing GFW
    Int_t pows[] = {3,0,2,0,3};
    Int_t powsFull[] = {5,0,4,0,3};
    Int_t powsPOI[] = {3,0,2,0,3};
    fGFW = new AliGFW();
    fGFW->AddRegion("refN",5,pows,-0.8,-fEtaV2Sep,1,1);
    fGFW->AddRegion("refP",5,pows,fEtaV2Sep,0.8,1,1);
    fGFW->AddRegion("mid",5,powsFull,-0.8,0.8,1,2);
    //No need to do full-blown PID, limit only with charged flow
    /*
    fGFW->AddRegion("refN",5,pows,-0.8,-0.4,1,1);
    fGFW->AddRegion("refP",5,pows,0.4,0.8,1,1);
    fGFW->AddRegion("chN",3,powsPOI,-0.8,-0.4,1,2);
    fGFW->AddRegion("chP",3,powsPOI,0.4,0.8,1,2);
    fGFW->AddRegion("piN",3,powsPOI,-0.8,-0.4,1,4);
    fGFW->AddRegion("piP",3,powsPOI,0.4,0.8,1,4);
    fGFW->AddRegion("kaN",3,powsPOI,-0.8,-0.4,1,8);
    fGFW->AddRegion("kaP",3,powsPOI,0.4,0.8,1,8);
    fGFW->AddRegion("prN",3,powsPOI,-0.8,-0.4,1,16);
    fGFW->AddRegion("prP",3,powsPOI,0.4,0.8,1,16);
    fGFW->AddRegion("OLchN",5,pows,-0.8,-0.4,1,32);
    fGFW->AddRegion("OLchP",5,pows,0.4,0.8,1,32);
    fGFW->AddRegion("OLpiN",5,pows,-0.8,-0.4,1,64);
    fGFW->AddRegion("OLpiP",5,pows,0.4,0.8,1,64);
    fGFW->AddRegion("OLkaN",5,pows,-0.8,-0.4,1,128);
    fGFW->AddRegion("OLkaP",5,pows,0.4,0.8,1,128);
    fGFW->AddRegion("OLprN",5,pows,-0.8,-0.4,1,256);
    fGFW->AddRegion("OLprP",5,pows,0.4,0.8,1,256);*/
    CreateCorrConfigs();
    //Covariance
    fCovList = new TList();
    fCovList->SetOwner(kTRUE);
    fCovariance = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fCovList->Add(new TProfile(Form("cov_%s",spNames[i].Data()),Form("cov_%s",spNames[i].Data()),fNMultiBins,fMultiBins));
      fCovariance[i] = (TProfile*)fCovList->At(i);
    };
    PostData(3,fCovList);
    fV2dPtList = new TList();
    // fV2dPtList->SetName(Form("MPtV2_%i",fSystFlag));
    fV2dPtList->SetOwner(kTRUE);
    fV2dPtList->Add(fV2dPtMulti);
    // delete oba;
    oba = new TObjArray();
    oba->Add(new TNamed("ChGap22","ChGap22"));
    if(!(20%fNmptBinModifier==0)) { AliFatal("dmpt bin modifier not divisable by 20"); return; }
    const int NmptBins = 20/fNmptBinModifier;
    for(Int_t j=0;j<fV2dPtMulti->GetNbinsX();j++) {
      AliGFWFlowContainer *fPV = new AliGFWFlowContainer();
      fPV->SetName(Form("v2dpt_%i",j));
      Double_t mptbins[NmptBins+1];
      for(Int_t i=0;i<=NmptBins;i++) mptbins[i] = (fNmptBinModifier*i - 10.)/100;
      fPV->Initialize(oba,NmptBins,mptbins);
      fV2dPtList->Add(fPV);
    };
    delete oba;
    PostData(4,fV2dPtList);

    fQAList = new TList();
    fQAList->SetOwner(kTRUE);
    fhCentvsNch = new TH2D("hCentVsCharged","Charged tracks vs Centrality",100,0,100,100,0,3000);
    fQAList->Add(fhCentvsNch);
    PostData(5,fQAList);
  }
  if(fStageSwitch==4) {
    fRequireReloadOnRunChange = kFALSE;
    fMPTList = new TList();
    fMPTList->SetOwner(kTRUE);
    fmPT = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fmPT[i] = new TProfile(Form("MeanPt_%s",spNames[i].Data()),Form("MeanPt_%s",spNames[i].Data()),fNMultiBins,fMultiBins);
      fMPTList->Add(fmPT[i]);
    }
    Double_t lV0Mbins[] = {0,5,10,20,30,40,50,60,70,80,90};
    fNchVsMulti = new TProfile("nChVsMulti","nChVsMulti",10,lV0Mbins);
    fNchInBins  = new TProfile("nChInBins" ,"nChInBins",fNMultiBins,fMultiBins);
    fMPTList->Add(fNchVsMulti);
    fMPTList->Add(fNchInBins);
    PostData(1,fMPTList);
  };
  if(fStageSwitch==5) {
    fRequireReloadOnRunChange = kFALSE;
    fMPTList = (TList*)GetInputData(1);
    if(!fMPTList) AliFatal("Could not fetch input mean pT list!\n");
    fmPT = new TProfile*[4];
    for(Int_t i=0;i<4;i++)
      fmPT[i] = (TProfile*)fMPTList->At(i);
    fptVarList = new TList();
    fptVarList->SetOwner(kTRUE);
    fptvar = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fptVarList->Add(new TProfile(Form("ptvar_%s",spNames[i].Data()),Form("ptvar_%s",spNames[i].Data()),fNMultiBins,fMultiBins));
      fptvar[i] = (TProfile*)fptVarList->At(i);
    };
    PostData(1,fptVarList);
  };
  if(fStageSwitch==6) {
    fRequireReloadOnRunChange = kFALSE;
    fSpectraList = new TList();
    fSpectraList->SetOwner(kTRUE);
    fSpectra = new TH2D*[4];
    TString lNames[] = {"ch","pi","ka","pr"};
    for(Int_t i=0;i<4;i++) {
      lNames[i].Prepend("Spectra_");
      // fSpectra[i] = new TH2D(lNames[i].Data(),lNames[i].Data(),nPtBins,PtBins,nV0MBins,lV0MBins);
      // fSpectra[i] = new TH2D(lNames[i].Data(),lNames[i].Data(),nNchPtBins,lNchPtBins,nV0MBins,lV0MBins);
      fSpectra[i] = new TH2D(lNames[i].Data(),lNames[i].Data(),fNPtBins,fPtBins,fNMultiBins,fMultiBins);
      fSpectraList->Add(fSpectra[i]);
    }
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",fNMultiBins,fMultiBins);
    fSpectraList->Add(fV0MMulti);
    PostData(1,fSpectraList);
  }
  if(fStageSwitch==7) {
    fRequireReloadOnRunChange = kFALSE;
    fEfficiencyList = new TList();
    fEfficiencyList->SetOwner(kTRUE);
    fEfficiency = new TH2D*[12];
    TString lNames[] = {"ch","pi","ka","pr"};
    for(Int_t i=0;i<4;i++) {
      lNames[i].Prepend("Spectra_");
      fEfficiency[i] = new TH2D(lNames[i].Data(),lNames[i].Data(),fNPtBins,fPtBins,fNMultiBins,fMultiBins);
      lNames[i].Append("_Gen");
      fEfficiency[4+i] = new TH2D(lNames[i].Data(),lNames[i].Data(),fNPtBins,fPtBins,fNMultiBins,fMultiBins);
      lNames[i].Append("_Sec");
      fEfficiency[8+i] = new TH2D(lNames[i].Data(),lNames[i].Data(),fNPtBins,fPtBins,fNMultiBins,fMultiBins);
      fEfficiencyList->Add(fEfficiency[i]);
      fEfficiencyList->Add(fEfficiency[i+4]);
      fEfficiencyList->Add(fEfficiency[i+8]);
    }
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",fNMultiBins,fMultiBins);
    fEfficiencyList->Add(fV0MMulti);
    PostData(1,fEfficiencyList);
  }
  if(fExtendV0MAcceptance) {
    fEventCuts.OverrideCentralityFramework(1);
    fEventCuts.SetCentralityEstimators("V0M","CL0");
    fEventCuts.SetCentralityRange(0.f,101.f);
  }
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fBayesPID = new AliPIDCombined();
  fBayesPID->SetDefaultTPCPriors();
  fBayesPID->SetSelectedSpecies(AliPID::kSPECIES);
  fBayesPID->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  LoadWeightAndMPT();
};
void AliAnalysisTaskDeform::UserExec(Option_t*) {
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCEvent) return;
  }
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t l_Cent = lMultSel->GetMultiplicityPercentile(fCentEst->Data());
  if(!CheckTrigger(l_Cent)) return;
  Double_t vtxXYZ[] = {0.,0.,0.};
  if(!AcceptAOD(fAOD, vtxXYZ)) return;
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  fGFWSelection->AcceptVertex(fAOD);
  if(fStageSwitch==1)
    FillWeights(fAOD, vz,l_Cent);
  if(fStageSwitch==2)
    fIsMC?FillMeanPtMC(fAOD,vz,l_Cent):FillMeanPt(fAOD, vz, l_Cent);
  if(fStageSwitch==3)
    FillCK(fAOD,vz,l_Cent);
  if(fStageSwitch==4)
    ProduceALICEPublished_MptProd(fAOD,vz,l_Cent);
  if(fStageSwitch==5)
    ProduceALICEPublished_CovProd(fAOD,vz,l_Cent);
  if(fStageSwitch==6)
    ProduceFBSpectra(fAOD,vz,l_Cent);
  if(fStageSwitch==7)
    ProduceEfficiencies(fAOD,vz,l_Cent);
};
void AliAnalysisTaskDeform::NotifyRun() {
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  //Reinitialize AliEventCuts (done automatically on check):
  Bool_t dummy = fEventCuts.AcceptEvent(fAOD);
  //Then override PU cut if required:
  if(fGFWSelection->GetSystFlagIndex()==15)
    fEventCuts.fESDvsTPConlyLinearCut[0] = 1500.;
}
void AliAnalysisTaskDeform::Terminate(Option_t*) {
  // fSpectraList->ls();
  // delete fSpectraList;
  // delete fSpectra;
  // delete fV0MMulti;
  // fGFWSelection->PrintSetup();
  // printf("TPC linear cut: %f\n",fEventCuts.fESDvsTPConlyLinearCut[0]);
};
Bool_t AliAnalysisTaskDeform::CheckTrigger(Double_t lCent) {
  fTriggerType = AliVEvent::kCentral;
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(fSelMask&AliVEvent::kMB) return kTRUE;
  if(fSelMask&AliVEvent::kINT7) return kTRUE;
  if((fSelMask&AliVEvent::kCentral) && lCent<10) return kTRUE;
  if((fSelMask&AliVEvent::kSemiCentral) && lCent<50 && lCent>30) return kTRUE;
  return kFALSE;
};
Bool_t AliAnalysisTaskDeform::AcceptAOD(AliAODEvent *inEv, Double_t *lvtxXYZ) {
  AliMultSelection* multSelection = nullptr;
  if(fIsPP)
  {
    if(inEv->IsPileupFromSPDInMultBins() ) { return kFALSE; }
    multSelection = (AliMultSelection*) inEv->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return kFALSE; }
    if(!multSelection->GetThisEventIsNotPileup() || !multSelection->GetThisEventIsNotPileupInMultBins() || !multSelection->GetThisEventHasNoInconsistentVertices() || !multSelection->GetThisEventPassesTrackletVsCluster()) { return kFALSE; }
    Int_t nTracksPrim = inEv->GetPrimaryVertex()->GetNContributors();
    if(nTracksPrim < 0.5) { return kFALSE; }
  }
  else
  {
    if(!fEventCuts.AcceptEvent(inEv)) return 0;
  }
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
Bool_t AliAnalysisTaskDeform::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, const Double_t &ptMin, const Double_t &ptMax, const Int_t &FilterBit) {
  // if(TMath::Abs(mtr->Eta())>0.8) return kFALSE; //0.8 is ALICE acceptance, so it's fine here if its hardcoded. Need to check explicitly in other places though
  if(mtr->Pt()<ptMin) return kFALSE;
  if(mtr->Pt()>ptMax) return kFALSE;
  // if(!mtr->TestFilterBit(FilterBit)) return kFALSE;
  // if(mtr->GetTPCNclsF()<70) return kFALSE;
  if(ltrackXYZ)
    mtr->GetXYZ(ltrackXYZ);
  return fGFWSelection->AcceptTrack(mtr,ltrackXYZ,0,kTRUE);
  // return kTRUE;
};
Bool_t AliAnalysisTaskDeform::AcceptParticle(AliVParticle *mpa) {
  if(!mpa->IsPhysicalPrimary()) return kFALSE;
  if(mpa->Charge()==0) return kFALSE;
  if(TMath::Abs(mpa->Eta())>fEta) return kFALSE;
  // if(mpa->Pt()<0.5) return kFALSE;
  // if(mpa->Pt()>2) return kFALSE;
  return kTRUE;
};
Int_t AliAnalysisTaskDeform::GetStageSwitch(TString instr) {
  if(instr.Contains("weights")) return 1;
  if(instr.Contains("meanpt")) return 2;
  if(instr.Contains("full")) return 3;
  if(instr.Contains("ALICEMpt")) return 4;
  if(instr.Contains("ALICECov")) return 5;
  if(instr.Contains("FBSpectra")) return 6;
  if(instr.Contains("Efficiency")) return 7;
  return 0;
}
void AliAnalysisTaskDeform::FillWeights(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  //MC generated
  AliVParticle *lPart;
  AliAODTrack *lTrack;
  Double_t trackXYZ[3];
  Double_t dummyDouble[] = {0.,0.};
  TClonesArray *tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
  Double_t ptMin = fPtBins[0];
  Double_t ptMax = fPtBins[fNPtBins];
  for(Int_t i=0;i<tca->GetEntries();i++) {
    lPart = (AliAODMCParticle*)tca->At(i);
    if(!AcceptParticle(lPart)) continue;
    if(!fGFWSelection->AcceptParticle(lPart,0,ptMin,ptMax)) continue;
    fWeights[0]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    Int_t pdgCode = TMath::Abs(lPart->PdgCode());
    if(pdgCode==211) fWeights[1]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    if(pdgCode==321) fWeights[2]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    if(pdgCode==2212) fWeights[3]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);

  };
  //MC reconstructed
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    lPart = (AliAODMCParticle*)tca->At(TMath::Abs(lTrack->GetLabel()));
    if(!AcceptAODTrack(lTrack,trackXYZ,ptMin,ptMax)) continue;
    if(TMath::Abs(lTrack->Eta())>fEta) continue;
    if(!fGFWSelection->AcceptTrack(lTrack,dummyDouble)) continue;
    fWeights[0]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
    if(fDisablePID) continue;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    if(PIDIndex) fWeights[PIDIndex]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
  };
  PostData(1,fWeightList);
}
void AliAnalysisTaskDeform::FillMeanPtCounter(Double_t pt, Double_t &l_sum, Double_t &l_count, AliGFWWeights *inWeight) {
  Double_t w = inWeight?inWeight->GetIntegratedEfficiency(pt):1;
  if(w==0) return;
  l_sum+=pt/w;
  l_count+=1./w;
}
void AliAnalysisTaskDeform::FillMeanPtCounterWW(const Double_t &pt, Double_t &l_sum, Double_t &l_count, const Double_t &w) {
  if(w==0) return;
  l_sum+=pt/w;
  l_count+=1./w;
}
void AliAnalysisTaskDeform::FillMeanPt(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  // LoadWeightAndMPT(fAOD);
  AliAODTrack *lTrack;
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  Double_t nTotNoTracks=0;
  Int_t iCent = fV0MMulti->FindBin(l_Cent);
  if(!iCent || iCent>fV0MMulti->GetNbinsX()) return;
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!AcceptAODTrack(lTrack,trackXYZ,0.2,3,fFilterBit)) continue;
    Double_t leta = lTrack->Eta();
    if(TMath::Abs(leta)<fEtaNch) nTotNoTracks+=1; //Nch calculated in EtaNch region
    if(TMath::Abs(leta)>fEta) continue; //<pt> calculated in fEta region
    // if(TMath::Abs(lTrack->Eta())<0.8 && lTrack->Pt()>0.2 && lTrack->Pt()<3)  nTotNoTracks++;
    Double_t lpt = lTrack->Pt();
    Double_t l_weight = fEfficiencies[iCent-1]->GetBinContent(fEfficiencies[iCent-1]->FindBin(lpt));
    FillMeanPtCounterWW(lpt,l_ptsum[0],l_ptCount[0],l_weight);
    if(fDisablePID) continue;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    if(PIDIndex) FillMeanPtCounter(lpt,l_ptsum[PIDIndex],l_ptCount[PIDIndex],fWeights[PIDIndex]);
  };
  if(l_ptCount[0]==0) return;
  Double_t lMulti  = fUseNch?nTotNoTracks:l_Cent; //Whatever the multiplicity is
  for(Int_t i=0;i<4;i++) {
    if(!l_ptCount[i]) continue;
    Double_t fillWeight = fUseWeightsOne?1:l_ptCount[i];
    fmPT[i]->Fill(lMulti,l_ptsum[i]/l_ptCount[i],fillWeight);
  }
  fMultiDist->Fill(lMulti);
  fV0MMulti->Fill(l_Cent);
  PostData(1,fMPTList);
};
void AliAnalysisTaskDeform::FillMeanPtMC(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  Double_t nTotNoTracks=0;
  Int_t iCent = fV0MMulti->FindBin(l_Cent);
  if(!iCent || iCent>fV0MMulti->GetNbinsX()) return;
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
    if(TMath::Abs(leta)<fEtaNch) nTotNoTracks+=1; //Nch calculated in EtaNch region
    if(TMath::Abs(leta)>fEta) continue; //<pt> calculated in fEta region
    Double_t lpt = lPart->Pt();
    FillMeanPtCounterWW(lpt,l_ptsum[0],l_ptCount[0],1); //MC truth, so weight = 1
  };
  if(l_ptCount[0]==0) return;
  Double_t lMulti  = fUseNch?nTotNoTracks:l_Cent; //Whatever the multiplicity is
  for(Int_t i=0;i<1;i++) { //No PID = index is only 1
    if(!l_ptCount[i]) continue;
    Double_t fillWeight = fUseWeightsOne?1:l_ptCount[i];
    fmPT[i]->Fill(lMulti,l_ptsum[i]/l_ptCount[i],fillWeight);
  }
  fMultiDist->Fill(lMulti);
  fV0MMulti->Fill(l_Cent);
  PostData(1,fMPTList);
};
void AliAnalysisTaskDeform::FillWPCounter(Double_t inArr[5], Double_t w, Double_t p) {
  inArr[0] += w;       // = w1p0
  inArr[1] += w*p;     // = w1p1
  inArr[2] += w*w*p*p; // = w2p2
  inArr[3] += w*w*p;   // = w2p1
  inArr[4] += w*w;     // = w2p0
}
void AliAnalysisTaskDeform::CalculateMptValues(Double_t outArr[4], Double_t inArr[5]) {
  //Input:
  //inArr[0] = w1p0
  //inArr[1] = w1p1
  //inArr[2] = w2p2
  //inArr[3] = w2p1
  //inArr[4] = w2p0
  //outArr[0] = <pT> (avg. over all events), has to be preset when calling the function
  //Output:
  //outArr[1] = variance
  //outAtt[2] = norm
  //outArr[3] = [pT] in this event (M(pt))
  //Assuming outArr[0] is preset to meanPt; outArr[1] = variance; outArr[2] = norm; outArr[3] = mpt in this event
  outArr[1] = TMath::Power(inArr[1] - outArr[0]*inArr[0], 2) //(w1p1 - l_meanPt*w1p0) * (w1p1 - l_meanPt*w1p0)
              - inArr[2] + 2*outArr[0]*inArr[3] - outArr[0]*outArr[0]*inArr[4]; //- w2p2 + 2*l_meanPt*w2p1 - l_meanPt*l_meanPt*w2p0;
  outArr[2] = inArr[0]*inArr[0] - inArr[4];
  outArr[3] = inArr[1]/inArr[0];
}
void AliAnalysisTaskDeform::FillCK(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  AliAODTrack *lTrack;
  Double_t wp[4][5] = {{0,0,0,0,0}, {0,0,0,0,0},
                       {0,0,0,0,0}, {0,0,0,0,0}}; //Initial values, [species][w*p]
  Double_t outVals[4][4] = {{0,0,0,0}, {0,0,0,0},
                            {0,0,0,0}, {0,0,0,0}};
  Double_t trackXYZ[3];
  fGFW->Clear();
  Double_t nTotNoTracks=0;
  Double_t ptmins[] = {0.2,0.2,0.3,0.5};
  Double_t ptmaxs[] = {10.,10.,6.0,6.0};
  Int_t iCent = fV0MMulti->FindBin(l_Cent);
  if(!iCent || iCent>fV0MMulti->GetNbinsX()) return;
  iCent--;
  if(fIsMC) {
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
      if(TMath::Abs(leta)<fEtaNch) nTotNoTracks+=1; //Nch calculated in EtaNch region
      Double_t lpt = lPart->Pt();
      if(TMath::Abs(leta)<fEta)  { //for mean pt, only consider -0.4-0.4 region
        FillWPCounter(wp[0],1,lpt); //weight = 1, naturally
      }  //Actually, no need for if() statememnt now since GFW knows about eta's, so I can fill it all the time
      fGFW->Fill(leta,1,lPart->Phi(),1,3); //filling both gap (bit mask 1) and full (bit mas 2). Since this is MC, weight is 1.
      // FillMeanPtCounterWW(lpt,l_ptsum[0],l_ptCount[0],1); //MC truth, so weight = 1
    };
  } else {
    if(!LoadMyWeights(fAOD->GetRunNumber())) return; //Only load wieghts for data
    for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
      lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
      if(!lTrack) continue;
      Double_t leta = lTrack->Eta();
      Double_t trackXYZ[] = {0.,0.,0.};
      if(!AcceptAODTrack(lTrack,trackXYZ,0.2,3,fFilterBit)) continue;
      if(TMath::Abs(leta)<fEtaNch) nTotNoTracks+=1;
      Double_t p1 = lTrack->Pt();
      Double_t weff = fEfficiencies[iCent]->GetBinContent(fEfficiencies[iCent]->FindBin(p1));
      if(weff==0) continue;
      Double_t wacc = fWeights[0]->GetNUA(lTrack->Phi(),lTrack->Eta(),vz);
      weff = 1./weff;
      if(TMath::Abs(lTrack->Eta())<fEta)  { //for mean pt, only consider -0.4-0.4 region
        FillWPCounter(wp[0],weff,p1);
      }  //Actually, no need for if() statememnt now since GFW knows about eta's, so I can fill it all the time
      fGFW->Fill(lTrack->Eta(),1,lTrack->Phi(),wacc*weff,3); //filling both gap (bit mask 1) and full (bit mas 2)
    };
  };
  if(wp[0][0]==0) return; //if no single charged particles, then surely no PID either, no sense to continue
  //Filling pT variance
  Double_t l_Multi = fUseNch?nTotNoTracks:l_Cent;
  fhCentvsNch->Fill(l_Cent,nTotNoTracks);
  //A check in case l_Multi is completely off the charts (in MC, sometimes it ends up being... -Xe-310???)
  if(fUseNch && l_Multi<1) return;
  for(Int_t i=0;i<1;i++) {
    if(!wp[i][0]) continue;
    outVals[i][0] = fmPT[i]->GetBinContent(fmPT[i]->FindBin(l_Multi));
    CalculateMptValues(outVals[i],wp[i]);
    Double_t ptvarw = fUseWeightsOne?1:outVals[i][2];
    if(outVals[i][2]!=0)
      fptvar[i]->Fill(l_Multi,outVals[i][1]/outVals[i][2],ptvarw);
  };
  fV0MMulti->Fill(l_Cent);
  fMultiDist->Fill(l_Multi);
  PostData(1,fptVarList);
  //Filling FCs
  for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
    Bool_t filled = FillFCs(corrconfigs.at(l_ind),l_Multi,0);
  };
  PostData(2,fFC);
  for(Int_t i=0;i<1;i++) {
    FillCovariance(fCovariance[i],corrconfigs.at(i*4),l_Multi,outVals[i][3]-outVals[i][0],wp[i][0]);
    //following is not necessary since we don't have any POIs
  };
  PostData(3,fCovList);
  if(outVals[0][0]==0) return;
  Int_t indx =   fV2dPtMulti->FindBin(l_Multi);
  //To avoid filling out of boundaries -- aparently, important for MC
  if(indx<1 || indx>fV2dPtMulti->GetNbinsX()) return;
  fV2dPtMulti->Fill(l_Multi);
  // printf("Will use dpt v2 profile index %i (out of %i-1), multiplicity is %f\n",indx,fV2dPtList->GetEntries(),l_Multi);
  Fillv2dPtFCs(corrconfigs.at(0),outVals[0][3]/outVals[0][0]-1,0,indx);
  PostData(4,fV2dPtList);
}
void AliAnalysisTaskDeform::ProduceALICEPublished_MptProd(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  AliAODTrack *lTrack;
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  Int_t nTotNoTracks=0;
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    Double_t lpt = lTrack->Pt();
    if(!AcceptAODTrack(lTrack,trackXYZ,0.5,2,fFilterBit)) continue;
    nTotNoTracks++;
    FillMeanPtCounter(lpt,l_ptsum[0],l_ptCount[0],0);
    if(fDisablePID) continue;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    if(PIDIndex) FillMeanPtCounter(lpt,l_ptsum[PIDIndex],l_ptCount[PIDIndex],0);
  };
  if(l_ptCount[0]==0) return;
  for(Int_t i=0;i<4;i++) {
    if(!l_ptCount[i]) continue;
    fmPT[i]->Fill(nTotNoTracks,l_ptsum[i]/l_ptCount[i],l_ptCount[i]);
  }
  fNchVsMulti->Fill(l_Cent,nTotNoTracks);
  fNchInBins->Fill(nTotNoTracks, nTotNoTracks);
  PostData(1,fMPTList);
}
void AliAnalysisTaskDeform::ProduceALICEPublished_CovProd(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  AliAODTrack *lTrack;
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  Int_t nTotNoTracks=0;
  Double_t wp[4][5] = {{0,0,0,0,0}, {0,0,0,0,0},
                       {0,0,0,0,0}, {0,0,0,0,0}}; //Initial values, [species][w*p]
  Double_t outVals[4][4] = {{0,0,0,0}, {0,0,0,0},
                            {0,0,0,0}, {0,0,0,0}};
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!AcceptAODTrack(lTrack,trackXYZ,0.5,2,fFilterBit)) continue;
    nTotNoTracks++;
    Double_t p1 = lTrack->Pt();
    FillWPCounter(wp[0],1,p1);
    if(fDisablePID) continue;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    if(PIDIndex) FillWPCounter(wp[PIDIndex],1,p1); //should be different weight here
  };
  if(wp[0][0]==0) return; //if no single charged particles, then surely no PID either, no sense to continue
  //Filling pT variance
  for(Int_t i=0;i<4;i++) {
    if(!wp[i][0]) continue;
    outVals[i][0] = fmPT[i]->GetBinContent(fmPT[i]->FindBin(nTotNoTracks));
    CalculateMptValues(outVals[i],wp[i]);
    if(outVals[i][2]!=0)
      fptvar[i]->Fill(nTotNoTracks,outVals[i][1]/outVals[i][2],outVals[i][2]);
  };
  PostData(1,fptVarList);
}
void AliAnalysisTaskDeform::ProduceFBSpectra(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  AliAODTrack *lTrack;
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    Double_t lpt = lTrack->Pt();
    if(!AcceptAODTrack(lTrack,trackXYZ,0.15,20)) continue;
    fSpectra[0]->Fill(lpt,l_Cent);
    if(fDisablePID) continue;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    if(PIDIndex) fSpectra[PIDIndex]->Fill(lpt,l_Cent);
  };
  fV0MMulti->Fill(l_Cent);//Do not care about nTracks here
  PostData(1,fSpectraList);
}
void AliAnalysisTaskDeform::ProduceEfficiencies(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  AliAODTrack *lTrack;
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  Double_t ptMin = fPtBins[0];
  Double_t ptMax = fPtBins[fNPtBins];
  Float_t dcaxy, dcaz;
  fV0MMulti->Fill(l_Cent);
  TClonesArray *tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
  Int_t nPrim = tca->GetEntries();
  AliAODMCParticle *lPart;
  Int_t partNotFetched=0;
  for (Int_t ipart = 0; ipart < nPrim; ipart++) {
    lPart = (AliAODMCParticle*)tca->At(ipart);
    if (!lPart) { partNotFetched++; continue; };
    /* get particlePDG */
    Int_t pdgcode = TMath::Abs(lPart->GetPdgCode());
    if (!lPart->IsPhysicalPrimary()) continue;
    if (lPart->Charge()==0.) continue;
    if (TMath::Abs(lPart->Eta()) > fEta) continue;
    Double_t pt = lPart->Pt();
    if (pt<0.15 || pt>50.) continue;
    fEfficiency[4]->Fill(pt,l_Cent);
    Int_t pidind = GetPIDIndex(pdgcode);
    if(pidind) fEfficiency[4+pidind]->Fill(pt,l_Cent);
  };
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    if(!AcceptAODTrack(lTrack,trackXYZ,ptMin,ptMax,fFilterBit)) continue;
    Int_t fLabel = lTrack->GetLabel();
    Int_t index = TMath::Abs(fLabel);
    if (index < 0) continue;
    lPart = (AliAODMCParticle*)tca->At(index);//fMCEvent->Particle(index);
    if(!lPart) continue;
    Int_t pdgcode = lPart->GetPdgCode();
    Int_t pidind = GetPIDIndex(pdgcode);
    Double_t lpt = lTrack->Pt();
    if(lPart->IsPhysicalPrimary()) {
        fEfficiency[0]->Fill(lPart->Pt(),l_Cent);
        if(pidind)
          fEfficiency[pidind]->Fill(lPart->Pt(),l_Cent);
    }
    if(lPart->IsSecondaryFromWeakDecay() || lPart->IsSecondaryFromMaterial()) {
        fEfficiency[8]->Fill(lPart->Pt(),l_Cent);
        if(pidind)
          fEfficiency[pidind+8]->Fill(lPart->Pt(),l_Cent);
    };
  };
  PostData(1,fEfficiencyList);
}

Bool_t AliAnalysisTaskDeform::FillFCs(const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &rndmn) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(corconf.Head.Data(),cent,val,fUseWeightsOne?1:dnx,rndmn);
    return kTRUE;
  };
  return kTRUE;
};
Bool_t AliAnalysisTaskDeform::Fillv2dPtFCs(const AliGFW::CorrConfig &corconf, const Double_t &dpt, const Double_t &rndmn, const Int_t index) {
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

Bool_t AliAnalysisTaskDeform::FillCovariance(TProfile *target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      target->Fill(cent,val*d_mpt,fUseWeightsOne?1:dnx*dw_mpt);
    return kTRUE;
  };
  return kTRUE;
};
void AliAnalysisTaskDeform::CreateCorrConfigs() {

  corrconfigs.push_back(GetConf("ChGap22","refP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("ChGap24","refP {2 2} refN {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("ChFull22","mid {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("ChFull24","mid {2 2 -2 -2}", kFALSE));
  return;

  //ditch the last code for now, since we don't need PID
  corrconfigs.push_back(GetConf("ChPos22","chP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("ChNeg22","chN {2} refP {-2}", kFALSE));
  corrconfigs.push_back(GetConf("ChPos24","chP refP | OLchP {2 2} refN {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("ChNeg24","chN refN | OLchN {2 2} refP {-2 -2}", kFALSE));
  return;
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
void AliAnalysisTaskDeform::GetSingleWeightFromList(AliGFWWeights **inWeights, TString pf) {
  (*inWeights) = (AliGFWWeights*)fWeightList->FindObject(Form("weight_%s",pf.Data()));
  if(!(*inWeights)) AliFatal(Form("Could not find weight %s in weight list\n", pf.Data()));
  if(!(*inWeights)->CalculateIntegratedEff()) AliFatal("Could not calculate integrated efficiency!\n");
  (*inWeights)->CreateNUA();
};
void AliAnalysisTaskDeform::LoadWeightAndMPT() {//AliAODEvent *inEv) {
  if(!fRequireReloadOnRunChange) return;
  if(!fWeightList) AliFatal("Weight list not set!\n");

  // Int_t l_RunNo = inEv->GetRunNumber();
  TString spNames[] = {"ch","pi","ka","pr"};
  fWeights = new AliGFWWeights*[4];
  for(Int_t i=0;i<4;i++) GetSingleWeightFromList(&fWeights[i],spNames[i]);
  if(fStageSwitch==3) { //if on switch 3 (full), then also check if need to preload dif. weight
    if(fmPT) delete [] fmPT;
    fmPT = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fmPT[i] = (TProfile*)fMPTList->FindObject(Form("MeanPt_%s",spNames[i].Data()));
      if(!fmPT[i]) AliFatal(Form("Could not find mean pT for %s in the list\n",spNames[i].Data()));
    };
  }
}
Bool_t AliAnalysisTaskDeform::WithinSigma(Double_t SigmaCut, AliAODTrack *inTrack, AliPID::EParticleType partType) {
  if(!fPIDResponse) return kFALSE;
  Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(inTrack,partType);
  Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(inTrack,partType);
  return (TMath::Sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF) < SigmaCut);
}
Int_t AliAnalysisTaskDeform::GetBayesPIDIndex(AliAODTrack *l_track) {
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
Bool_t AliAnalysisTaskDeform::LoadMyWeights(const Int_t &lRunNo) {
  if(!fWeightList) AliFatal("NUA list not set or does not exist!\n");
  if(lRunNo && lRunNo == fRunNo) return kTRUE;
  // if(!fWeights) { fWeights = new AliGFWWeights*[1]; };
  // if(fWeights[0]) delete fWeights[0];
  fWeights[0] = (AliGFWWeights*)fWeightList->FindObject(Form("w%i",lRunNo));
  if(!fWeights[0]) AliFatal(Form("Weights w%i not not found in the list provided!\n",lRunNo));
  fWeights[0]->CreateNUA();
  return kTRUE;
}
Double_t AliAnalysisTaskDeform::GetMyWeight(Double_t eta, Double_t phi, Int_t pidind) {
  Int_t etaind = fNUAHist[pidind]->GetXaxis()->FindBin(eta);
  Int_t phiind = fNUAHist[pidind]->GetYaxis()->FindBin(phi);
  return fNUAHist[pidind]->GetBinContent(etaind,phiind);
}
void AliAnalysisTaskDeform::SetPtBins(Int_t nPtBins, Double_t *PtBins) {
  if(fPtAxis) delete fPtAxis;
  fPtAxis = new TAxis(nPtBins, PtBins);
}
void AliAnalysisTaskDeform::SetMultiBins(Int_t nMultiBins, Double_t *multibins) {
  if(fMultiAxis) delete fMultiAxis;
  fMultiAxis = new TAxis(nMultiBins, multibins);
}
void AliAnalysisTaskDeform::SetV2dPtMultiBins(Int_t nMultiBins, Double_t *multibins) {
  if(fV2dPtMulti) delete fV2dPtMulti;
  fV2dPtMulti = new TH1D("v2_vs_mpt_mbins","v2_vs_mpt_mbins",nMultiBins, multibins);
}
Double_t *AliAnalysisTaskDeform::GetBinsFromAxis(TAxis *inax) {
  Int_t lBins = inax->GetNbins();
  Double_t *retBins = new Double_t[lBins+1];
  for(Int_t i=0;i<lBins;i++)
    retBins[i] = inax->GetBinLowEdge(i+1);
  retBins[lBins] = inax->GetBinUpEdge(lBins);
  return retBins;
}
Int_t AliAnalysisTaskDeform::GetPIDIndex(const Int_t &pdgcode) {
  if(TMath::Abs(pdgcode)==211) return 1;
  if(TMath::Abs(pdgcode)==321) return 2;
  if(TMath::Abs(pdgcode)==2212) return 3;
  return 0;
}
