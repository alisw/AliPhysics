#include "AliAnalysisTaskGammaSoft.h"
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
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliGenHijingEventHeader.h"

ClassImp(AliAnalysisTaskGammaSoft);

AliAnalysisTaskGammaSoft::AliAnalysisTaskGammaSoft():
  AliAnalysisTaskSE(),
  fSystFlag(0),
  fEventCutFlag(0),
  fContSubfix(0),
  fCentEst(0),
  fExtendV0MAcceptance(kTRUE),
  fIsMC(kFALSE),
  fBypassTriggerAndEventCuts(kFALSE),
  fDisablePileup(kFALSE),
  fUseOldPileup(kFALSE),
  fFillStdMethod(kTRUE),
  fDCAxyFunctionalForm(0),
  fOnTheFly(false),
  fGenerator("AMPT"),
  fMCEvent(0),
  fUseRecoNchForMC(kFALSE),
  fRndm(0),
  fNBootstrapProfiles(10),
  fFillAdditionalQA(kFALSE),
  fPtAxis(0),
  fEtaAxis(0),
  fMultiAxis(0),
  fV0MMultiAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fEtaBins(0),
  fNEtaBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fV0MBinsDefault(0),
  fNV0MBinsDefault(0),
  fUseNch(kFALSE),
  fUseNUAOne(kFALSE),
  fUseNUEOne(kFALSE),
  fUseEventWeightOne(kFALSE),
  fPtMpar(8),
  fEtaMpt(0.4),
  fEtaAcceptance(0.8),
  fEtaV2Sep(0.4),
  fQAList(0),
  fEventCount(0),
  fMultiDist(0),
  fNchTrueVsReco(0),
  fptList(0),
  fPtCont(0),
  fCovList(0),
  fCovariance(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fWeightList(0),
  fWeights(0),
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
  fIP(0),
  fSPDCutPU(0),
  fV0CutPU(0),
  fCenCutLowPU(0),
  fCenCutHighPU(0),
  fMultCutPU(0),
  fCentralPU(1500),
  fPhiEtaVz(0),
  fPt(0),
  fDCAxy(0),
  fDCAz(0),
  fChi2TPCcls(0),
  fTPCcls(0),
  fEtaMptAcceptance(0),
  fPtMptAcceptance(0),
  fAcceptedNch(0),
  fImpactParameterMC(-1.0),
  fStdTPCITS2011(0),
  fEventWeight(PtPtSpace::kUnity),
  fConsistencyFlag(3),
  fRequireReloadOnRunChange(kFALSE),
  fEnableFB768dcaxy(kFALSE),
  wp(0),
  abcd(0),
  wabcd(0)
{
};
AliAnalysisTaskGammaSoft::AliAnalysisTaskGammaSoft(const char *name, Bool_t IsMC, TString ContSubfix):
  AliAnalysisTaskSE(name),
  fSystFlag(0),
  fEventCutFlag(0),
  fContSubfix(0),
  fCentEst(0),
  fExtendV0MAcceptance(kTRUE),
  fIsMC(IsMC),
  fBypassTriggerAndEventCuts(kFALSE),
  fDisablePileup(kFALSE),
  fUseOldPileup(kFALSE),
  fFillStdMethod(kTRUE),
  fDCAxyFunctionalForm(0),
  fOnTheFly(false),
  fGenerator("AMPT"),
  fMCEvent(0),
  fUseRecoNchForMC(kFALSE),
  fRndm(0),
  fNBootstrapProfiles(10),
  fFillAdditionalQA(kFALSE),
  fPtAxis(0),
  fEtaAxis(0),
  fMultiAxis(0),
  fV0MMultiAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fEtaBins(0),
  fNEtaBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fV0MBinsDefault(0),
  fNV0MBinsDefault(0),
  fUseNch(kFALSE),
  fUseNUAOne(kFALSE),
  fUseNUEOne(kFALSE),
  fUseEventWeightOne(kFALSE),
  fPtMpar(8),
  fEtaMpt(0.4),
  fEtaAcceptance(0.8),
  fEtaV2Sep(0.4),
  fQAList(0),
  fEventCount(0),
  fMultiDist(0),
  fNchTrueVsReco(0),
  fPtCont(0),
  fCovList(0),
  fptList(0),
  fCovariance(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fWeightList(0),
  fWeights(0),
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
  fIP(0),
  fSPDCutPU(0),
  fV0CutPU(0),
  fCenCutLowPU(0),
  fCenCutHighPU(0),
  fMultCutPU(0),
  fCentralPU(1500),
  fPhiEtaVz(0),
  fPt(0),
  fDCAxy(0),
  fDCAz(0),
  fChi2TPCcls(0),
  fTPCcls(0),
  fEtaMptAcceptance(0),
  fPtMptAcceptance(0),
  fAcceptedNch(0),
  fImpactParameterMC(-1.0),
  fStdTPCITS2011(0),
  fEventWeight(PtPtSpace::kUnity),
  fConsistencyFlag(3),
  fRequireReloadOnRunChange(kFALSE),
  fEnableFB768dcaxy(kFALSE),
  wp(0),
  abcd(0),
  wabcd(0)
{
  SetContSubfix(ContSubfix);
  fCentEst = new TString("V0M");
  if(!fIsMC) { //Efficiency and NUA only important for data
    DefineInput(1,TList::Class()); //NUA
    DefineInput(2,TList::Class());  //NUE
  };
  DefineOutput(1,TList::Class());
  DefineOutput(2,AliGFWFlowContainer::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,TList::Class());
};
AliAnalysisTaskGammaSoft::~AliAnalysisTaskGammaSoft() {
};
void AliAnalysisTaskGammaSoft::UserCreateOutputObjects(){
  if(!fGFWSelection) SetSystFlag(0);
  fGFWSelection->PrintSetup();
  fSystFlag = fGFWSelection->GetSystFlagIndex();
  if(fGFWSelection->GetSystFlagIndex() == 20) SetCentralityEstimator("CL0");
  else if(fGFWSelection->GetSystFlagIndex() == 21) SetCentralityEstimator("CL1");
  if(!fDCAxyFunctionalForm.IsNull()) { fGFWSelection->SetPtDepDCAXY(fDCAxyFunctionalForm); }
  OpenFile(1);
  SetupAxes();
  CreateVnMptOutputObjects();
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
  fGFWNtotSelection->SetEta(fEtaAcceptance);
};
void AliAnalysisTaskGammaSoft::CreateVnMptOutputObjects(){
    fRndm = new TRandom(0);
    fRequireReloadOnRunChange = kFALSE;
    if(!fIsMC) LoadCorrectionsFromLists(); //Efficiencies and NUA are only for the data or if specified for pseudoefficiencies
    if(fOnTheFly)
    {
      printf("Creating OTF objects\n");
      printf("Generator is %s\n",fGenerator.Data());
      if(centralitymap.empty() && fGenerator.EqualTo("AMPT")) {
        vector<double> b = {0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,13.50,14.51,100.0};
        vector<double> cent = {0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0};
        for(size_t i(0); i<b.size(); ++i) centralitymap[b[i]]=cent[i];
      }
      if(centralitymap.empty() && fGenerator.EqualTo("HIJING")) {
        vector<double> b = {0.0,1.60,2.27,2.79,3.22,3.60,5.09,7.20,8.83,10.20,11.40,12.49,13.49,14.44,15.46,100.0};
        vector<double> cent = {0.0,1.0,2.0,3.0,4.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
        for(size_t i(0); i<b.size(); ++i) centralitymap[b[i]]=cent[i];
      }

      fIP = new TH1D("fIP","Impact parameter",1000,0.0,30.0);
      printf("OTF objects created\n");
    }
    // if(!LoadMyWeights(0)) return; //Loading run-avg NUA weights
    printf("Creating pt-correlation objects\n");
    fptList = new TList();
    fptList->SetOwner(kTRUE);
    fPtCont = new AliPtPtContainer("ptcont","ptcont",fNMultiBins,fMultiBins,fPtMpar);
    fPtCont->SetEventWeight(fEventWeight);
    fptList->Add(fPtCont);
    if(fNBootstrapProfiles) fPtCont->InitializeSubsamples(fNBootstrapProfiles);
    printf("pt-correlation objects created\n");
    printf("Creating multiplicity objects\n");
    fMultiDist = new TH1D("MultiDistribution","Multiplicity distribution; #it{N}_{ch}; N(events)",fNMultiBins,fMultiBins);
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",fNV0MBinsDefault,fV0MBinsDefault);
    fptList->Add(fMultiDist);
    fptList->Add(fV0MMulti);
    if(fIsMC) {
      fNchTrueVsReco = new TH2D("NchTrueVsReco",";Nch (MC-true); Nch (MC-reco)",fNMultiBins,fMultiBins,fNMultiBins,fMultiBins);
      fptList->Add(fNchTrueVsReco);
    }
    printf("Multiplicity objects created\n");

    PostData(1,fptList);
    //Setting up the FlowContainer
    printf("Creating flow container\n");
    TObjArray *oba = new TObjArray();
    oba->Add(new TNamed("ChGap22","ChGap22")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24","ChGap24")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap32","ChGap32")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap42","ChGap42")); //gap case
    oba->Add(new TNamed("ChFull22","ChFull22")); //no-gap case
    oba->Add(new TNamed("ChFull24","ChFull24")); //no-gap case
    oba->Add(new TNamed("ChFull26","ChFull26")); //no-gap case
    oba->Add(new TNamed("ChFull28","ChFull28"));
    oba->Add(new TNamed("ChFull32","ChFull32")); //no-gap case
    oba->Add(new TNamed("ChFull34","ChFull34")); //no-gap case

    fFC = new AliGFWFlowContainer();
    TString fcname("FlowContainer");
    if(!fContSubfix->IsNull()) fcname.Append(fContSubfix->Data());
    fFC->SetName(fcname.Data());
    fFC->Initialize(oba,fNMultiBins,fMultiBins,fNBootstrapProfiles);
    delete oba;
    PostData(2,fFC);
    fGFW = new AliGFW();
    fGFW->AddRegion("refN",-fEtaAcceptance,-fEtaV2Sep,1,1);
    fGFW->AddRegion("refP",fEtaV2Sep,fEtaAcceptance,1,1);
    fGFW->AddRegion("mid",-fEtaAcceptance,fEtaAcceptance,1,2);
    CreateCorrConfigs();
    fGFW->CreateRegions();
    printf("Flow container created\n");
    //Covariance
    printf("Creating covariance objects\n");
    fCovList = new TList();
    fCovList->SetOwner(kTRUE);
    const int ncovprofiles = 9;
    fCovariance = new AliProfileBS*[ncovprofiles];
    //Default terms, standard method
    fCovariance[0] = new AliProfileBS("v24pt2","v24pt2",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[0]);
    fCovariance[1] = new AliProfileBS("v24pt","v24pt",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[1]);
    fCovariance[2] = new AliProfileBS("v22pt2","v22pt2",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[2]);
    fCovariance[3] = new AliProfileBS("v22pt","v22pt",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[3]);
    //Extra terms for central moment method
    fCovariance[4] = new AliProfileBS("v24pt_6w","v24pt_6w",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[4]);
    fCovariance[5] = new AliProfileBS("v24_6w","v24_6w",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[5]);
    fCovariance[6] = new AliProfileBS("v22pt_4w","v22pt_4w",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[6]);
    fCovariance[7] = new AliProfileBS("v22_4w","v22_4w",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[7]);
    fCovariance[8] = new AliProfileBS("v22_3w","v22_3w",fNMultiBins,fMultiBins);
    fCovList->Add(fCovariance[8]);
    if(fNBootstrapProfiles) for(Int_t i=0;i<ncovprofiles;i++) fCovariance[i]->InitializeSubsamples(fNBootstrapProfiles);
    printf("Covariance objects created\n");
    PostData(3,fCovList);
    printf("Creating QA objects\n");
    fQAList = new TList();
    fQAList->SetOwner(kTRUE);
    fEventCuts.AddQAplotsToList(fQAList,kTRUE);
    int nEventCutLabel = 7;
    fEventCount = new TH1D("fEventCount","Event counter",nEventCutLabel,0,nEventCutLabel);
    TString eventCutLabel[7]={"Input","Centrality","Trigger","AliEventCuts","Vertex","Pileup","Tracks"};
    for(int i=0;i<nEventCutLabel;++i) fEventCount->GetXaxis()->SetBinLabel(i+1,eventCutLabel[i].Data());
    fQAList->Add(fEventCount);
    int NNchBins = 3000;
    double* NchBins = new double[NNchBins+1];
    for(int i(0);i<=NNchBins;++i) NchBins[i] = i+0.5;
    int NdummyCentBins = 5;
    double dummyCentBins[] = {0,2,4,6,8,10};
    if(fFillAdditionalQA) {
      fQAList->Add(fIP);
      fPhiEtaVz = new TH3D*[2];
      fPt = new TH2D*[2];
      fDCAxy = new TH2D*[2];
      fDCAz = new TH2D*[2];
      fChi2TPCcls = new TH1D*[2];
      fTPCcls = new TH1D*[2];
      TString str_cut[] = {"beforeCuts","afterCuts"};
      for(int i(0);i<2;++i){
        fPhiEtaVz[i] = new TH3D(Form("hPhiEtaVz_%s",str_cut[i].Data()),Form("#phi,#eta,v_{z} %s;#varphi;#eta;v_{z};Counts",str_cut[i].Data()),60,0,TMath::TwoPi(),64,-1.6,1.6,40,-10,10);
        fQAList->Add(fPhiEtaVz[i]);
        fPt[i] = new TH2D(Form("hPt_%s",str_cut[i].Data()),Form("#it{p}_{T} %s;#it{p}_{T};Counts",str_cut[i].Data()),fNPtBins,fPtBins,fNMultiBins,fMultiBins);
        fQAList->Add(fPt[i]);
        fDCAxy[i] = new TH2D(Form("hDCAxy_%s",str_cut[i].Data()),Form("DCAxy vs pt %s;#it{p}_{T};DCA_{xy}",str_cut[i].Data()),100,0.2,3.0,250,0,2.5);
        fQAList->Add(fDCAxy[i]);
        fDCAz[i] = new TH2D(Form("hDCAz_%s",str_cut[i].Data()),Form("DCAz vs pt %s;#it{p}_{T};DCA_{z}",str_cut[i].Data()),100,0.2,3.0,200,0,4);
        fQAList->Add(fDCAz[i]);
        fChi2TPCcls[i] = new TH1D(Form("chi2prTPCcls_%s",str_cut[i].Data()),Form("Chi2TPCcls %s;#chi^{2} pr. TPC cluster;Counts",str_cut[i].Data()),100,0,6);
        fQAList->Add(fChi2TPCcls[i]);
        fTPCcls[i] = new TH1D(Form("TPCcls_%s",str_cut[i].Data()),Form("TPCcls %s;#TPC cluster;Counts",str_cut[i].Data()),100,0,6);
        fQAList->Add(fTPCcls[i]);
      }
      fEtaMptAcceptance = new TH1D("hEtaMptAcceptance","#eta in [#it{p}_{T}] acceptance;#eta;Counts",100,-1.1,1.1);
      fQAList->Add(fEtaMptAcceptance);
      fPtMptAcceptance = new TH1D("hPtMptAcceptance","#it{p}_{T} in [#it{p}_{T}] acceptance;#it{p}_{T};Counts",100,0,5);
      fQAList->Add(fPtMptAcceptance);
      fAcceptedNch = new TH1D("AcceptedNch","Accepted N_{ch};N_{ch};Counts",4000,0.5,4000.5);
      fQAList->Add(fAcceptedNch);
    }
    fhQAEventsfMult32vsCentr = new TH2D("fhQAEventsfMult32vsCentr", "; centrality V0M; TPC multiplicity (FB32)", 100, 0, 100, 100, 0, 3000);
    fQAList->Add(fhQAEventsfMult32vsCentr);
    fhQAEventsMult128vsCentr = new TH2D("fhQAEventsfMult128vsCentr", "; centrality V0M; TPC multiplicity (FB128)", 100, 0, 100, 100, 0, 5000);
    fQAList->Add(fhQAEventsMult128vsCentr);
    fhQAEventsfMultTPCvsTOF = new TH2D("fhQAEventsfMultTPCvsTOF", "; TPC FB32 multiplicity; TOF multiplicity", 200, 0, 4000, 200, 0, 2000);
    fQAList->Add(fhQAEventsfMultTPCvsTOF);
    fhQAEventsfMultTPCvsESD = new TH2D("fhQAEventsfMultTPCvsESD", "; TPC FB128 multiplicity; ESD multiplicity", 200, 0, 7000, 300, -1000, 35000);
    fQAList->Add(fhQAEventsfMultTPCvsESD);
    printf("QA objects created!\n");
    PostData(4,fQAList);
}
void AliAnalysisTaskGammaSoft::UserExec(Option_t*) {
  if(fOnTheFly) { ProcessOnTheFly(); return; }
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCEvent) return;
  }
  fEventCount->Fill("Input",1);
  AliMultSelection *l_MultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  if(!l_MultSel) { AliFatal("MultSelection not found\n"); return; }
  Double_t l_Cent  = l_MultSel->GetMultiplicityPercentile(fCentEst->Data());
  if(l_Cent<0) return;
  fEventCount->Fill("Centrality",1);
  if(!fBypassTriggerAndEventCuts)
    if(!CheckTrigger(l_Cent)) return;
  fEventCount->Fill("Trigger",1);
  Double_t vtxXYZ[] = {0.,0.,0.};
  if(!AcceptAOD(fAOD, vtxXYZ)) return;
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  if(!fGFWSelection->AcceptVertex(fAOD)) return;
  fEventCount->Fill("Vertex",1);
  if(fUseOldPileup && IsPileupEvent(fAOD,l_Cent)) return;
  if(l_Cent < 10) fEventCuts.fESDvsTPConlyLinearCut[0] = fCentralPU;
  else fEventCuts.fESDvsTPConlyLinearCut[0] = 15000.;
  fEventCount->Fill("Pileup",1);
  ProcessTracks(fAOD,vz,l_Cent,vtxXYZ);
  return;
};
void AliAnalysisTaskGammaSoft::NotifyRun() {
  if(!fIsMC) LoadWeights(fInputEvent->GetRunNumber());
  if(!fEventCutFlag || fEventCutFlag>100) { //Only relevant if we're using the standard AliEventCuts
    //Reinitialize AliEventCuts (done automatically on check):
    Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
    if(!fDisablePileup) fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);

    //Then override PU cut if required:
    if(fGFWSelection->GetSystFlagIndex()==22)
      fEventCuts.fESDvsTPConlyLinearCut[0] = 1500.;
  };
}
void AliAnalysisTaskGammaSoft::Terminate(Option_t*) {
};
Bool_t AliAnalysisTaskGammaSoft::CheckTrigger(Double_t lCent) {
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  //Apparently, MB trigger can also mark special triggers, leaving depleted regions in multi. To avoid this, pass true, if MB has been triggered.
  //This would fail if spec. triggers would also flag MB trigger, which seems to NOT be the case.
  if(!(fTriggerType&fSelMask)) { return kFALSE; }; //printf("Returning from the generic check\n");
  if(fSelMask&(fTriggerType&(AliVEvent::kINT7+AliVEvent::kMB))) {return kTRUE; }; //printf("Passed by MB trigger!\n");
  if((fSelMask&fTriggerType&AliVEvent::kCentral) && lCent>10) {return kFALSE; }; //printf("Returnning from kCent case\n");
  if((fSelMask&fTriggerType&AliVEvent::kSemiCentral) && (lCent<30 || lCent>50)) {return kFALSE; }; //printf("Returning from kSC case\n");
  return kTRUE;
};
AliMCEvent *AliAnalysisTaskGammaSoft::getMCEvent() {
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!ev) { AliFatal("MC event not found!"); return 0; }
  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!"); return 0; }
  AliCollisionGeometry* headerH;
  TString genName;
  TList *ltgen = (TList*)ev->GetCocktailList();
  if (ltgen) {
  for(auto&& listObject: *ltgen){
    genName = Form("%s",listObject->GetName());
    if (genName.Contains("Hijing")) {
      headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
      break;
      }
    }
  }
  else headerH = dynamic_cast<AliCollisionGeometry*>(ev->GenEventHeader());
  if(headerH){
      fImpactParameterMC = headerH->ImpactParameter();
  }
  return ev;
}
double AliAnalysisTaskGammaSoft::getGeneratorCentrality()
{
  vector<double> b;
  if(centralitymap.empty()) AliFatal("Centralitymap is empty!");
  for (auto const& element : centralitymap) b.push_back(element.first);
  vector<double>::iterator it = upper_bound(b.begin(),b.end(),fImpactParameterMC);
  double l_cent = (fImpactParameterMC<0)?-1.0:(centralitymap[b[it-b.begin()]]+centralitymap[b[it-b.begin()-1]])/2.0;
  return l_cent;
}
Bool_t AliAnalysisTaskGammaSoft::IsPileupEvent(AliAODEvent* ev, double centrality){
  // Check for additional pile-up rejection in Run 2 Pb-Pb collisions (15o, 17n)
  // based on multiplicity correlations
  // ***************************************************************************

  Bool_t bIs17n = kFALSE;
  Bool_t bIs15o = kFALSE;
  Bool_t bIs18qr = kFALSE;

  Int_t iRunNumber = ev->GetRunNumber();
  if(iRunNumber >= 244824 && iRunNumber <= 246994) { bIs15o = kTRUE; }
  else if(iRunNumber == 280235 || iRunNumber == 20234) { bIs17n = kTRUE; }
  else if(iRunNumber >= 295585 && iRunNumber <= 297595 ) { bIs18qr = kTRUE; }
  else { return kFALSE; }

  // recounting multiplcities
  const Int_t multESD = ((AliAODHeader*) ev->GetHeader())->GetNumberOfESDTracks();
  const Int_t nTracks = ev->GetNumberOfTracks();
  Int_t multTPC32 = 0;
  Int_t multTPC128 = 0;
  Int_t multTOF = 0;
  Int_t multTrk = 0;
  Double_t multESDTPCdif = 0.0;
  Double_t v0Centr = 0.0;

  for(Int_t it(0); it < nTracks; it++)
  {
    AliAODTrack* track = (AliAODTrack*) ev->GetTrack(it);
    if(!track) { continue; }

    if(track->TestFilterBit(32))
    {
      multTPC32++;
      if(TMath::Abs(track->GetTOFsignalDz()) <= 10.0 && track->GetTOFsignal() >= 12000.0 && track->GetTOFsignal() <= 25000.0) { multTOF++; }
      if((TMath::Abs(track->Eta())) < 0.8 && (track->GetTPCNcls() >= 70) && (track->Pt() >= 0.2) && (track->Pt() < 3)) { multTrk++; }
    }

    if(track->TestFilterBit(128)) { multTPC128++; }
  }
  int fPileupCut = 15000;
  if(centrality < 10) fPileupCut = fCentralPU;
  if(bIs17n)
  {
    multESDTPCdif = multESD - (6.6164 + 3.64583*multTPC128 + 0.000126397*multTPC128*multTPC128);
    if(multESDTPCdif > 1000) { return kTRUE; }
    if( ((AliAODHeader*) ev->GetHeader())->GetRefMultiplicityComb08() < 0) { return kTRUE; }
  }

  if(bIs15o)
  {
    multESDTPCdif = multESD - 3.38*multTPC128;
    if(multESDTPCdif > fPileupCut) { return kTRUE; }

    TF1 fMultTOFLowCut = TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFLowCut.SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) < fMultTOFLowCut.Eval(Double_t (multTPC32))) { return kTRUE; }

    TF1 fMultTOFHighCut = TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFHighCut.SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) > fMultTOFHighCut.Eval(Double_t (multTPC32))) { return kTRUE; }

    AliMultSelection* multSelection = (AliMultSelection*) ev->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }
    v0Centr = multSelection->GetMultiplicityPercentile("V0M");

    TF1 fMultCentLowCut = TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCentLowCut.SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
    if(Double_t(multTrk) < fMultCentLowCut.Eval(v0Centr)) { return kTRUE; }
  }

  if(bIs18qr)
  {
    multESDTPCdif = multESD - 3.38*multTPC128;
    if(multESDTPCdif > fPileupCut) { return kTRUE; }

  }

  // QA Plots
  fhQAEventsfMult32vsCentr->Fill(v0Centr, multTrk);
  fhQAEventsMult128vsCentr->Fill(v0Centr, multTPC128);
  fhQAEventsfMultTPCvsTOF->Fill(multTPC32, multTOF);
  fhQAEventsfMultTPCvsESD->Fill(multTPC128, multESD);

  return kFALSE;
}
Bool_t AliAnalysisTaskGammaSoft::AcceptAOD(AliAODEvent *inEv, Double_t *lvtxXYZ) {
  if(!fBypassTriggerAndEventCuts) {
    if(!fEventCutFlag) { if(!fEventCuts.AcceptEvent(inEv)) return 0; } //Don't perform AcceptEvent if not relevant
    else if(!AcceptCustomEvent(inEv)) return 0;
    if(fEventCutFlag>100) Bool_t dummy = fEventCuts.AcceptEvent(inEv); //if flag > 100, then also store QA output from AcceptEvent
  };
  fEventCount->Fill("AliEventCuts",1);
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
Bool_t AliAnalysisTaskGammaSoft::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp) {
  if(mtr->Pt()<ptMin) return kFALSE;
  if(mtr->Pt()>ptMax) return kFALSE;
  if(ltrackXYZ && vtxp) {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2]-vtxp[2];
  } else return kFALSE; //DCA cut is a must for now
  return fGFWSelection->AcceptTrack(mtr,(fSystFlag==1&&!fEnableFB768dcaxy)?0:ltrackXYZ,0,kFALSE); //All complementary DCA track cuts for FB768 are disabled
};
Bool_t AliAnalysisTaskGammaSoft::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp, Int_t &nTot) {
  if(mtr->Pt()<ptMin) return kFALSE;
  if(mtr->Pt()>ptMax) return kFALSE;
  if(ltrackXYZ && vtxp) {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2]-vtxp[2];
  } else return kFALSE; //DCA cut is a must for now
  if(fGFWNtotSelection->AcceptTrack(mtr,ltrackXYZ,0,kFALSE)) nTot++;
  return fGFWSelection->AcceptTrack(mtr,(fSystFlag==1&&!fEnableFB768dcaxy)?0:ltrackXYZ,0,kFALSE); //All complementary DCA track cuts for FB768 are disabled
};
Int_t AliAnalysisTaskGammaSoft::GetNtotTracks(AliAODEvent* lAOD, const Double_t &ptmin, const Double_t &ptmax, Double_t *vtxp) {
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
void AliAnalysisTaskGammaSoft::FillWPCounter(vector<vector<double>> &inarr, double w, double p)
{
  for(int i=0;i<=fPtMpar;++i)
  {
    for(int j=0;j<=fPtMpar;++j)
    {
      inarr[i][j] += pow(w,i)*pow(p,j);
    }
  }
  return;
}
template<typename T>
void AliAnalysisTaskGammaSoft::FillABCDCounter(vector<vector<vector<vector<T>>>> &inarr, T a, T b, T c, T d)
{
  for(int i = 0; i<3; ++i)
    for(int j = 0; j < 3; ++j)
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < 3; ++l)
          inarr[i][j][k][l] += pow(a,i)*pow(b,j)*pow(c,k)*pow(d,l);
  return;
}
void AliAnalysisTaskGammaSoft::ProcessTracks(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp) {
  AliAODTrack *lTrack;
  wp.clear(); wp.resize(fPtMpar+1,vector<double>(fPtMpar+1));
  abcd.clear(); abcd.resize(3, vector<vector<vector<std::complex<double>>>>(3,vector<vector<std::complex<double>>>(3,vector<std::complex<double>>(3))));
  wabcd.clear(); wabcd.resize(3, vector<vector<vector<std::complex<double>>>>(3,vector<vector<std::complex<double>>>(3,vector<std::complex<double>>(3))));
  Int_t iCent = fV0MMulti->FindBin(l_Cent);
  if(!iCent || iCent>fV0MMulti->GetNbinsX()) return;
  iCent--;
  Int_t lPosCount=0, lNegCount=0, lMidCount=0;
  Double_t ptMin = fPtBins[0];
  Double_t ptMax = fPtBins[fNPtBins];
  Int_t nTotNoTracks=0;
  Int_t nTotTracksFB128=0;
  fGFW->Clear();
  if(fIsMC) {
    Int_t nTotNoTracksMC=0;
    Int_t nTotNoTracksReco=0;
    if(fUseRecoNchForMC) nTotNoTracksReco = GetNtotTracks(fAOD,ptMin,ptMax,vtxp);
    TClonesArray *tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
    Int_t nPrim = tca->GetEntries();
    if(nPrim<1) return;
    AliAODMCParticle *lPart;
    for(Int_t ipart = 0; ipart < nPrim; ipart++) {
      lPart = (AliAODMCParticle*)tca->At(ipart);
      if (!lPart->IsPhysicalPrimary()) continue;
      if (lPart->Charge()==0.) continue;
      //Hardcoded cuts to inhereted from AcceptAODTrack
      Double_t leta = lPart->Eta();
      if(TMath::Abs(leta) > 0.8) continue;
      Double_t pt = lPart->Pt();
      if(pt<ptMin || pt>ptMax) continue;
      if(leta<-fEtaV2Sep) lNegCount++;
      if(leta>fEtaV2Sep) lPosCount++;
      if(TMath::Abs(leta)<fEtaAcceptance) nTotNoTracksMC++; //Nch calculated in EtaNch region
      if(TMath::Abs(leta)<fEtaMpt)
        FillWPCounter(wp,1,pt);
      fGFW->Fill(leta,1,lPart->Phi(),1,3); //filling both gap (bit mask 1) and full (bit maks 2). Since this is MC, weight is 1.
      if(fFillAdditionalQA) {
        fPhiEtaVz[1]->Fill(lPart->Phi(),lPart->Eta(),vz);
        fPt[1]->Fill(lPart->Pt(),l_Cent);
        if(TMath::Abs(leta)<fEtaMpt){
          fEtaMptAcceptance->Fill(lPart->Eta());
          fPtMptAcceptance->Fill(lPart->Pt());
        }
      }
      if(fFillStdMethod){
        FillABCDCounter(abcd,Q(1.,2*lPart->Phi()),Q(1.,-2*lPart->Phi()),std::complex<double>(pt,0.),std::complex<double>(1.,0.));
        FillABCDCounter(wabcd,std::complex<double>(1.,0.),std::complex<double>(1.,0.),std::complex<double>(1.,0.),std::complex<double>(1.,0.));
      }


    };
    nTotNoTracks = fUseRecoNchForMC?nTotNoTracksReco:nTotNoTracksMC;
    if(fUseRecoNchForMC) fNchTrueVsReco->Fill(nTotNoTracksMC,nTotNoTracksReco);
  } else {
    //if(!LoadMyWeights(fAOD->GetRunNumber())) return; //Only load wieghts for data
    Bool_t usingPseudoEff = (fPseudoEfficiency<1);
    for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
      if(usingPseudoEff) if(fRndm->Uniform()>fPseudoEfficiency) continue;
      lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
      if(!lTrack) continue;
      if(fFillAdditionalQA) FillAdditionalTrackQAPlots(*lTrack,l_Cent,1,1,vz,vtxp,kTRUE);
      Double_t leta = lTrack->Eta();
      Double_t trackXYZ[] = {0.,0.,0.};
      //Counting FB128 for QA:
      if(lTrack->TestFilterBit(128)) nTotTracksFB128++;
      if(!AcceptAODTrack(lTrack,trackXYZ,ptMin,ptMax,vtxp)) continue;
      nTotNoTracks++;
      if(leta<-fEtaV2Sep) lNegCount++;
      if(leta>fEtaV2Sep) lPosCount++;
      if(fEtaV2Sep>0 && TMath::Abs(leta)<fEtaV2Sep) lMidCount++;
      Double_t lpt = lTrack->Pt();
      Double_t weff = fEfficiencies[iCent]->GetBinContent(fEfficiencies[iCent]->FindBin(lpt));
      if(weff==0.0) continue;
      weff = 1./weff;
      if(fUseNUEOne) weff = 1.;
      if(TMath::Abs(lTrack->Eta())<fEtaMpt) FillWPCounter(wp,weff,lpt);
      Double_t wacc = fWeights[0]->GetNUA(lTrack->Phi(),leta,vz);
      if(fUseNUAOne) wacc = 1.;
      fGFW->Fill(leta,1,lTrack->Phi(),wacc*weff,3); //filling both gap (bit mask 1) and full (bit mask 2)
      if(fFillAdditionalQA) FillAdditionalTrackQAPlots(*lTrack,l_Cent,wacc,weff,vz,vtxp,kFALSE);
      if(fFillStdMethod){
        FillABCDCounter(abcd,Q(weff*wacc,2*lTrack->Phi()),Q(weff*wacc,-2*lTrack->Phi()),std::complex<double>(weff*lpt,0.),std::complex<double>(weff,0.));
        FillABCDCounter(wabcd,std::complex<double>(weff*wacc,0.),std::complex<double>(weff*wacc,0.),std::complex<double>(weff,0.),std::complex<double>(weff,0.));
      }
    };
  };
  if(wp[1][0]==0) return; //if no single charged particles, then surely no PID either, no sense to continue
  fEventCount->Fill("Tracks",1);
  //here in principle one could use the GFW output to check if the values are calculated, but this is more efficient
  if(fConsistencyFlag&1) if(!lPosCount || !lNegCount) return; // only events where v2{2, gap} could be calculated
  if(fConsistencyFlag&2) if(nTotNoTracks<4) return; //only events where v2{4} can be calculated (assuming same region as nch)
  if(fConsistencyFlag&4) if(lPosCount < 2 && lNegCount < 2) return; //only events where v2{4} gap can be calculated
  //Filling pT variance
  Double_t l_Multi = fUseNch?(1.0*nTotNoTracks):l_Cent;
  //A check in case l_Multi is completely off the charts (in MC, sometimes it ends up being... -Xe-310???)
  if(fUseNch && l_Multi<1) return;
  Double_t l_Random = fRndm->Rndm();
  fPtCont->CalculateCorrelations(wp);
  fPtCont->FillProfiles(l_Multi,l_Random);
  fPtCont->FillCMProfiles(wp,l_Multi,l_Random);
  fV0MMulti->Fill(l_Cent);
  fMultiDist->Fill(l_Multi);
  PostData(1,fptList);
  //Filling FCs
  for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
    FillFCs(corrconfigs.at(l_ind),l_Multi,l_Random);
  };
  PostData(2,fFC);
  if(!fFillStdMethod){
    double wmpt2 = wp[1][0]*wp[1][0]-wp[2][0];
    if(wmpt2!=0.){
      double mpt2 = (wp[1][1]*wp[1][1]-wp[2][2])/wmpt2;
      double mpt_2w = (wp[1][0]*wp[1][1]-wp[2][1])/wmpt2;
      FillCovariance(fCovariance[0],corrconfigs.at(1),l_Multi,mpt2,wmpt2,l_Random);
      FillCovariance(fCovariance[2],corrconfigs.at(0),l_Multi,mpt2,wmpt2,l_Random);
      FillCovariance(fCovariance[4],corrconfigs.at(1),l_Multi,mpt_2w,wmpt2,l_Random);
      FillCovariance(fCovariance[5],corrconfigs.at(1),l_Multi,1.,wmpt2,l_Random);
      FillCovariance(fCovariance[6],corrconfigs.at(0),l_Multi,mpt_2w,wmpt2,l_Random);
      FillCovariance(fCovariance[7],corrconfigs.at(0),l_Multi,1.,wmpt2,l_Random);
    }
    double mpt = wp[1][1]/wp[1][0];
    FillCovariance(fCovariance[1],corrconfigs.at(1),l_Multi,mpt,wp[1][0],l_Random);
    FillCovariance(fCovariance[3],corrconfigs.at(0),l_Multi,mpt,wp[1][0],l_Random);
    FillCovariance(fCovariance[8],corrconfigs.at(0),l_Multi,1.,wp[1][0],l_Random);
  }
  else {
    double wAABBCC = getStdAABBCC(wabcd);
    if(wAABBCC!=0.) fCovariance[0]->FillProfile(l_Multi,getStdAABBCC(abcd)/wAABBCC,(fUseEventWeightOne)?1.:wAABBCC,l_Random);
    double wAABBC = getStdAABBC(wabcd);
    if(wAABBC!=0.) fCovariance[1]->FillProfile(l_Multi,getStdAABBC(abcd)/wAABBC,(fUseEventWeightOne)?1.:wAABBC,l_Random);
    double wABCC = getStdABCC(wabcd);
    if(wABCC!=0.) fCovariance[2]->FillProfile(l_Multi,getStdABCC(abcd)/wABCC,(fUseEventWeightOne)?1.:wABCC,l_Random);
    double wABC = getStdABC(wabcd);
    if(wABC!=0.) fCovariance[3]->FillProfile(l_Multi,getStdABC(abcd)/wABC,(fUseEventWeightOne)?1.:wABC,l_Random);
    double wAABBCD = getStdAABBCD(wabcd);
    if(wAABBCD!=0.) fCovariance[4]->FillProfile(l_Multi,getStdAABBCD(abcd)/wAABBCD,(fUseEventWeightOne)?1.:wAABBCD,l_Random);
    double wAABBDD = getStdAABBDD(wabcd);
    if(wAABBDD!=0.) fCovariance[5]->FillProfile(l_Multi,getStdAABBDD(abcd)/wAABBDD,(fUseEventWeightOne)?1.:wAABBDD,l_Random);
    double wABCD = getStdABCD(wabcd);
    if(wABCD!=0.) fCovariance[6]->FillProfile(l_Multi,getStdABCD(abcd)/wABCD,(fUseEventWeightOne)?1.:wABCD,l_Random);
    double wABDD = getStdABDD(wabcd);
    if(wABDD!=0.) fCovariance[7]->FillProfile(l_Multi,getStdABDD(abcd)/wABDD,(fUseEventWeightOne)?1.:wABDD,l_Random);
    double wABD = getStdABD(wabcd);
    if(wABD!=0.) fCovariance[8]->FillProfile(l_Multi,getStdABD(abcd)/wABD,(fUseEventWeightOne)?1.:wABD,l_Random);
  }

  PostData(3,fCovList);
}
Bool_t AliAnalysisTaskGammaSoft::FillFCs(const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &rndmn, const Bool_t debug) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).real();
  if(debug) printf("FillFCs: dnx = %f\n",dnx);
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).real()/dnx;
    if(debug) printf("FillFCs: val = %f\n",val);
    if(TMath::Abs(val)<1)
      fFC->FillProfile(corconf.Head.c_str(),cent,val,(fUseEventWeightOne)?1.0:dnx,rndmn);
    return kTRUE;
  };
  return kTRUE;
};
void AliAnalysisTaskGammaSoft::ProcessOnTheFly() {
  fMCEvent = getMCEvent();
  fIP->Fill(fImpactParameterMC);
  Double_t l_Cent = getGeneratorCentrality();
  Int_t nTracks = fMCEvent->GetNumberOfPrimaries();
  if(nTracks < 1) { return; }
  wp.clear(); wp.resize(fPtMpar+1,vector<double>(fPtMpar+1));
  fGFW->Clear();
  Double_t ptMin = fPtBins[0];
  Double_t ptMax = fPtBins[fNPtBins];
  for(Int_t i=0;i<nTracks;i++) {
    AliMCParticle* lPart = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(i));
    if(!lPart) { continue; };
    if(!lPart->IsPhysicalPrimary()) continue;
    Double_t l_pt=lPart->Pt();
    Double_t l_phi=lPart->Phi();
    Double_t l_eta=lPart->Eta();
    if (TMath::Abs(l_eta) > fEtaAcceptance) continue;
    if (l_pt<ptMin || l_pt>ptMax) continue;
    if(TMath::Abs(l_eta)<fEtaMpt) //for mean pt, only consider -0.4-0.4 region
      FillWPCounter(wp,1,l_pt);
    fGFW->Fill(l_eta,1,l_phi,1,3); //filling both gap (bit mask 1) and full (bit mas 2). Since this is MC, weight is 1.
  };
  if(wp[0][0]==0) return; //if no single charged particles, then surely no PID either, no sense to continue
  Double_t l_Random = fRndm->Rndm();
  fPtCont->CalculateCorrelations(wp);
  fPtCont->FillProfiles(l_Cent,l_Random);
  fPtCont->FillCMProfiles(wp,l_Cent,l_Random);
  fV0MMulti->Fill(l_Cent);
  fMultiDist->Fill(l_Cent);
  PostData(1,fptList);
  //Filling FCs
  for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
    FillFCs(corrconfigs.at(l_ind),l_Cent,l_Random);
  };
  PostData(2,fFC);
  PostData(3,fCovList);
  return;
}
Bool_t AliAnalysisTaskGammaSoft::FillCovariance(AliProfileBS *target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt, const Double_t &l_rndm) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).real();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).real()/dnx;
    if(TMath::Abs(val)<1)
      target->FillProfile(cent,val*d_mpt,(fUseEventWeightOne)?1.0:dnx*dw_mpt,l_rndm);
    return kTRUE;
  };
  return kTRUE;
};
void AliAnalysisTaskGammaSoft::CreateCorrConfigs() {

  corrconfigs.push_back(GetConf("ChGap22","refP {2} refN {-2}", kFALSE));     //ChGap22 0
  corrconfigs.push_back(GetConf("ChGap24","refP {2 2} refN {-2 -2}", kFALSE));  //ChGap24 1
  corrconfigs.push_back(GetConf("ChGap32","refP {3} refN {-3}", kFALSE));   //ChGap32 2
  corrconfigs.push_back(GetConf("ChGap42","refP {4} refN {-4}", kFALSE));   //ChGap42 3
  corrconfigs.push_back(GetConf("ChFull22","mid {2 -2}", kFALSE));  //ChFull22 4
  corrconfigs.push_back(GetConf("ChFull24","mid {2 2 -2 -2}", kFALSE));   //ChFull24 5
  corrconfigs.push_back(GetConf("ChFull26","mid {2 2 2 -2 -2 -2}",kFALSE));  //  ChFull26 6
  corrconfigs.push_back(GetConf("ChFull28","mid {2 2 2 2 -2 -2 -2 -2}",kFALSE));  //  ChFull28 7
  corrconfigs.push_back(GetConf("ChFull32","mid {3 -3}", kFALSE));   //ChFull32 8
  corrconfigs.push_back(GetConf("ChFull34","mid {3 3 -3 -3}", kFALSE));    //ChFull34 9
  return;
};
void AliAnalysisTaskGammaSoft::FillAdditionalTrackQAPlots(AliAODTrack &track, const Double_t &cent, Double_t weff, Double_t wacc, const Double_t &vz, Double_t* vtxp, Bool_t beforeCuts){
  Double_t trackXYZ[] = {0.,0.,0.};
  track.GetXYZ(trackXYZ);
  trackXYZ[0] = trackXYZ[0]-vtxp[0];
  trackXYZ[1] = trackXYZ[1]-vtxp[1];
  trackXYZ[2] = trackXYZ[2]-vtxp[2];
  if(beforeCuts){
    fDCAxy[0]->Fill(track.Pt(),TMath::Sqrt(trackXYZ[0]*trackXYZ[0]+trackXYZ[1]*trackXYZ[1]));
    fDCAz[0]->Fill(track.Pt(),TMath::Abs(trackXYZ[2]));
    fChi2TPCcls[0]->Fill(track.GetTPCchi2perCluster());
    fTPCcls[0]->Fill(track.GetTPCncls());
  }
  else {
    fPhiEtaVz[0]->Fill(track.Phi(),track.Eta(),vz);
    fPt[0]->Fill(track.Pt(),cent);
    fPhiEtaVz[1]->Fill(track.Phi(),track.Eta(),vz,wacc);
    fPt[1]->Fill(track.Pt(),cent,weff);
    if(TMath::Abs(track.Eta())<fEtaMpt) {
      fEtaMptAcceptance->Fill(track.Eta());
      fPtMptAcceptance->Fill(track.Pt(),weff);
    }
    fDCAxy[1]->Fill(track.Pt(),TMath::Sqrt(trackXYZ[0]*trackXYZ[0]+trackXYZ[1]*trackXYZ[1]));
    fDCAz[1]->Fill(track.Pt(),TMath::Abs(trackXYZ[2]));
    fChi2TPCcls[1]->Fill(track.GetTPCchi2perCluster());
    fTPCcls[1]->Fill(track.GetTPCncls());
  }
}
double AliAnalysisTaskGammaSoft::getStdAABBCC(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
  std::complex<double> a = abcdvec[1][0][0][0];
  std::complex<double> b = abcdvec[0][1][0][0];
  std::complex<double> c = abcdvec[0][0][1][0];
  std::complex<double> aa = abcdvec[2][0][0][0];
  std::complex<double> bb = abcdvec[0][2][0][0];
  std::complex<double> cc = abcdvec[0][0][2][0];
  std::complex<double> ab = abcdvec[1][1][0][0];
  std::complex<double> ac = abcdvec[1][0][1][0];
  std::complex<double> bc = abcdvec[0][1][1][0];
  std::complex<double> aab = abcdvec[2][1][0][0];
  std::complex<double> aac = abcdvec[2][0][1][0];
  std::complex<double> abb = abcdvec[1][2][0][0];
  std::complex<double> acc = abcdvec[1][0][2][0];
  std::complex<double> abc = abcdvec[1][1][1][0];
  std::complex<double> bbc = abcdvec[0][2][1][0];
  std::complex<double> bcc = abcdvec[0][1][2][0];
  std::complex<double> aabb = abcdvec[2][2][0][0];
  std::complex<double> aacc = abcdvec[2][0][2][0];
  std::complex<double> aabc = abcdvec[2][1][1][0];
  std::complex<double> abbc = abcdvec[1][2][1][0];
  std::complex<double> abcc = abcdvec[1][1][2][0];
  std::complex<double> bbcc = abcdvec[0][2][2][0];
  std::complex<double> aabbc = abcdvec[2][2][1][0];
  std::complex<double> aabcc = abcdvec[2][1][2][0];
  std::complex<double> abbcc = abcdvec[0][0][0][0];
  std::complex<double> aabbcc = abcdvec[2][2][2][0];
  return (a*a*b*b*c*c - aa*b*b*c*c - a*a*bb*c*c - a*a*b*b*cc - 4.*a*ab*b*c*c -
 4.*a*ac*b*b*c - 4.*a*a*b*bc*c + 4.*aab*b*c*c + 4.*aac*b*b*c +
 4.*a*abb*c*c + 4.*a*acc*b*b + 4.*a*a*bbc*c + 4.*a*a*b*bcc +
 16.*a*abc*b*c + aa*bb*c*c + aa*b*b*cc + a*a*bb*cc + 2.*ab*ab*c*c +
 2.*ac*ac*b*b + 2.*a*a*bc*bc + 4.*aa*b*bc*c + 4.*a*ac*bb*c +
 4.*a*ab*b*cc + 8.*ab*ac*b*c + 8.*a*ab*bc*c + 8.*a*ac*b*bc - 6.*aabb*c*c -
 24.*aabc*b*c - 6.*aacc*b*b - 24.*abbc*a*c - 24.*abcc*a*b - 6.*bbcc*a*a -
 8.*aab*bc*c - 8.*aac*b*bc - 4.*aac*bb*c - 4.*aab*b*cc - 8.*abb*ac*c -
 4.*abb*a*cc - 8.*acc*ab*b - 4.*acc*a*bb - 8.*bbc*a*ac - 4.*bbc*aa*c -
 8.*bcc*a*ab - 4.*bcc*aa*b - 16.*abc*ab*c - 16.*abc*ac*b - 16.*abc*a*bc -
 aa*bb*cc - 2.*ab*ab*cc - 2.*ac*ac*bb - 2.*bc*bc*aa - 8.*ab*ac*bc +
 48.*aabbc*c + 48.*aabcc*b + 48.*abbcc*a + 6.*aabb*cc + 6.*aacc*bb +
 6.*bbcc*aa + 24.*aabc*bc + 24.*abbc*ac + 24.*abcc*ab + 8.*aab*bcc +
 8.*aac*bbc + 8.*abb*acc + 16.*abc*abc - 120.*aabbcc).real();
}
double AliAnalysisTaskGammaSoft::getStdAABBCD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
  std::complex<double> a = abcdvec[1][0][0][0];
  std::complex<double> b = abcdvec[0][1][0][0];
  std::complex<double> c = abcdvec[0][0][1][0];
  std::complex<double> d = abcdvec[0][0][0][1];
  std::complex<double> aa = abcdvec[2][0][0][0];
  std::complex<double> bb = abcdvec[0][2][0][0];
  std::complex<double> ab = abcdvec[1][1][0][0];
  std::complex<double> ac = abcdvec[1][0][1][0];
  std::complex<double> ad = abcdvec[1][0][0][1];
  std::complex<double> bc = abcdvec[0][1][1][0];
  std::complex<double> bd = abcdvec[0][1][0][1];
  std::complex<double> cd = abcdvec[0][0][1][1];
  std::complex<double> aab = abcdvec[2][1][0][0];
  std::complex<double> aac = abcdvec[2][0][1][0];
  std::complex<double> aad = abcdvec[2][0][0][1];
  std::complex<double> abb = abcdvec[1][2][0][0];
  std::complex<double> abc = abcdvec[1][1][1][0];
  std::complex<double> abd = abcdvec[1][1][0][1];
  std::complex<double> acd = abcdvec[1][0][1][1];
  std::complex<double> bbc = abcdvec[0][2][1][0];
  std::complex<double> bbd = abcdvec[0][2][0][1];
  std::complex<double> bcd = abcdvec[0][1][1][1];
  std::complex<double> aabb = abcdvec[2][2][0][0];
  std::complex<double> aabc = abcdvec[2][1][1][0];
  std::complex<double> aabd = abcdvec[2][1][0][1];
  std::complex<double> aacd = abcdvec[2][0][1][1];
  std::complex<double> abbc = abcdvec[1][2][1][0];
  std::complex<double> abbd = abcdvec[1][2][0][1];
  std::complex<double> abcd = abcdvec[0][1][1][1];
  std::complex<double> bbcd = abcdvec[0][2][1][1];
  std::complex<double> aabbc = abcdvec[2][2][1][0];
  std::complex<double> aabbd = abcdvec[2][2][0][1];
  std::complex<double> aabcd = abcdvec[2][1][1][1];
  std::complex<double> abbcd = abcdvec[1][2][1][1];
  std::complex<double> aabbcd = abcdvec[2][2][1][1];
  return (-120.*aabbcd + 48.*a*abbcd + 24.*ab*abcd + 16.*abc*abd + 12.*abbd*ac +
 8.*abb*acd + 12.*abbc*ad + 48.*aabcd*b - 24.*a*abcd*b - 8.*abd*ac*b -
 8.*ab*acd*b - 8.*abc*ad*b - 6.*aacd*b*b + 4.*a*acd*b*b + 2.*ac*ad*b*b +
 6.*aacd*bb - 4.*a*acd*bb - 2.*ac*ad*bb + 4.*aad*bbc - 4.*a*ad*bbc -
 6.*a*a*bbcd + 6.*aa*bbcd + 4.*aac*bbd - 4.*a*ac*bbd + 12.*aabd*bc -
 8.*a*abd*bc - 4.*ab*ad*bc - 4.*aad*b*bc + 4.*a*ad*b*bc + 8.*aab*bcd -
 8.*a*ab*bcd + 4.*a*a*b*bcd - 4.*aa*b*bcd + 12.*aabc*bd - 8.*a*abc*bd -
 4.*ab*ac*bd - 4.*aac*b*bd + 4.*a*ac*b*bd + 2.*a*a*bc*bd - 2.*aa*bc*bd +
 24.*aabbd*c - 12.*a*abbd*c - 8.*ab*abd*c - 4.*abb*ad*c - 12.*aabd*b*c +
 8.*a*abd*b*c + 4.*ab*ad*b*c + 2.*aad*b*b*c - 2.*a*ad*b*b*c -
 2.*aad*bb*c + 2.*a*ad*bb*c + 2.*a*a*bbd*c - 2.*aa*bbd*c - 4.*aab*bd*c +
 4.*a*ab*bd*c - 2.*a*a*b*bd*c + 2.*aa*b*bd*c + 6.*aabb*cd - 2.*ab*ab*cd -
 4.*a*abb*cd - 4.*aab*b*cd + 4.*a*ab*b*cd - a*a*b*b*cd + aa*b*b*cd +
 a*a*bb*cd - aa*bb*cd + 24.*aabbc*d - 12.*a*abbc*d - 8.*ab*abc*d -
 4.*abb*ac*d - 12.*aabc*b*d + 8.*a*abc*b*d + 4.*ab*ac*b*d + 2.*aac*b*b*d -
 2.*a*ac*b*b*d - 2.*aac*bb*d + 2.*a*ac*bb*d + 2.*a*a*bbc*d - 2.*aa*bbc*d -
 4.*aab*bc*d + 4.*a*ab*bc*d - 2.*a*a*b*bc*d + 2.*aa*b*bc*d - 6.*aabb*c*d +
 2.*ab*ab*c*d + 4.*a*abb*c*d + 4.*aab*b*c*d - 4.*a*ab*b*c*d +
 a*a*b*b*c*d - aa*b*b*c*d - a*a*bb*c*d + aa*bb*c*d).real();
}
double AliAnalysisTaskGammaSoft::getStdAABBDD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
  std::complex<double> a = abcdvec[1][0][0][0];
  std::complex<double> b = abcdvec[0][1][0][0];
  std::complex<double> d = abcdvec[0][0][1][1];
  std::complex<double> aa = abcdvec[2][0][0][0];
  std::complex<double> bb = abcdvec[0][2][0][0];
  std::complex<double> dd = abcdvec[0][0][0][2];
  std::complex<double> ab = abcdvec[1][1][0][0];
  std::complex<double> ad = abcdvec[1][0][0][1];
  std::complex<double> bd = abcdvec[0][1][0][1];
  std::complex<double> aab = abcdvec[2][1][0][0];
  std::complex<double> aad = abcdvec[2][0][0][1];
  std::complex<double> abb = abcdvec[1][2][0][0];
  std::complex<double> add = abcdvec[1][0][0][2];
  std::complex<double> abd = abcdvec[1][1][0][1];
  std::complex<double> bbd = abcdvec[0][2][0][1];
  std::complex<double> bdd = abcdvec[0][1][0][2];
  std::complex<double> aabb = abcdvec[2][2][0][0];
  std::complex<double> aadd = abcdvec[2][0][0][2];
  std::complex<double> aabd = abcdvec[2][1][0][1];
  std::complex<double> abbd = abcdvec[1][2][0][1];
  std::complex<double> abdd = abcdvec[1][1][0][2];
  std::complex<double> bbdd = abcdvec[0][2][0][2];
  std::complex<double> aabbd = abcdvec[2][2][0][1];
  std::complex<double> aabdd = abcdvec[2][1][0][2];
  std::complex<double> abbdd = abcdvec[0][0][0][2];
  std::complex<double> aabbdd  = abcdvec[2][2][0][2];
  return (-120.*aabbdd + 48.*a*abbdd + 16.*abd*abd + 24.*ab*abdd + 24.*abbd*ad +
 8.*abb*add + 48.*aabdd*b - 24.*a*abdd*b - 16.*abd*ad*b - 8.*ab*add*b -
 6.*aadd*b*b + 2.*ad*ad*b*b + 4.*a*add*b*b + 6.*aadd*bb - 2.*ad*ad*bb -
 4.*a*add*bb + 8.*aad*bbd - 8.*a*ad*bbd - 6.*a*a*bbdd + 6.*aa*bbdd +
 24.*aabd*bd - 16.*a*abd*bd - 8.*ab*ad*bd - 8.*aad*b*bd + 8.*a*ad*b*bd +
 2.*a*a*bd*bd - 2.*aa*bd*bd + 8.*aab*bdd - 8.*a*ab*bdd + 4.*a*a*b*bdd -
 4.*aa*b*bdd + 48.*aabbd*d - 24.*a*abbd*d - 16.*ab*abd*d - 8.*abb*ad*d -
 24.*aabd*b*d + 16.*a*abd*b*d + 8.*ab*ad*b*d + 4.*aad*b*b*d -
 4.*a*ad*b*b*d - 4.*aad*bb*d + 4.*a*ad*bb*d + 4.*a*a*bbd*d - 4.*aa*bbd*d -
 8.*aab*bd*d + 8.*a*ab*bd*d - 4.*a*a*b*bd*d + 4.*aa*b*bd*d - 6.*aabb*d*d +
 2.*ab*ab*d*d + 4.*a*abb*d*d + 4.*aab*b*d*d - 4.*a*ab*b*d*d +
 a*a*b*b*d*d - aa*b*b*d*d - a*a*bb*d*d + aa*bb*d*d + 6.*aabb*dd -
 2.*ab*ab*dd - 4.*a*abb*dd - 4.*aab*b*dd + 4.*a*ab*b*dd - a*a*b*b*dd +
 aa*b*b*dd + a*a*bb*dd - aa*bb*dd).real();
}
double AliAnalysisTaskGammaSoft::getStdAABBC(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
  std::complex<double> a = abcdvec[1][0][0][0];
  std::complex<double> b = abcdvec[0][1][0][0];
  std::complex<double> c = abcdvec[0][0][1][0];
  std::complex<double> aa = abcdvec[2][0][0][0];
  std::complex<double> ab = abcdvec[1][1][0][0];
  std::complex<double> ac = abcdvec[1][0][1][0];
  std::complex<double> bb = abcdvec[0][2][0][0];
  std::complex<double> bc = abcdvec[0][1][1][0];
  std::complex<double> aab = abcdvec[2][1][0][0];
  std::complex<double> aac =  abcdvec[2][0][1][0];
  std::complex<double> abb = abcdvec[1][2][0][0];
  std::complex<double> abc = abcdvec[1][1][1][0];
  std::complex<double> bbc = abcdvec[0][2][1][0];
  std::complex<double> aabb = abcdvec[2][2][0][0];
  std::complex<double> aabc = abcdvec[2][1][1][0];
  std::complex<double> abbc = abcdvec[1][2][1][0];
  std::complex<double> aabbc = abcdvec[2][2][1][0];
  return (a*a*b*b*c - aa*b*b*c - a*a*bb*c - 4.*ab*a*b*c - 2.*a*ac*b*b - 2.*a*a*bc*b
  + 2.*ab*ab*c + 4.*ab*ac*b + 4.*ab*bc*a + 8.*abc*a*b + 4.*aab*b*c + 2.*aac*b*b + 4.*abb*a*c + 2.*bbc*a*a + aa*bb*c + 2.*aa*b*bc + 2.*bb*a*ac
  - 12.*aabc*b - 12.*abbc*a - 6.*aabb*c - 8.*abc*ab - 2.*bbc*aa - 2.*aac*bb - 4.*aab*bc - 4.*abb*ac
  + 24.*aabbc).real(); }
double AliAnalysisTaskGammaSoft::getStdABCC(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
  std::complex<double> a = abcdvec[1][0][0][0];
  std::complex<double> b = abcdvec[0][1][0][0];
  std::complex<double> c = abcdvec[0][0][1][0];
  std::complex<double> ab  = abcdvec[1][1][0][0];
  std::complex<double> ac = abcdvec[1][0][1][0];
  std::complex<double> bc = abcdvec[0][1][1][0];
  std::complex<double> cc = abcdvec[0][0][2][0];
  std::complex<double> abc = abcdvec[1][1][1][0];
  std::complex<double> acc = abcdvec[1][0][2][0];
  std::complex<double> bcc = abcdvec[0][1][2][0];
  std::complex<double> abcc = abcdvec[1][1][2][0];
  return (a*b*c*c - a*b*cc - 2.*a*bc*c - 2.*ac*b*c - ab*c*c
  + 2.*acc*b + 2.*a*bcc + 4.*abc*c + ab*cc + 2.*ac*bc
  - 6.*abcc).real(); }
double AliAnalysisTaskGammaSoft::getStdABCD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
  std::complex<double> a = abcdvec[1][0][0][0];
  std::complex<double> b = abcdvec[0][1][0][0];
  std::complex<double> c = abcdvec[0][0][1][0];
  std::complex<double> d = abcdvec[0][0][0][1];
  std::complex<double> ab  = abcdvec[1][1][0][0];
  std::complex<double> ac = abcdvec[1][0][1][0];
  std::complex<double> ad = abcdvec[1][0][0][1];
  std::complex<double> bc = abcdvec[0][1][1][0];
  std::complex<double> bd = abcdvec[0][1][0][1];
  std::complex<double> cd = abcdvec[0][0][1][1];
  std::complex<double> abc = abcdvec[1][1][1][0];
  std::complex<double> abd = abcdvec[1][1][0][1];
  std::complex<double> acd = abcdvec[1][0][1][1];
  std::complex<double> bcd = abcdvec[0][1][1][1];
  std::complex<double> abcd = abcdvec[1][1][0][1];
  return (-6.*abcd + 2.*acd*b + ad*bc + 2.*a*bcd + ac*bd + 2.*abd*c - ad*b*c -
 a*bd*c + ab*cd - a*b*cd + 2.*abc*d - ac*b*d - a*bc*d - ab*c*d +
 a*b*c*d).real();
 }
double AliAnalysisTaskGammaSoft::getStdABDD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
  std::complex<double> a = abcdvec[1][0][0][0];
  std::complex<double> b = abcdvec[0][1][0][0];
  std::complex<double> d = abcdvec[0][0][0][1];
  std::complex<double> ab  = abcdvec[1][1][0][0];
  std::complex<double> ad = abcdvec[1][0][0][1];
  std::complex<double> bd = abcdvec[0][1][0][1];
  std::complex<double> dd = abcdvec[0][0][0][2];
  std::complex<double> abd = abcdvec[1][1][0][1];
  std::complex<double> add = abcdvec[1][0][0][2];
  std::complex<double> bdd = abcdvec[0][1][0][2];
  std::complex<double> abdd = abcdvec[1][1][0][2];
  return (a*b*d*d - a*b*dd - 2.*a*bd*d - 2.*ad*b*d - ab*d*d
  + 2.*add*b + 2.*a*bdd + 4.*abd*d + ab*dd + 2.*ad*bd
  - 6.*abdd).real(); }
double AliAnalysisTaskGammaSoft::getStdABC(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
std::complex<double> a = abcdvec[1][0][0][0];
std::complex<double> b = abcdvec[0][1][0][0];
std::complex<double> c = abcdvec[0][0][1][0];
std::complex<double> ab = abcdvec[1][1][0][0];
std::complex<double> ac = abcdvec[1][0][1][0];
std::complex<double> bc = abcdvec[0][1][1][0];
std::complex<double> abc = abcdvec[1][1][1][0];
  return (a*b*c - ab*c - ac*b - a*bc + 2.*abc).real(); }
double AliAnalysisTaskGammaSoft::getStdABD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec){
std::complex<double> a = abcdvec[1][0][0][0];
std::complex<double> b = abcdvec[0][1][0][0];
std::complex<double> d = abcdvec[0][0][0][1];
std::complex<double> ab = abcdvec[1][1][0][0];
std::complex<double> ad = abcdvec[1][0][0][1];
std::complex<double> bd = abcdvec[0][1][0][1];
std::complex<double> abd = abcdvec[1][1][0][1];
  return (a*b*d - ab*d - ad*b - a*bd + 2.*abd).real(); }
template <typename... args>
std::complex<double> AliAnalysisTaskGammaSoft::Q(double w, double nphi, args... wnphi){
  std::complex<double> q = w*(TMath::Cos(nphi)+1i*TMath::Sin(nphi));
  return q*Q(wnphi...);
}
std::complex<double> AliAnalysisTaskGammaSoft::Q(double w, double nphi){
  std::complex<double> q = w*(TMath::Cos(nphi)+1i*TMath::Sin(nphi));
  return q;
}
template <typename... args>
double AliAnalysisTaskGammaSoft::P(double w, double pt, args... wpt){
  return w*pt*P(wpt...);
}
double AliAnalysisTaskGammaSoft::P(double w, double pt){
  return w*pt;
}
Bool_t AliAnalysisTaskGammaSoft::LoadWeights(const Int_t &runno) { //Cannot be used when running on the trains
  if(fWeightList) {
    // fWeights[0] = (AliGFWWeights*)fWeightList->FindObject(Form("w%i",runno));
    fWeights[0] = (AliGFWWeights*)fWeightList->FindObject(Form("w%i%s",runno,fGFWSelection->GetSystPF()));
    if(!fWeights) {
      AliFatal("Weights could not be found in the list!\n");
      return kFALSE;
    };
    fWeights[0]->CreateNUA();
    return kTRUE;
  } else {
    AliFatal("Weight list (for some reason) not set!\n");
    return kFALSE;
  };
};
void AliAnalysisTaskGammaSoft::LoadCorrectionsFromLists(){
  fWeightList = (TList*)GetInputData(1);
  fWeights = new AliGFWWeights*[4];
  fEfficiencyList = (TList*)GetInputData(2); //Efficiencies start from input slot 2
  fEfficiencies = new TH1D*[fNV0MBinsDefault];
  for(Int_t i=0;i<fNV0MBinsDefault;i++) {
      fEfficiencies[i] = (TH1D*)fEfficiencyList->FindObject(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
      if(fEfficiencies[i] && fPseudoEfficiency<1) fEfficiencies[i]->Scale(fPseudoEfficiency);
      if(!fEfficiencies[i]) {
        if(!i) AliFatal("Could not fetch efficiency!\n");
        printf("Could not find efficiency for V0M bin no. %i! Cloning the previous efficiency instead...\n",i);
        fEfficiencies[i] = (TH1D*)fEfficiencies[i-1]->Clone(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
      };
  }
  return;
}
void AliAnalysisTaskGammaSoft::SetupAxes() {
 const Int_t temp_NV0MBinsDefault=fExtendV0MAcceptance?11:10;
  Double_t temp_V0MBinsDefault[12] = {0,5,10,20,30,40,50,60,70,80,90,101}; //Last bin to include V0M beyond anchor point
  if(!fV0MMultiAxis) SetV0MBins(temp_NV0MBinsDefault,temp_V0MBinsDefault);
  fV0MBinsDefault=GetBinsFromAxis(fV0MMultiAxis);
  fNV0MBinsDefault=fV0MMultiAxis->GetNbins();
  if(fV0MBinsDefault[fNV0MBinsDefault]>90) fExtendV0MAcceptance = kTRUE; //If V0M is beyond 90, then we need to extend the V0M acceptance!
  if(!fMultiAxis) SetMultiBins(fNV0MBinsDefault,fV0MBinsDefault);
  fMultiBins = GetBinsFromAxis(fMultiAxis);
  fNMultiBins = fMultiAxis->GetNbins();
  const Int_t l_NPtBinsDefault = 25;
  Double_t l_PtBinsDefault[l_NPtBinsDefault+1] = {0.20, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
                     1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90,
                     2.00, 2.20, 2.40, 2.60, 2.80, 3.00};
  if(!fPtAxis) SetPtBins(l_NPtBinsDefault,l_PtBinsDefault);
  fPtBins = GetBinsFromAxis(fPtAxis);
  fNPtBins = fPtAxis->GetNbins();
  int Neta_Default = 1;
  double l_eta_Default[] = {-0.8,0.8};
  if(!fEtaAxis) { printf("Setting default eta bins\n"); SetEtaBins(Neta_Default,l_eta_Default);}
  fEtaBins=GetBinsFromAxis(fEtaAxis);
  fNEtaBins=fEtaAxis->GetNbins();
  return;
}
void AliAnalysisTaskGammaSoft::SetPtBins(Int_t nPtBins, Double_t *PtBins) {
  if(fPtAxis) delete fPtAxis;
  fPtAxis = new TAxis(nPtBins, PtBins);
}
void AliAnalysisTaskGammaSoft::SetEtaBins(Int_t nbins, Double_t *etabins) {
  if(fEtaAxis) delete fEtaAxis;
  fEtaAxis = new TAxis(nbins,etabins);
}
void AliAnalysisTaskGammaSoft::SetMultiBins(Int_t nMultiBins, Double_t *multibins) {
  if(fMultiAxis) delete fMultiAxis;
  fMultiAxis = new TAxis(nMultiBins, multibins);
}
void AliAnalysisTaskGammaSoft::SetV0MBins(Int_t nMultiBins, Double_t *multibins) {
  if(fV0MMultiAxis) delete fV0MMultiAxis;
  fV0MMultiAxis = new TAxis(nMultiBins, multibins);
}
Double_t *AliAnalysisTaskGammaSoft::GetBinsFromAxis(TAxis *inax) {
  Int_t lBins = inax->GetNbins();
  Double_t *retBins = new Double_t[lBins+1];
  for(Int_t i=0;i<lBins;i++)
    retBins[i] = inax->GetBinLowEdge(i+1);
  retBins[lBins] = inax->GetBinUpEdge(lBins);
  return retBins;
}
Bool_t AliAnalysisTaskGammaSoft::AcceptCustomEvent(AliAODEvent* fAOD) { //From Alex
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


