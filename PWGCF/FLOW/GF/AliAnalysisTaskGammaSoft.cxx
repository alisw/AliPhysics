#include "AliAnalysisTaskGammaSoft.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVEvent.h"
#include "AliAODVertex.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TNamed.h"
#include "AliMCEvent.h"
#include "AliExternalTrackParam.h"

ClassImp(AliAnalysisTaskGammaSoft);

AliAnalysisTaskGammaSoft::AliAnalysisTaskGammaSoft():
  AliAnalysisTaskSE(),
  fSystFlag(0),
  fEventCutFlag(0),
  fContSubfix(0),
  fCentEst(0),
  fIsMC(kFALSE),
  fBypassTriggerAndEventCuts(kFALSE),
  fDisablePileup(kFALSE),
  fUseOldPileup(kFALSE),
  fDCAxyFunctionalForm(0),
  fUseRecoNchForMC(kFALSE),
  fRndm(0),
  fNBootstrapProfiles(10),
  fPtAxis(0),
  fEtaAxis(0),
  fMultiAxis(0),
  fV0MMultiAxis(0),
  fDCAAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fEtaBins(0),
  fNEtaBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fDCABins(0),
  fNDCABins(0),
  fV0MBinsDefault(0),
  fNV0MBinsDefault(0),
  fUseNch(kFALSE),
  fUseUnityParticleWeights(kFALSE),
  fUseUnityEventWeights(kFALSE),
  fEtaMpt(0.4),
  fEtaAcceptance(0.8),
  fEtaV2Sep(0.4),
  fQAList(0),
  fEventCount(0),
  fMultiDist(0),
  fMultiVsV0MCorr(0),
  fNchTrueVsReco(0),
  fESDvsFB128(0),
  fPtList(0),
  fCovList(0),
  fv24deltapt2_v24mpt2(0),
  fv24deltapt2_v24mpt(0),
  fv24(0),
  fdeltapt2_mpt2(0),
  fdeltapt2_mpt(0),
  fmmpt(0),
  fv2deltapt2_v2mpt2(0),
  fv2deltapt2_v2mpt(0),
  fv2(0),
  fv2deltapt_v2mpt(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fWeightList(0),
  fWeights(0),
  fRunNo(0),
  fGFWSelection(0),
  fFC(0),
  fGFW(0),
  corrconfigs(0),
  fEfficiencyList(0),
  fEfficiency(0),
  fEfficiencies(0),
  fPseudoEfficiency(2),
  fDCAxyVsPt_noChi2(0),
  fWithinDCAvsPt_withChi2(0),
  fDCAxyVsPt_withChi2(0),
  fWithinDCAvsPt_noChi2(0),
  fV0MMulti(0),
  fITSvsTPCMulti(0),
  fSPDCutPU(0),
  fV0CutPU(0),
  fCenCutLowPU(0),
  fCenCutHighPU(0),
  fMultCutPU(0),
  fCentralPU(0),
  fPhiEtaVz(0),
  fptDCAxyDCAz(0),
  fChi2TPCcls(0),
  fTPCcls(0),
  fTPCcrsrws(0),
  fhQAEventsfMult32vsCentr(0),
  fhQAEventsMult128vsCentr(0),
  fhQAEventsfMultTPCvsTOF(0),
  fhQAEventsfMultTPCvsESD(0),
  fEventWeight(0),
  wp(0)
{
};
AliAnalysisTaskGammaSoft::AliAnalysisTaskGammaSoft(const char *name, Bool_t IsMC, TString ContSubfix):
  AliAnalysisTaskSE(name),
  fSystFlag(0),
  fEventCutFlag(0),
  fContSubfix(0),
  fCentEst(0),
  fIsMC(kFALSE),
  fBypassTriggerAndEventCuts(kFALSE),
  fDisablePileup(kFALSE),
  fUseOldPileup(kFALSE),
  fDCAxyFunctionalForm(0),
  fUseRecoNchForMC(kFALSE),
  fRndm(0),
  fNBootstrapProfiles(10),
  fPtAxis(0),
  fEtaAxis(0),
  fMultiAxis(0),
  fV0MMultiAxis(0),
  fDCAAxis(0),
  fPtBins(0),
  fNPtBins(0),
  fEtaBins(0),
  fNEtaBins(0),
  fMultiBins(0),
  fNMultiBins(0),
  fDCABins(0),
  fNDCABins(0),
  fV0MBinsDefault(0),
  fNV0MBinsDefault(0),
  fUseNch(kFALSE),
  fUseUnityParticleWeights(kFALSE),
  fUseUnityEventWeights(kFALSE),
  fEtaMpt(0.4),
  fEtaAcceptance(0.8),
  fEtaV2Sep(0.4),
  fQAList(0),
  fEventCount(0),
  fMultiDist(0),
  fMultiVsV0MCorr(0),
  fNchTrueVsReco(0),
  fESDvsFB128(0),
  fPtList(0),
  fCovList(0),
  fv24deltapt2_v24mpt2(0),
  fv24deltapt2_v24mpt(0),
  fv24(0),
  fdeltapt2_mpt2(0),
  fdeltapt2_mpt(0),
  fmmpt(0),
  fv2deltapt2_v2mpt2(0),
  fv2deltapt2_v2mpt(0),
  fv2(0),
  fv2deltapt_v2mpt(0),
  fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
  fWeightList(0),
  fWeights(0),
  fRunNo(0),
  fGFWSelection(0),
  fFC(0),
  fGFW(0),
  corrconfigs(0),
  fEfficiencyList(0),
  fEfficiency(0),
  fEfficiencies(0),
  fPseudoEfficiency(2),
  fDCAxyVsPt_noChi2(0),
  fWithinDCAvsPt_withChi2(0),
  fDCAxyVsPt_withChi2(0),
  fWithinDCAvsPt_noChi2(0),
  fV0MMulti(0),
  fITSvsTPCMulti(0),
  fSPDCutPU(0),
  fV0CutPU(0),
  fCenCutLowPU(0),
  fCenCutHighPU(0),
  fMultCutPU(0),
  fCentralPU(0),
  fPhiEtaVz(0),
  fptDCAxyDCAz(0),
  fChi2TPCcls(0),
  fTPCcls(0),
  fTPCcrsrws(0),
  fhQAEventsfMult32vsCentr(0),
  fhQAEventsMult128vsCentr(0),
  fhQAEventsfMultTPCvsTOF(0),
  fhQAEventsfMultTPCvsESD(0),
  fEventWeight(0),
  wp(0)
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

    fRndm = new TRandom(0);
    fRequireReloadOnRunChange = kFALSE;
    if(!fIsMC) LoadCorrectionsFromLists(); //Efficiencies and NUA are only for the data or if specified for pseudoefficiencies
    printf("Creating pt-correlation objects\n");
    fPtList = new TList();
    fPtList->SetOwner(kTRUE);

    printf("pt-correlation objects created\n");
    printf("Creating multiplicity objects\n");
    fMultiDist = new TH1D("MultiDistribution","Multiplicity distribution; #it{N}_{ch}; N(events)",fNMultiBins,fMultiBins);
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",fNV0MBinsDefault,fV0MBinsDefault);
    fPtList->Add(fMultiDist);
    fPtList->Add(fV0MMulti);
    fMultiVsV0MCorr = new TH2D*[2];
    fMultiVsV0MCorr[0] = new TH2D("MultVsV0M_BeforeConsistency","MultVsV0M_BeforeConsistency",103,0,103,fNMultiBins,fMultiBins[0],fMultiBins[fNMultiBins]);
    fMultiVsV0MCorr[1] = new TH2D("MultVsV0M_AfterConsistency","MultVsV0M_AfterConsistency",103,0,103,fNMultiBins,fMultiBins[0],fMultiBins[fNMultiBins]);
    fESDvsFB128 = new TH2D("ESDvsFB128","; N(FB128); N(ESD)",500,-0.5,4999.5,1500,-0.5,14999.5);
    fPtList->Add(fMultiVsV0MCorr[0]);
    fPtList->Add(fMultiVsV0MCorr[1]);
    fPtList->Add(fESDvsFB128);
    //ITS vs TPC tracklets cut for PU
    fITSvsTPCMulti = new TH2I("TPCvsITSclusters",";TPC clusters; ITS clusters",1000,0,10000,5000,0,50000);
    fPtList->Add(fITSvsTPCMulti);
    if(fIsMC) {
      fNchTrueVsReco = new TH2D("NchTrueVsReco",";Nch (MC-true); Nch (MC-reco)",fNMultiBins,fMultiBins,fNMultiBins,fMultiBins);
      fPtList->Add(fNchTrueVsReco);
    }
    printf("Multiplicity objects created\n");
    PostData(1,fPtList);

    //Setting up the FlowContainer
    printf("Creating flow container\n");
    TObjArray *oba = new TObjArray();
    oba->Add(new TNamed("ChGap22","ChGap22")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap24","ChGap24")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChFull22","ChFull22")); //no-gap case
    oba->Add(new TNamed("ChFull24","ChFull24")); //no-gap case

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

    fv24deltapt2_v24mpt2 = new AliProfileBS("fv24deltapt2_v24mpt2","fv24deltapt2_v24mpt2",fNMultiBins,fMultiBins);
    fv24deltapt2_v24mpt2->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fv24deltapt2_v24mpt2);
    fv24deltapt2_v24mpt = new AliProfileBS("fv24deltapt2_v24mpt","fv24deltapt2_v24mpt",fNMultiBins,fMultiBins);
    fv24deltapt2_v24mpt->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fv24deltapt2_v24mpt);

    fv24 = new AliProfileBS("fv24","fv24",fNMultiBins,fMultiBins);
    fv24->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fv24);

    fdeltapt2_mpt2 = new AliProfileBS("fdeltapt2_mpt2","fdeltapt2_mpt2",fNMultiBins,fMultiBins);
    fdeltapt2_mpt2->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fdeltapt2_mpt2);
    fdeltapt2_mpt = new AliProfileBS("fdeltapt2_mpt","fdeltapt2_mpt",fNMultiBins,fMultiBins);
    fdeltapt2_mpt->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fdeltapt2_mpt);

    fmmpt = new AliProfileBS("fmmpt","fmmpt",fNMultiBins,fMultiBins);
    fmmpt->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fmmpt);

    fv2deltapt2_v2mpt2 = new AliProfileBS("fv2deltapt2_v2mpt2","fv2deltapt2_v2mpt2",fNMultiBins,fMultiBins);
    fv2deltapt2_v2mpt2->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fv2deltapt2_v2mpt2);
    fv2deltapt2_v2mpt = new AliProfileBS("fv2deltapt2_v2mpt","fv2deltapt2_v2mpt",fNMultiBins,fMultiBins);
    fv2deltapt2_v2mpt->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fv2deltapt2_v2mpt);
    fv2 = new AliProfileBS("fv2","fv2",fNMultiBins,fMultiBins);
    fv2->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fv2);

    fv2deltapt_v2mpt = new AliProfileBS("fv2deltapt_v2mpt","fv2deltapt_v2mpt",fNMultiBins,fMultiBins);
    fv2deltapt_v2mpt->InitializeSubsamples(fNBootstrapProfiles);
    fCovList->Add(fv2deltapt_v2mpt);

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

    fPhiEtaVz = new TH3D*[2];
    fChi2TPCcls = new TH1D*[2];
    fTPCcls = new TH1D*[2];
    fTPCcrsrws = new TH1D*[2];
    fptDCAxyDCAz = new TH3D*[2];
    TString str_cut[] = {"beforeCuts","afterCuts"};
    for(int i(0);i<2;++i){
      fPhiEtaVz[i] = new TH3D(Form("hPhiEtaVz_%s",str_cut[i].Data()),Form("#phi,#eta,v_{z} %s;#varphi;#eta;v_{z};Counts",str_cut[i].Data()),60,0,TMath::TwoPi(),64,-1.6,1.6,40,-10,10);
      fQAList->Add(fPhiEtaVz[i]);
      fptDCAxyDCAz[i] = new TH3D(Form("fptDCAxyDCA_%s",str_cut[i].Data()),Form("#it{p}_{T},DCA_{xy},DCA_{z} %s;#it{p_{T}};DCA_{xy};DCA_{z}",str_cut[i].Data()),fNPtBins,fPtBins,fNDCABins,fDCABins,fNDCABins,fDCABins);
      fQAList->Add(fptDCAxyDCAz[i]);
      fChi2TPCcls[i] = new TH1D(Form("chi2prTPCcls_%s",str_cut[i].Data()),Form("Chi2TPCcls %s;#chi^{2} pr. TPC cluster;Counts",str_cut[i].Data()),100,0,6);
      fQAList->Add(fChi2TPCcls[i]);
      fTPCcls[i] = new TH1D(Form("TPCcls_%s",str_cut[i].Data()),Form("TPCcls %s;#TPC cluster;Counts",str_cut[i].Data()),159,0,159);
      fQAList->Add(fTPCcls[i]);
      fTPCcrsrws[i] = new TH1D(Form("TPCcrsrws_%s",str_cut[i].Data()),Form("TPCcrsrws %s;#TPC crossed rows;Counts",str_cut[i].Data()),159,0,159);
      fQAList->Add(fTPCcrsrws[i]);
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

    fEventCuts.OverrideAutomaticTriggerSelection(fTriggerType,true);
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
};
void AliAnalysisTaskGammaSoft::UserExec(Option_t*) {
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

  AliAODTrack *lTrack;
  wp.clear(); wp.resize(3,vector<double>(3));
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
    if(fUseRecoNchForMC) nTotNoTracksReco = GetNtotTracks(fAOD,ptMin,ptMax,vtxXYZ);
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
      if(TMath::Abs(leta) > fEtaAcceptance) continue;
      Double_t pt = lPart->Pt();
      if(pt<ptMin || pt>ptMax) continue;
      if(TMath::Abs(leta)<fEtaMpt)  {
        FillWPCounter(wp,1,pt);
      }
      fGFW->Fill(leta,1,lPart->Phi(),1,3); //filling both gap (bit mask 1) and full (bit maks 2). Since this is MC, weight is 1.
    };
    nTotNoTracks = fUseRecoNchForMC?nTotNoTracksReco:nTotNoTracksMC;
    if(fUseRecoNchForMC) fNchTrueVsReco->Fill(nTotNoTracksMC,nTotNoTracksReco);
  } else {
    Bool_t usingPseudoEff = (fPseudoEfficiency<1);
    for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
      if(usingPseudoEff) if(fRndm->Uniform()>fPseudoEfficiency) continue;
      lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
      if(!lTrack) continue;
      Double_t leta = lTrack->Eta();
      Double_t trackXYZ[] = {0.,0.,0.};
      double dcaxyz[2];
      DCAxyz(lTrack,fAOD,dcaxyz);
      fptDCAxyDCAz[0]->Fill(lTrack->Pt(),dcaxyz[0],dcaxyz[1]);
      fChi2TPCcls[0]->Fill(lTrack->GetTPCchi2perCluster());
      fTPCcls[0]->Fill(lTrack->GetTPCncls());
      fTPCcrsrws[0]->Fill(lTrack->GetTPCNCrossedRows());
      fPhiEtaVz[0]->Fill(lTrack->Phi(),leta,vz);
      //Counting FB128 for QA:
      if(lTrack->TestFilterBit(128)) nTotTracksFB128++;
      if(!AcceptAODTrack(lTrack,trackXYZ,ptMin,ptMax,vtxXYZ)) continue;
      nTotNoTracks++;
      if(fEtaV2Sep>0 && leta < -fEtaV2Sep) lNegCount++;
      if(fEtaV2Sep>0 && leta > fEtaV2Sep) lPosCount++;
      Double_t lpt = lTrack->Pt();
      Double_t weff = fEfficiencies[iCent]->GetBinContent(fEfficiencies[iCent]->FindBin(lpt));
      if(weff==0.0) continue;
      weff = 1./weff;
      if(TMath::Abs(leta<fEtaMpt)) {
        lMidCount++;
        FillWPCounter(wp,(fUseUnityParticleWeights)?1.0:weff,lpt);
      }
      Double_t wacc = fWeights[0]->GetNUA(lTrack->Phi(),leta,vz);
      fGFW->Fill(leta,1,lTrack->Phi(),(fUseUnityParticleWeights)?1.:wacc*weff,3); //filling both gap (bit mask 1) and full (bit mask 2)

      fptDCAxyDCAz[1]->Fill(lTrack->Pt(),dcaxyz[0],dcaxyz[1]);
      fChi2TPCcls[1]->Fill(lTrack->GetTPCchi2perCluster());
      fTPCcls[1]->Fill(lTrack->GetTPCncls());
      fTPCcrsrws[1]->Fill(lTrack->GetTPCNCrossedRows());
      fPhiEtaVz[1]->Fill(lTrack->Phi(),leta,vz);
    };
  };
  if(wp[1][0]==0) return; //if no single charged particles, then surely no PID either, no sense to continue
  fEventCount->Fill("Tracks",1);
  fMultiVsV0MCorr[0]->Fill(l_Cent,nTotNoTracks);
  //here in principle one could use the GFW output to check if the values are calculated, but this is more efficient
  if(fConsistencyFlag&1) if(!lPosCount || !lNegCount) return; // only events where v2{2, gap} could be calculated
  if(fConsistencyFlag&2) if(nTotNoTracks<4) return; //only events where v2{4} can be calculated (assuming same region as nch)
  if(fConsistencyFlag&4) if(lPosCount < 2 || lNegCount < 2) return; //only events where v2{4, gap} could be calculated
  if(fConsistencyFlag&8) if(lMidCount < 2) return;//only events where variance of pt can be calculated in mid
  fMultiVsV0MCorr[1]->Fill(l_Cent,nTotNoTracks);

  fITSvsTPCMulti->Fill(fAOD->GetNumberOfITSClusters(0)+fAOD->GetNumberOfITSClusters(1),fAOD->GetNumberOfTPCClusters());
  //Filling pT variance
  Double_t l_Multi = fUseNch?(1.0*nTotNoTracks):l_Cent;
  //A check in case l_Multi is completely off the charts (in MC, sometimes it ends up being... -Xe-310???)
  if(fUseNch && l_Multi<1) return;
  //Fetching number of ESD tracks -> for QA. Only after all the events are/were rejected
  AliAODHeader *head = (AliAODHeader*)fAOD->GetHeader();
  Int_t nESD = head->GetNumberOfESDTracks();
  fESDvsFB128->Fill(nTotTracksFB128,nESD);
  Double_t l_Random = fRndm->Rndm();
  fV0MMulti->Fill(l_Cent);
  fMultiDist->Fill(l_Multi);
  PostData(1,fPtList);
  //Filling FCs
  for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
    FillFCs(corrconfigs.at(l_ind),l_Multi,l_Random);
  };
  PostData(2,fFC);

  double wpt2 = wp[1][0]*wp[1][0] - wp[2][0];
  FillCovariance(fv24deltapt2_v24mpt2,corrconfigs.at(1),l_Multi,(wp[1][1]*wp[1][1] - wp[2][2])/wpt2,1,l_Random);
  FillCovariance(fv24deltapt2_v24mpt,corrconfigs.at(1),l_Multi,(wp[1][1]*wp[1][0]-wp[2][1])/wpt2,1,l_Random);

  fdeltapt2_mpt2->FillProfile(l_Multi,(wp[1][1]*wp[1][1] - wp[2][2])/wpt2,1,l_Random);
  fdeltapt2_mpt->FillProfile(l_Multi,(wp[1][1]*wp[1][0]-wp[2][1])/wpt2,1,l_Random);
  fmmpt->FillProfile(l_Multi,wp[1][1]/wp[1][0],1,l_Random);

  FillCovariance(fv2deltapt2_v2mpt2,corrconfigs.at(0),l_Multi,(wp[1][1]*wp[1][1] - wp[2][2])/wpt2,1,l_Random);
  FillCovariance(fv2deltapt2_v2mpt,corrconfigs.at(0),l_Multi,(wp[1][1]*wp[1][0]-wp[2][1])/wpt2,1,l_Random);
  FillCovariance(fv2deltapt_v2mpt,corrconfigs.at(0),l_Multi,wp[1][1]/wp[1][0],1,l_Random);

  FillCovariance(fv24,corrconfigs.at(0),l_Multi,1,1,l_Random);
  FillCovariance(fv2,corrconfigs.at(0),l_Multi,1,1,l_Random);

  PostData(3,fCovList);
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
  if(TMath::Abs(mtr->Eta()) > fEtaAcceptance) return kFALSE;
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
void AliAnalysisTaskGammaSoft::FillWPCounter(vector<vector<double>> &inarr, double w, double p)
{
  for(int i=0;i<=2;++i)
  {
    for(int j=0;j<=2;++j)
    {
      inarr[i][j] += pow(w,i)*pow(p,j);
    }
  }
  return;
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
      fFC->FillProfile(corconf.Head.c_str(),cent,val,(fUseUnityEventWeights)?1.0:dnx,rndmn);
    return kTRUE;
  };
  return kTRUE;
};
Bool_t AliAnalysisTaskGammaSoft::FillCovariance(AliProfileBS *target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt, const Double_t &l_rndm) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).real();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).real()/dnx;
    if(TMath::Abs(val)<1)
      target->FillProfile(cent,val*d_mpt,(fUseUnityEventWeights)?1.0:dnx*dw_mpt,l_rndm);
    return kTRUE;
  };
  return kTRUE;
};
void AliAnalysisTaskGammaSoft::CreateCorrConfigs() {

  corrconfigs.push_back(GetConf("ChGap22","refP {2} refN {-2}", kFALSE));     //ChGap22 0
  corrconfigs.push_back(GetConf("ChGap24","refP {2 2} refN {-2 -2}", kFALSE));  //ChGap24 1
  corrconfigs.push_back(GetConf("ChFull22","mid {2 -2}", kFALSE));  //ChFull22 2
  corrconfigs.push_back(GetConf("ChFull24","mid {2 2 -2 -2}", kFALSE));   //ChFull24 3

  return;
};
void AliAnalysisTaskGammaSoft::DCAxyz(const AliAODTrack *track, const AliVEvent *evt, Double_t (&dcaxyz)[2])
{
  //standard "error" values:
  dcaxyz[0] = -9999.;
  dcaxyz[1] = -9999.;

  if(!track) return;

  // Create an external parameter from the AODtrack
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);

  Double_t covar[3]={0.,0.,0.};
  if(!etp.PropagateToDCA(evt->GetPrimaryVertex(),evt->GetMagneticField(),10.,dcaxyz,covar))
    {
      dcaxyz[0] = -9999.;
      dcaxyz[1] = -9999.;
      return;
    }
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
 const Int_t temp_NV0MBinsDefault= 10 ;
  Double_t temp_V0MBinsDefault[12] = {0,5,10,20,30,40,50,60,70,80,90}; //Last bin to include V0M beyond anchor point
  if(!fV0MMultiAxis) SetV0MBins(temp_NV0MBinsDefault,temp_V0MBinsDefault);
  fV0MBinsDefault=GetBinsFromAxis(fV0MMultiAxis);
  fNV0MBinsDefault=fV0MMultiAxis->GetNbins();
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
  if(!fDCAAxis) SetDCABins(400,-2,2);
  fDCABins = GetBinsFromAxis(fDCAAxis);
  fNDCABins = fDCAAxis->GetNbins();
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
void AliAnalysisTaskGammaSoft::SetDCABins(Int_t nDCABins, Double_t *DCAbins) {
  if(fDCAAxis) delete fDCAAxis;
  fDCAAxis = new TAxis(nDCABins, DCAbins);
}
void AliAnalysisTaskGammaSoft::SetDCABins(Int_t nDCABins, Double_t low, Double_t high) {
  if(fDCAAxis) delete fDCAAxis;
  fDCAAxis = new TAxis(nDCABins, low, high);
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