/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
*/
#include "AliAnalysisTaskGFWPIDFlow.h"
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

ClassImp(AliAnalysisTaskGFWPIDFlow);

AliAnalysisTaskGFWPIDFlow::AliAnalysisTaskGFWPIDFlow():
  AliAnalysisTaskSE(),
  fStageSwitch(0),
  fIsMC(kFALSE),
  fPIDResponse(0),
  fMPTList(0),
  fmPT(0),
  fMultiDist(0),
  fptvar(0),
  fCovariance(0),
  fTriggerType(AliVEvent::kMB),
  fUseRunAveragedWeights(kFALSE),
  fWeightList(0),
  fWeights(0),
  fWeights_pi(0),
  fWeights_ka(0),
  fWeights_pr(0),
  fRunNo(0),
  fMidSelection(0),
  fFWSelection(0),
  fFC(0),
  fGFW(0),
  fGFWMode(0),
  fWeightArray(0),
  fBayesPID(0),
  fPtAxis(0),
  fRndm(0)
{
};
AliAnalysisTaskGFWPIDFlow::AliAnalysisTaskGFWPIDFlow(const char *name, Bool_t IsMC, TString stageSwitch):
  AliAnalysisTaskSE(name),
  fStageSwitch(0),
  fIsMC(IsMC),
  fPIDResponse(0),
  fMPTList(0),
  fmPT(0),
  fMultiDist(0),
  fptvar(0),
  fCovariance(0),
  fTriggerType(AliVEvent::kMB),
  fUseRunAveragedWeights(kFALSE),
  fWeightList(0),
  fWeights(0),
  fWeights_pi(0),
  fWeights_ka(0),
  fWeights_pr(0),
  fRunNo(0),
  fMidSelection(0),
  fFWSelection(0),
  fFC(0),
  fGFW(0),
  fGFWMode(0),
  fWeightArray(0),
  fBayesPID(0),
  fPtAxis(0),
  fRndm(0)
{
  fStageSwitch = GetStageSwitch(stageSwitch);
  if(!fStageSwitch) AliFatal("Stage switch is 0, not sure what should be done!\n");
  if(fStageSwitch==1)
    DefineOutput(1,TList::Class());
  if(fStageSwitch==2) {
    DefineInput(1,TList::Class());
    DefineOutput(1,TList::Class());
    DefineOutput(2,TH1D::Class());
  };
  if(fStageSwitch==3) {
    DefineInput(1,TList::Class());
    DefineInput(2,TList::Class());
    DefineOutput(1,TProfile::Class());
    DefineOutput(2,AliGFWFlowContainer::Class());
    DefineOutput(3,TProfile::Class());
    DefineOutput(4,TProfile::Class());
  }
  if(fStageSwitch==4) {
    DefineInput(1,TList::Class());
    DefineInput(2,TList::Class());
    DefineOutput(1,AliGFWFlowContainer::Class());
  }
  if(fStageSwitch==5) {
    DefineOutput(1,TList::Class());
  }

};
AliAnalysisTaskGFWPIDFlow::~AliAnalysisTaskGFWPIDFlow() {
  if(fStageSwitch==4 || fStageSwitch==5) {
    delete [] fWeightArray;
  }
};
void AliAnalysisTaskGFWPIDFlow::UserCreateOutputObjects(){
  OpenFile(1);
  const Int_t nMultiBins = 200;
  Double_t lMultiBins[nMultiBins+1];
  for(Int_t i=0;i<=nMultiBins;i++) lMultiBins[i] = i*10;
  if(fStageSwitch==1) {
    fWeightList = new TList();
    fWeightList->SetOwner(kTRUE);
      const Int_t NbinsPtForV2=31;
      Double_t binsPtForV2[NbinsPtForV2+1] = {
      0.3, 0.4, 0.5, 0.6,
      0.7, 0.8, 0.9, 1.0, 1.25,
      1.5, 1.75, 2., 2.25, 2.5,
      2.75, 3.0, 3.25, 3.50, 3.75,
      4.0, 4.5, 5.0, 5.5, 6.0,
      7.0, 8.0, 9.0, 10.0, 12.0,
      14.0, 16.0, 20.0};
      fWeights = new AliGFWWeights();
      fWeights->SetPtBins(NbinsPtForV2,binsPtForV2);
      fWeights->SetName("SomeWeight");
      fWeights->Init(kFALSE,kTRUE);
      fWeightList->Add(fWeights);
      fWeights = new AliGFWWeights();
      fWeights->SetPtBins(NbinsPtForV2,binsPtForV2);
      fWeights->SetName("SomeWeight_pi");
      fWeights->Init(kFALSE,kTRUE);
      fWeightList->Add(fWeights);
      fWeights = new AliGFWWeights();
      fWeights->SetPtBins(NbinsPtForV2,binsPtForV2);
      fWeights->SetName("SomeWeight_ka");
      fWeights->Init(kFALSE,kTRUE);
      fWeightList->Add(fWeights);
      fWeights = new AliGFWWeights();
      fWeights->SetPtBins(NbinsPtForV2,binsPtForV2);
      fWeights->SetName("SomeWeight_pr");
      fWeights->Init(kFALSE,kTRUE);
      fWeightList->Add(fWeights);
      PostData(1,fWeightList);
  };
  if(fStageSwitch==2) {
    fWeightList = (TList*)GetInputData(1);
    // fWeights = (AliGFWWeights*)GetInputData(1);
    // if(!fWeights) AliFatal("Could not fetch input weights!\n");
    // if(!fWeights->CalculateIntegratedEff()) AliFatal("Could not calculate integrated efficiency!\n");
    fMPTList = new TList();
    fMPTList->SetOwner(kTRUE);
    fmPT = new TProfile("MeanPt","MeanPt",nMultiBins,lMultiBins);
    fmPT_pi = new TProfile("MeanPt_pi","MeanPt_pi",nMultiBins,lMultiBins);
    fmPT_ka = new TProfile("MeanPt_ka","MeanPt_ka",nMultiBins,lMultiBins);
    fmPT_pr = new TProfile("MeanPt_pr","MeanPt_pr",nMultiBins,lMultiBins);
    fMPTList->Add(fmPT);
    fMPTList->Add(fmPT_pi);
    fMPTList->Add(fmPT_ka);
    fMPTList->Add(fmPT_pr);
    PostData(1,fMPTList);
    fMultiDist = new TH1D("MultiDistribution","Multiplicity distribution; #it{N}_{ch}; N(events)",nMultiBins,lMultiBins);
    PostData(2,fMultiDist);
  };
  if(fStageSwitch==3) {
    fWeightList = (TList*)GetInputData(1);
    if(!fWeightList) AliFatal("Could not fetch weight list!\n");
    fMPTList = (TList*)GetInputData(2);
    if(!fMPTList) AliFatal("Could not fetch input mean pT list!\n");
    fptvar = new TProfile("varpt","varpt",nMultiBins,lMultiBins);
    PostData(1,fptvar);
    //Setting up the FlowContainer
    TObjArray *oba = new TObjArray();
    oba->Add(new TNamed("MidV22","MidV22"));
    oba->Add(new TNamed("MidV24","MidV24"));
    fFC = new AliGFWFlowContainer();
    fFC->SetName("FlowContainer");
    fFC->Initialize(oba,nMultiBins,lMultiBins);
    delete oba;
    PostData(2,fFC);
    //Initializing GFW
    Int_t pows[] = {3,0,2,0,3};
    fGFW = new AliGFW();
    fGFW->AddRegion("refN",5,pows,-0.8,-0.4,1,1);
    fGFW->AddRegion("refP",5,pows,0.4,0.8,1,1);
    CreateCorrConfigs();
    //Covariance
    fCovariance = new TProfile("cov","Covariance",nMultiBins,lMultiBins);
    PostData(3,fCovariance);
  }
  if(fStageSwitch==4) {
    const Int_t NbinsPtForV2=28;
    Double_t binsPtForV2[NbinsPtForV2+1] = {
    0.2, 0.3, 0.4, 0.5, 0.6,
    0.7, 0.8, 0.9, 1.0, 1.25,
    1.5, 1.75, 2., 2.25, 2.5,
    2.75, 3.0, 3.25, 3.50, 3.75,
    4.0, 4.5, 5.0, 5.5, 6.0,
    7.0, 8.0, 9.0, 10.0};
    fPtAxis = new TAxis(NbinsPtForV2,binsPtForV2);
    if(fGFWMode==0 || fGFWMode==1) fWeightList = (TList*)GetInputData(1);
    if(fGFWMode==2 || fGFWMode==3) fWeightList = (TList*)GetInputData(2);// Load ZM weights instead, if required
    if(!fWeightList) AliFatal("Could not fetch weight list!\n");
    if(fGFWMode==2 || fGFWMode==3) LoadZMWeights();
    TObjArray *oba = new TObjArray();
    AddToOBA(oba,"MidV22");
    AddToOBA(oba,"MidV24");
    AddToOBA(oba,"MidV26");
    AddToOBA(oba,"MidV28");
    AddToOBA(oba,"MidVCh22",NbinsPtForV2);
    AddToOBA(oba,"MidVCh24",NbinsPtForV2);
    AddToOBA(oba,"MidVCh26",NbinsPtForV2);
    AddToOBA(oba,"MidVCh28",NbinsPtForV2);
    AddToOBA(oba,"MidVPi22",NbinsPtForV2);
    AddToOBA(oba,"MidVPi24",NbinsPtForV2);
    AddToOBA(oba,"MidVPi26",NbinsPtForV2);
    AddToOBA(oba,"MidVPi28",NbinsPtForV2);
    AddToOBA(oba,"MidVKa22",NbinsPtForV2);
    AddToOBA(oba,"MidVKa24",NbinsPtForV2);
    AddToOBA(oba,"MidVKa26",NbinsPtForV2);
    AddToOBA(oba,"MidVKa28",NbinsPtForV2);
    AddToOBA(oba,"MidVPr22",NbinsPtForV2);
    AddToOBA(oba,"MidVPr24",NbinsPtForV2);
    AddToOBA(oba,"MidVPr26",NbinsPtForV2);
    AddToOBA(oba,"MidVPr28",NbinsPtForV2);

    fFC = new AliGFWFlowContainer();
    fFC->SetName(Form("FlowContainer_%i",fGFWMode));
    fFC->SetXAxis(fPtAxis);
    Double_t l_MultiBins[] = {5,10,20,30,40,50,60};
    fFC->Initialize(oba,6,l_MultiBins,10);
    delete oba;
    PostData(1,fFC);
    //Initializing GFW

    Int_t pows[] = {9,0,8,0,7,0,6,0,5};
    fGFW = new AliGFW();
    fGFW->AddRegion("ref",9,pows,-0.8,0.8,1,1);
    fGFW->AddRegion("poiCh",9,pows,-0.8,0.8,fPtAxis->GetNbins()+1,2);
    fGFW->AddRegion("poiPi",9,pows,-0.8,0.8,fPtAxis->GetNbins()+1,4);
    fGFW->AddRegion("poiKa",9,pows,-0.8,0.8,fPtAxis->GetNbins()+1,8);
    fGFW->AddRegion("poiPr",9,pows,-0.8,0.8,fPtAxis->GetNbins()+1,16);

    //Overlap
    fGFW->AddRegion("olCh",9,pows,-0.8,0.8,fPtAxis->GetNbins()+1,32);
    fGFW->AddRegion("olPi",9,pows,-0.8,0.8,fPtAxis->GetNbins()+1,64);
    fGFW->AddRegion("olKa",9,pows,-0.8,0.8,fPtAxis->GetNbins()+1,128);
    fGFW->AddRegion("olPr",9,pows,-0.8,0.8,fPtAxis->GetNbins()+1,256);

    // Int_t pows[] = {5,0,4,0,3};
    // fGFW = new AliGFW();
    // fGFW->AddRegion("ref",5,pows,-0.8,0.8,1,1);
    // fGFW->AddRegion("poiCh",5,pows,-0.8,0.8,fPtAxis->GetNbins()+1,2);
    // fGFW->AddRegion("poiPi",5,pows,-0.8,0.8,fPtAxis->GetNbins()+1,4);
    // fGFW->AddRegion("poiKa",5,pows,-0.8,0.8,fPtAxis->GetNbins()+1,8);
    // fGFW->AddRegion("poiPr",5,pows,-0.8,0.8,fPtAxis->GetNbins()+1,16);

    CreateCorrConfigs();

    fBayesPID = new AliPIDCombined();
    fBayesPID->SetDefaultTPCPriors();
    fBayesPID->SetSelectedSpecies(AliPID::kSPECIES);
    fBayesPID->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask

    fRndm = new TRandom(0);
  }
  if(fStageSwitch==5) {
    fWeightList = new TList();
    fWeightList->SetOwner(kTRUE);
    fWeightList->Add(new TH2D("Refs","ChIncl",16,-0.8,0.8,100,0,TMath::TwoPi()));
    fWeightList->Add(new TH2D("Charged","ChExcl",16,-0.8,0.8,100,0,TMath::TwoPi()));
    fWeightList->Add(new TH2D("Pion","ChPi",16,-0.8,0.8,100,0,TMath::TwoPi()));
    fWeightList->Add(new TH2D("Kaon","Kacl",16,-0.8,0.8,100,0,TMath::TwoPi()));
    fWeightList->Add(new TH2D("Proton","PrIncl",16,-0.8,0.8,100,0,TMath::TwoPi()));
    PostData(1,fWeightList);
    fBayesPID = new AliPIDCombined();
    fBayesPID->SetDefaultTPCPriors();
    fBayesPID->SetSelectedSpecies(AliPID::kSPECIES);
    fBayesPID->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask
  }

  fMidSelection = new AliGFWCuts();
  fMidSelection->SetupCuts(0);
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
};
void AliAnalysisTaskGFWPIDFlow::UserExec(Option_t*) {
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t l_Cent = lMultSel->GetMultiplicityPercentile("V0M");
  if(l_Cent<5 || l_Cent>70) return; //Lets only consider 5-70%
  if(!CheckTrigger(l_Cent)) return;
  Double_t vtxXYZ[] = {0.,0.,0.};
  if(!AcceptAOD(fAOD, vtxXYZ)) return;
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  if(fStageSwitch==1)
    FillWeights(fAOD, vz,l_Cent);
  if(fStageSwitch==2)
    FillMeanPt(fAOD, vz, l_Cent);
  if(fStageSwitch==3)
    FillCK(fAOD,vz,l_Cent);
  if(fStageSwitch==4)
    DevFunction(fAOD,vz,l_Cent);
  if(fStageSwitch==5)
    FillCustomWeights(fAOD,vz,l_Cent);
};
void AliAnalysisTaskGFWPIDFlow::Terminate(Option_t*) {
};
Bool_t AliAnalysisTaskGFWPIDFlow::CheckTrigger(Double_t lCent) {
  fTriggerType = AliVEvent::kCentral;
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(fSelMask&AliVEvent::kMB) return kTRUE;
  if(fSelMask&AliVEvent::kINT7) return kTRUE;
  if((fSelMask&AliVEvent::kCentral) && lCent<10) return kTRUE;
  if((fSelMask&AliVEvent::kSemiCentral) && lCent<50 && lCent>30) return kTRUE;
  return kFALSE;
};
Bool_t AliAnalysisTaskGFWPIDFlow::AcceptAOD(AliAODEvent *inEv, Double_t *lvtxXYZ) {
  if(!fEventCuts.AcceptEvent(inEv)) return 0;
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
Bool_t AliAnalysisTaskGFWPIDFlow::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ) {
  if(TMath::Abs(mtr->Eta())>0.8) return kFALSE;
  if(mtr->Pt()<0.5) return kFALSE;
  if(mtr->Pt()>2) return kFALSE;
  if(!mtr->TestFilterBit(96)) return kFALSE;
  if(mtr->GetTPCNclsF()<70) return kFALSE;
  if(ltrackXYZ)
    mtr->GetXYZ(ltrackXYZ);
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWPIDFlow::AcceptParticle(AliVParticle *mpa) {
  if(!mpa->IsPhysicalPrimary()) return kFALSE;
  if(mpa->Charge()==0) return kFALSE;
  if(TMath::Abs(mpa->Eta())>0.4) return kFALSE;
  if(mpa->Pt()<0.5) return kFALSE;
  if(mpa->Pt()>2) return kFALSE;
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWPIDFlow::devAcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ) {
  if(TMath::Abs(mtr->Eta())>0.8) return kFALSE;
  if(mtr->Pt()<0.2) return kFALSE;
  if(mtr->Pt()>10) return kFALSE;
  if(!mtr->TestFilterBit(96)) return kFALSE;
  if(mtr->GetTPCNclsF()<70) return kFALSE;
  if(ltrackXYZ)
    mtr->GetXYZ(ltrackXYZ);
  return kTRUE;
};

Int_t AliAnalysisTaskGFWPIDFlow::GetStageSwitch(TString instr) {
  if(instr.Contains("weights")) return 1;
  if(instr.Contains("meanpt")) return 2;
  if(instr.Contains("full")) return 3;
  if(instr.Contains("dev")) return 4;
  if(instr.Contains("CustomWeights")) return 5;
  return 0;
}
void AliAnalysisTaskGFWPIDFlow::FillWeights(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  //MC generated
  AliVParticle *lPart;
  AliAODTrack *lTrack;
  Double_t trackXYZ[3];
  Double_t dummyDouble[] = {0.,0.};
  TClonesArray *tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
  for(Int_t i=0;i<tca->GetEntries();i++) {
    lPart = (AliAODMCParticle*)tca->At(i);
    if(!AcceptParticle(lPart)) continue;
    if(!fMidSelection->AcceptParticle(lPart,0)) continue;
    ((AliGFWWeights*)fWeightList->At(0))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    Int_t pdgCode = TMath::Abs(lPart->PdgCode());
    if(pdgCode==211) ((AliGFWWeights*)fWeightList->At(1))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    if(pdgCode==321) ((AliGFWWeights*)fWeightList->At(2))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);
    if(pdgCode==2212) ((AliGFWWeights*)fWeightList->At(3))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,2);

  };
  //MC reconstructed
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    lPart = (AliAODMCParticle*)tca->At(TMath::Abs(lTrack->GetLabel()));
    if(!AcceptAODTrack(lTrack,trackXYZ)) continue;
    if(!fMidSelection->AcceptTrack(lTrack,dummyDouble)) continue;
    ((AliGFWWeights*)fWeightList->At(0))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
    if(WithinSigma(2.5,lTrack,AliPID::kPion)) ((AliGFWWeights*)fWeightList->At(1))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
    if(WithinSigma(2.5,lTrack,AliPID::kKaon)) ((AliGFWWeights*)fWeightList->At(2))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
    if(WithinSigma(2.5,lTrack,AliPID::kProton)) ((AliGFWWeights*)fWeightList->At(3))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);

  };
  PostData(1,fWeightList);
}
void AliAnalysisTaskGFWPIDFlow::FillMeanPtCounter(Double_t pt, Double_t &l_sum, Double_t &l_count, AliGFWWeights *inWeight) {
  Double_t w = inWeight->GetIntegratedEfficiency(pt);
  if(w==0) return;
  l_sum+=pt/w;
  l_count+=1./w;
}
void AliAnalysisTaskGFWPIDFlow::FillMeanPt(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  LoadWeightAndMPT(fAOD);
  AliAODTrack *lTrack;
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!AcceptAODTrack(lTrack,trackXYZ)) continue;
    if(TMath::Abs(lTrack->Eta())>0.4) continue; //for mean pt, only consider -0.4-0.4 region
    Double_t lpt = lTrack->Pt();
    FillMeanPtCounter(lpt,l_ptsum[0],l_ptCount[0],fWeights);

    // Double_t w = fWeights->GetIntegratedEfficiency(lpt);
    // if(w==0) continue;
    // l_ptsum[0]+=lpt/w;
    // l_ptCount[0]+=1./w;
    if(WithinSigma(2.5,lTrack,AliPID::kPion)) FillMeanPtCounter(lpt,l_ptsum[1],l_ptCount[1],fWeights_pi);
    if(WithinSigma(2.5,lTrack,AliPID::kKaon)) FillMeanPtCounter(lpt,l_ptsum[2],l_ptCount[2],fWeights_ka);
    if(WithinSigma(2.5,lTrack,AliPID::kProton)) FillMeanPtCounter(lpt,l_ptsum[3],l_ptCount[3],fWeights_pr);
  };
  if(l_ptCount[0]==0) return;
  l_ptsum[0]=l_ptsum[0]/l_ptCount[0];
  fmPT->Fill(l_ptCount[0],l_ptsum[0],l_ptCount[0]);
  if(l_ptCount[1]!=0) fmPT_pi->Fill(l_ptCount[0],l_ptsum[1]/l_ptCount[1],l_ptCount[1]);
  if(l_ptCount[2]!=0) fmPT_ka->Fill(l_ptCount[0],l_ptsum[2]/l_ptCount[2],l_ptCount[2]);
  if(l_ptCount[3]!=0) fmPT_pr->Fill(l_ptCount[0],l_ptsum[3]/l_ptCount[3],l_ptCount[3]);
  PostData(1,fMPTList);
  fMultiDist->Fill(l_ptCount[0]);
  PostData(2,fMultiDist);
};
void AliAnalysisTaskGFWPIDFlow::FillCK(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  LoadWeightAndMPT(fAOD);
  AliAODTrack *lTrack;
  Double_t w1p0=0;
  Double_t w1p1=0;
  Double_t w2p2=0;
  Double_t w2p1=0;
  Double_t w2p0=0;
  Double_t trackXYZ[3];
  fGFW->Clear();
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!AcceptAODTrack(lTrack,trackXYZ)) continue;
    Double_t p1 = lTrack->Pt();
    Double_t weff = fWeights->GetIntegratedEfficiency(p1);
    Double_t wacc = fWeights->GetNUA(lTrack->Phi(),lTrack->Eta(),vz);
    if(TMath::Abs(lTrack->Eta())<0.4)  { //for mean pt, only consider -0.4-0.4 region
      if(weff==0) continue;
      Double_t w = 1./weff;
      w1p0 += w;
      w1p1 += w*p1;
      w2p2 += w*w*p1*p1;
      w2p1 += w*w*p1;
      w2p0 += w*w;
    } else { //Otherwise, we consider it for vn calculations
      fGFW->Fill(lTrack->Eta(),1,lTrack->Phi(),wacc,1);
    };
  };
  if(w1p0==0) return;
  Double_t l_meanPt = fmPT->GetBinContent(fmPT->FindBin(w1p0)); //l_Cent should be replaced with w1p0 here (->weighted Nch)
  Double_t l_val = (w1p1 - l_meanPt*w1p0) * (w1p1 - l_meanPt*w1p0)
                   - w2p2 + 2*l_meanPt*w2p1 - l_meanPt*l_meanPt*w2p0;
  Double_t l_norm= w1p0*w1p0-w2p0;
  Double_t mpt_local = w1p1/w1p0;
  if(l_norm!=0) {
    fptvar->Fill(w1p0,l_val/l_norm, l_norm);
    PostData(1,fptvar);
  };
  //Filling FCs
  for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
    Bool_t filled = FillFCs(corrconfigs.at(l_ind),w1p0,0);
  };
  PostData(2,fFC);
  FillCovariance(corrconfigs.at(0),w1p0,mpt_local-l_meanPt,w1p0);
  PostData(3,fCovariance);
}
void AliAnalysisTaskGFWPIDFlow::DevFunction(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  if(fGFWMode==0 || fGFWMode==1) LoadMyWeights(fAOD); //My mode, my weights, run-by-run
  AliAODTrack *lTrack;
  Double_t trackXYZ[3];
  fGFW->Clear();
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!devAcceptAODTrack(lTrack,trackXYZ)) continue;
    Double_t p1 = lTrack->Pt();
    Int_t ptind = fPtAxis->FindBin(p1)-1;
    Double_t l_eta = lTrack->Eta();
    Double_t l_phi = lTrack->Phi();
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    // Int_t spIndex = PIDIndex+1;
    Double_t ptmins[] = {0.2,0.2,0.3,0.5};
    Double_t ptmaxs[] = {10.,10.,6.0,6.0};
    Bool_t WithinRef=(p1>0.2 && p1<5);
    Bool_t WithinPOI=(p1>ptmins[PIDIndex] && p1<ptmaxs[PIDIndex]);
    Bool_t WithinNch=(p1>ptmins[0] && p1<ptmaxs[0]); //Within Ncharged (important for e.g. protons)
    if(!WithinRef && !WithinPOI) continue;
    Bool_t REFnotPOI = (WithinRef && !WithinPOI); // Particles that could be POI, but fall out of pT range
    if(fGFWMode==0 || fGFWMode==1) { //My modes, to be tested against morphed weights
      Double_t wRef = GetMyWeight(l_eta,l_phi,0);
      Double_t wPOI = WithinPOI?GetMyWeight(l_eta,l_phi,PIDIndex+1):GetMyWeight(l_eta,l_phi,1); //If not within POI (e.g. low pT protons), then use Nch weights
      if(fGFWMode == 0) { //This should be the default mode
        if(WithinRef && WithinPOI)
          if(PIDIndex) wRef = wPOI; //All there is. If particle is both (then it's overlap), override ref with POI
      }
      if(fGFWMode == 1) {//This is to check whether we can use Nch weights for ref particles. Same as before, but no PIDIndex check
        if(WithinRef && WithinPOI)
          wRef = wPOI;
      }
      if(WithinRef) fGFW->Fill(l_eta,ptind,l_phi,wRef,1); //Filling in ref flow
      if(WithinPOI && PIDIndex) fGFW->Fill(l_eta,ptind,l_phi,wPOI,(1<<(PIDIndex+1))); //Filling POI flow for ID'ed
      if(WithinNch) fGFW->Fill(l_eta,ptind,l_phi,wPOI,2); //Filling POI flow for ID'ed
      //Filling overlaps:
      if(WithinPOI && PIDIndex && WithinRef) fGFW->Fill(l_eta,ptind,l_phi,wPOI,1<<(PIDIndex+5)); //Filling POI flow for ID'ed
      if(WithinNch && WithinRef) fGFW->Fill(l_eta,ptind,l_phi,wPOI,32); //Filling POI flow for ID'ed
    };
    if(fGFWMode==2) { //Approach #2 with morphed weights
      Double_t wRef = GetZMWeight(l_eta,l_phi,0);
      Double_t wPOI = GetZMWeight(l_eta,l_phi,PIDIndex+1);
      Double_t wCha = GetZMWeight(l_eta,l_phi,1);//Need NCh weight anyways
      if(WithinRef) fGFW->Fill(l_eta,0,l_phi,wRef,1);
      if(WithinPOI) {
        fGFW->Fill(l_eta,ptind,l_phi,wCha,2); //Fill all charged
        if(PIDIndex) fGFW->Fill(l_eta,ptind,l_phi,wPOI,(1<<(PIDIndex+1))); //Explicitly treating PID to avoid double-counting for Nch
      };
      if(WithinRef&&WithinPOI) {
        if(PIDIndex) fGFW->Fill(l_eta,ptind,l_phi,wPOI,(1<<(PIDIndex+5)),wRef); //Explicitly treat PID to avoid double-counting for Nch
        fGFW->Fill(l_eta,ptind,l_phi,wCha,32,wRef); //Filling all charged
      }
    }
    if(fGFWMode==3) { //The old method
      Double_t wRef = GetZMWeight(l_eta,l_phi,0);
      Double_t wPOI = GetZMWeight(l_eta,l_phi,PIDIndex+1);
      Double_t wCha = GetZMWeight(l_eta,l_phi,1);//Need NCh weight anyways
      if(WithinRef) fGFW->Fill(l_eta,0,l_phi,wRef,1);
      if(WithinPOI) {
        fGFW->Fill(l_eta,ptind,l_phi,wCha,2); //Fill all charged
        if(PIDIndex) fGFW->Fill(l_eta,ptind,l_phi,wPOI,(1<<(PIDIndex+1))); //Explicitly treating PID to avoid double-counting for Nch
      };
      if(WithinRef&&WithinPOI) {
        if(PIDIndex) fGFW->Fill(l_eta,ptind,l_phi,wPOI,(1<<(PIDIndex+5))); //Explicitly treat PID to avoid double-counting for Nch
        fGFW->Fill(l_eta,ptind,l_phi,wCha,32); //Filling all charged
      }
    }

    /* trying to figure out whether it's a POI or a ref or both.
    This is in particular stupid, because e.g. a proton with pT < 500 MeV should go with... ref flow weight?
    ------------------------------
    Under construction for now. But first need to quickly commit a weight production
    */

    // Double_t wref = GetMyWeight(l_eta,l_phi,0); //ref weight
    // Double_t wpoi = GetMyWeight(l_eta,l_phi,PIDIndex+1); //poi weight
    //This is stupid. If e.g.

    // Double_t wref = GetMyWeight(l_eta,l_phi,0);
    // fGFW->Fill(lTrack->Eta(),ptind,lTrack->Phi(),wref,3); // masks 1+2 (ref + ncharged)
    // Int_t PIDIndex = GetBayesPIDIndex(lTrack);
    // if(PIDIndex>=0) {
    //   wref = GetMyWeight(l_eta,l_phi,1+PIDIndex);
    //   fGFW->Fill(l_eta,ptind,l_phi,wref,(1<<(PIDIndex+2)));
    //   // if(PIDIndex==2) printf("Filling a proton with eta (%f), ptInd (%i), phi (%f), weight (%f) and mask (%i)\n",
    //   //                         l_eta,ptind,l_phi,wref,(1<<(PIDIndex+2)));
    // }
  };
  Double_t rndm = fRndm->Rndm();
  for(Int_t i=0;i<corrconfigs.size();i++)  Bool_t dm = FillFCs(corrconfigs.at(i),l_Cent,rndm);
  PostData(1,fFC);
}
void AliAnalysisTaskGFWPIDFlow::FillCustomWeights(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  AliAODTrack *lTrack;
  Double_t trackXYZ[3];
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!devAcceptAODTrack(lTrack,trackXYZ)) continue;
    Double_t p1 = lTrack->Pt();
    Double_t l_eta = lTrack->Eta();
    Double_t l_phi = lTrack->Phi();
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1; //-1 for unID, 0 for pi, etc.; shift everything by 1
    Double_t ptmins[] = {0.2,0.2,0.3,0.5};
    Double_t ptmaxs[] = {10.,10.,6.0,6.0};
    Bool_t WithinRef=(p1>0.2 && p1<5);
    Bool_t WithinPOI=(p1>ptmins[PIDIndex] && p1<ptmaxs[PIDIndex]);
    if(WithinPOI) ((TH2D*)fWeightList->At(PIDIndex+1))->Fill(l_eta, l_phi);
    if(WithinRef && !PIDIndex) ((TH2D*)fWeightList->At(0))->Fill(l_eta, l_phi); //This becomes a pt-subset of Nch
  };
  PostData(1,fWeightList);
}

Bool_t AliAnalysisTaskGFWPIDFlow::GetIntValAndDNX(AliGFW::CorrConfig corconf, Double_t &l_val, Double_t &l_dnx) {
  l_dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(l_dnx==0) return kFALSE;
  l_val = fGFW->Calculate(corconf,0,kFALSE).Re()/l_dnx;
  if(TMath::Abs(l_val)>1) return kFALSE;
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWPIDFlow::FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndmn, Bool_t EnableDebug) {
  Double_t dnx, val;
  if(!corconf.pTDif) {
    dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
    if(dnx==0) return kFALSE;
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(corconf.Head.Data(),cent,val,dnx,rndmn);
    return kTRUE;
  } else {
    for(Int_t i=1; i<=fPtAxis->GetNbins();i++) {
      dnx = fGFW->Calculate(corconf,i-1,kTRUE).Re();
      if(dnx==0) continue;
      val = fGFW->Calculate(corconf,i-1,kFALSE).Re()/dnx;
      if(EnableDebug) printf("dnx cut passed. Dnx = %f\t val = %f\n",dnx,val);
      if(TMath::Abs(val)<1) {
        fFC->FillProfile(Form("%s_pt_%i",corconf.Head.Data(),i),cent,val,dnx,rndmn);
      if(EnableDebug) printf("Just filled %s with %f and %f\n",Form("%s_pt_%i",corconf.Head.Data(),i),val,dnx);
      }
    }
  }
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWPIDFlow::FillCovariance(AliGFW::CorrConfig corconf, Double_t cent, Double_t d_mpt, Double_t dw_mpt) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      fCovariance->Fill(cent,val*d_mpt,dnx*dw_mpt);
    return kTRUE;
  };
  return kTRUE;
};
void AliAnalysisTaskGFWPIDFlow::CreateCorrConfigs() {
  // corrconfigs.push_back(GetConf("MidV22","refP {2} refN {-2}", kFALSE));
  // corrconfigs.push_back(GetConf("MidV24","refP {2 2} refN {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV22","ref {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV24","ref {2 2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV26","ref {2 2 2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV28","ref {2 2 2 2 -2 -2 -2 -2}", kFALSE));
  //2-part corr
  corrconfigs.push_back(GetConf("MidVCh22","ref {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVPi22","ref {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVKa22","ref {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVPr22","ref {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVCh22","poiCh ref | olCh {2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVPi22","poiPi ref | olPi {2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVKa22","poiKa ref | olKa {2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVPr22","poiPr ref | olPr {2 -2}", kTRUE));
  //4-part corr
  corrconfigs.push_back(GetConf("MidVCh24","ref {2 2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVPi24","ref {2 2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVKa24","ref {2 2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVPr24","ref {2 2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVCh24","poiCh ref | olCh{2 2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVPi24","poiPi ref | olPi{2 2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVKa24","poiKa ref | olKa{2 2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVPr24","poiPr ref | olPr{2 2 -2 -2}", kTRUE));
  //6-part corr
  corrconfigs.push_back(GetConf("MidVCh26","ref {2 2 2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVPi26","ref {2 2 2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVKa26","ref {2 2 2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVPr26","ref {2 2 2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVCh26","poiCh ref | olCh {2 2 2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVPi26","poiPi ref | olPi {2 2 2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVKa26","poiKa ref | olKa {2 2 2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVPr26","poiPr ref | olPr {2 2 2 -2 -2 -2}", kTRUE));
  //8-part corr
  corrconfigs.push_back(GetConf("MidVCh28","ref {2 2 2 2 -2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVPi28","ref {2 2 2 2 -2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVKa28","ref {2 2 2 2 -2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVPr28","ref {2 2 2 2 -2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidVCh28","poiCh ref | olCh {2 2 2 2 -2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVPi28","poiPi ref | olPi {2 2 2 2 -2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVKa28","poiKa ref | olKa {2 2 2 2 -2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidVPr28","poiPr ref | olPr {2 2 2 2 -2 -2 -2 -2}", kTRUE));

};
void AliAnalysisTaskGFWPIDFlow::GetSingleWeightFromList(AliGFWWeights **inWeights, Int_t l_RunNo, TString pf) {
  if((*inWeights)) { delete (*inWeights); (*inWeights)=0; };
  (*inWeights) = (AliGFWWeights*)fWeightList->FindObject(Form("w%i%s",l_RunNo,pf.Data()));
  if(!(*inWeights)) AliFatal(Form("Could not find weight %i in weight list\n",l_RunNo));
  if(!(*inWeights)->CalculateIntegratedEff()) AliFatal("Could not calculate integrated efficiency!\n");
  (*inWeights)->CreateNUA();
};
void AliAnalysisTaskGFWPIDFlow::LoadWeightAndMPT(AliAODEvent *inEv) {
  if(!fWeightList) AliFatal("Weight list not set!\n");
  Int_t l_RunNo = inEv->GetRunNumber();
  if(!fRunNo || fRunNo != l_RunNo) {
    GetSingleWeightFromList(&fWeights,l_RunNo);
    GetSingleWeightFromList(&fWeights_pi,l_RunNo,"_pi");
    GetSingleWeightFromList(&fWeights_ka,l_RunNo,"_ka");
    GetSingleWeightFromList(&fWeights_pr,l_RunNo,"_pr");
    if(fStageSwitch==3) { //if on switch 3 (full), then also check if need to preload dif. weight
      if(fmPT) delete fmPT;
      fmPT = (TProfile*)fMPTList->FindObject(Form("mpt%i",l_RunNo));
      if(!fmPT) AliFatal(Form("Could not find mean pT for run %i in the list\n",l_RunNo));
    }
    fRunNo = l_RunNo;
  };
}
void AliAnalysisTaskGFWPIDFlow::LoadMyWeights(AliAODEvent* mev) {
  Int_t lRunNo = mev->GetRunNumber(); //246048 for testing purposes
  if(fRunNo == lRunNo) return;
  fRunNo = lRunNo;
  if(!fWeightList) AliFatal("Weight list not set!\n");
  TString wNames[] = {"Refs","Charged","Pion","Kaon","Proton"};
  if(!fWeightArray) fWeightArray = new TH2D*[5];
  for(Int_t i=0;i<5;i++) {
    wNames[i].Prepend(Form("w%i_",fRunNo));
    fWeightArray[i] = (TH2D*)fWeightList->FindObject(wNames[i].Data());
    if(!fWeightArray[i]) AliFatal(Form("Could not get %s weights!\n",wNames[i].Data()));
  };
}
void AliAnalysisTaskGFWPIDFlow::LoadZMWeights() {
  if(!fWeightList) AliFatal("Weight list not set!\n");
  TString wNames[] = {"Refs","Charged","Pion","Kaon","Proton"};
  if(!fWeightArray) fWeightArray = new TH2D*[5];
  for(Int_t i=0;i<5;i++) {
    fWeightArray[i] = (TH2D*)fWeightList->FindObject(wNames[i].Data());
    if(!fWeightArray[i]) AliFatal(Form("Could not get %s weights!\n",wNames[i].Data()));
  };
}
Bool_t AliAnalysisTaskGFWPIDFlow::WithinSigma(Double_t SigmaCut, AliAODTrack *inTrack, AliPID::EParticleType partType) {
  if(!fPIDResponse) return kFALSE;
  Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(inTrack,partType);
  Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(inTrack,partType);
  return (TMath::Sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF) < SigmaCut);
}
Double_t AliAnalysisTaskGFWPIDFlow::GetMyWeight(Double_t eta, Double_t phi, Int_t pidind) {
  Int_t etaind = fWeightArray[pidind]->GetXaxis()->FindBin(eta);
  Int_t phiind = fWeightArray[pidind]->GetYaxis()->FindBin(phi);
  return fWeightArray[pidind]->GetBinContent(etaind,phiind);
}
Double_t AliAnalysisTaskGFWPIDFlow::GetZMWeight(Double_t eta, Double_t phi, Int_t pidind) {
  Int_t phiind = fWeightArray[pidind]->GetXaxis()->FindBin(phi);
  Int_t etaind = fWeightArray[pidind]->GetYaxis()->FindBin(eta);
  return fWeightArray[pidind]->GetBinContent(phiind,etaind);
};
Bool_t AliAnalysisTaskGFWPIDFlow::HasTPCPID(AliAODTrack* l_track) {
  if(!l_track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus l_TPCExists = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, l_track);
  return (l_TPCExists == AliPIDResponse::kDetPidOk);
}
// ============================================================================
Bool_t AliAnalysisTaskGFWPIDFlow::HasTOFPID(AliAODTrack* l_track) {
  if(!l_track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus l_TOFExists = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, l_track);
  return ((l_TOFExists == AliPIDResponse::kDetPidOk) && (l_track->GetStatus()& AliVTrack::kTOFout) && (l_track->GetStatus()& AliVTrack::kTIME));
}
Int_t AliAnalysisTaskGFWPIDFlow::GetBayesPIDIndex(AliAODTrack *l_track) {
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
void AliAnalysisTaskGFWPIDFlow::AddToOBA(TObjArray *oba, TString l_name, Int_t nPT) {
  oba->Add(new TNamed(l_name.Data(),l_name.Data()));
  TString skel(l_name);
  TString skelTitle(l_name);
  skel.Append("_pt_%i");
  skelTitle.Append("_pTDiff");
  if(nPT>0) for(Int_t i=0;i<nPT;i++) oba->Add(new TNamed(Form(skel.Data(),i+1),skelTitle.Data()));
}
