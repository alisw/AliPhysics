#include "AliAnalysisTaskMeanPtV2Corr.h"
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

ClassImp(AliAnalysisTaskMeanPtV2Corr);

AliAnalysisTaskMeanPtV2Corr::AliAnalysisTaskMeanPtV2Corr():
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
  fWeightList(0),
  fWeights(0),
  fWeights_pi(0),
  fWeights_ka(0),
  fWeights_pr(0),
  fRunNo(0),
  fMidSelection(0),
  fFWSelection(0),
  fFC(0),
  fGFW(0)
{
};
AliAnalysisTaskMeanPtV2Corr::AliAnalysisTaskMeanPtV2Corr(const char *name, Bool_t IsMC, TString stageSwitch):
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
  fWeightList(0),
  fWeights(0),
  fWeights_pi(0),
  fWeights_ka(0),
  fWeights_pr(0),
  fRunNo(0),
  fMidSelection(0),
  fFWSelection(0),
  fFC(0),
  fGFW(0)
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
};
AliAnalysisTaskMeanPtV2Corr::~AliAnalysisTaskMeanPtV2Corr() {
};
void AliAnalysisTaskMeanPtV2Corr::UserCreateOutputObjects(){
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
  fMidSelection = new AliGFWCuts();
  fMidSelection->SetupCuts(0);
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
};
void AliAnalysisTaskMeanPtV2Corr::UserExec(Option_t*) {
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t l_Cent = lMultSel->GetMultiplicityPercentile("V0M");
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
};
void AliAnalysisTaskMeanPtV2Corr::Terminate(Option_t*) {
};
Bool_t AliAnalysisTaskMeanPtV2Corr::CheckTrigger(Double_t lCent) {
  fTriggerType = AliVEvent::kCentral;
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(fSelMask&AliVEvent::kMB) return kTRUE;
  if(fSelMask&AliVEvent::kINT7) return kTRUE;
  if((fSelMask&AliVEvent::kCentral) && lCent<10) return kTRUE;
  if((fSelMask&AliVEvent::kSemiCentral) && lCent<50 && lCent>30) return kTRUE;
  return kFALSE;
};
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptAOD(AliAODEvent *inEv, Double_t *lvtxXYZ) {
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
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ) {
  if(TMath::Abs(mtr->Eta())>0.8) return kFALSE;
  if(mtr->Pt()<0.5) return kFALSE;
  if(mtr->Pt()>2) return kFALSE;
  if(!mtr->TestFilterBit(96)) return kFALSE;
  if(mtr->GetTPCNclsF()<70) return kFALSE;
  if(ltrackXYZ)
    mtr->GetXYZ(ltrackXYZ);
  return kTRUE;
};
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptParticle(AliVParticle *mpa) {
  if(!mpa->IsPhysicalPrimary()) return kFALSE;
  if(mpa->Charge()==0) return kFALSE;
  if(TMath::Abs(mpa->Eta())>0.4) return kFALSE;
  if(mpa->Pt()<0.5) return kFALSE;
  if(mpa->Pt()>2) return kFALSE;
  return kTRUE;
};
Int_t AliAnalysisTaskMeanPtV2Corr::GetStageSwitch(TString instr) {
  if(instr.Contains("weights")) return 1;
  if(instr.Contains("meanpt")) return 2;
  if(instr.Contains("full")) return 3;
  return 0;
}
void AliAnalysisTaskMeanPtV2Corr::FillWeights(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
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
void AliAnalysisTaskMeanPtV2Corr::FillMeanPtCounter(Double_t pt, Double_t &l_sum, Double_t &l_count, AliGFWWeights *inWeight) {
  Double_t w = inWeight->GetIntegratedEfficiency(pt);
  if(w==0) return;
  l_sum+=pt/w;
  l_count+=1./w;
}
void AliAnalysisTaskMeanPtV2Corr::FillMeanPt(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
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
void AliAnalysisTaskMeanPtV2Corr::FillCK(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
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
Bool_t AliAnalysisTaskMeanPtV2Corr::FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndmn) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(corconf.Head.Data(),cent,val,dnx,rndmn);
    return kTRUE;
  };
  return kTRUE;
};
Bool_t AliAnalysisTaskMeanPtV2Corr::FillCovariance(AliGFW::CorrConfig corconf, Double_t cent, Double_t d_mpt, Double_t dw_mpt) {
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
void AliAnalysisTaskMeanPtV2Corr::CreateCorrConfigs() {
  corrconfigs.push_back(GetConf("MidV22","refP {2} refN {-2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV24","refP {2 2} refN {-2 -2}", kFALSE));
};
void AliAnalysisTaskMeanPtV2Corr::GetSingleWeightFromList(AliGFWWeights **inWeights, Int_t l_RunNo, TString pf) {
  if((*inWeights)) { delete (*inWeights); (*inWeights)=0; };
  (*inWeights) = (AliGFWWeights*)fWeightList->FindObject(Form("w%i%s",l_RunNo,pf.Data()));
  if(!(*inWeights)) AliFatal(Form("Could not find weight %i in weight list\n",l_RunNo));
  if(!(*inWeights)->CalculateIntegratedEff()) AliFatal("Could not calculate integrated efficiency!\n");
  (*inWeights)->CreateNUA();
};
void AliAnalysisTaskMeanPtV2Corr::LoadWeightAndMPT(AliAODEvent *inEv) {
  if(!fWeightList) AliFatal("Weight list not set!\n");
  Int_t l_RunNo = inEv->GetRunNumber();
  if(!fRunNo || fRunNo != l_RunNo) {
    GetSingleWeightFromList(&fWeights,l_RunNo);
    GetSingleWeightFromList(&fWeights_pi,l_RunNo,"_pi");
    GetSingleWeightFromList(&fWeights_ka,l_RunNo,"_ka");
    GetSingleWeightFromList(&fWeights_pr,l_RunNo,"_pr");
    // if(fWeights) delete fWeights;
    // fWeights = (AliGFWWeights*)fWeightList->FindObject(Form("w%i",l_RunNo));
    // if(!fWeights) AliFatal(Form("Could not find weight %i in weight list\n",l_RunNo));
    // if(!fWeights->CalculateIntegratedEff()) AliFatal("Could not calculate integrated efficiency!\n");
    // fWeights->CreateNUA();
    if(fStageSwitch==3) { //if on switch 3 (full), then also check if need to preload dif. weight
      if(fmPT) delete fmPT;
      fmPT = (TProfile*)fMPTList->FindObject(Form("mpt%i",l_RunNo));
      if(!fmPT) AliFatal(Form("Could not find mean pT for run %i in the list\n",l_RunNo));
    }
    fRunNo = l_RunNo;
  };
}
Bool_t AliAnalysisTaskMeanPtV2Corr::WithinSigma(Double_t SigmaCut, AliAODTrack *inTrack, AliPID::EParticleType partType) {
  if(!fPIDResponse) return kFALSE;
  Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(inTrack,partType);
  Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(inTrack,partType);
  return (TMath::Sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF) < SigmaCut);
}
