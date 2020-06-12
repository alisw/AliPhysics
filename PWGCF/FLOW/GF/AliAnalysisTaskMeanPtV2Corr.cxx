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
#include "AliPIDCombined.h"

ClassImp(AliAnalysisTaskMeanPtV2Corr);

AliAnalysisTaskMeanPtV2Corr::AliAnalysisTaskMeanPtV2Corr():
  AliAnalysisTaskSE(),
  fStageSwitch(0),
  fIsMC(kFALSE),
  fPIDResponse(0),
  fBayesPID(0),
  fMPTList(0),
  fmPT(0),
  fMultiDist(0),
  fNchVsMulti(0),
  fptVarList(0),
  fptvar(0),
  fCovList(0),
  fCovariance(0),
  fmptSet(kFALSE),
  fTriggerType(AliVEvent::kMB),
  fWeightList(0),
  fWeights(0),
  fNUAList(0),
  fNUAHist(0),
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
  fBayesPID(0),
  fMPTList(0),
  fmPT(0),
  fmptSet(kFALSE),
  fMultiDist(0),
  fNchVsMulti(0),
  fptVarList(0),
  fptvar(0),
  fCovList(0),
  fCovariance(0),
  fTriggerType(AliVEvent::kMB),
  fWeightList(0),
  fWeights(0),
  fNUAList(0),
  fNUAHist(0),
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
    DefineInput(1,TList::Class()); //NUE weights; ultimately, should be combined with NUA, but don't want to rerun now
    DefineInput(2,TList::Class()); //Mean Pt, should be rerun with Bayes PID
    DefineInput(3,TList::Class()); //NUA weights from other analysis; quickfix
    DefineOutput(1,TList::Class());
    DefineOutput(2,AliGFWFlowContainer::Class());
    DefineOutput(3,TList::Class());
  }
  if(fStageSwitch==4) {
    DefineOutput(1,TList::Class());
  }
  if(fStageSwitch==5) {
    DefineInput(1,TList::Class());
    DefineOutput(1,TList::Class());
  }
};
AliAnalysisTaskMeanPtV2Corr::~AliAnalysisTaskMeanPtV2Corr() {
};
void AliAnalysisTaskMeanPtV2Corr::UserCreateOutputObjects(){
  OpenFile(1);
  const Int_t nMultiBins = 200;
  Double_t lMultiBins[nMultiBins+1];
  for(Int_t i=0;i<=nMultiBins;i++) lMultiBins[i] = i*10;
  TString spNames[] = {"ch","pi","ka","pr"};
  printf("Stage switch is %i\n\n\n",fStageSwitch);
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
      TString wNames[] = {"ch","pi","ka","pr"};
      fWeights = new AliGFWWeights*[4];
      for(Int_t i=0; i<4;i++) {
        fWeights[i] = new AliGFWWeights();
        fWeights[i]->SetPtBins(NbinsPtForV2,binsPtForV2);
        fWeights[i]->SetName(Form("weight_%s",wNames[i].Data()));
        fWeights[i]->Init(kFALSE,kTRUE);
        fWeightList->Add(fWeights[i]);
      }
      PostData(1,fWeightList);
  };
  if(fStageSwitch==2) {
    fWeightList = (TList*)GetInputData(1);
    fWeights = new AliGFWWeights*[4];
    fMPTList = new TList();
    fMPTList->SetOwner(kTRUE);
    fmPT = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fmPT[i] = new TProfile(Form("MeanPt_%s",spNames[i].Data()),Form("MeanPt_%s",spNames[i].Data()),nMultiBins,lMultiBins);
      fMPTList->Add(fmPT[i]);
    }
    PostData(1,fMPTList);
    fMultiDist = new TH1D("MultiDistribution","Multiplicity distribution; #it{N}_{ch}; N(events)",nMultiBins,lMultiBins);
    PostData(2,fMultiDist);
  };
  if(fStageSwitch==3) {
    fWeightList = (TList*)GetInputData(1);
    if(!fWeightList) AliFatal("Could not fetch weight list!\n");
    fMPTList = (TList*)GetInputData(2);
    if(!fMPTList) AliFatal("Could not fetch input mean pT list!\n");
    fNUAList = (TList*)GetInputData(3);
    if(!LoadMyWeights(0)) return; //Loading run-avg NUA weights
    fptVarList = new TList();
    fptVarList->SetOwner(kTRUE);
    fptvar = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fptVarList->Add(new TProfile(Form("varpt_%s",spNames[i].Data()),Form("varpt_%s",spNames[i].Data()),nMultiBins,lMultiBins));
      fptvar[i] = (TProfile*)fptVarList->At(i);
    }
    PostData(1,fptVarList);
    //Setting up the FlowContainer
    TObjArray *oba = new TObjArray();
    oba->Add(new TNamed("ChPos22","ChPos22"));
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
    oba->Add(new TNamed("PrNeg24","PrNeg24"));
    fFC = new AliGFWFlowContainer();
    fFC->SetName("FlowContainer");
    fFC->Initialize(oba,nMultiBins,lMultiBins);
    delete oba;
    PostData(2,fFC);
    //Initializing GFW
    Int_t pows[] = {3,0,2,0,3};
    Int_t powsPOI[] = {3,0,2,0,3};
    fGFW = new AliGFW();
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
    fGFW->AddRegion("OLprP",5,pows,0.4,0.8,1,256);
    CreateCorrConfigs();
    //Covariance
    fCovList = new TList();
    fCovList->SetOwner(kTRUE);
    fCovariance = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fCovList->Add(new TProfile(Form("cov_%s",spNames[i].Data()),Form("cov_%s",spNames[i].Data()),nMultiBins,lMultiBins));
      fCovariance[i] = (TProfile*)fCovList->At(i);
    };
    PostData(3,fCovList);
  }
  if(fStageSwitch==4) {
    fMPTList = new TList();
    fMPTList->SetOwner(kTRUE);
    fmPT = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fmPT[i] = new TProfile(Form("MeanPt_%s",spNames[i].Data()),Form("MeanPt_%s",spNames[i].Data()),nMultiBins,lMultiBins);
      fMPTList->Add(fmPT[i]);
    }
    Double_t lV0Mbins[] = {0,5,10,20,30,40,50,60,70,80,90};
    fNchVsMulti = new TProfile("nChVsMulti","nChVsMulti",10,lV0Mbins);
    fMPTList->Add(fNchVsMulti);
    PostData(1,fMPTList);
  };
  if(fStageSwitch==5) {
    fMPTList = (TList*)GetInputData(1);
    if(!fMPTList) AliFatal("Could not fetch input mean pT list!\n");
    fmPT = new TProfile*[4];
    for(Int_t i=0;i<4;i++)
      fmPT[i] = (TProfile*)fMPTList->At(i);
    fptVarList = new TList();
    fptVarList->SetOwner(kTRUE);
    fptvar = new TProfile*[4];
    for(Int_t i=0;i<4;i++) {
      fptVarList->Add(new TProfile(Form("ptvar_%s",spNames[i].Data()),Form("ptvar_%s",spNames[i].Data()),nMultiBins,lMultiBins));
      fptvar[i] = (TProfile*)fptVarList->At(i);
    };
    PostData(1,fptVarList);
  };

  fMidSelection = new AliGFWCuts();
  fMidSelection->SetupCuts(0);
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
  if(fStageSwitch==4)
    ProduceALICEPublished_MptProd(fAOD,vz,l_Cent);
  if(fStageSwitch==5)
    ProduceALICEPublished_CovProd(fAOD,vz,l_Cent);
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
Bool_t AliAnalysisTaskMeanPtV2Corr::AcceptAODTrackALICEPublished(AliAODTrack *mtr, Double_t *ltrackXYZ) {
  if(TMath::Abs(mtr->Eta())>0.8) return kFALSE;
  if(mtr->Pt()<0.15) return kFALSE;
  if(mtr->Pt()>2) return kFALSE;
  if(!mtr->TestFilterBit(128)) return kFALSE;
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
  if(instr.Contains("ALICEMpt")) return 4;
  if(instr.Contains("ALICECov")) return 5;
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
    if(!AcceptAODTrack(lTrack,trackXYZ)) continue;
    if(!fMidSelection->AcceptTrack(lTrack,dummyDouble)) continue;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    fWeights[0]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
    if(PIDIndex) fWeights[PIDIndex]->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),l_Cent,1);
  };
  PostData(1,fWeightList);
}
void AliAnalysisTaskMeanPtV2Corr::FillMeanPtCounter(Double_t pt, Double_t &l_sum, Double_t &l_count, AliGFWWeights *inWeight) {
  Double_t w = inWeight?inWeight->GetIntegratedEfficiency(pt):1;
  if(w==0) return;
  l_sum+=pt/w;
  l_count+=1./w;
}
void AliAnalysisTaskMeanPtV2Corr::FillMeanPt(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  // LoadWeightAndMPT(fAOD);
  AliAODTrack *lTrack;
  Double_t l_ptsum[]={0,0,0,0};
  Double_t l_ptCount[]={0,0,0,0};
  Double_t trackXYZ[3];
  Int_t nTotNoTracks=0;
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!AcceptAODTrack(lTrack,trackXYZ)) continue;
    if(TMath::Abs(lTrack->Eta())<0.8 && lTrack->Pt()>0.2 && lTrack->Pt()<3)  nTotNoTracks++;
    if(TMath::Abs(lTrack->Eta())>0.4) continue; //for mean pt, only consider -0.4-0.4 region
    Double_t lpt = lTrack->Pt();
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    FillMeanPtCounter(lpt,l_ptsum[0],l_ptCount[0],fWeights[0]);
    if(PIDIndex) FillMeanPtCounter(lpt,l_ptsum[PIDIndex],l_ptCount[PIDIndex],fWeights[PIDIndex]);
  };
  if(l_ptCount[0]==0) return;
  for(Int_t i=0;i<4;i++) {
    if(!l_ptCount[i]) continue;
    fmPT[i]->Fill(nTotNoTracks,l_ptsum[i]/l_ptCount[i],l_ptCount[i]);
  }
  PostData(1,fMPTList);
  fMultiDist->Fill(l_ptCount[0]);
  PostData(2,fMultiDist);
};
void AliAnalysisTaskMeanPtV2Corr::FillWPCounter(Double_t inArr[5], Double_t w, Double_t p) {
  inArr[0] += w;       // = w1p0
  inArr[1] += w*p;     // = w1p1
  inArr[2] += w*w*p*p; // = w2p2
  inArr[3] += w*w*p;   // = w2p1
  inArr[4] += w*w;     // = w2p0
}
void AliAnalysisTaskMeanPtV2Corr::CalculateMptValues(Double_t outArr[4], Double_t inArr[5]) {
  //Assuming outArr[0] is preset to meanPt; outArr[1] = variance; outArr[2] = norm; outArr[3] = mpt in this event
  outArr[1] = TMath::Power(inArr[1] - outArr[0]*inArr[0], 2) //(w1p1 - l_meanPt*w1p0) * (w1p1 - l_meanPt*w1p0)
              - inArr[2] + 2*outArr[0]*inArr[3] - outArr[0]*outArr[0]*inArr[4]; //- w2p2 + 2*l_meanPt*w2p1 - l_meanPt*l_meanPt*w2p0;
  outArr[2] = inArr[0]*inArr[0] - inArr[4];
  outArr[3] = inArr[1]/inArr[0];
}
void AliAnalysisTaskMeanPtV2Corr::FillCK(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
  AliAODTrack *lTrack;
  Double_t wp[4][5] = {{0,0,0,0,0}, {0,0,0,0,0},
                       {0,0,0,0,0}, {0,0,0,0,0}}; //Initial values, [species][w*p]
  Double_t outVals[4][4] = {{0,0,0,0}, {0,0,0,0},
                            {0,0,0,0}, {0,0,0,0}};
  Double_t trackXYZ[3];
  fGFW->Clear();
  Int_t nTotNoTracks=0;
  Double_t ptmins[] = {0.2,0.2,0.3,0.5};
  Double_t ptmaxs[] = {10.,10.,6.0,6.0};
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    Double_t trackXYZ[] = {0.,0.,0.};
    if(!AcceptAODTrack(lTrack,trackXYZ)) continue;
    Double_t p1 = lTrack->Pt();
    if(TMath::Abs(lTrack->Eta())<0.8 && lTrack->Pt()>0.2 && p1<3)  nTotNoTracks++;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    Double_t weff = fWeights[PIDIndex]->GetIntegratedEfficiency(p1);
    Double_t wacc = GetMyWeight(lTrack->Eta(),lTrack->Phi(),PIDIndex);//POI weight
    Double_t waccRef = GetMyWeight(lTrack->Eta(),lTrack->Phi(),0);//Nch weight
    Bool_t WithinRef=(p1>0.2 && p1<5);
    Bool_t WithinPOI=(p1>ptmins[PIDIndex] && p1<ptmaxs[PIDIndex]);
    Bool_t WithinNch=(p1>ptmins[0] && p1<ptmaxs[0]); //Within Ncharged (important for e.g. protons)
    if(TMath::Abs(lTrack->Eta())<0.4)  { //for mean pt, only consider -0.4-0.4 region
      if(weff==0) continue;
      Double_t w = 1./weff;
      FillWPCounter(wp[0],w,p1);
      if(PIDIndex) FillWPCounter(wp[PIDIndex],w,p1); //should be different weight here
    } else { //Otherwise, we consider it for vn calculations
      if(!WithinPOI && !WithinRef) continue;
      if(WithinPOI && WithinRef) waccRef = wacc; //If overlapping, override ref weight
      if(WithinRef) fGFW->Fill(lTrack->Eta(),1,lTrack->Phi(),waccRef,1); //Filling ref flow
      if(WithinPOI && PIDIndex) fGFW->Fill(lTrack->Eta(),1,lTrack->Phi(),wacc,(1<<(1+PIDIndex))); //Filling POI/only identified
      if(WithinNch) fGFW->Fill(lTrack->Eta(),1,lTrack->Phi(),wacc,2); //always filling for Nch
      if(WithinPOI && PIDIndex && WithinRef) fGFW->Fill(lTrack->Eta(),1,lTrack->Phi(),wacc,1<<(PIDIndex+5));
      if(WithinNch && WithinRef) fGFW->Fill(lTrack->Eta(),1,lTrack->Phi(),wacc,32); //Filling POI flow for ID'ed
    };
  };
  if(wp[0][0]==0) return; //if no single charged particles, then surely no PID either, no sense to continue
  //Filling pT varianve
  for(Int_t i=0;i<4;i++) {
    if(!wp[i][0]) continue;
    outVals[i][0] = fmPT[i]->GetBinContent(fmPT[i]->FindBin(nTotNoTracks));
    CalculateMptValues(outVals[i],wp[i]);
    if(outVals[i][2]!=0)
      fptvar[i]->Fill(nTotNoTracks,outVals[i][1]/outVals[i][2],outVals[i][2]);
  };
  PostData(1,fptVarList);
  //Filling FCs
  for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
    Bool_t filled = FillFCs(corrconfigs.at(l_ind),nTotNoTracks,0);
  };
  PostData(2,fFC);
  for(Int_t i=0;i<4;i++) {
    FillCovariance(fCovariance[i],corrconfigs.at(i*4),nTotNoTracks,outVals[i][3]-outVals[i][0],wp[i][0]);
    FillCovariance(fCovariance[i],corrconfigs.at(i*4+1),nTotNoTracks,outVals[i][3]-outVals[i][0],wp[i][0]);
  };
  PostData(3,fCovList);
  //Assuming outArr[0] is preset to meanPt; outArr[1] = variance; outArr[2] = norm; outArr[3] = mpt in this event
}
void AliAnalysisTaskMeanPtV2Corr::ProduceALICEPublished_MptProd(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
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
    if(!AcceptAODTrackALICEPublished(lTrack,trackXYZ)) continue;
    nTotNoTracks++;
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    FillMeanPtCounter(lpt,l_ptsum[0],l_ptCount[0],0);
    if(PIDIndex) FillMeanPtCounter(lpt,l_ptsum[PIDIndex],l_ptCount[PIDIndex],0);
  };
  if(l_ptCount[0]==0) return;
  for(Int_t i=0;i<4;i++) {
    if(!l_ptCount[i]) continue;
    fmPT[i]->Fill(nTotNoTracks,l_ptsum[i]/l_ptCount[i],l_ptCount[i]);
  }
  fNchVsMulti->Fill(l_Cent,nTotNoTracks);
  PostData(1,fMPTList);
}
void AliAnalysisTaskMeanPtV2Corr::ProduceALICEPublished_CovProd(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent) {
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
    if(!AcceptAODTrackALICEPublished(lTrack,trackXYZ)) continue;
    nTotNoTracks++;
    Double_t p1 = lTrack->Pt();
    Int_t PIDIndex = GetBayesPIDIndex(lTrack)+1;
    FillWPCounter(wp[0],1,p1);
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
Bool_t AliAnalysisTaskMeanPtV2Corr::FillCovariance(TProfile *target, AliGFW::CorrConfig corconf, Double_t cent, Double_t d_mpt, Double_t dw_mpt) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      target->Fill(cent,val*d_mpt,dnx*dw_mpt);
    return kTRUE;
  };
  return kTRUE;
};
void AliAnalysisTaskMeanPtV2Corr::CreateCorrConfigs() {
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
  if(fStageSwitch==1 || fStageSwitch==4 || fStageSwitch==5) return;
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
Bool_t AliAnalysisTaskMeanPtV2Corr::WithinSigma(Double_t SigmaCut, AliAODTrack *inTrack, AliPID::EParticleType partType) {
  if(!fPIDResponse) return kFALSE;
  Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(inTrack,partType);
  Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(inTrack,partType);
  return (TMath::Sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF) < SigmaCut);
}
Int_t AliAnalysisTaskMeanPtV2Corr::GetBayesPIDIndex(AliAODTrack *l_track) {
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
Bool_t AliAnalysisTaskMeanPtV2Corr::LoadMyWeights(Int_t lRunNo) {
  if(!fNUAList) AliFatal("NUA list not set or does not exist!\n");
  if(lRunNo && lRunNo == fRunNo) return kTRUE;
  if(fNUAHist) { delete fNUAHist; };
  fNUAHist = new TH2D*[4];
  TString nuaNames[] = {"Charged","Pion","Kaon","Proton"};
  for(Int_t i=0; i<4;i++) {
    if(lRunNo) nuaNames[i].Prepend(Form("w%i_",lRunNo));
    fNUAHist[i] = (TH2D*)fNUAList->FindObject(nuaNames[i].Data());
    if(!fNUAHist[i]) AliFatal(Form("%s could not be found in the list!\n",nuaNames[i].Data()));
  }
  return kTRUE;
}
Double_t AliAnalysisTaskMeanPtV2Corr::GetMyWeight(Double_t eta, Double_t phi, Int_t pidind) {
  Int_t etaind = fNUAHist[pidind]->GetXaxis()->FindBin(eta);
  Int_t phiind = fNUAHist[pidind]->GetYaxis()->FindBin(phi);
  return fNUAHist[pidind]->GetBinContent(etaind,phiind);
}
