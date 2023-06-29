/*
Author: Vytautas Vislavicius
Contains the additional event and track selection used within the <AliGFW> framework.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#include "AliGFWCuts.h"
const Int_t AliGFWCuts::fNTrackFlags=16;
const Int_t AliGFWCuts::fNEventFlags=10;
AliESDtrackCuts *AliGFWCuts::fTCFB32=0;
AliESDtrackCuts *AliGFWCuts::fTCFB64=0;
AliESDtrackCuts *AliGFWCuts::fTCFB256=0;
AliESDtrackCuts *AliGFWCuts::fTCFB512=0;
AliESDtrackCuts *AliGFWCuts::fTCFB16=0;
AliGFWCuts::AliGFWCuts():
  fSystFlag(0),
  fFilterBit(96),
  fDCAxyCut(7.),
  fDCAzCut(2),
  fTPCNcls(70),
  fTPCChi2PerCluster(2.5),
  fVtxZ(10),
  fEta(0.8),
  fPtDepXYCut(0x0),
  fRequiresExtraWeight(kTRUE)
{
};
AliGFWCuts::~AliGFWCuts() {
};
Int_t AliGFWCuts::AcceptTrack(AliAODTrack *&l_Tr, Double_t *l_DCA, const Int_t &BitShift, const Bool_t &lDisableDCAxyCheck) {
  if(TMath::Abs(l_Tr->Eta())>fEta) return 0;
  if(!l_Tr->TestFilterBit(fFilterBit)) return 0;
  if(fFilterBit!=2) {//Check is not valid for ITSsa tracks
    if(l_Tr->GetTPCNclsF()<fTPCNcls) return 0;
  } else {
    UInt_t status = l_Tr->GetStatus();
    if ((status&AliESDtrack::kITSrefit)==0) return 0;
    if ((status & AliESDtrack::kITSin) == 0 || (status & AliESDtrack::kTPCin)) return 0;
  };
  if(l_Tr->GetTPCchi2perCluster()>fTPCChi2PerCluster) return 0;
  if(!l_DCA) return 1<<BitShift;
  if(TMath::Abs(l_DCA[2])>fDCAzCut) return 0;
  if(lDisableDCAxyCheck) return 1<<BitShift;
  Double_t DCAxyValue = TMath::Sqrt(l_DCA[0]*l_DCA[0] + l_DCA[1]*l_DCA[1]);
  if(DCAxyValue > fPtDepXYCut->Eval(l_Tr->Pt())) return 0;
  return 1<<BitShift;
  // Double_t DCAxycut;
  // if(fFilterBit!=2) DCAxycut = 0.0105+0.0350/TMath::Power(l_Tr->Pt(),1.1);//*fDCAxyCut/7.; //TPC tracks and my ITS cuts
  // else DCAxycut = 0.0231+0.0315/TMath::Power(l_Tr->Pt(),1.3);
  // if(DCAxyValue>DCAxycut*(fDCAxyCut/7.))
  //   return 0;
  // return 1<<BitShift;
};
Int_t AliGFWCuts::AcceptTrack(AliESDtrack *&l_Tr, Double_t *l_DCA, const Int_t &BitShift, UInt_t &PrimFlags) {
  //Initialize ESDtrackCuts if needed:
  if(!fTCFB32) fTCFB32 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  if(!fTCFB64) {
    fTCFB64 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    fTCFB64->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
    fTCFB64->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kFirst);
  };
  //For FB768, DCA cut is little bit more complicated and primary estimation is not _that_ straightforward. For now, we'll fit just the xy distribution
  if(!fTCFB256) {
    fTCFB256 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    fTCFB256->SetMaxDCAToVertexXY(2.4);
    fTCFB256->SetMaxDCAToVertexZ(3.2);
    fTCFB256->SetDCAToVertex2D(kTRUE);
    // fTCFB256->SetMaxChi2TPCConstrainedGlobal(36);
    fTCFB256->SetMaxFractionSharedTPCClusters(0.4);
  }
  if(!fTCFB512) {
    fTCFB512 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    fTCFB512->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    fTCFB512->SetRequireITSRefit(kTRUE);
    //Not sure if the following are relevant or not. But they should be? Otherwise, this is left without any DCA cuts, etc.
    fTCFB512->SetMaxDCAToVertexXY(2.4);
    fTCFB512->SetMaxDCAToVertexZ(3.2);
    fTCFB512->SetDCAToVertex2D(kTRUE);
    // fTCFB512->SetMaxChi2TPCConstrainedGlobal(36);
    fTCFB512->SetMaxFractionSharedTPCClusters(0.4);
  }
  if(!fTCFB16) {
    fTCFB16 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    fTCFB16->SetMaxDCAToVertexXY(2.4);
    fTCFB16->SetMaxDCAToVertexZ(3.2);
    fTCFB16->SetDCAToVertex2D(kTRUE);
  }

  Double_t l_Pt = l_Tr->Pt();
  if(l_Pt<0.01) return 0; //If below 10 MeV, throw away the track immediatelly -> Otherwise, problems calculating DCAxy value
  PrimFlags=0; //Initial value for primary selection -- initially, not a primary
  if(TMath::Abs(l_Tr->Eta())>fEta) return 0;
  if(fFilterBit==96) {
    if(!fTCFB32->AcceptTrack(l_Tr) && !fTCFB64->AcceptTrack(l_Tr)) return 0;
    if(l_DCA) if(TMath::Abs(l_DCA[0]) < fPtDepXYCut->Eval(l_Pt)) PrimFlags+=1;
  }
  if(fFilterBit==768) {
    if(!fTCFB256->AcceptTrack(l_Tr) && !fTCFB512->AcceptTrack(l_Tr)) return 0;
    if(l_DCA) if(TMath::Abs(l_DCA[0]) < fPtDepXYCut->Eval(l_Pt)) PrimFlags+=1; //Here assume DCA cut is passed by default -- it's implemented in ESDtrackCuts
  }
  if(fFilterBit==16) {
    if(!fTCFB16->AcceptTrack(l_Tr)) return 0;
    if(l_DCA) if(TMath::Abs(l_DCA[0]) < fPtDepXYCut->Eval(l_Pt)) PrimFlags+=1; //Here assume DCA cut is passed by default -- it's implemented in ESDtrackCuts
  }
  //chi2 TPC vs global constrained could be a custom number, but it's 36 everywhere anyways
  Double_t chi2TPCConstrained = GetChi2TPCConstrained(l_Tr);
  if(!(chi2TPCConstrained<0) && chi2TPCConstrained<36.) PrimFlags+=2;
  Double_t TPCChi2PerCluster = l_Tr->GetTPCchi2()/Double_t(l_Tr->GetTPCclusters(0));
  if(TPCChi2PerCluster>fTPCChi2PerCluster || TPCChi2PerCluster<0) return 0;
  if(!l_DCA) return 1<<BitShift;
  if(l_DCA[1]>fDCAzCut) return 0;
  // if(lDisableDCAxyCheck) return 1<<BitShift;
  // Double_t DCAxycut;
  // if(fFilterBit!=2) DCAxycut = 0.0105+0.0350/TMath::Power(l_Tr->Pt(),1.1);//*fDCAxyCut/7.; //TPC tracks and my ITS cuts
  // else DCAxycut = 0.0231+0.0315/TMath::Power(l_Tr->Pt(),1.3);
  // Double_t DCAxyValue = l_DCA[0];//TMath::Sqrt(l_DCA[0]*l_DCA[0] + l_DCA[1]*l_DCA[1]);
  // if(DCAxyValue>DCAxycut*(fDCAxyCut/7.))
  //   return 0;
  return 1<<BitShift;
};


Int_t AliGFWCuts::AcceptParticle(AliVParticle *l_Pa, Int_t BitShift, Double_t ptLow, Double_t ptHigh) {
  if(TMath::Abs(l_Pa->Eta())>fEta) return 0;
  if(ptLow>0) if(l_Pa->Pt()<ptLow) return 0;
  if(ptHigh>0) if(l_Pa->Pt()>ptHigh) return 0;
  return 1<<BitShift;
};
Int_t AliGFWCuts::AcceptVertex(AliAODEvent *l_Ev, Int_t BitShift) {
  Double_t lvtxz = TMath::Abs(l_Ev->GetPrimaryVertex()->GetZ());
  if(lvtxz>fVtxZ) return 0;
  return 1<<BitShift;
};
Int_t AliGFWCuts::AcceptVertex(AliESDEvent *l_Ev, Int_t BitShift) {
  Double_t lvtxz = TMath::Abs(l_Ev->GetPrimaryVertex()->GetZ());
  if(lvtxz>fVtxZ) return 0;
  return 1<<BitShift;
};

void AliGFWCuts::ResetCuts() {
  fSystFlag=0;
  fFilterBit=96;
  fDCAxyCut=7;
  fDCAzCut=2;
  fTPCNcls=70;
  fVtxZ=10;
  fEta=0.8;
  fTPCChi2PerCluster=2.5;
  fRequiresExtraWeight=kTRUE;
  SetPtDepDCAXY("[0]*(0.0026+0.005/(x^1.01))");//[0]*(0.0015 + 0.005/(x^1.1))");//"0.0105+0.0350/(x^1.1)");
};
void AliGFWCuts::PrintSetup() {
  printf("**********\n");
  printf("Syst. flag: %i\n",fSystFlag);
  printf("Eta: %f\n",fEta);
  printf("(Flag 1, 10-16) Filter bit: %i\n",fFilterBit);
  printf("(Flag 2,3,12,16) DCAxy cut: %2.0f sigma \n",fDCAxyCut);
  printf("(Flag 4,5) DCAz cut: %.1f\n",fDCAzCut);
  printf("(Flag 6-8) TPC Ncls: %i\n",fTPCNcls);
  printf("(Flag 10-16) TPC chi2/Ncls: %.1f \n",fTPCChi2PerCluster);
  printf("Rest of the flags are global per event. Total flag = %i + vtx/ev flag\n",fNTrackFlags);
  printf("(Flag 1-3) Vertex selection: |z|<%2.1f\n",fVtxZ);
  printf("(Flag 4-5) CL1, CL2 multi. estimator (no weights)\n");
  printf("(Flag 6) pile-up 15000 (-> 1500) cut\n");
  // printf("Extra cuts (not following the prev. nomenclature for consistency):\n");
  // printf("(NUA weight reload flag set to 0)\n");
  printf("(Flag 7-8): MF++/--\n");
  printf("(Flag 9-10): Extreme efficiency, +-4\%\n");
  //printf("(Flag 12, disabled) ITS tracks (filter bit %i, TPC Ncls = %i)\n",fFilterBit,fTPCNcls);
  printf("**********\n");
};
void AliGFWCuts::SetupTrackCuts(Int_t sysflag) {
  switch(sysflag) {
  case 1:
    fFilterBit=768;
    fRequiresExtraWeight=kTRUE;
    break;
  case 2:
    fDCAxyCut=6;
    fRequiresExtraWeight=kTRUE;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  case 3:
    fDCAxyCut=4;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    fRequiresExtraWeight=kTRUE;
    break;
  case 4:
    fDCAzCut=1;
    fRequiresExtraWeight=kTRUE;
    break;
  case 5:
    fDCAzCut=0.5;
    fRequiresExtraWeight=kTRUE;
    break;
  case 6:
    fTPCNcls=80;
    fRequiresExtraWeight=kTRUE;
    break;
  case 7:
    fTPCNcls=90;
    fRequiresExtraWeight=kTRUE;
    break;
  case 8:
    fTPCNcls=65;
    fRequiresExtraWeight=kTRUE;
    break;
  case 9:
    fRequiresExtraWeight=kTRUE;
    fFilterBit=768;
    fDCAzCut=10;//Arbitrary large value to disable the cut -> implemented as 2.4 in the AliESDtrackCuts (2D cut)
    fDCAxyCut=7;
    fTPCChi2PerCluster=2.5;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  case 10:
    fRequiresExtraWeight=kTRUE;
    fFilterBit=768;
    fDCAzCut=10;//Arbitrary large value to disable the cut -> implemented as 2.4 in the AliESDtrackCuts (2D cut)
    fDCAxyCut=7;
    fTPCChi2PerCluster=2;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  case 11:
    fRequiresExtraWeight=kTRUE;
    fFilterBit=768;
    fDCAzCut=10;//Arbitrary large value to disable the cut -> implemented as 2.4 in the AliESDtrackCuts (2D cut)
    fDCAxyCut=7;
    fTPCChi2PerCluster=3;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  case 12:
    fRequiresExtraWeight=kTRUE;
    fFilterBit=768;
    fDCAzCut=10;//Arbitrary large value to disable the cut -> implemented as 2.4 in the AliESDtrackCuts (2D cut)
    fDCAxyCut=8;
    fTPCChi2PerCluster=2.5;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  case 13:
    fRequiresExtraWeight=kTRUE;
    fFilterBit=16;
    fDCAzCut=10;//Arbitrary large value to disable the cut -> implemented as 2.4 in the AliESDtrackCuts (2D cut)
    fDCAxyCut=7;
    fTPCChi2PerCluster=2.5;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  case 14:
    fRequiresExtraWeight=kTRUE;
    fFilterBit=16;
    fDCAzCut=10;//Arbitrary large value to disable the cut -> implemented as 2.4 in the AliESDtrackCuts (2D cut)
    fDCAxyCut=7;
    fTPCChi2PerCluster=2;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  case 15:
    fRequiresExtraWeight=kTRUE;
    fFilterBit=16;
    fDCAzCut=10;//Arbitrary large value to disable the cut -> implemented as 2.4 in the AliESDtrackCuts (2D cut)
    fDCAxyCut=7;
    fTPCChi2PerCluster=3;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  case 16:
    fRequiresExtraWeight=kTRUE;
    fFilterBit=16;
    fDCAzCut=10;//Arbitrary large value to disable the cut -> implemented as 2.4 in the AliESDtrackCuts (2D cut)
    fDCAxyCut=8;
    fTPCChi2PerCluster=2.5;
    fPtDepXYCut->SetParameter(0,fDCAxyCut);
    break;
  default:
    break;
  };
};
void AliGFWCuts::SetupEventCuts(Int_t sysflag) {
  switch(sysflag) {
  case 1:
    fVtxZ = 8;
    fRequiresExtraWeight=kTRUE;
    break;
  case 2:
    fVtxZ = 7;
    fRequiresExtraWeight=kTRUE;
    break;
  case 3:
    fVtxZ = 9;
    fRequiresExtraWeight=kTRUE;
    break;
  case 4:
    printf("Warning! Event flag %i (syst. flag %i), CL1 estmator: make sure the proper estimator is used in the task!\n",sysflag, sysflag+fNTrackFlags);
    fRequiresExtraWeight=kFALSE;
    break;
  case 5:
    printf("Warning! Event flag %i (syst. flag %i), CL0 estmator: make sure the proper estimator is used in the task!\n",sysflag, sysflag+fNTrackFlags);
    fRequiresExtraWeight=kFALSE;
    break;
  case 6:
    printf("Warning! Event flag %i (syst. flag %i), PU cuts: make sure proper PU cuts are used in the task!\n",sysflag, sysflag+fNTrackFlags);
    fRequiresExtraWeight=kTRUE;
    break;
  case 7:
    printf("Warning! Event flag %i (syst. flag %i), magnetic field configuration ++: no cuts here, please make sure the proper runlist is used!\n",sysflag,sysflag+fNTrackFlags);
    fRequiresExtraWeight=kFALSE;
    break;
  case 8:
    printf("Warning! Event flag %i (syst. flag %i), magnetic field configuration --: no cuts here, please make sure the proper runlist is used!\n",sysflag,sysflag+fNTrackFlags);
    fRequiresExtraWeight=kFALSE;
    break;
  case 9:
    printf("Warning! Make sure the correct efficiency is provided (+4\%)\n",sysflag,sysflag+fNTrackFlags);
    fRequiresExtraWeight=kFALSE;
    break;
  case 10:
    printf("Warning! Make sure the correct efficiency is provided (-4\%)\n",sysflag,sysflag+fNTrackFlags);
    fRequiresExtraWeight=kFALSE;
    break;
  default:
    break;
  };
};
void AliGFWCuts::SetupCuts(Int_t sysflag) {
  ResetCuts();
  fSystFlag=sysflag;
  if(sysflag==0 || sysflag>fNTrackFlags+fNEventFlags) return;
  if(sysflag<=fNTrackFlags) SetupTrackCuts(sysflag);
  else SetupEventCuts(sysflag-fNTrackFlags);
};

TString *AliGFWCuts::GetTrackFlagDescriptor(Int_t sysflag) {
  TString *retstr = new TString("");
  switch(sysflag) {
  case 1:
    retstr->Append("Filter bit 768");
    break;
  case 2:
    retstr->Append("DCA_{xy} < 10 sigma");
    break;
  case 3:
    retstr->Append("DCA_{xy} < 4 sigma");
    break;
  case 4:
    retstr->Append("DCA_{z} < 1 cm");
    break;
  case 5:
    retstr->Append("DCA_{z} < 0.5 cm");
    break;
  case 6:
    retstr->Append("TPC N_{Cls} > 80");
    break;
  case 7:
    retstr->Append("TPC N_{Cls} > 90");
    break;
  case 8:
    retstr->Append("TPC N_{Cls} > 100");
    break;
  case 9:
    retstr->Append("FB768, #chi^{2}/N_{TPCcls} < 2.5");
    break;
  case 10:
    retstr->Append("FB768, #chi^{2}/N_{TPCcls} < 2");
    break;
  case 11:
    retstr->Append("FB768, #chi^{2}/N_{TPCcls} < 3");
    break;
  case 12:
    retstr->Append("FB768, DCA_{xy} 8#sigma");
    break;
  case 13:
    retstr->Append("FB16, #chi^{2}/N_{TPCcls} < 2.5");
    break;
  case 14:
    retstr->Append("FB16, #chi^{2}/N_{TPCcls} < 2");
    break;
  case 15:
    retstr->Append("FB16, #chi^{2}/N_{TPCcls} < 3");
    break;
  case 16:
    retstr->Append("FB16, DCA_{xy} 8#sigma");
    break;
  default:
    break;
  };
  return retstr;
};
TString* AliGFWCuts::GetEventFlagDescriptor(Int_t sysflag) {
  TString *retstr = new TString("");
  switch(sysflag) {
  case 1:
    retstr->Append("|z_{vtx}| < 5 cm");
    break;
  case 2:
    retstr->Append("|z_{vtx}| < 7 cm");
    break;
  case 3:
    retstr->Append("|z_{vtx}| < 9 cm");
    break;
  case 4:
    retstr->Append("CL1 multi.");
    break;
  case 5:
    retstr->Append("CL0 multi.");
    break;
  case 6:
    retstr->Append("PU cut 1500");
    break;
  case 7:
    retstr->Append("MF ++");
    break;
  case 8:
    retstr->Append("MF --");
    break;
  case 9:
    retstr->Append("Efficiency +4\%");
    break;
  case 10:
    retstr->Append("Efficiency -4\%");
    break;
  default:
    break;
  };
  return retstr;
};
TString *AliGFWCuts::GetFlagDescription(Int_t sysflag) {
  if(sysflag>0 && sysflag <= fNTrackFlags + fNEventFlags) {
    if(sysflag>fNTrackFlags) return GetEventFlagDescriptor(sysflag-fNTrackFlags);
    return GetTrackFlagDescriptor(sysflag);
  };
  TString *retst = new TString("");
  if(sysflag==0) retst->Append("Nominal");
  else retst->Append("Unknown_%i",sysflag);
  return retst;
};
Double_t AliGFWCuts::GetChi2TPCConstrained(const AliESDtrack *l_Tr) {
  const AliESDEvent* esdEvent = l_Tr->GetESDEvent();
  // get vertex
  const AliESDVertex* vertex = 0;
  Double_t chi2TPCConstrainedVsGlobal=-2; //Initial value
  vertex = esdEvent->GetPrimaryVertexTracks();
  if (!vertex || !vertex->GetStatus()) vertex = esdEvent->GetPrimaryVertexSPD();
  //TPC vertex not considered
  // if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & kVertexTPC) vertex = esdEvent->GetPrimaryVertexTPC();
  if (vertex->GetStatus()) chi2TPCConstrainedVsGlobal = l_Tr->GetChi2TPCConstrainedVsGlobal(vertex);
  return chi2TPCConstrainedVsGlobal;
}
