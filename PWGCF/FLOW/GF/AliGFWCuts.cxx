#include "AliGFWCuts.h"
const Int_t AliGFWCuts::fNTrackFlags=9;
const Int_t AliGFWCuts::fNEventFlags=8;
AliGFWCuts::AliGFWCuts():
  fSystFlag(0),
  fFilterBit(96),
  fDCAxyCut(7),
  fDCAzCut(2),
  fTPCNcls(70),
  fVtxZ(10),
  fEta(0.8),
  fRequiresExtraWeight(kTRUE)
{
};
AliGFWCuts::~AliGFWCuts() {
};
Int_t AliGFWCuts::AcceptTrack(AliAODTrack* l_Tr, Double_t* l_DCA, Int_t BitShift) {
  if(TMath::Abs(l_Tr->Eta())>fEta) return 0;
  if(!l_Tr->TestFilterBit(fFilterBit)) return 0;
  if(fFilterBit!=2) {//Check is not valid for ITSsa tracks
    if(l_Tr->GetTPCNclsF()<fTPCNcls) return 0;
  } else {
    UInt_t status = l_Tr->GetStatus();
    if ((status&AliESDtrack::kITSrefit)==0) return 0;
    if ((status & AliESDtrack::kITSin) == 0 || (status & AliESDtrack::kTPCin)) return 0;
  };
  if(l_DCA[0]>fDCAzCut) return 0;
  Double_t DCAxycut;
  if(fFilterBit!=2) DCAxycut = 0.0105+0.0350/TMath::Power(l_Tr->Pt(),1.1);//*fDCAxyCut/7.; //TPC tracks and my ITS cuts
  else DCAxycut = 0.0231+0.0315/TMath::Power(l_Tr->Pt(),1.3);
  if(l_DCA[1]>DCAxycut*(fDCAxyCut/7.))
    return 0;
  return 1<<BitShift;
};

Int_t AliGFWCuts::AcceptParticle(AliVParticle *l_Pa, Int_t BitShift, Double_t ptLow, Double_t ptHigh) {
  if(TMath::Abs(l_Pa->Eta())>fEta) return 0;
  if(ptLow>0) if(l_Pa->Pt()<ptLow) return 0;
  if(ptHigh>0) if(l_Pa->Pt()>ptHigh) return 0;
  // if(l_Pa->Pt()<0.3) return 0;
  // if(l_Pa->Pt()>20) return 0;
  // if(!l_Pa->IsMCPrimary()) return 0; //Not sure if I need this one here?
  return 1<<BitShift;
};
Int_t AliGFWCuts::AcceptVertex(AliAODEvent *l_Ev, Int_t BitShift) {
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
  fRequiresExtraWeight=kTRUE;
};
void AliGFWCuts::PrintSetup() {
  printf("**********\n");
  printf("Syst. flag: %i\n",fSystFlag);
  printf("Eta: %f\n",fEta);
  printf("(Flag 1) Filter bit: %i\n",fFilterBit);
  printf("(Flag 2,3) DCAxy cut: %2.0f sigma\n",fDCAxyCut);
  printf("(Flag 4,5) DCAz cut: %f\n",fDCAzCut);
  printf("(Flag 6-8) TPC Ncls: %i\n",fTPCNcls);
  printf("(Flag 9) ITS tracks\n");
  printf("Rest of the flags are global per event. Total flag = %i + vtx/ev flag\n",fNTrackFlags);
  printf("(Flag 1-3) Vertex selection: |z|<%2.1f\n",fVtxZ);
  printf("(Flag 4-5) CL1, CL2 multi. estimator (no weights)\n");
  printf("(Flag 6) pile-up 1500 cut\n");
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
    fDCAxyCut=10.;
    fRequiresExtraWeight=kTRUE;
    break;
  case 3:
    fDCAxyCut=4.;
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
    fTPCNcls=100;
    fRequiresExtraWeight=kTRUE;
    break;
  case 9:
    fFilterBit=2;
    fEta=1.6;
    fRequiresExtraWeight=kTRUE;
  default:
    break;
  };
};
void AliGFWCuts::SetupEventCuts(Int_t sysflag) {
  switch(sysflag) {
  case 1:
    fVtxZ = 5;
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
    retstr->Append("DCA_{xy} < 10 (old:8) sigma");
    break;
  case 3:
    retstr->Append("DCA_{xy} < 4 (old:6) sigma");
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
    retstr->Append("ITS tracklets");
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
