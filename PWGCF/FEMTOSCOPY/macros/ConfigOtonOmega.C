//ConfigOtonOmega called by AddTaskOtonOmega to set additional cuts or overwrite the standards

void ConfigOtonOmega(
AliFemtoDreamCascadeCuts *CascadeCuts,
//AliOtonOmegaCascadeCuts *CascadeCuts,
AliFemtoDreamTrackCuts *XiPosCuts,
AliFemtoDreamTrackCuts *XiNegCuts,
AliFemtoDreamTrackCuts *XiBachCuts,
AliFemtoDreamCascadeCuts *AntiCascadeCuts,
//AliOtonOmegaCascadeCuts *AntiCascadeCuts,
AliFemtoDreamTrackCuts *AntiXiPosCuts,
AliFemtoDreamTrackCuts *AntiXiNegCuts,
AliFemtoDreamTrackCuts *AntiXiBachCuts,
AliFemtoDreamCascadeCuts * CascadeOmegaCuts,
//AliOtonOmegaCascadeCuts * CascadeOmegaCuts,
AliFemtoDreamTrackCuts *OmegaPosCuts,
AliFemtoDreamTrackCuts *OmegaNegCuts,
AliFemtoDreamTrackCuts *OmegaBachCuts,
AliFemtoDreamCascadeCuts *AntiCascadeOmegaCuts,
//AliOtonOmegaCascadeCuts *AntiCascadeOmegaCuts,
AliFemtoDreamTrackCuts *AntiOmegaPosCuts,
AliFemtoDreamTrackCuts *AntiOmegaNegCuts,
AliFemtoDreamTrackCuts *AntiOmegaBachCuts,
AliFemtoDreamTrackCuts *TrackCuts,
AliFemtoDreamTrackCuts *AntiTrackCuts,
int iCuts = 0
){

//cut values: (iCuts 0 -> std, iCuts 1 -> open, for syst errors, iCuts 2 -> std+proton filterbit 96)
//-----------
Float_t XiDaughterDCA = .8; //std
if(iCuts==1) XiDaughterDCA = 1.1;

Float_t v0MaxDaughterDCA = 1.2; //std
if(iCuts==1) v0MaxDaughterDCA = 1.6;

Float_t v0MinDistToPrimVtx = 0.06; //std
if(iCuts==1) v0MinDistToPrimVtx = 0.03;

Float_t XiCPA = 0.995; //std
if(iCuts==1) XiCPA = 0.993;

Float_t v0CPA = 0.97; //std
if(iCuts==1) v0CPA = 0.94;

Float_t XiMinDistBachToPrimVtx = 0.04; //std
if(iCuts==1) XiMinDistBachToPrimVtx = 0.02;

Float_t v0MinDaugDistToPrimVtx = 0.04; //std
if(iCuts==1) v0MinDaugDistToPrimVtx = 0.02;

Float_t Eta = 0.8; //std
if(iCuts==1) Eta = 0.9;

Float_t nSigma = 4.; //std
if(iCuts==1) nSigma = 5.;

Int_t TPCcrR = 60; //std
if(iCuts==1) TPCcrR = 50;

Float_t CrF = 0.75; //std
if(iCuts==1) CrF = 0.65;

Float_t v0TransverseRadius = .000001; //std, no cut (cut in Vertexer?)

Float_t XiTransverseRadius = .000001; //std, no cut (cut in Vertexer?)

Float_t v0MassRange = 0.006; //std

Float_t RejectXis = 0.008; //std

Float_t PtRangeXi = 0.000001; //std, no cut (cut in Vertexer?)

bool PionCheckPileUp = false; //std


  //Background selection
  //--------------------
  CascadeCuts->SetXiMassRange(1.67245, 0.0575, 0.005);
  CascadeCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  CascadeCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  CascadeCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  CascadeCuts->SetCutXiCPA(XiCPA);
  CascadeCuts->SetCutv0CPA(v0CPA);
  CascadeCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  CascadeCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  XiPosCuts->SetEtaRange(-1.*Eta,Eta);
  XiNegCuts->SetEtaRange(-1.*Eta,Eta);
  XiBachCuts->SetEtaRange(-1.*Eta,Eta);
  XiPosCuts->SetPID(AliPID::kProton, 999., nSigma);
  XiNegCuts->SetPID(AliPID::kPion, 999., nSigma);
  //XiBachCuts->SetPID(AliPID::kKaon, 999., 4, true, 3); //old, for allowing ITS only
  XiBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  XiPosCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  XiNegCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  XiBachCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  CascadeCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  CascadeCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  CascadeCuts->Setv0MassRange(1.116, v0MassRange);
  CascadeCuts->SetRejectXis(1.322, RejectXis);
  CascadeCuts->SetPtRangeXi(PtRangeXi, 999.9);
  XiNegCuts->SetCheckPileUp(PionCheckPileUp);

//nsigma proton
//fix in mT.C DCA xi standard 20 y change 0.3!!!!!

  AntiCascadeCuts->SetXiMassRange(1.67245, 0.0575, 0.005);
  AntiCascadeCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  AntiCascadeCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  AntiCascadeCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  AntiCascadeCuts->SetCutXiCPA(XiCPA);
  AntiCascadeCuts->SetCutv0CPA(v0CPA);
  AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  AntiXiPosCuts->SetEtaRange(-1.*Eta,Eta);
  AntiXiNegCuts->SetEtaRange(-1.*Eta,Eta);
  AntiXiBachCuts->SetEtaRange(-1.*Eta,Eta);
  AntiXiPosCuts->SetPID(AliPID::kPion, 999., nSigma);
  AntiXiNegCuts->SetPID(AliPID::kProton, 999., nSigma);
  //AntiXiBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);
  AntiXiBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  AntiXiPosCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiXiNegCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiXiBachCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiCascadeCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  AntiCascadeCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  AntiCascadeCuts->Setv0MassRange(1.116, v0MassRange);
  AntiCascadeCuts->SetRejectXis(1.322, RejectXis);
  AntiCascadeCuts->SetPtRangeXi(PtRangeXi, 999.9); 
  AntiXiPosCuts->SetCheckPileUp(PionCheckPileUp);

  //Signal selection
  //----------------
  CascadeOmegaCuts->SetXiMassRange(1.67245, 0.005);
  CascadeOmegaCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  CascadeOmegaCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  CascadeOmegaCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  CascadeOmegaCuts->SetCutXiCPA(XiCPA);
  CascadeOmegaCuts->SetCutv0CPA(v0CPA);
  CascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  CascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  OmegaPosCuts->SetPID(AliPID::kProton, 999., nSigma);
  OmegaNegCuts->SetPID(AliPID::kPion, 999., nSigma);
  //OmegaBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);
  OmegaBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  OmegaPosCuts->SetEtaRange(-1.*Eta,Eta);
  OmegaNegCuts->SetEtaRange(-1.*Eta,Eta);
  OmegaBachCuts->SetEtaRange(-1.*Eta,Eta);
  OmegaPosCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  OmegaNegCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  OmegaBachCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  CascadeOmegaCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  CascadeOmegaCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  CascadeOmegaCuts->Setv0MassRange(1.116, v0MassRange);
  CascadeOmegaCuts->SetRejectXis(1.322, RejectXis);
  CascadeOmegaCuts->SetPtRangeXi(PtRangeXi, 999.9); 
  OmegaNegCuts->SetCheckPileUp(PionCheckPileUp);

  AntiCascadeOmegaCuts->SetXiMassRange(1.67245, 0.005);
  AntiCascadeOmegaCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  AntiCascadeOmegaCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  AntiCascadeOmegaCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  AntiCascadeOmegaCuts->SetCutXiCPA(XiCPA);
  AntiCascadeOmegaCuts->SetCutv0CPA(v0CPA);
  AntiCascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  AntiCascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  AntiOmegaPosCuts->SetEtaRange(-1.*Eta,Eta);
  AntiOmegaNegCuts->SetEtaRange(-1.*Eta,Eta);
  AntiOmegaBachCuts->SetEtaRange(-1.*Eta,Eta);
  AntiOmegaPosCuts->SetPID(AliPID::kPion, 999., nSigma);
  AntiOmegaNegCuts->SetPID(AliPID::kProton, 999., nSigma);
  //AntiOmegaBachCuts->SetPID(AliPID::kKaon, 999., 4,true,3);
  AntiOmegaBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  AntiOmegaPosCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiOmegaNegCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiOmegaBachCuts->SetCutTPCCrossedRows(true, TPCcrR, CrF);
  AntiCascadeOmegaCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  AntiCascadeOmegaCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  AntiCascadeOmegaCuts->Setv0MassRange(1.116, v0MassRange);
  AntiCascadeOmegaCuts->SetRejectXis(1.322, RejectXis);
  AntiCascadeOmegaCuts->SetPtRangeXi(PtRangeXi, 999.9);    
  AntiOmegaPosCuts->SetCheckPileUp(PionCheckPileUp);

  //Proton cuts:
  //------------
  if(iCuts==1){
   TrackCuts->SetPtRange(0.4, 4.05); //std 0.5
   TrackCuts->SetEtaRange(-0.9,0.9); //std 0.8
   TrackCuts->SetNClsTPC(70); //std 80
   TrackCuts->SetDCAVtxXY(2.5); //std 0.1
   TrackCuts->SetPID(AliPID::kProton, 0.75, 4.); //std 3

   AntiTrackCuts->SetPtRange(0.4, 4.05);
   AntiTrackCuts->SetEtaRange(-0.9, 0.9);
   AntiTrackCuts->SetNClsTPC(70);
   AntiTrackCuts->SetDCAVtxXY(2.5);
   AntiTrackCuts->SetPID(AliPID::kProton, 0.75, 4.); 
  }
  if(iCuts==2) {
   TrackCuts->SetFilterBit(96);
   AntiTrackCuts->SetFilterBit(96);
  }

}
