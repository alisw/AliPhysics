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
AliFemtoDreamTrackCuts *AntiTrackCuts
){

  //Background selection
  //--------------------


  CascadeCuts->SetXiMassRange(1.6725, 0.0575, 0.005);
//  CascadeCuts->SetXiMassRange(1.6725, 0.0575, 1.6725, 0.005);
  CascadeCuts->SetCutXiDaughterDCA(.8);
  CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.001);
  CascadeCuts->SetCutXiCPA(0.995);
  CascadeCuts->SetCutXiTransverseRadius(0.001, 200);
  CascadeCuts->Setv0MassRange(1.116, 0.006);
  CascadeCuts->SetCutv0MaxDaughterDCA(1.2);
  CascadeCuts->SetCutv0CPA(0.97);
  CascadeCuts->SetCutv0TransverseRadius(.001, 200);
  CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
  CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.001);
  //CascadeCuts->SetRejectMass(1.322, 0.005, 3312);
  CascadeCuts->SetRejectOmegas(1.322, 0.005);
  CascadeCuts->SetPtRangeXi(0.001, 999.9);
XiNegCuts->SetCheckPileUp(false);
  XiPosCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  XiNegCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  XiBachCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  //XiBachCuts->SetPID(AliPID::kKaon, 999., 4, true, 3);

  AntiCascadeCuts->SetXiMassRange(1.6725, 0.0575, 0.005);
  //AntiCascadeCuts->SetXiMassRange(1.6725, 0.0575, 1.6725, 0.005);
  AntiCascadeCuts->SetCutXiDaughterDCA(.8);
  AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.001);
  AntiCascadeCuts->SetCutXiCPA(0.995);
  AntiCascadeCuts->SetCutXiTransverseRadius(0.001, 200);
  AntiCascadeCuts->Setv0MassRange(1.116, 0.006);
  AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.2);
  AntiCascadeCuts->SetCutv0CPA(0.97);
  AntiCascadeCuts->SetCutv0TransverseRadius(.001, 200);
  AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);
  AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.001);
  //AntiCascadeCuts->SetRejectMass(1.322, 0.005, 3312);
  AntiCascadeCuts->SetRejectOmegas(1.322, 0.005);
  AntiCascadeCuts->SetPtRangeXi(0.001, 999.9); 
AntiXiNegCuts->SetCheckPileUp(false);
  AntiXiPosCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  AntiXiNegCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  AntiXiBachCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  AntiXiBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);

  CascadeOmegaCuts->SetXiMassRange(1.6725, 0.005);
//  CascadeOmegaCuts->SetXiMassRange(1.6725, 0.005, 0., 0.);
  CascadeOmegaCuts->SetCutXiDaughterDCA(.8);
  CascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(0.001);
  CascadeOmegaCuts->SetCutXiCPA(0.995);
  CascadeOmegaCuts->SetCutXiTransverseRadius(0.001, 200);
  CascadeOmegaCuts->Setv0MassRange(1.116, 0.006);
  CascadeOmegaCuts->SetCutv0MaxDaughterDCA(1.2);
  CascadeOmegaCuts->SetCutv0CPA(0.97);
  CascadeOmegaCuts->SetCutv0TransverseRadius(.001, 200);
  CascadeOmegaCuts->SetCutv0MinDistToPrimVtx(0.06);
  CascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.001);
  //CascadeOmegaCuts->SetRejectMass(1.322, 0.005, 3312);
  CascadeOmegaCuts->SetRejectOmegas(1.322, 0.005);
  CascadeOmegaCuts->SetPtRangeXi(0.001, 999.9); 
OmegaNegCuts->SetCheckPileUp(false);
  OmegaPosCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  OmegaNegCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  OmegaBachCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  OmegaBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);

  AntiCascadeOmegaCuts->SetXiMassRange(1.6725, 0.005);
//  AntiCascadeOmegaCuts->SetXiMassRange(1.6725, 0.005, 0., 0.);
  AntiCascadeOmegaCuts->SetCutXiDaughterDCA(.8);
  AntiCascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(0.001);
  AntiCascadeOmegaCuts->SetCutXiCPA(0.995);
  AntiCascadeOmegaCuts->SetCutXiTransverseRadius(0.001, 200);
  AntiCascadeOmegaCuts->Setv0MassRange(1.116, 0.006);
  AntiCascadeOmegaCuts->SetCutv0MaxDaughterDCA(1.2);
  AntiCascadeOmegaCuts->SetCutv0CPA(0.97);
  AntiCascadeOmegaCuts->SetCutv0TransverseRadius(.001, 200);
  AntiCascadeOmegaCuts->SetCutv0MinDistToPrimVtx(0.06);
  AntiCascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.001);
  //AntiCascadeOmegaCuts->SetRejectMass(1.322, 0.005, 3312);
  AntiCascadeOmegaCuts->SetRejectOmegas(1.322, 0.005);
  AntiCascadeOmegaCuts->SetPtRangeXi(0.001, 999.9);    
AntiOmegaNegCuts->SetCheckPileUp(false);
  AntiOmegaPosCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  AntiOmegaNegCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  AntiOmegaBachCuts->SetCutTPCCrossedRows(true, 60, 0.75);
  AntiOmegaBachCuts->SetPID(AliPID::kKaon, 999., 4,true,3);

//Proton cuts:
}
