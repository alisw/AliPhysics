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


//cut values:

Float_t XiDaughterDCA = .8;
//Float_t XiDaughterDCA = 1.1;

Float_t XiMinDistBachToPrimVtx = 0.001;

Float_t XiCPA = 0.995;
//Float_t XiCPA = 0.99;

Float_t XiTransverseRadius = .001;

Float_t v0MassRange = 0.006;

//Float_t v0MaxDaughterDCA = 1.2;
Float_t v0MaxDaughterDCA = 1.6;

//Float_t v0CPA = 0.97;
Float_t v0CPA = 0.94;

Float_t v0TransverseRadius = .001;

Float_t v0MinDistToPrimVtx = 0.06;

Float_t v0MinDaugDistToPrimVtx = 0.001;

Float_t RejectOmegas = 0.005;

Float_t PtRangeXi = 0.001;

bool CheckPileUp = false; 

//Int_t TPCclusters = 60;
Int_t TPCclusters = 40;

Float_t CrF = 0.75;

//Int_t nSigma = 4;
Int_t nSigma = 5;


  //Background selection
  //--------------------

  CascadeCuts->SetXiMassRange(1.6725, 0.0575, 0.005);
//  CascadeCuts->SetXiMassRange(1.6725, 0.0575, 1.6725, 0.005);
  CascadeCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  CascadeCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  CascadeCuts->SetCutXiCPA(XiCPA);
  CascadeCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  CascadeCuts->Setv0MassRange(1.116, v0MassRange);
  CascadeCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  CascadeCuts->SetCutv0CPA(v0CPA);
  CascadeCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  CascadeCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  CascadeCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  //CascadeCuts->SetRejectMass(1.322, 0.005, 3312);
  CascadeCuts->SetRejectOmegas(1.322, RejectOmegas);
  CascadeCuts->SetPtRangeXi(PtRangeXi, 999.9);
XiNegCuts->SetCheckPileUp(CheckPileUp);
  XiPosCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  XiNegCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  XiBachCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  ////XiBachCuts->SetPID(AliPID::kKaon, 999., 4, true, 3);
  //XiBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  XiPosCuts->SetPID(AliPID::kProton, 999., nSigma);
  XiNegCuts->SetPID(AliPID::kPion, 999., nSigma);
  XiBachCuts->SetPID(AliPID::kKaon, 999., nSigma);


  AntiCascadeCuts->SetXiMassRange(1.6725, 0.0575, 0.005);
  //AntiCascadeCuts->SetXiMassRange(1.6725, 0.0575, 1.6725, 0.005);
  AntiCascadeCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  AntiCascadeCuts->SetCutXiCPA(XiCPA);
  AntiCascadeCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  AntiCascadeCuts->Setv0MassRange(1.116, v0MassRange);
  AntiCascadeCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  AntiCascadeCuts->SetCutv0CPA(v0CPA);
  AntiCascadeCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  AntiCascadeCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  //AntiCascadeCuts->SetRejectMass(1.322, 0.005, 3312);
  AntiCascadeCuts->SetRejectOmegas(1.322, RejectOmegas);
  AntiCascadeCuts->SetPtRangeXi(PtRangeXi, 999.9); 
AntiXiNegCuts->SetCheckPileUp(CheckPileUp);
  AntiXiPosCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  AntiXiNegCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  AntiXiBachCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  ////AntiXiBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);
  //AntiXiBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  AntiXiPosCuts->SetPID(AliPID::kPion, 999., nSigma);
  AntiXiNegCuts->SetPID(AliPID::kProton, 999., nSigma);
  AntiXiBachCuts->SetPID(AliPID::kKaon, 999., nSigma);


  //Signal selection
  //--------------------

  CascadeOmegaCuts->SetXiMassRange(1.6725, 0.005);
//  CascadeOmegaCuts->SetXiMassRange(1.6725, 0.005, 0., 0.);
  CascadeOmegaCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  CascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  CascadeOmegaCuts->SetCutXiCPA(XiCPA);
  CascadeOmegaCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  CascadeOmegaCuts->Setv0MassRange(1.116, v0MassRange);
  CascadeOmegaCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  CascadeOmegaCuts->SetCutv0CPA(v0CPA);
  CascadeOmegaCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  CascadeOmegaCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  CascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  //CascadeOmegaCuts->SetRejectMass(1.322, 0.005, 3312);
  CascadeOmegaCuts->SetRejectOmegas(1.322, RejectOmegas);
  CascadeOmegaCuts->SetPtRangeXi(PtRangeXi, 999.9); 
OmegaNegCuts->SetCheckPileUp(CheckPileUp);
  OmegaPosCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  OmegaNegCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  OmegaBachCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  ////OmegaBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);
  //OmegaBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  OmegaPosCuts->SetPID(AliPID::kProton, 999., nSigma);
  OmegaNegCuts->SetPID(AliPID::kPion, 999., nSigma);
  OmegaBachCuts->SetPID(AliPID::kKaon, 999., nSigma);

  AntiCascadeOmegaCuts->SetXiMassRange(1.6725, 0.005);
//  AntiCascadeOmegaCuts->SetXiMassRange(1.6725, 0.005, 0., 0.);
  AntiCascadeOmegaCuts->SetCutXiDaughterDCA(XiDaughterDCA);
  AntiCascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(XiMinDistBachToPrimVtx);
  AntiCascadeOmegaCuts->SetCutXiCPA(XiCPA);
  AntiCascadeOmegaCuts->SetCutXiTransverseRadius(XiTransverseRadius, 200);
  AntiCascadeOmegaCuts->Setv0MassRange(1.116, v0MassRange);
  AntiCascadeOmegaCuts->SetCutv0MaxDaughterDCA(v0MaxDaughterDCA);
  AntiCascadeOmegaCuts->SetCutv0CPA(v0CPA);
  AntiCascadeOmegaCuts->SetCutv0TransverseRadius(v0TransverseRadius, 200);
  AntiCascadeOmegaCuts->SetCutv0MinDistToPrimVtx(v0MinDistToPrimVtx);
  AntiCascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(v0MinDaugDistToPrimVtx);
  //AntiCascadeOmegaCuts->SetRejectMass(1.322, 0.005, 3312);
  AntiCascadeOmegaCuts->SetRejectOmegas(1.322, RejectOmegas);
  AntiCascadeOmegaCuts->SetPtRangeXi(PtRangeXi, 999.9);    
AntiOmegaNegCuts->SetCheckPileUp(CheckPileUp);
  AntiOmegaPosCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  AntiOmegaNegCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  AntiOmegaBachCuts->SetCutTPCCrossedRows(true, TPCclusters, CrF);
  ////AntiOmegaBachCuts->SetPID(AliPID::kKaon, 999., 4,true,3);
  //AntiOmegaBachCuts->SetPID(AliPID::kKaon, 999., nSigma);
  AntiOmegaPosCuts->SetPID(AliPID::kPion, 999., nSigma);
  AntiOmegaNegCuts->SetPID(AliPID::kProton, 999., nSigma);
  AntiOmegaBachCuts->SetPID(AliPID::kKaon, 999., nSigma);

//Proton cuts:
}
