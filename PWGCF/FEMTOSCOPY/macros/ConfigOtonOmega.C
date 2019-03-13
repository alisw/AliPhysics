//ConfigOtonOmega called by AddTaskOtonOmega to set additional cuts or overwrite the standards

void ConfigOtonOmega(
AliOtonOmegaCascadeCuts *CascadeCuts,
AliFemtoDreamTrackCuts *XiPosCuts,
AliFemtoDreamTrackCuts *XiNegCuts,
AliFemtoDreamTrackCuts *XiBachCuts,
AliOtonOmegaCascadeCuts *AntiCascadeCuts,
AliFemtoDreamTrackCuts *AntiXiPosCuts,
AliFemtoDreamTrackCuts *AntiXiNegCuts,
AliFemtoDreamTrackCuts *AntiXiBachCuts,
AliOtonOmegaCascadeCuts * CascadeOmegaCuts,
AliFemtoDreamTrackCuts *OmegaPosCuts,
AliFemtoDreamTrackCuts *OmegaNegCuts,
AliFemtoDreamTrackCuts *OmegaBachCuts,
AliOtonOmegaCascadeCuts *AntiCascadeOmegaCuts,
AliFemtoDreamTrackCuts *AntiOmegaPosCuts,
AliFemtoDreamTrackCuts *AntiOmegaNegCuts,
AliFemtoDreamTrackCuts *AntiOmegaBachCuts
){


  CascadeCuts->SetXiMassRange(1.672, 0.040, 1.672,0.005);
  CascadeCuts->SetCutXiDaughterDCA(2.);                
  CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.001);    
  CascadeCuts->SetCutXiCPA(0.98);                  
  CascadeCuts->SetCutXiTransverseRadius(0.001, 200);   
  CascadeCuts->Setv0MassRange(1.116, 0.006);         
  CascadeCuts->SetCutv0MaxDaughterDCA(1.5);          
  CascadeCuts->SetCutv0CPA(0.97);                    
  CascadeCuts->SetCutv0TransverseRadius(.001, 200);    
  CascadeCuts->SetCutv0MinDistToPrimVtx(0.06);       
  CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.001);   
  CascadeCuts->SetRejectMass(1.322, 0.005, 3312);    
  CascadeCuts->SetPtRangeXi(0.001, 999.9);             
//XiPosCuts->SetCheckPileUp(false);
XiNegCuts->SetCheckPileUp(false);
//XiBachCuts->SetCheckPileUp(false);
  XiPosCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  XiNegCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  XiBachCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  XiBachCuts->SetPID(AliPID::kKaon, 999., 4, true, 3);



  AntiCascadeCuts->SetXiMassRange(1.672, 0.040, 1.672,0.005);
  AntiCascadeCuts->SetCutXiDaughterDCA(2.);                
  AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.001);    
  AntiCascadeCuts->SetCutXiCPA(0.98);                  
  AntiCascadeCuts->SetCutXiTransverseRadius(0.001, 200);   
  AntiCascadeCuts->Setv0MassRange(1.116, 0.006);         
  AntiCascadeCuts->SetCutv0MaxDaughterDCA(1.5);          
  AntiCascadeCuts->SetCutv0CPA(0.97);                    
  AntiCascadeCuts->SetCutv0TransverseRadius(.001, 200);    
  AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.06);       
  AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.001);   
  AntiCascadeCuts->SetRejectMass(1.322, 0.005, 3312);    
  AntiCascadeCuts->SetPtRangeXi(0.001, 999.9);            //AntiXiPosCuts->SetCheckPileUp(false);
AntiXiNegCuts->SetCheckPileUp(false);
//AntiXiBachCuts->SetCheckPileUp(false);
  AntiXiPosCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  AntiXiNegCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  AntiXiBachCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  AntiXiBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);


  CascadeOmegaCuts->SetXiMassRange(1.672, 0.005, 0., 0.);
  CascadeOmegaCuts->SetCutXiDaughterDCA(2.);                
  CascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(0.001);    
  CascadeOmegaCuts->SetCutXiCPA(0.98);                  
  CascadeOmegaCuts->SetCutXiTransverseRadius(0.001, 200);   
  CascadeOmegaCuts->Setv0MassRange(1.116, 0.006);         
  CascadeOmegaCuts->SetCutv0MaxDaughterDCA(1.5);          
  CascadeOmegaCuts->SetCutv0CPA(0.97);                    
  CascadeOmegaCuts->SetCutv0TransverseRadius(.001, 200);    
  CascadeOmegaCuts->SetCutv0MinDistToPrimVtx(0.06);       
  CascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.001);   
  CascadeOmegaCuts->SetRejectMass(1.322, 0.005, 3312);    
  CascadeOmegaCuts->SetPtRangeXi(0.001, 999.9);            //OmegaPosCuts->SetCheckPileUp(false);
OmegaNegCuts->SetCheckPileUp(false);
//OmegaBachCuts->SetCheckPileUp(false);
  OmegaPosCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  OmegaNegCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  OmegaBachCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  OmegaBachCuts->SetPID(AliPID::kKaon, 999., 4,true, 3);



  AntiCascadeOmegaCuts->SetXiMassRange(1.672, 0.005, 0., 0.);
  AntiCascadeOmegaCuts->SetCutXiDaughterDCA(2.);                
  AntiCascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(0.001);    
  AntiCascadeOmegaCuts->SetCutXiCPA(0.98);                  
  AntiCascadeOmegaCuts->SetCutXiTransverseRadius(0.001, 200);   
  AntiCascadeOmegaCuts->Setv0MassRange(1.116, 0.006);         
  AntiCascadeOmegaCuts->SetCutv0MaxDaughterDCA(1.5);          
  AntiCascadeOmegaCuts->SetCutv0CPA(0.97);                    
  AntiCascadeOmegaCuts->SetCutv0TransverseRadius(.001, 200);    
  AntiCascadeOmegaCuts->SetCutv0MinDistToPrimVtx(0.06);       
  AntiCascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.001);   
  AntiCascadeOmegaCuts->SetRejectMass(1.322, 0.005, 3312);    
  AntiCascadeOmegaCuts->SetPtRangeXi(0.001, 999.9);            //AntiOmegaPosCuts->SetCheckPileUp(false);
AntiOmegaNegCuts->SetCheckPileUp(false);
//AntiOmegaBachCuts->SetCheckPileUp(false);
  AntiOmegaPosCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  AntiOmegaNegCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  AntiOmegaBachCuts->SetCutTPCCrossedRows(false, 50, 0.5);
  AntiOmegaBachCuts->SetPID(AliPID::kKaon, 999., 4,true,3);

}
