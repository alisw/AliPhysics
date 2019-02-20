//ConfigOtonOmega called by AddTaskOtonOmega to set additional cuts or overwrite the standards

void ConfigOtonOmega(AliOtonOmegaCascadeCuts *CascadeCuts,AliFemtoDreamTrackCuts *XiPosCuts,AliFemtoDreamTrackCuts *XiNegCuts,AliFemtoDreamTrackCuts *XiBachCuts,AliOtonOmegaCascadeCuts *AntiCascadeCuts,AliFemtoDreamTrackCuts *AntiXiPosCuts,AliFemtoDreamTrackCuts *AntiXiNegCuts,AliFemtoDreamTrackCuts *AntiXiBachCuts,AliOtonOmegaCascadeCuts * CascadeOmegaCuts,AliFemtoDreamTrackCuts *OmegaPosCuts,AliFemtoDreamTrackCuts *OmegaNegCuts,AliFemtoDreamTrackCuts *OmegaBachCuts,AliOtonOmegaCascadeCuts *AntiCascadeOmegaCuts,AliFemtoDreamTrackCuts *AntiOmegaPosCuts,AliFemtoDreamTrackCuts *AntiOmegaNegCuts,AliFemtoDreamTrackCuts *AntiOmegaBachCuts){


  CascadeCuts->SetXiMassRange(1.322, 0.005, 0., 0.);
  CascadeCuts->SetCutXiDaughterDCA(.9);             //YES
  CascadeCuts->SetCutXiMinDistBachToPrimVtx(0.03);  //as ESD new vertexer cut
  CascadeCuts->SetCutXiCPA(0.991);                  //YES
  CascadeCuts->SetCutXiTransverseRadius(0.4, 200);  //as ESD new vertexer cut
  CascadeCuts->Setv0MassRange(1.116, 0.006);        //YES
  CascadeCuts->SetCutv0MaxDaughterDCA(2.);          //ESD new vertexer cut
  CascadeCuts->SetCutv0CPA(0.99);                   //YES
  CascadeCuts->SetCutv0TransverseRadius(1., 200);   //as ESD new vertexer cut
  CascadeCuts->SetCutv0MinDistToPrimVtx(0.05);      //as ESD new vertexer cut
  CascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.02);  //as ESD new vertexer cut
  CascadeCuts->SetRejectMass(1.672, 0.005, 3334);   //YES
  CascadeCuts->SetPtRangeXi(0., 999.9);
//XiPosCuts->SetCheckPileUp(false);
XiNegCuts->SetCheckPileUp(false);
//XiBachCuts->SetCheckPileUp(false)

  AntiCascadeCuts->SetXiMassRange(1.322, 0.005, 0., 0.);
  AntiCascadeCuts->SetCutXiDaughterDCA(.9);             //YES
  AntiCascadeCuts->SetCutXiMinDistBachToPrimVtx(0.03);  //as ESD new vertexer cut
  AntiCascadeCuts->SetCutXiCPA(0.991);                  //YES
  AntiCascadeCuts->SetCutXiTransverseRadius(0.4, 200);  //as ESD new vertexer cut
  AntiCascadeCuts->Setv0MassRange(1.116, 0.006);        //YES
  AntiCascadeCuts->SetCutv0MaxDaughterDCA(2.);          //ESD new vertexer cut
  AntiCascadeCuts->SetCutv0CPA(0.99);                   //YES
  AntiCascadeCuts->SetCutv0TransverseRadius(1., 200);   //as ESD new vertexer cut
  AntiCascadeCuts->SetCutv0MinDistToPrimVtx(0.05);      //as ESD new vertexer cut
  AntiCascadeCuts->SetCutv0MinDaugDistToPrimVtx(0.02);  //as ESD new vertexer cut
  AntiCascadeCuts->SetRejectMass(1.672, 0.005, 3334);   //YES
  AntiCascadeCuts->SetPtRangeXi(0., 999.9);
//AntiXiPosCuts->SetCheckPileUp(false);
AntiXiNegCuts->SetCheckPileUp(false);
//AntiXiBachCuts->SetCheckPileUp(false);


  CascadeOmegaCuts->SetXiMassRange(1.672, 0.005, 0., 0.);
  CascadeOmegaCuts->SetCutXiDaughterDCA(.8);                  //YES
  CascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(0.03);       //as ESD new vertexer cut
  CascadeOmegaCuts->SetCutXiCPA(0.995);                       //YES
  CascadeOmegaCuts->SetCutXiTransverseRadius(0.4, 200);       //as ESD new vertexer cut
  CascadeOmegaCuts->Setv0MassRange(1.116, 0.006);             //YES
  CascadeOmegaCuts->SetCutv0MaxDaughterDCA(1.2);              //YES
  CascadeOmegaCuts->SetCutv0CPA(0.99);                        //YES
  CascadeOmegaCuts->SetCutv0TransverseRadius(1., 200);        //as ESD new vertexer cut
  CascadeOmegaCuts->SetCutv0MinDistToPrimVtx(0.06);           //YES
  CascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.02);       //as ESD new vertexer cut
  CascadeOmegaCuts->SetRejectMass(1.322, 0.013, 3312);        //YES...but should be changed to asymmetric from 1.322+0.013
  CascadeOmegaCuts->SetPtRangeXi(0., 999.9);                  //as ESD new vertexer cut
//OmegaPosCuts->SetCheckPileUp(false);
OmegaNegCuts->SetCheckPileUp(false);
//OmegaBachCuts->SetCheckPileUp(false);



  AntiCascadeOmegaCuts->SetXiMassRange(1.672, 0.005, 0., 0.);
  AntiCascadeOmegaCuts->SetCutXiDaughterDCA(.8);                  //YES
  AntiCascadeOmegaCuts->SetCutXiMinDistBachToPrimVtx(0.03);       //as ESD new vertexer cut
  AntiCascadeOmegaCuts->SetCutXiCPA(0.995);                       //YES
  AntiCascadeOmegaCuts->SetCutXiTransverseRadius(0.4, 200);       //as ESD new vertexer cut
  AntiCascadeOmegaCuts->Setv0MassRange(1.116, 0.006);             //YES
  AntiCascadeOmegaCuts->SetCutv0MaxDaughterDCA(1.2);              //YES
  AntiCascadeOmegaCuts->SetCutv0CPA(0.99);                        //YES
  AntiCascadeOmegaCuts->SetCutv0TransverseRadius(1., 200);        //as ESD new vertexer cut
  AntiCascadeOmegaCuts->SetCutv0MinDistToPrimVtx(0.06);           //YES
  AntiCascadeOmegaCuts->SetCutv0MinDaugDistToPrimVtx(0.02);       //as ESD new vertexer cut
  AntiCascadeOmegaCuts->SetRejectMass(1.322, 0.013, 3312);        //YES...but should be changed to asymmetric from 1.322+0.013
//AntiOmegaPosCuts->SetCheckPileUp(false);
AntiOmegaNegCuts->SetCheckPileUp(false);
//AntiOmegaBachCuts->SetCheckPileUp(false);



}
