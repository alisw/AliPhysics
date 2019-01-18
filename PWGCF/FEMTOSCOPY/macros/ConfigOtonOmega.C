//ConfigOtonOmega called by AddTaskOtonOmega to set additional cuts or overwrite the standards

void ConfigOtonOmega(AliOtonOmegaCascadeCuts *CascadeCuts,AliFemtoDreamTrackCuts *XiPosCuts,AliFemtoDreamTrackCuts *XiNegCuts,AliFemtoDreamTrackCuts *XiBachCuts,AliOtonOmegaCascadeCuts *AntiCascadeCuts,AliFemtoDreamTrackCuts *AntiXiPosCuts,AliFemtoDreamTrackCuts *AntiXiNegCuts,AliFemtoDreamTrackCuts *AntiXiBachCuts,AliOtonOmegaCascadeCuts * CascadeOmegaCuts,AliFemtoDreamTrackCuts *OmegaPosCuts,AliFemtoDreamTrackCuts *OmegaNegCuts,AliFemtoDreamTrackCuts *OmegaBachCuts,AliOtonOmegaCascadeCuts *AntiCascadeOmegaCuts,AliFemtoDreamTrackCuts *AntiOmegaPosCuts,AliFemtoDreamTrackCuts *AntiOmegaNegCuts,AliFemtoDreamTrackCuts *AntiOmegaBachCuts){


  CascadeCuts->SetXiMassRange(1.322, 0.020);
//  XiPosCuts->SetCheckPileUp(false);
//  XiNegCuts->SetCheckPileUp(false);
//  XiBachCuts->SetCheckPileUp(false);

  AntiCascadeCuts->SetXiMassRange(1.322, 0.020);
//  AntiXiPosCuts->SetCheckPileUp(false);
//  AntiXiNegCuts->SetCheckPileUp(false);
//  AntiXiBachCuts->SetCheckPileUp(false);

  CascadeOmegaCuts->SetXiMassRange(1.672, 0.020);
//  OmegaPosCuts->SetCheckPileUp(false);
//  OmegaNegCuts->SetCheckPileUp(false);
//  OmegaBachCuts->SetCheckPileUp(false);

  AntiCascadeOmegaCuts->SetXiMassRange(1.672, 0.020);
//  AntiOmegaPosCuts->SetCheckPileUp(false);
//  AntiOmegaNegCuts->SetCheckPileUp(false);
//  AntiOmegaBachCuts->SetCheckPileUp(false);

}
