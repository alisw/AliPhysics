AliAnalysisTaskFlowQnSPCascade *AddTaskFlowEPCascade(Int_t harmonic,Float_t centrMin, Float_t centrMax,Bool_t WeightCorr)
{


  TString suffixName="mySuffix";
  Bool_t fStrange = kTRUE;     //
  AliAnalysisTaskFlowQnSPCascade * task
    = new AliAnalysisTaskFlowQnSPCascade(Form("v2EP_%s",
                                            suffixName.Data()),
                                       centrMin, centrMax,WeightCorr);
  if (fStrange)
  {  
      task->SetAnaObjectMultiStrange();
  }
  else
  {
      task->SetAnaObjectStrange();  
  }
   task->SetHarmonicOrder(harmonic);
   task->SetUseHybridGlobalTrack();
   task->SetSepAnalysis(kFALSE);//Distinguish the charge of Xi and Omega.Also distinguish lambda and anti-lambda
//-------------------Xi------------------------------------------------
   task->SetPrimaryTrackEta(0.8);
   task->SetXiEtaMin(-0.5);
   task->SetXiEtaMax(0.5);
   task->SetV0RadiusXiMin(0.9);
   task->SetV0RadiusXiMax(100);
   task->SetXiRadiusMin(0.9);
   task->SetXiRadiusMax(100);
   task->SetDCAXiDaughtersMax(0.3);
   task->SetXiCosOfPointingAngleMin(0.999);
   task->SetDCAV0ToPrimaryVtxXiMin(0.05);
   task->SetDCABachToPrimaryVtxXiMin(0.03);
   task->SetLambdaMassWindow(0.008);
   task->SetDCAV0DaughtersXi(1);
   task->SetV0CosOfPointingAngleXiMin(0.998);
   task->SetDCAPosToPrimaryVtxXiMin(0.1);
   task->SetDCANegToPrimaryVtxXiMin(0.1);
//--------------------V0---------------------------------------------------
//   task->SetV0Eta(0.8);
//   task->SetV0PtMin(0.2);

//   task->SetV0DecayLength(3);
   task->SetV0Rapidity(0.5);
   task->SetV0DecayRadius(5.0);
   task->SetV0DCADaughtersMax(1);
   task->SetV0CosinePointingAngleMin(0.998);
   task->SetV0DCAToPrimVertexMin(0.1);
   task->SetV0LifeTimeMax(15);
   task->SetV0DecayLength(3);

//-------------------track-------------------------------------------------  
   task->SetRPTrackFromTPC(2);
   task->SetTPCNcls(70);
   task->SetTrackEta(0.8);
   task->SetTrackPtMin(0.15);
//-------------------PID-------------------------------------------------        
   task->SetXiPIDSigma(3);
   task->SetV0PIDSigma(3);

  return task;
}

