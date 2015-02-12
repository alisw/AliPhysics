/***************************************************************************
              fbellini@cern.ch - last modified on 06/08/2012

// *** Configuration script for K*, anti-K* analysis with 2010 PbPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigTPCanalysisKStarTest
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Float_t                nsigmaPi = 2.0,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    Int_t                  Pdg = 313
   )
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  TString opt = "PbPb";
  TString schemePi="";  
  TString schemeK="";  
  TString cutnameK = "K_KS";
  TString cutnamePi = "Pi_KS";
  if (!opt.IsNull()) cutnameK += Form("_%s",opt.Data());
  if (!opt.IsNull()) cutnamePi += Form("_%s",opt.Data());

  AliRsnCutSet * cutSetKaon = new AliRsnCutSet(cutnameK.Data(), AliRsnTarget::kDaughter);
  AliRsnCutSet * cutSetPion = new AliRsnCutSet(cutnamePi.Data(), AliRsnTarget::kDaughter);

  AliRsnCutTrackQuality *fQualityTrackCut = new AliRsnCutTrackQuality("AliRsnCutTrackQuality");

  //
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;

  Int_t clusterCut = 1;
  Bool_t selPrimaries = kTRUE;  
  if(clusterCut == 0)  esdTrackCuts->SetMinNClustersTPC(50);
  else if (clusterCut == 1) {
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  }
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  if(selPrimaries) {
    // 7*(0.0015+0.0050/pt^1.1)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
    esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  }
  esdTrackCuts->SetMaxDCAToVertexZ(2);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMaxChi2PerClusterITS(36);

  fQualityTrackCut->SetESDtrackCuts(esdTrackCuts);

  //PID selection
  cutSetKaon->AddCut(fQualityTrackCut);
  if (!schemeK.IsNull()) schemeK += "&";
  schemeK += fQualityTrackCut->GetName(); 
  
  cutSetPion->AddCut(fQualityTrackCut);
  if (!schemePi.IsNull()) schemePi += "&";
  schemePi += fQualityTrackCut->GetName(); 



  AliRsnCutPIDNSigma *cutPiTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCPi",AliPID::kPion,AliRsnCutPIDNSigma::kTPC);
  cutPiTPC->SinglePIDRange(nsigmaPi);
  cutSetPion->AddCut(cutPiTPC);
  if (!schemePi.IsNull()) schemePi += "&";
  schemePi += cutPiTPC->GetName();

  AliRsnCutPIDNSigma *cutKTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
  cutKTPC->SinglePIDRange(nsigmaKa);
  cutSetKaon->AddCut(cutKTPC);
  if (!schemeK.IsNull()) schemeK += "&";
  schemeK += cutKTPC->GetName();

  Printf ("CUT Scheme for KAON is '%s'",schemeK.Data());
  Printf ("CUT Scheme for PION is '%s'",schemePi.Data());

  cutSetPion->SetCutScheme(schemePi.Data());
  cutSetKaon->SetCutScheme(schemeK.Data());

  Int_t iCutPi = task->AddTrackCuts(cutSetPion);
  Int_t iCutK = task->AddTrackCuts(cutSetKaon);

  
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPion->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetKaon->GetMonitorOutput());
  }  
  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --
  Bool_t  use     [10] = { !IsMcTrueOnly,  !IsMcTrueOnly,  !IsMcTrueOnly,  !IsMcTrueOnly ,  !IsMcTrueOnly, !IsMcTrueOnly,  isMC   ,   isMC   ,  isMC   ,   isMC   };
  Bool_t  useIM   [10] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  1      ,   1      ,  0      ,   0      };
  TString name    [10] = {"UnlikePM", "UnlikeMP", "MixingPM", "MixingMP", "LikePP", "LikeMM", "TruesPM",  "TruesMP", "ResPM"  ,  "ResMP"  };
  TString comp    [10] = {"PAIR"   , "PAIR"   , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  ,  "TRUE"  , "TRUE"  ,  "TRUE"  };
  //TString output  [10] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };
  TString output  [10] = {"SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"  , "SPARSE"  , "SPARSE"  ,  "SPARSE"  , "SPARSE"  ,  "SPARSE"  };
  Char_t  charge1 [10] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-'     };
  Char_t  charge2 [10] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+'     };
  Int_t   cutID1  [10] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK  };
  Int_t   cutID2  [10] = { iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi ,  iCutPi ,  iCutPi ,   iCutPi ,  iCutPi ,   iCutPi };
  
  for (Int_t i = 0; i < 10; i++) {
    if (!use[i]) continue;
    if(Pdg > 0) AliRsnMiniOutput *out = task->CreateOutput(Form("kstar1_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    if(Pdg < 0) AliRsnMiniOutput *out = task->CreateOutput(Form("kstar2_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(Pdg);//313
    out->SetMotherMass(0.89594);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 90, 0.6, 1.5);
    //else
    //out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID, 300, 0.0, 30.0);
    
    // axis Z: centrality-multiplicity
    if (!isPP)
      out->AddAxis(centID, 100, 0.0, 100.0);
    else 
      out->AddAxis(centID, 400, 0.0, 400.0);
    
    // axis W: pseudorapidity
    //    out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    //out->AddAxis(yID, 32, -0.8, 0.8);
    
  }
  
  if (isMC){
    // create output
    
    if(Pdg > 0) {AliRsnMiniOutput *outm = task->CreateOutput(Form("kstar_Mother1%s", suffix), "SPARSE", "MOTHER");}
    if(Pdg < 0) {AliRsnMiniOutput *outm = task->CreateOutput(Form("kstar_Mother2%s", suffix), "SPARSE", "MOTHER");}
    outm->SetDaughter(0, AliRsnDaughter::kKaon);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(Pdg);//313
    outm->SetMotherMass(0.89594);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 90, 0.6, 1.5);
    outm->AddAxis(ptID, 300, 0.0, 30.0);
    if (!isPP){
    outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else { 
     outm->AddAxis(centID, 400, 0.0, 400.0);
    }
    //outm->AddAxis(yID, 32, -0.8, 0.8);
  }
  return kTRUE;
}
