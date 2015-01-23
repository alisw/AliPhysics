/***************************************************************************
              fbellini@cern.ch - last modified on 06/08/2012

// *** Configuration script for Lambda*, anti-Lambda* analysis with 2011 PbPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigLambdaStarPbPbMC
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidKstarPP2010,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidKstarPP2010,
    Float_t                nsigmaPr = 2.0,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    Int_t                  Pdg = 3124,
    Int_t                  aodN = 0
)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetPr;
  AliRsnCutSetDaughterParticle * cutSetK;
  
  cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kKaon, -1.0, aodFilterBit);
  cutSetPr = new AliRsnCutSetDaughterParticle(Form("cutProtonTPCpp2010_%2.1fsigma",nsigmaPr), cutPrCandidate, AliPID::kProton, nsigmaPr, aodFilterBit);
  cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaonTPCpp2010_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit);
  
  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutPr = task->AddTrackCuts(cutSetPr);
  Int_t iCutK = task->AddTrackCuts(cutSetK);
  
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetPr->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput());
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
  Int_t   cutID2  [10] = { iCutPr  ,  iCutPr  ,  iCutPr  ,  iCutPr  ,  iCutPr ,  iCutPr ,  iCutPr ,   iCutPr ,  iCutPr ,   iCutPr };
  
  for (Int_t i = 0; i < 10; i++) {
    if (!use[i]) continue;
    if(Pdg > 0) AliRsnMiniOutput *out = task->CreateOutput(Form("lstar1_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    if(Pdg < 0) AliRsnMiniOutput *out = task->CreateOutput(Form("lstar2_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kProton);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(Pdg);//2122
    out->SetMotherMass(1.520);//GeV
    out->SetPairCuts(cutsPair);
    
    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 900, 1.3, 2.2);
    //else
    //out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID, 500, 0.0, 50.0);
    
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
    
    if(Pdg > 0) {AliRsnMiniOutput *outm = task->CreateOutput(Form("lstar_Mother1%s", suffix), "SPARSE", "MOTHER");}
    if(Pdg < 0) {AliRsnMiniOutput *outm = task->CreateOutput(Form("lstar_Mother2%s", suffix), "SPARSE", "MOTHER");}
    outm->SetDaughter(0, AliRsnDaughter::kKaon);
    outm->SetDaughter(1, AliRsnDaughter::kProton);
    outm->SetMotherPDG(Pdg);//313
    outm->SetMotherMass(1.520);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 900, 1.3, 2.2);
    outm->AddAxis(ptID, 500, 0.0, 50.0);
    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
    //outm->AddAxis(yID, 32, -0.8, 0.8);
  }
  return kTRUE;
}
