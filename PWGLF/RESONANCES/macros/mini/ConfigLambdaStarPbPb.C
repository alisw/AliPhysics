/***************************************************************************
            
/
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigLambdaStarPbPb
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC , 
    Bool_t                 isPP ,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstarPbPb2011,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstarPbPb2011,
    Float_t                nsigmaPrTPC = 3.0,
    Float_t                nsigmaKaTPC = 3.0,
    Float_t                nsigmaPrTOF = 3.0,
    Float_t                nsigmaKaTOF = 3.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetPr;
  AliRsnCutSetDaughterParticle * cutSetK;

  cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kKaon, -1.0, -1.0, aodFilterBit, kTRUE);
  cutSetQ->SetPtRange(0.2, 1.e3);
  cutSetPr = new AliRsnCutSetDaughterParticle(Form("cutPro_%2.1fsTPC_%2.1fsTOF",nsigmaPrTPC, nsigmaPrTOF), cutPrCandidate, AliPID::kProton, nsigmaPrTPC, nsigmaPrTOF, aodFilterBit, kTRUE); 
  cutSetPr->SetUse2011StdQualityCuts(kTRUE); 
  cutSetPr->SetPtRange(0.2, 1.e3);

  cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKa_%2.1fsTPC_%2.1fsTOF",nsigmaKaTPC, nsigmaKaTOF), cutKaCandidate, AliPID::kKaon, nsigmaKaTPC, nsigmaKaTOF, aodFilterBit, kTRUE);
  cutSetK->SetUse2011StdQualityCuts(kTRUE); 
  cutSetK->SetPtRange(0.2, 1.e3);



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
  // /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  // /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --
  Bool_t  use     [10] = { !IsMcTrueOnly,  !IsMcTrueOnly,  !IsMcTrueOnly,  !IsMcTrueOnly ,  !IsMcTrueOnly, !IsMcTrueOnly,  isMC   ,   isMC   ,  isMC   ,   isMC  };
  Bool_t  useIM   [10] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  1      ,   1      ,  0      ,   0      };
  TString name    [10] = {"UnlikePM", "UnlikeMP", "MixingPM", "MixingMP", "LikePP", "LikeMM", "TruesPM",  "TruesMP", "ResPM"  ,  "ResMP"  };
  TString comp    [10] = {"PAIR"   , "PAIR"   , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  ,  "TRUE"  , "TRUE"  ,  "TRUE"  };
  //TString output  [10] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };
  TString output  [10] = {"SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"  , "SPARSE"  , "SPARSE"  ,  "SPARSE"  , "SPARSE"  ,  "SPARSE"  };
  Char_t  charge1 [10] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-'     };
  Char_t  charge2 [10] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+'     };
  Int_t   cutID1  [10] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK  };
  Int_t cutID2    [10] = { iCutPr  ,  iCutPr  ,  iCutPr  ,  iCutPr  ,iCutPr  ,  iCutPr  , iCutPr  ,  iCutPr  ,  iCutPr  ,  iCutPr };
  
  for (Int_t i = 0; i < 10; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("lstar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kProton);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(3124);
    out->SetMotherMass(1.520);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 400, 1.4, 1.8);
    // axis Y: transverse momentum
    out->AddAxis(ptID, 300, 0.0, 30.0);    
    // axis Z: centrality-multiplicity
    if (!isPP)
      out->AddAxis(centID, 100, 0.0, 100.0);
    else 
      out->AddAxis(centID, 400, 0.0, 400.0); 
  }
  
  if (isMC){   
    //get mothers for L* PDG = 3124
    AliRsnMiniOutput *outm = task->CreateOutput(Form("Ls_Mother%s", suffix), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kProton);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(3124);
    outm->SetMotherMass(1.520);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID, 400, 1.4, 1.8);
    outm->AddAxis(ptID, 300, 0.0, 30.0);
    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
    
    //get mothers for antiL* PDG = -3124
    AliRsnMiniOutput *outam = task->CreateOutput(Form("antiLs_Mother%s", suffix), "SPARSE", "MOTHER");
    outam->SetDaughter(0, AliRsnDaughter::kProton);
    outam->SetDaughter(1, AliRsnDaughter::kKaon);
    outam->SetMotherPDG(-3124);
    outam->SetMotherMass(1.520);
    outam->SetPairCuts(cutsPair);
    outam->AddAxis(imID, 400, 1.4, 1.8);
    outam->AddAxis(ptID, 300, 0.0, 30.0);
    if (!isPP){
      outam->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outam->AddAxis(centID, 400, 0.0, 400.0);
    }
  }
   
  return kTRUE;
}
