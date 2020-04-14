/***************************************************************************
// rama.chandra.baral@cern.ch & sarita.sahoo@cern.ch - last modified on 12/06/2014
// adapted for Pb-Pb analysis: Neelima Agrawal, Roberto Preghenella on 11/04/2016
//
// *** Configuration script for L*, anti-L*, syst. analysis for pp and p-Pb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigLStar_PbPb2018(
			    AliRsnMiniAnalysisTask *task, 
			    Float_t                 yCut              = 0.5,
			    Int_t                   aodFilterBit      = 5,
			    Bool_t                  useTPCCrossedRows = kTRUE,
			    Int_t                   qualityCut        = AliRsnCutSetDaughterParticle::kQualityStd2011,
			    Int_t                   pidCut            = AliRsnCutSetDaughterParticle::kTPCTOFpidTunedPbPbTOFneed_2018,
			    Float_t                 nsPr              = 1.0, // factor wrt. default n-sigma
			    Float_t                 nsKa              = 1.0, // factor wrt. default n-sigma
			    Bool_t                  isMC              = kFALSE, 
			    const char             *suffix            = ""
			    )
{

  //
  printf("******************************** \n");
  printf(" ConfigLStar_PbPb2018 \n");
  printf("******************************** \n");
  printf(" yCut              = %g \n", yCut);
  printf(" aodFilterBit      = %d \n", aodFilterBit);
  printf(" useTPCCrossedRows = %d \n", useTPCCrossedRows);
  printf(" qualityCut        = %d \n", qualityCut);
  printf(" pidCut            = %d \n", pidCut);
  printf(" nsPr              = %g \n", nsPr);
  printf(" nsKa              = %g \n", nsKa);
  printf(" isMC              = %d \n", isMC);
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n");    
  
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);

  // -- Cuts ---------------------------------------------------------------------------------------

  // set mother cuts

  AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutY", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-yCut, yCut);      
  AliRsnCutSet *cutM = new AliRsnCutSet("cutMother", AliRsnTarget::kMother);
  cutM->AddCut(cutY);
  cutM->SetCutScheme(cutY->GetName());

  // set daughter cuts
  
  // quality cut
  AliRsnCutSetDaughterParticle *cutQ = new AliRsnCutSetDaughterParticle("cutQuality", 
									qualityCut, 
									AliPID::kPion, 
									-1.0, 
									aodFilterBit, 
									useTPCCrossedRows);
  task->AddTrackCuts(cutQ);
  // proton cut
  AliRsnCutSetDaughterParticle *cutP = new AliRsnCutSetDaughterParticle("cutProton", 
									pidCut, 
									AliPID::kProton, 
									nsPr,
									aodFilterBit, 
									useTPCCrossedRows);
  Int_t icutP = task->AddTrackCuts(cutP);
  // kaon cut
  AliRsnCutSetDaughterParticle *cutK = new AliRsnCutSetDaughterParticle("cutKaon", 
									pidCut, 
									AliPID::kKaon, 
									nsKa,
									aodFilterBit, 
									useTPCCrossedRows);
  Int_t icutK = task->AddTrackCuts(cutK);

  // -- Monitoring --------------------------------------------------------------------------------

  Printf("======== Cut monitoring enabled");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
  AddMonitorOutput(isMC, cutQ->GetMonitorOutput(), "");
  AddMonitorOutput(isMC, cutP->GetMonitorOutput(), "");
  AddMonitorOutput(isMC, cutK->GetMonitorOutput(), "");

  // -- Values ------------------------------------------------------------------------------------

  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  //  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  //  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  //  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  //  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  //  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  //  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);
  //  /* Dip Angle        */ Int_t OpAn   = task->CreateValue(AliRsnMiniValue::kDipAngle, kFALSE);
  //  /* kPtRatio         */ Int_t PtRat  = task->CreateValue(AliRsnMiniValue::kPtRatio, kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------

  Bool_t  use       [14] = { !isMC      ,  !isMC     ,  !isMC   , !isMC    , !isMC      , !isMC      , !isMC      , !isMC      , isMC     , isMC     , isMC     , isMC     , isMC     , isMC     };
  TString name      [14] = { "UnlikePM" , "UnlikeMP" , "LikePP" , "LikeMM" , "MixingPM" , "MixingMP" , "MixingPP" , "MixingMM" , "RecoPM" , "RecoMP" , "ResoPM" , "ResoMP" , "TruePM" , "TrueMP" };
  TString comp      [14] = { "PAIR"     , "PAIR"     , "PAIR"   , "PAIR"   , "MIX"      , "MIX"      , "MIX"      , "MIX"      , "TRUE"   , "TRUE"   , "TRUE"   , "TRUE"   , "MOTHER" , "MOTHER" };
  Char_t  charge1   [14] = { '+'        , '-'        , '+'      , '-'      , '+'        , '-'        , '+'        , '-'        , '+'      , '-'      , '+'      , '-'      , 0        , 0        };
  Char_t  charge2   [14] = { '-'        , '+'        , '+'      , '-'      , '-'        , '+'        , '+'        , '-'        , '-'      , '+'      , '-'      , '+'      , 0        , 0        };
  Int_t   motherPDG [14] = { 0          , 0          , 0        , 0        , 0          , 0          , 0          , 0          , 3124     , -3124    , 3124     , -3124    , 3124     , -3124    };
  Int_t   xaxis     [14] = { imID       , imID       , imID     , imID     , imID       , imID       , imID       , imID       , imID     , imID     , resID    , resID    , imID     , imID     };
  Int_t   cutID1         =   icutP   ;
  Int_t   cutID2         =   icutK   ;
  TString output         =  "SPARSE" ;
  
  for (Int_t i = 0; i < 14; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("Lstar_PbPb_%s%s", name[i].Data(), suffix), output.Data(), comp[i].Data());
    //
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetMotherPDG(motherPDG[i]);
    out->SetMotherMass(1.51953);
    out->SetPairCuts(cutM);
    //
    if (!name[i].Contains("True")) {
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetCutID(0, cutID1);
      out->SetCutID(1, cutID2);
    }
    //
    // axis X: invariant mass / invariant-mass resolution
    if (name[i].Contains("Reso")) out->AddAxis(resID, 400, -0.01, 0.01);
    else                          out->AddAxis(imID, 400, 1.4, 2.4);
    // axis Y: transverse momentum of pair
    out->AddAxis(ptID, 40, 0.0, 10.0);
    // axis Z: centrality
    out->AddAxis(centID, 100, 0.0, 100.0);
  }   

  return kTRUE;
}
