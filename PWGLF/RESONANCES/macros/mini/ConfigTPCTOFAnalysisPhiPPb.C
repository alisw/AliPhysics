/***************************************************************************
 * 	   edgar.perez.lezama@cern.ch - last modified on 28/06/2015
 * 
// *** Configuration script for Phi-Meson analysis with 2013 PPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigTPCTOFAnalysisPhiPPb
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    AliPID::EParticleType  type1,
    Float_t                nsigmaKa = 3.0,
    Bool_t                 enableMonitor = kTRUE
)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // retrieve mass from PDG database
   Int_t         pdg  = 333;
   TDatabasePDG *db   = TDatabasePDG::Instance();
   TParticlePDG *part = db->GetParticle(pdg);
   Double_t mass      = part->Mass();

  TString opt = "PPb";
  TString scheme="";  
  TString cutname = "K_Phi";
  if (!opt.IsNull()) cutname += Form("_%s",opt.Data());


  AliRsnCutSet * cutSetKp = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

  AliRsnCutSet * cutSetTPC = new AliRsnCutSet("set_TPC", AliRsnTarget::kDaughter);

  AliRsnCutSet * cutSetTPCTOF = new AliRsnCutSet("set_TPCTOF", AliRsnTarget::kDaughter);
  

  
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts"); 
  esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
  esdTrackCuts->SetMinNCrossedRowsTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.01");   //"0.0105+0.0350/pt^1.1" = 7sigma
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  
  AliRsnCutTrackQuality *fQualityTrackCut = new AliRsnCutTrackQuality("AliRsnCutTrackQuality");
  fQualityTrackCut->SetDefaults2010();
  fQualityTrackCut->SetESDtrackCuts(esdTrackCuts);
  fQualityTrackCut->SetPtRange(0.15,10);
  fQualityTrackCut->SetTPCmaxChi2(4.0); 
  fQualityTrackCut->SetRejectKinkDaughters();
  fQualityTrackCut->SetDCAZmax(2.0);
  

  cutSetKp->AddCut(fQualityTrackCut);
  cutSetTPC->AddCut(fQualityTrackCut);
  cutSetTPCTOF->AddCut(fQualityTrackCut);


  AliRsnCutTOFMatch  *cutTOFMatch = new AliRsnCutTOFMatch("cutTOFMatch");


  
  AliRsnCutPIDNSigma *cutKTOF = new AliRsnCutPIDNSigma("cutNSigmaTOFK",type1,AliRsnCutPIDNSigma::kTOF);
  cutKTOF->SinglePIDRange(nsigmaKa);
  cutSetKp->AddCut(cutTOFMatch);	
  cutSetKp->AddCut(cutKTOF);
  if (!scheme.IsNull()) scheme += "&";
  scheme += cutKTOF->GetName();
  scheme = "AliRsnCutTrackQuality & cutTOFMatch & cutNSigmaTOFK";//&cutNSigmaTOFK";
  

  AliRsnCutPIDNSigma *cutKTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCK",type1,AliRsnCutPIDNSigma::kTPC);
  cutKTPC->SinglePIDRange(nsigmaKa);
  cutSetTPC->AddCut(cutKTPC);	


  
  AliRsnCutPIDNSigma *cutKNoTOF = new AliRsnCutPIDNSigma("cutNSigmaNoTOFK",type1,AliRsnCutPIDNSigma::kTOF);
  cutKNoTOF->SinglePIDRange(1000);
  cutSetKp->AddCut(cutKNoTOF);


  cutSetTPCTOF->AddCut(cutKTPC);
  cutSetTPCTOF->AddCut(cutKTOF);
  cutSetTPCTOF->AddCut(cutKNoTOF);





  if (scheme.IsNull()) Fatal("ConfigTPCTOFAnalysisPhiPPb.C","No Scheme");
  Printf ("CUT Scheme for KAON is '%s'",scheme.Data());
  cutSetKp->SetCutScheme(scheme.Data());
//   cutSetTPC->SetCutScheme("AliRsnCutTrackQuality & cutNSigmaTPCK ");
  cutSetTPCTOF->SetCutScheme("AliRsnCutTrackQuality & (cutNSigmaTOFK | (cutNSigmaTPCK & !cutNSigmaNoTOFK)) ");


  
  Int_t iCutKp = task->AddTrackCuts(cutSetKp);
//   Int_t iCutTPC= task->AddTrackCuts(cutSetTPC);
  Int_t iCutTPCTOF = task->AddTrackCuts(cutSetTPCTOF);

  
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetKp->GetMonitorOutput());
  }  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [7] = { 1      , 1       , 1       , 1       , isMC    , isMC    , isMC    };
  Bool_t  useIM   [7] = { 1      , 1       , 1       , 1       , 1       , 1       , 0       };
  TString name    [7] = {"Unlike", "Mixing", "LikePP", "LikeMM", "MCGen" , "Trues" , "Res"};
  TString comp    [7] = {"PAIR"  , "MIX"   , "PAIR"  , "PAIR"  , "MOTHER", "TRUE"  , "TRUE"  };
  TString output  [7] = {"SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE"}; 
  Char_t  charge1 [7] = {'+'     , '+'     , '+'     , '-'     , '+'     , '+'       , '+'   };
  Char_t  charge2 [7] = {'-'     , '-'     , '+'     , '-'     , '-'     , '-'       , '-'   };


  Int_t   cut_TPCTOF   [7] = {iCutTPCTOF  , iCutTPCTOF  , iCutTPCTOF  , iCutTPCTOF  , iCutTPCTOF  , iCutTPCTOF    , iCutTPCTOF};

  for (Int_t i = 0; i < 7; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cut_TPCTOF[i]);
    out->SetCutID(1, cut_TPCTOF[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdg);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) out->AddAxis(imID, 120, 0.98, 1.1);             //(i, nbins, min, max);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    out->AddAxis(ptID, 100, 0.0, 10.0);
    if (!isPP)
      out->AddAxis(centID, 100, 0.0, 100.0);    //Centrality
    else 
      out->AddAxis(centID, 400, 0.0, 400.0);
    
    out->AddAxis(yID, 400, -1.0, 1.0);
    
  }

  //
  // -- Create output for MC generated ------------------------------------------------------------
  //
  if (isMC) {
    // create ouput
    AliRsnMiniOutput *outMC = task->CreateOutput("TruePair", "SPARSE", "TRUE");
    // selection settings
    outMC->SetDaughter(0, AliRsnDaughter::kKaon);
    outMC->SetDaughter(1, AliRsnDaughter::kKaon);
    outMC->SetMotherPDG(pdg);
    outMC->SetMotherMass(mass);
    // pair cuts
    if (cutsPair) outMC->SetPairCuts(cutsPair);
    // axis X: invmass
    outMC->AddAxis(imID, 120, 0.98, 1.1);
    
    // axis Y: transverse momentum
    outMC->AddAxis(ptID, 100, 0.0, 10.0);
    
    // axis Z: centrality-multiplicity
    if (!isPP)
      outMC->AddAxis(centID, 100, 0.0, 100.0);
    else 
      outMC->AddAxis(centID, 400, 0.0, 400.0);
    
    outMC->AddAxis(yID, 400, -1.0, 1.0);
  }


  return kTRUE;
}

