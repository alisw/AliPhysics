/***************************************************************************
              adash@cern.ch - last modified on 03/04/2013

// *** Configuration script for Phi-Meson analysis with 2013 PPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigPhiPPb_TPC(  
		       AliRsnMiniAnalysisTask *task, 
		       Bool_t                 isMC, 
		       Bool_t                 isPP,
		       const char             *suffix,
		       AliRsnCutSet           *cutsPair,
		       AliPID::EParticleType  type,
		       Float_t                nsigmaKp = 2.0,
		       Bool_t                 enableMonitor = kTRUE,
		       TString                optSyt="Default"
			 )
{

  //These are the Default values for 2011 ESD track cuts
  
  Double_t  dcaxymax                  = 2.0;
  TString   PtDcaFormula              = "0.0105+0.0350/pt^1.1";//7sigma
  Double_t  dcazmax                   = 2.0;
  Double_t  minNcls                   = 50;
  Double_t  maxX2TPCcls               = 4.0;
  Double_t  maxX2ITScls               = 36.0;
  Double_t  minCrossedRows            = 70.0;
  Double_t  minRatioClsCrRowsOverFCls = 0.8;
  
  if(optSyt.Contains("PtDCAXY5s")) {PtDcaFormula = "0.0075+0.025/pt^1.1";}
  if(optSyt.Contains("PtDCAXY6s")) {PtDcaFormula = "0.0090+0.030/pt^1.1";}
  if(optSyt.Contains("PtDCAXY7s")) {PtDcaFormula = "0.0105+0.035/pt^1.1";}
  if(optSyt.Contains("PtDCAXY8s")) {PtDcaFormula = "0.0120+0.040/pt^1.1";}
  if(optSyt.Contains("PtDCAXY9s")) {PtDcaFormula = "0.0135+0.045/pt^1.1";}

  if(optSyt.Contains("FixDCAZ1")) {dcazmax = 1.;}
  if(optSyt.Contains("FixDCAZ2")) {dcazmax = 2.;}
  if(optSyt.Contains("FixDCAZ3")) {dcazmax = 3.;}
  
  if(optSyt.Contains("NCrRows60")){minCrossedRows = 60;}
  if(optSyt.Contains("NCrRows70")){minCrossedRows = 70;}
  if(optSyt.Contains("NCrRows80")){minCrossedRows = 80;}
  if(optSyt.Contains("NCrRows90")){minCrossedRows = 90;}

  if(optSyt.Contains("RClsCrRowsOvFCls0.7")){minRatioClsCrRowsOverFCls = 0.7;}
  if(optSyt.Contains("RClsCrRowsOvFCls0.8")){minRatioClsCrRowsOverFCls = 0.8;}
  if(optSyt.Contains("RClsCrRowsOvFCls0.9")){minRatioClsCrRowsOverFCls = 0.9;}

  if(optSyt.Contains("ChiSqrPerTPCCls3")) {maxX2TPCcls = 3.0;}
  if(optSyt.Contains("ChiSqrPerTPCCls4")) {maxX2TPCcls = 4.0;}
  if(optSyt.Contains("ChiSqrPerTPCCls5")) {maxX2TPCcls = 5.0;}

  if(optSyt.Contains("ChiSqrPerITSCls30")) {maxX2ITScls = 30.0;}
  if(optSyt.Contains("ChiSqrPerITSCls36")) {maxX2ITScls = 36.0;}
  if(optSyt.Contains("ChiSqrPerITSCls45")) {maxX2ITScls = 45.0;}
  
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

  AliRsnCutSet * cutSetKaon = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);
  AliRsnCutTrackQuality *fQualityTrackCut = new AliRsnCutTrackQuality("AliRsnCutTrackQuality");
  //
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  
  if(optSyt.Contains("Default")){
    ::Info("ConfigPhipPbTPC", Form("Default 2011 ESD track cuts used : %s\n",optSyt.Data()));
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
  }
  else{
    ::Info("ConfigPhipPbTPC", Form("User Defined ESD track cuts used for Sys : %s\n",optSyt.Data()));
    Int_t clusterCut = 1;
    Bool_t selPrimaries = kTRUE;  
    if(clusterCut == 0)  esdTrackCuts->SetMinNClustersTPC(50);
    else if (clusterCut == 1) {
      esdTrackCuts->SetMinNCrossedRowsTPC(minCrossedRows);
      esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioClsCrRowsOverFCls);
    }
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxX2TPCcls);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    // ITS
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					   AliESDtrackCuts::kAny);
    if(selPrimaries) {
      // 7*(0.0015+0.0050/pt^1.1)
      esdTrackCuts->SetMaxDCAToVertexXYPtDep(PtDcaFormula.Data());
      esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
    }
    //esdTrackCuts->SetMaxDCAToVertexXY(2)
    esdTrackCuts->SetMaxDCAToVertexZ(dcazmax);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxChi2PerClusterITS(maxX2ITScls);
  }
  
    fQualityTrackCut->SetESDtrackCuts(esdTrackCuts);
    fQualityTrackCut->SetPtRange(0.15,10);
    fQualityTrackCut->SetEtaRange(-0.8,0.8); 
    
  //fQualityTrackCut->Print();
  cutSetKaon->AddCut(fQualityTrackCut);
  if (!scheme.IsNull()) scheme += "&";
  scheme += fQualityTrackCut->GetName(); 

  //AliRsnCutPIDNSigma *cutKTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
  AliRsnCutPIDNSigma *cutKTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCK",type,AliRsnCutPIDNSigma::kTPC);
  cutKTPC->SinglePIDRange(nsigmaKp);
  cutSetKaon->AddCut(cutKTPC);
  if (!scheme.IsNull()) scheme += "&";
  scheme += cutKTPC->GetName();
  //Printf ("CUT Scheme for KAON is '%s'",scheme.Data());
  cutSetKaon->SetCutScheme(scheme.Data());
  
  Int_t iCutK = task->AddTrackCuts(cutSetKaon);
  
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetKaon->GetMonitorOutput());
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
  TString name    [7] = {"Unlike", "Mixing", "LikePP", "LikeMM", "MCGen" , "Trues" , "Res"   };
  TString comp    [7] = {"PAIR"  , "MIX"   , "PAIR"  , "PAIR"  , "MOTHER", "TRUE"  , "TRUE"  };
  TString output  [7] = {"SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE"};
  Char_t  charge1 [7] = {'+'     , '+'     , '+'     , '-'     , '+'     , '+'       , '+'   };
  Char_t  charge2 [7] = {'-'     , '-'     , '+'     , '-'     , '-'     , '-'       , '-'   };
  Int_t   cutID   [7] = {iCutK   , iCutK   , iCutK   , iCutK   , iCutK   , iCutK     , iCutK };
  
  for (Int_t i = 0; i < 7; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("PhiMeson_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID[i]);
    out->SetCutID(1, cutID[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdg);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 120, 0.98, 1.1);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID, 100, 0.0, 10.0);
    
    // axis Z: centrality-multiplicity
    if (!isPP)
      out->AddAxis(centID, 100, 0.0, 100.0);
    else 
      out->AddAxis(centID, 400, 0.0, 400.0);
      
    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    //out->AddAxis(yID, 90, -4.5, 4.5);
    
  }
  return kTRUE;
}
