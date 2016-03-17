/***************************************************************************
               adash@cern.ch - last modified on 03/04/2013
// Modified by skundu@cern.ch 
// *** Configuration script for Kstar-Meson analysis with pp@13TeV runs ***
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigKstarPP13TeV(  
			AliRsnMiniAnalysisTask *task, 
			Bool_t                 isMC, 
			Bool_t                 isPP,
			const char             *suffix,
			AliRsnCutSet           *cutsPair,
			Float_t                nsigmaPi = 3.0,
			Float_t                nsigmaK  = 3.0,
			Bool_t                 enableMonitor = kTRUE,
			TString                optSys = "DefaultITSTPC2011"
			  )
{
  
  //These are the Default values for 2011 ESD track cuts
  
  AliPID::EParticleType  type1   = AliPID::kPion;
  AliPID::EParticleType  type2   = AliPID::kKaon;

  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // retrieve mass from PDG database
  Int_t         pdg  = 313;
  TDatabasePDG *db   = TDatabasePDG::Instance();
  TParticlePDG *part = db->GetParticle(pdg);
  Double_t mass      = part->Mass();
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle* cutSetPi;
  AliRsnCutSetDaughterParticle* cutSetK;

  /* AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(optSyt.Contains("DefaultITSTPC2011")){ 
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    Printf(Form("::::: SetCustomQualityCut:: using standard 2011 track quality cuts"));
  }
  else if(optSyt.Contains("DefaultTPCOnly")){
    trkQualityCut->SetDefaultsTPCOnly(kTRUE);
    Printf(Form("::::: SetCustomQualityCut:: using TPC-only track quality cuts"));
  }
  trkQualityCut->Print();*/
  
   AliRsnCutTrackQuality *fQualityTrackCut = new AliRsnCutTrackQuality("AliRsnCutTrackQuality");

  //Analysis Track cuts are implemented here
  AliESDtrackCuts * esdTrackCuts = MyTrackCuts(1, kTRUE,optSys.Data());
  fQualityTrackCut->SetESDtrackCuts(esdTrackCuts);
  fQualityTrackCut->SetPtRange(0.15,30);
  fQualityTrackCut->SetEtaRange(-0.8,0.8);   
  //kTPCpidTOFveto4s
  //kTPCpidTOFveto3s
  //kFastTPCpidNsigma
  cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,nsigmaPi),fQualityTrackCut,AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,AliPID::kPion,nsigmaPi);
  cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, nsigmaK),fQualityTrackCut,AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,AliPID::kKaon,nsigmaK);
  
  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutK  = task->AddTrackCuts(cutSetK);
  
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput());
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

  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC    };
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,0       ,0       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  TString output  [12] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE"};
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };
  Int_t   cutIDPi [12] = {iCutPi    ,iCutPi    ,iCutPi    ,iCutPi    ,iCutPi  ,iCutPi  ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi  ,iCutPi  };
  Int_t   cutIDK  [12] = {iCutK     ,iCutK     ,iCutK     ,iCutK     ,iCutK   ,iCutK   ,iCutK    ,iCutK    ,iCutK    ,iCutK    ,iCutK   ,iCutK   };
  Int_t   PDGCode [12] = {313       ,313       ,313       ,313       ,313     ,313     ,313      ,-313     ,313      ,313      ,313     ,-313    };

  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("KSTAR_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCutID(0, cutIDK[i]);
    out->SetCutID(1, cutIDPi[i]);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 90, 0.6, 1.5);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID, 200, 0.0, 20.0);
    
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

AliESDtrackCuts* MyTrackCuts(Int_t clusterCut = 1,  Bool_t selPrimaries = kTRUE, TString optSyt="DefaultITSTPC2011")
{
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
  if(optSyt.Contains("PtDCAXY7s")) {PtDcaFormula = "0.0105+0.035/pt^1.1";}//Defult
  if(optSyt.Contains("PtDCAXY8s")) {PtDcaFormula = "0.0120+0.040/pt^1.1";}
  if(optSyt.Contains("PtDCAXY9s")) {PtDcaFormula = "0.0135+0.045/pt^1.1";}

  if(optSyt.Contains("FixDCAZ1")) {dcazmax = 1.0;}
  if(optSyt.Contains("FixDCAZ2")) {dcazmax = 2.0;}//Defult
  if(optSyt.Contains("FixDCAZ3")) {dcazmax = 3.0;}
  
  if(optSyt.Contains("NCrRows60")){minCrossedRows = 60;}
  if(optSyt.Contains("NCrRows70")){minCrossedRows = 70;}//Defult
  if(optSyt.Contains("NCrRows80")){minCrossedRows = 80;}
  if(optSyt.Contains("NCrRows90")){minCrossedRows = 90;}

  if(optSyt.Contains("RClsCrRowsOvFCls0.7")){minRatioClsCrRowsOverFCls = 0.7;}
  if(optSyt.Contains("RClsCrRowsOvFCls0.8")){minRatioClsCrRowsOverFCls = 0.8;}//Defult
  if(optSyt.Contains("RClsCrRowsOvFCls0.9")){minRatioClsCrRowsOverFCls = 0.9;}

  if(optSyt.Contains("ChiSqrPerTPCCls3")) {maxX2TPCcls = 3.0;}
  if(optSyt.Contains("ChiSqrPerTPCCls4")) {maxX2TPCcls = 4.0;}//Defult
  if(optSyt.Contains("ChiSqrPerTPCCls5")) {maxX2TPCcls = 5.0;}

  if(optSyt.Contains("ChiSqrPerITSCls30")) {maxX2ITScls = 30.0;}
  if(optSyt.Contains("ChiSqrPerITSCls36")) {maxX2ITScls = 36.0;}//Defult
  if(optSyt.Contains("ChiSqrPerITSCls45")) {maxX2ITScls = 45.0;}
  


  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  if(optSyt.Contains("DefaultITSTPC2011")){
    ::Info("Config KSTAR ", Form("Default 2011 ESD track cuts used : %s\n",optSyt.Data()));
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
    ::Info("Config KSTAR ", Form("User Defined ESD track cuts used for Sys : %s\n",optSyt.Data()));
    if(clusterCut == 0)  esdTrackCuts->SetMinNClustersTPC(minNcls);
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
    esdTrackCuts->SetMaxDCAToVertexZ(dcazmax);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxChi2PerClusterITS(maxX2ITScls);
  }
  return esdTrackCuts;
}


