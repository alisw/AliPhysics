/***************************************************************************
              adash@cern.ch - last modified on 05/09/2016

// *** Configuration script for KStar-Meson analysis with 2016 PbPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigLStarpPbRunII(AliRsnMiniAnalysisTask *task, 
			    Bool_t                 isMC, 
			  
                            Bool_t                 isPP,
			    AliRsnCutSet           *cutsPair,
			    Int_t                  Strcut = 2011,
			    Int_t                  customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
			   // AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCTOFpidphikstarpPb2016,
                            AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar,
                            AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar,

			    Float_t                nsigmaP = 2.0,
			    Float_t                nsigmaK  = 2.0,
                            Float_t                nsigmatofP = 3.0,
			    Float_t                nsigmatofK  = 3.0,
                            Bool_t                 enableMonitor = kTRUE

			    )
{
  // retrieve mass from PDG database
  Int_t         pdg  = 3124;
  TDatabasePDG *db   = TDatabasePDG::Instance();
  TParticlePDG *part = db->GetParticle(pdg);
  Double_t mass      = part->Mass();

  // set daughter cuts
  AliRsnCutSetDaughterParticle* cutSetP;
  AliRsnCutSetDaughterParticle* cutSetK;
  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(!trkQualityCut) return kFALSE;

  if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,Strcut)){
    cutSetP =new AliRsnCutSetDaughterParticle(Form("cutProton%i_%2.1fsigma",cutPrCandidate,nsigmaP,nsigmatofP),trkQualityCut,cutPrCandidate,AliPID::kProton,nsigmaP,nsigmatofP);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaK,nsigmatofK),trkQualityCut,cutKaCandidate,AliPID::kKaon,nsigmaK,nsigmatofK);
  }
  else{
    printf("Doughter Track cuts has been selected =================\n");
    return kFALSE;
  }
  Int_t iCutP = task->AddTrackCuts(cutSetP);
  Int_t iCutK  = task->AddTrackCuts(cutSetK);
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
   //  gROOT->LoadMacro("/home/vikas/alice/AliPhysics/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");

    // gROOT->LoadMacro("AddMonitorOutput.C");
    //gROOT->LoadMacro("/Users/dukhishyammallick/Job_Submit/LHC16rData/WSDD_DATA/job_macro_wsdd/PIDcut/MC_6s3s2s/AddMonitorOutput.C");
      //gROOT->LoadMacro("/Users/dukhishyammallick/Job_Submit/LHC16rData/WSDD_DATA/job_macro_wsdd/PIDcut/Data_6s3s2s/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetP->GetMonitorOutput());
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

  /*
  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC    };
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,0       ,0       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  TString output  [12] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE"};
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };
  Int_t   cutIDK  [12] = {iCutK     ,iCutK     ,iCutK     ,iCutK     ,iCutK   ,iCutK   ,iCutK    ,iCutK    ,iCutK    ,iCutK    ,iCutK   ,iCutK   };
  Int_t   cutIDPi [12] = {iCutPi    ,iCutPi    ,iCutPi    ,iCutPi    ,iCutPi  ,iCutPi  ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi  ,iCutPi  };
  Int_t   PDGCode [12] = {313       ,-313      ,313       ,313       ,313     ,313     ,313      ,-313     ,313      ,-313      ,313     ,-313   };
  */



  Bool_t  use     [10] = { 1, 1, 1,1,1,1, isMC    ,  isMC   , isMC   ,  isMC  };
  Bool_t  useIM   [10] = { 1            ,  1           ,  1           ,  1          ,  1          , 1           ,  1      ,   1     , 0      ,   0  };
  TString name    [10] = {"UnlikePM"    , "UnlikeMP"   , "MixingPM"   , "MixingMP"  , "LikePP"    , "LikeMM"    ,"TruesPM","TruesMP","ResPM" ,"ResMP"};
  TString comp    [10] = {"PAIR"        , "PAIR"       , "MIX"        , "MIX"       , "PAIR"      , "PAIR"      , "TRUE"  , "TRUE"  , "TRUE" , "TRUE"  };
  TString output  [10] = {"SPARSE"      , "SPARSE"     , "SPARSE"     , "SPARSE"    , "SPARSE"    , "SPARSE"    , "SPARSE","SPARSE" ,"SPARSE","SPARSE"};
  Char_t  charge1 [10] = {'+'           , '-'          , '+'          , '-'         , '+'         , '-'         , '+'     ,  '-'    , '+'    ,  '-'   };
  Char_t  charge2 [10] = {'-'           , '+'          , '-'          , '+'         , '+'         , '-'         , '-'     ,  '+'    , '-'    ,  '+' };
  Int_t   cutID1  [10] = { iCutP        ,  iCutP       ,  iCutP       ,  iCutP      ,  iCutP      ,  iCutP      ,  iCutP  ,  iCutP  ,  iCutP ,  iCutP  };
  Int_t   cutID2  [10] = { iCutK        ,  iCutK       ,  iCutK       ,  iCutK      ,  iCutK      ,  iCutK      ,  iCutK  ,  iCutK  ,  iCutK ,  iCutK};


  for (Int_t i = 0; i < 10; i++) {
    if (!use[i]) continue;
    //  AliRsnMiniOutput *out = task->CreateOutput(Form("Lstar_%s%s", name[i].Data(), opt.Data()), output[i].Data(), comp[i].Data());
    AliRsnMiniOutput *out = task->CreateOutput(Form("CustomId%d_%s", customQualityCutsID, name[i].Data()), output[i].Data(), comp[i].Data());

    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdg);
    out->SetMotherMass(1.51953);
    out->SetPairCuts(cutsPair);


    // axis X: invmass (or resolution)                                                                                                                                 
    if (useIM[i])
      out->AddAxis(imID,120, 1.4, 2.0);
    else
      out->AddAxis(resID, 500, -0.02, 0.02);

    // axis Y: transverse momentum of pair as default - else chosen value                                                                                             
/* 
   if (yaxisVar==AliRsnMiniValue::kFirstDaughterPt) {iCutP
      out->AddAxis(fdpt, 100, 0.0, 10.0);
      out->AddAxis(sdpt, 100, 0.0, 10.0);  }
      else
      if (yaxisVar==AliRsnMiniValue::kFirstDaughterP) {
      out->AddAxis(fdp, 100, 0.0, 10.0);
      out->AddAxis(sdp, 100, 0.0, 10.0);  }
      else
*/
    out->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt                                                                                                    

    // axis Z: centrality-multiplicity                                                                                                                                 
    if (!isPP)
      out->AddAxis(centID, 100, 0.0, 100.0);
    else
      out->AddAxis(centID, 400, 0.0, 400.0);







    //else 
    //out->AddAxis(centID, 400, 0.0, 400.0);
    //axis W: pseudorapidity
    //out->AddAxis(etaID, 20, -1.0, 1.0);
    //    axis J: rapidity
    //  out->AddAxis(yID, 18, -0.9, 0.9);
  }
  return kTRUE;
}
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID = 0,Int_t trCut = 2011)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.

  if ((!trkQualityCut)){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }

  if(customQualityCutsID>=1 && customQualityCutsID<100 && customQualityCutsID!=2){
    if(trCut == 2011){
      trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
      Printf(Form("::::: SetCustomQualityCut:: using standard 2011 track quality cuts"));
    }
    else if(trCut == 2015){
      trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
      trkQualityCut->GetESDtrackCuts()->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
      Printf(Form("::::: SetCustomQualityCut:: using standard 2015 track quality cuts"));
    }
    if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.015+0.05/pt^1.1");}//10Sig // D = 7*(0.0015+0.0050/pt^1.1)
    else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.02/pt^1.1");}//4Sig
    else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(3.);}// D = 2.
    else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);} 
    else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}
    else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(5.);}// D = 4
    else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3);}
    else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}// D = 70
    else if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);}
    else if(customQualityCutsID==12){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}
    else if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}// D = 8
    else if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
    else if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(49.);}// D = 36
    else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(25.);}
    else if(customQualityCutsID==17){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(49.);}// D = 36
    else if(customQualityCutsID==18){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(25.);}
    else if(customQualityCutsID==19){trkQualityCut->GetESDtrackCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}
    
    trkQualityCut->Print();
    return kTRUE;
  }else if(customQualityCutsID==2 || (customQualityCutsID>=100 && customQualityCutsID<200)){
    trkQualityCut->SetDefaultsTPCOnly(kTRUE);
    Printf(Form("::::: SetCustomQualityCut:: using TPC-only track quality cuts"));
    if(customQualityCutsID==103){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXY(3.);}
    else if(customQualityCutsID==104){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXY(1.);}
    else if(customQualityCutsID==105){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(4.);}
    else if(customQualityCutsID==106){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);}
    else if(customQualityCutsID==107){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(7.);}
    else if(customQualityCutsID==108){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.5);}
    else if(customQualityCutsID==109){trkQualityCut->GetESDtrackCuts()->SetMinNClustersTPC(30);}
    else if(customQualityCutsID==110){trkQualityCut->GetESDtrackCuts()->SetMinNClustersTPC(85);}

    trkQualityCut->Print();
    return kTRUE;
  }else{
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }
  trkQualityCut->SetPtRange(0.15, 100000.0);
  trkQualityCut->SetEtaRange(-0.8, 0.8);
  
  Printf(Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));
  trkQualityCut->Print();
  return kTRUE;
}
