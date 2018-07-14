/***************************************************************************
  priyanka.sett@cern.ch - last modified on 27/10/2016
  himani.bhatt@cern.ch  - last modified on 13/07/2018
// *** Configuration script for L*, anti-L*, syst. analysis for 13 TeV pp data  ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigureLstar13TeVpp_II
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    Int_t                  customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,
    Float_t                nsigmaPr = 2.0,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    Int_t                  signedPdg = 3124,
    TString                monitorOpt = "",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
    Bool_t                 useCrossedRows = kTRUE,
    const char*            yaxisVar = "",  //yaxisVar = "PtDaughter_PDaughter_cent"
    Bool_t                 useMixLS = 0
)
{
  // TString opt(yaxisVar);
  //opt.ToUpper();

  TString opt(suffix);
  if (strlen(suffix) > 0) opt = Form("_%s", suffix);
  
  Bool_t isPtDaughter  = opt.Contains("PTDAUGHTER") ;
  Bool_t isPDaughter   = opt.Contains("PDAUGHTER") ;
  Bool_t iscent        = opt.Contains("CENT") ;
  Bool_t iseta         = opt.Contains("ETA") ;
  Bool_t israpidity    = opt.Contains("RAPIDITY") ; 


  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetP;
  AliRsnCutSetDaughterParticle * cutSetK;


  AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("myQualityCut");
  if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit)){
    //Set custom quality cuts for systematic checks
  
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQuality_bit%i",aodFilterBit),trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2011,AliPID::kPion,-1.);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutKaon_%i_%2.1fsigma",cutKaCandidate, nsigmaKa),trkQualityCut,cutKaCandidate,AliPID::kKaon,nsigmaKa);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton_%i_%2.1fsigma",cutPrCandidate, nsigmaPr),trkQualityCut,cutPrCandidate,AliPID::kProton,nsigmaPr);
    
  }
  
  else{
    
    Bool_t useCrossedRows = 1;
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit),AliRsnCutSetDaughterParticle::kQualityStd2011,AliPID::kPion,-1.,aodFilterBit,useCrossedRows);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate,nsigmaKa),cutKaCandidate,AliPID::kKaon,nsigmaKa,aodFilterBit,useCrossedRows);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutPr%i_%2.1fsigma",cutPrCandidate,nsigmaPr),cutKaCandidate,AliPID::kProton,nsigmaPr,aodFilterBit,useCrossedRows);
  }

    

  
  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutP = task->AddTrackCuts(cutSetP);
  Int_t iCutK = task->AddTrackCuts(cutSetK);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetP->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput(), monitorOpt.Data());
  }  
  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);

  /* Dip Angle        */ Int_t OpAn   = task->CreateValue(AliRsnMiniValue::kDipAngle, kFALSE);
  /* kPtRatio         */ Int_t PtRat  = task->CreateValue(AliRsnMiniValue::kPtRatio, kFALSE);


  /* kptresolution    */ Int_t respT  = task->CreateValue(AliRsnMiniValue::kPairPtRes, kFALSE);
  /* kmassresolution  */ Int_t resmass  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kFALSE);
  /* IM diff          */ Int_t IMdif   = task->CreateValue(AliRsnMiniValue::kInvMassDiff, kTRUE);


  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,useMixLS    ,useMixLS};
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,1       ,1       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  //TString output  [12] = {"HIST"    ,"HIST"    ,"HIST"    ,"HIST"    ,"HIST"  ,"HIST"  ,"HIST"   ,"HIST"   ,"HIST"   ,"HIST"   ,"HIST"  ,"HIST"  };
  // TString output  [12] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE"};
  TString output  [12] = {"HIST"  ,"HIST"  ,"HIST"  ,"HIST"  ,"HIST","HIST","HIST" ,"HIST" ,"HIST" ,"HIST" ,"HIST","HIST"};
  // TString output [12] = {"HIST",  "HIST",   "HIST",  "HIST",   "HIST",   "HIST",  "HIST",  "HIST",  "HIST",  "HIST",  "HIST",  "HIST"}
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };
  Int_t   cutID1  [12] = {iCutP     ,iCutP     ,iCutP    ,iCutP      ,iCutP   ,iCutP   ,iCutP    ,iCutP    ,iCutP    ,iCutP    ,iCutP   ,iCutP   };
  Int_t   cutID2  [12] = {iCutK     ,iCutK     ,iCutK    ,iCutK      ,iCutK   ,iCutK   ,iCutK    ,iCutK    ,iCutK    ,iCutK    ,iCutK   ,iCutK   };
  Int_t   PDGCode [12] = {3124      ,3124      ,3124     ,3124       ,3124    ,3124    ,3124     ,-3124    ,3124     ,-3124    ,3124    ,-3124   };

  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("Lstar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(1.51953);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      //      out->AddAxis(imID, 300, 1.4, 1.7);
      //     out->AddAxis(imID, 300, 1.4, 2.2);
      out->AddAxis(imID, 800, 1.4, 2.2);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    

    // axis Y: transverse momentum of pair as default - else chosen value
    out->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt
    //out->AddAxis(OpAn,  50, 0.0, 25.0);
    //out->AddAxis(PtRat, 20, 0.0, 1.0);
    out->AddAxis(centID, 110, 0.0, 110.0); // adding multiplicity axis


    if(isPtDaughter) {
      out->AddAxis(fdpt, 100, 0.0, 10.0);
      out->AddAxis(sdpt, 100, 0.0, 10.0);     }
    
    if(isPDaughter) {
	out->AddAxis(fdp, 100, 0.0, 10.0);
	out->AddAxis(sdp, 100, 0.0, 10.0);    }
    
    // axis Z: centrality-multiplicity
    if(iscent) {
      if (!isPP)	out->AddAxis(centID, 100, 0.0, 100.0);
      else       	out->AddAxis(centID, 400, 0.0, 400.0);      }
    // axis W: pseudorapidity

    if(iseta)       out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    if(israpidity)  out->AddAxis(yID, 12, -0.6, 0.6);
    
  }
   
  if (isMC){
    cout<<"here###################################################################"<<endl;
    
    /* 
    AliRsnMiniOutput* outm=task->CreateOutput(Form("LStar_Resolution%s", suffix),"SPARSE","TRUE");
    outm->SetDaughter(0,AliRsnDaughter::kProton);
    outm->SetDaughter(1,AliRsnDaughter::kKaon);
    outm->SetMotherPDG(3124);
    outm->SetMotherMass(1.51953);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID, 800, 1.4, 2.2);
    outm->AddAxis(ptID, 100, 0., 10.);
    outm->AddAxis(IMdif, 1000, -0.5, 0.5);// mass diff (true -reco)

    AliRsnMiniOutput* outm_al=task->CreateOutput(Form("AntiLStar_Resolution%s", suffix),"SPARSE","TRUE");
    outm_al->SetDaughter(0,AliRsnDaughter::kProton);
    outm_al->SetDaughter(1,AliRsnDaughter::kKaon);
    outm_al->SetMotherPDG(-3124);
    outm_al->SetMotherMass(1.51953);
    outm_al->SetPairCuts(cutsPair);
    outm_al->AddAxis(imID, 800, 1.4, 2.2);
    outm_al->AddAxis(ptID, 100, 0., 10.);
    outm_al->AddAxis(IMdif, 1000, -0.5, 0.5);// mass diff (true -reco)
  

    AliRsnMiniOutput* outps=task->CreateOutput(Form("lambda_phaseSpace%s", suffix),"SPARSE","TRUE");
    outps->SetDaughter(0,AliRsnDaughter::kProton);
    outps->SetDaughter(1,AliRsnDaughter::kKaon);
    outps->SetMotherPDG(3124);
    outps->SetMotherMass(1.51953);
    outps->SetPairCuts(cutsPair);
    outps->AddAxis(ptID,200,0.,20.);
    outps->AddAxis(fdpt,100,0.,10.);
    outps->AddAxis(sdpt,100,0.,10.);
   

    AliRsnMiniOutput* outpsf=task->CreateOutput(Form("antilambda_phaseSpace%s", suffix),"SPARSE","TRUE");
    outpsf->SetDaughter(0,AliRsnDaughter::kProton);
    outpsf->SetDaughter(1,AliRsnDaughter::kKaon);
    outpsf->SetMotherPDG(-3124);
    outpsf->SetMotherMass(1.59153);
    outpsf->SetPairCuts(cutsPair);
    outpsf->AddAxis(ptID,200,0.,20.);
    outpsf->AddAxis(fdpt,100,0.,10.);
    outpsf->AddAxis(sdpt,100,0.,10.);
   

    AliRsnMiniOutput* outps_g=task->CreateOutput(Form("lambda_phaseSpace_gen%s", suffix),"SPARSE","MOTHER");
    outps_g->SetDaughter(0,AliRsnDaughter::kProton);
    outps_g->SetDaughter(1,AliRsnDaughter::kKaon);
    outps_g->SetMotherPDG(3124);
    outps_g->SetMotherMass(1.51953);
    outps_g->SetPairCuts(cutsPair);
    outps_g->AddAxis(ptID,200,0.,20.);
    outps_g->AddAxis(fdpt,100,0.,10.);
    outps_g->AddAxis(sdpt,100,0.,10.);
   

    AliRsnMiniOutput* outpsf_g=task->CreateOutput(Form("antilambda_phaseSpace_gen%s", suffix),"SPARSE","MOTHER");
    outpsf_g->SetDaughter(0,AliRsnDaughter::kProton);
    outpsf_g->SetDaughter(1,AliRsnDaughter::kKaon);
    outpsf_g->SetMotherPDG(-3124);
    outpsf_g->SetMotherMass(1.59153);
    outpsf_g->SetPairCuts(cutsPair);
    outpsf_g->AddAxis(ptID,200,0.,20.);
    outpsf_g->AddAxis(fdpt,100,0.,10.);
    outpsf_g->AddAxis(sdpt,100,0.,10.);
    
    */
  }  // isMC 
    
  return kTRUE;
}

/*
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID = 0, Int_t customFilterBit = 0)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.
  
  if ((!trkQualityCut)){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }
  
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
  trkQualityCut->SetPtRange(0.15, 100.0);
  trkQualityCut->SetEtaRange(-0.8, 0.8);
  
  Printf(Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));
  //  trkQualityCut->SetAODTestFilterBit(customFilterBit);
  //trkQualityCut->SetAODTestFilterBit(5); // hardcoded to 5 
  //trkQualityCut->SetCheckOnlyFilterBit(kFALSE);
  if(!customFilterBit){//ESD
    // if (customQualityCutsID==1) trkQualityCut->SetCheckOnlyFilterBit(kTRUE);
    // else if(customQualityCutsID==2) trkQualityCut->SetDCARPtFormula("0.006+0.0200/pt^1.1");
  else if(customQualityCutsID==2) trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.0200/pt^1.1");
  else if(customQualityCutsID==3) trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);
  else if(customQualityCutsID==4) trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.0);
  else if(customQualityCutsID==5) trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.3);
  else if(customQualityCutsID==6) trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3.5);
  else if(customQualityCutsID==7) trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);
  else if(customQualityCutsID==8) trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);
  else if(customQualityCutsID==9) trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  }else{ //AOD
    if (customQualityCutsID==1) trkQualityCut->SetCheckOnlyFilterBit(kTRUE);    
    else if(customQualityCutsID==2) trkQualityCut->SetDCARPtFormula("0.006+0.0200/pt^1.1");
    else if(customQualityCutsID==3) trkQualityCut->SetDCAZmax(0.2);
    else if(customQualityCutsID==4) trkQualityCut->SetDCAZmax(1.);
    else if(customQualityCutsID==5) trkQualityCut->SetTrackMaxChi2(2.3);
    else if(customQualityCutsID==6) trkQualityCut->SetTrackMaxChi2(3.5);
    else if(customQualityCutsID==7) trkQualityCut->SetMinNCrossedRowsTPC(100,kTRUE);
    else if(customQualityCutsID==8) trkQualityCut->SetMinNCrossedRowsTPC(80,kTRUE);
    else if(customQualityCutsID==9) trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.9,kTRUE);
  }
  
  
  trkQualityCut->Print();
  return kTRUE;
}
*/
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID = 0, Int_t customFilterBit = 0)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.

  if ((!trkQualityCut)){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }

  if(customQualityCutsID>=1 && customQualityCutsID<100 && customQualityCutsID!=12){
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    Printf(Form("::::: SetCustomQualityCut:: using standard 2011 track quality cuts"));

    if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.0150+0.0500/pt^1.1");}
    else if(customQualityCutsID==2){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.0200/pt^1.1");}
    else if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(5.);}
    else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}
    else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(5.);}
    else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.3);}
    else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}
    else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}
    else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
    else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
   
    else if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(49.);}
    else if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(4.);}
    else if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(49.);}
    else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(25.);}
    else if(customQualityCutsID==17){trkQualityCut->GetESDtrackCuts()->SetSPDminNClusters(0);}
    else if(customQualityCutsID==56){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);}
    else if(customQualityCutsID==58){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3.);}
    else if(customQualityCutsID==60){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);}
    else if(customQualityCutsID==64){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(25.);}

    trkQualityCut->Print();
    return kTRUE;
  }else if(customQualityCutsID==12 || (customQualityCutsID>=100 && customQualityCutsID<200)){
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

 
  trkQualityCut->SetPtRange(0.15, 20.0);
  trkQualityCut->SetEtaRange(-0.8, 0.8);
  
  Printf(Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));
  trkQualityCut->Print();
  return kTRUE;
}
