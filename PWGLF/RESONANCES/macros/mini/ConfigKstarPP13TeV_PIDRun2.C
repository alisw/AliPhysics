/***************************************************************************
             Modified by Kunal Garg - 2nd September 2016
             Anders Knospe - last modified on 9 July 2015

*** Configuration script for Kstar analysis of 2015 pp 13-TeV data ***
****************************************************************************/

Bool_t ConfigKstarPP13TeV_PIDRun2
(  
 AliRsnMiniAnalysisTask *task, 
 Bool_t                 isMC, 
 Bool_t                 isPP,
 const char             *suffix,
 AliRsnCutSet           *cutsPair,
 Int_t                  aodFilterBit=5,
 Int_t                  customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate=AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,    
 Float_t                nsigmaKa=2.0,
 Float_t                nsigmaPi = 2.0,
 Bool_t                 enableMonitor=kTRUE,
 Bool_t                 IsMcTrueOnly=kFALSE,
 TString                monitorOpt="",
 Bool_t                 useMixLS=0,
 Bool_t                 checkReflex=0,
 AliRsnMiniValue::EType yaxisVar=AliRsnMiniValue::kPt
)
{
  // manage suffix
  if(strlen(suffix)>0) suffix=Form("_%s",suffix);

  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetPi;
  AliRsnCutSetDaughterParticle * cutSetK;
  
  AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("myQualityCut");
    
    Bool_t check= SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit);
    
    cout<<"Checking which cut set is used   "<<check<<"     ENDING the printout"<<endl;
    
    
  if (SetCustomQualityCut(trkQualityCut, customQualityCutsID, aodFilterBit)) {
    //Set custom quality cuts for systematic checks
      
    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0);
      
      
    cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi), trkQualityCut, cutPiCandidate, AliPID::kPion, nsigmaPi);
    cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa), trkQualityCut, cutKaCandidate, AliPID::kKaon, nsigmaKa);
  } else {
    //use defult quality cuts (std 2010 or 2011)
      
    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, aodFilterBit, kTRUE);
    cutSetQ->SetUse2011StdQualityCuts(kTRUE);
      
    cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi), cutPiCandidate, AliPID::kPion, nsigmaPi, aodFilterBit,kTRUE);
    cutSetPi->SetUse2011StdQualityCuts(kTRUE);
      
    cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit,kTRUE);
    cutSetK->SetUse2011StdQualityCuts(kTRUE);
   }
  
  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutK = task->AddTrackCuts(cutSetK);
  
  
    if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput()), monitorOpt.Data();
  }  

  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes,kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP,kFALSE);
  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP,kFALSE);
  
  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

   Bool_t  use     [12] = {!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly ,isMC,isMC,isMC,isMC, useMixLS , useMixLS    };
  Bool_t  useIM   [12] = { 1        , 1         ,  1       ,  1       ,  1      ,  1      ,  1      ,  1      ,  0      , 0      , 1         , 1           };
  TString name    [12] = {"UnlikePM", "UnlikeMP","MixingPM","MixingMP", "LikePP", "LikeMM","TruesPM","TruesMP", "ResPM" ,"ResMP" , "MixingPP", "MixingMM"  };
  TString comp    [12] = {"PAIR"    , "PAIR"    , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  , "TRUE"  , "TRUE"  ,"TRUE"  , "MIX"     , "MIX"       };
  TString output  [12] = {"SPARSE"  , "SPARSE"  , "SPARSE" , "SPARSE" , "SPARSE", "SPARSE", "SPARSE","SPARSE" , "SPARSE","SPARSE", "SPARSE"  , "SPARSE"    };
  Int_t   pdgCode [12] = {313       , 313       , 313      , 313      , 313     , 313     , 313     , -313    ,  313    , -313   ,  313      , 313         };
  Char_t  charge1 [12] = {'+'       , '-'       , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'    , '+'     , '-'    , '+'       , '-'         };
  Char_t  charge2 [12] = {'-'       , '+'       , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'    , '-'     , '+'    , '+'       , '-'         };
  Int_t   cutID1  [12] = { iCutK    ,  iCutK    ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK ,  iCutK  , iCutK  , iCutK     , iCutK       };
  Int_t   cutID2  [12] = { iCutPi   ,  iCutPi   ,  iCutPi  ,  iCutPi  ,  iCutPi ,  iCutPi ,  iCutPi ,   iCutPi,  iCutPi , iCutPi , iCutPi    , iCutPi      };

  for(Int_t i=0;i<12;i++){
    if(!use[i]) continue;
    AliRsnMiniOutput* out=task->CreateOutput(Form("Kstar_%s%s",name[i].Data(),suffix),output[i].Data(),comp[i].Data());
    out->SetCutID(0,iCutK);
    out->SetCutID(1,iCutPi);
    out->SetDaughter(0,AliRsnDaughter::kKaon);
    out->SetDaughter(1,AliRsnDaughter::kPion);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(pdgCode[i]);
    out->SetMotherMass(0.89594);
    out->SetPairCuts(cutsPair);

    //axis X: invmass (or resolution)
    if(useIM[i]) out->AddAxis(imID,90,0.6,1.5);
    else out->AddAxis(resID,200,-0.02,0.02);

    //axis Y: transverse momentum of pair as default - else chosen value
    if(yaxisVar==AliRsnMiniValue::kFirstDaughterPt) out->AddAxis(fdpt,100,0.,10.);
    else if(yaxisVar==AliRsnMiniValue::kSecondDaughterPt) out->AddAxis(sdpt,100,0.,10.);
    else if(yaxisVar==AliRsnMiniValue::kFirstDaughterP) out->AddAxis(fdp,100,0.,10.);
    else if(yaxisVar==AliRsnMiniValue::kSecondDaughterP)  out->AddAxis(sdp,100,0.,10.);
    else out->AddAxis(ptID,200,0.,20.);//default use mother pt

    // axis Z: centrality-multiplicity
    if(!isPP) out->AddAxis(centID,100,0.,100.);
    else out->AddAxis(centID,400,0.5,400.5);
    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    // out->AddAxis(yID, 10, -0.5, 0.5);
  }   

  if(isMC){   
    //get mothers for Kstar PDG = 313
    AliRsnMiniOutput* outm=task->CreateOutput(Form("Kstar_Mother%s", suffix),"SPARSE","MOTHER");
    outm->SetDaughter(0,AliRsnDaughter::kKaon);
    outm->SetDaughter(1,AliRsnDaughter::kPion);
    outm->SetMotherPDG(313);
    outm->SetMotherMass(0.89594);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID,215,0.985,1.2);
    outm->AddAxis(ptID,200,0.,20.);
    if(!isPP) outm->AddAxis(centID,100,0.,100.);
    else outm->AddAxis(centID,400,0.5,400.5);

    //get phase space of the decay from mothers
    AliRsnMiniOutput* outps=task->CreateOutput(Form("Kstar_phaseSpace%s", suffix),"HIST","TRUE");
    outps->SetDaughter(0,AliRsnDaughter::kKaon);
    outps->SetDaughter(1,AliRsnDaughter::kPion);
    outps->SetCutID(0,iCutK);
    outps->SetCutID(1,iCutPi);
    outps->SetMotherPDG(313);
    outps->SetMotherMass(0.89594);
    outps->SetPairCuts(cutsPair);
    outps->AddAxis(fdpt,50,0.,5.);
    outps->AddAxis(sdpt,50,0.,5.);
    outps->AddAxis(ptID,200,0.,20.);

    //get reflections
    if(checkReflex){
      AliRsnMiniOutput* outreflex=task->CreateOutput(Form("Kstar_reflex%s", suffix),"SPARSE","TRUE");
      outreflex->SetDaughter(0,AliRsnDaughter::kKaon);
      outreflex->SetDaughter(1,AliRsnDaughter::kPion);
      outreflex->SetCutID(0,iCutK);
      outreflex->SetCutID(1,iCutM);
      outreflex->SetMotherPDG(313);
      outreflex->SetMotherMass(0.89594);
      outreflex->SetPairCuts(cutsPair);
      outreflex->AddAxis(imID,215,0.985,1.2);
      outreflex->AddAxis(ptID,200,0.,20.);
      if(!isPP) outreflex->AddAxis(centID,100,0.,100.);
      else outreflex->AddAxis(centID,400,0.5,400.5);
    }//end reflections
  }//end MC

  return kTRUE;
}

//-------------------------------------------------------  
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

  if(customQualityCutsID>=1 && customQualityCutsID<100 && customQualityCutsID!=2){
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    Printf(Form("::::: SetCustomQualityCut:: using standard 2011 track quality cuts"));

    if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.0150+0.0500/pt^1.1");}
    else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.0075+0.0250/pt^1.1");}
    else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(3.);}
    else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}
    else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(7.);}
    else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.5);}
    else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}
    else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}
    else if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
    else if(customQualityCutsID==12){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.99);}
    else if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(50);}
    else if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(4);}
    else if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(40);}
    else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(32);}

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

  //for pA 2013
  //trkQualityCut->SetDefaults2011();//with filter bit=10
  //reset filter bit to very loose cuts 
  trkQualityCut->SetAODTestFilterBit(customFilterBit); 
  //apply all other cuts "by hand"
  trkQualityCut->SetCheckOnlyFilterBit(kFALSE);
  trkQualityCut->SetMinNCrossedRowsTPC(70, kTRUE);
  trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.8, kTRUE);
  trkQualityCut->SetMaxChi2TPCConstrainedGlobal(36);//used for ESD only - for AOD does not correspond to any cut
  trkQualityCut->SetTPCmaxChi2(4.0); //already in filter bit 0
  trkQualityCut->SetRejectKinkDaughters(kTRUE); //already in filter bit 0
  trkQualityCut->SetSPDminNClusters(AliESDtrackCuts::kAny);
  trkQualityCut->SetITSmaxChi2(36);
  trkQualityCut->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);//already in defaults 2011
  trkQualityCut->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);//already in defaults 2011
  trkQualityCut->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);//already in defaults 2011

  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kFilterBitCustom) {
    trkQualityCut->SetCheckOnlyFilterBit(kTRUE);
  } 
  
  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdLooserDCAXY){
    trkQualityCut->SetDCARmax(2.4);
  } else {
    trkQualityCut->SetDCARPtFormula("0.0105+0.0350/pt^1.1");
  }
  
  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdLooserDCAZ){
    trkQualityCut->SetDCAZmax(3.2);
  } else {
    trkQualityCut->SetDCAZmax(2.0); 
  }
  
  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdCrossedRows60){
    trkQualityCut->SetMinNCrossedRowsTPC(60, kTRUE);
  }
  
  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdCrossedRows80){
    trkQualityCut->SetMinNCrossedRowsTPC(80, kTRUE);
  }
  
  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdRowsToCls075){
    trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.75, kTRUE);
  }
  
  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdRowsToCls085){
    trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.85, kTRUE);
  }
  
  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdCls70){
    trkQualityCut->SetAODTestFilterBit(10);
    trkQualityCut->SetTPCminNClusters(70);
  }
  
  if (customQualityCutsID==AliRsnCutSetDaughterParticle::kStdChi2TPCCls35){
    trkQualityCut->SetTPCmaxChi2(3.5);
  }
  
  trkQualityCut->SetPtRange(0.15, 20.0);
  trkQualityCut->SetEtaRange(-0.8, 0.8);
  
  Printf(Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));
  trkQualityCut->Print();
  return kTRUE;
}
