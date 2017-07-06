/***************************************************************************
              Anders Knospe - last modified on 26 March 2016

*** Configuration script for phi analysis of 2015 pp 13-TeV data ***
****************************************************************************/

Bool_t ConfigPhiPP13TeV_PID
(  
 AliRsnMiniAnalysisTask *task, 
 Bool_t                 isMC, 
 Bool_t                 isPP,
 const char             *suffix,
 AliRsnCutSet           *cutsPair,
 Int_t                  aodFilterBit=5,
 Int_t                  customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kQualityStd2010,
 Float_t                nsigmaKa=2.,
 Bool_t                 enableMonitor=kTRUE,
 Bool_t                 IsMcTrueOnly=kFALSE,
 TString                monitorOpt="",
 Bool_t                 useMixLS=0,
 Bool_t                 checkReflex=0,
 AliRsnMiniValue::EType yaxisVar=AliRsnMiniValue::kPt,
 TString                polarizationOpt="" /* J - Jackson,T - Transversity */
)
{
  // manage suffix
  if(strlen(suffix)>0) suffix=Form("_%s",suffix);

  // set daughter cuts
  AliRsnCutSetDaughterParticle* cutSetQ;
  AliRsnCutSetDaughterParticle* cutSetK;

  Int_t MultBins=aodFilterBit/100;
  aodFilterBit=aodFilterBit%100;

  Float_t nsigmaKaTPC=fmod(nsigmaKa,1000.);
  Float_t nsigmaKaTOF=(nsigmaKa-fmod(nsigmaKa,1000.))/1000.;
  if(nsigmaKaTOF<1.e-10) nsigmaKaTOF=-1.;

  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit)){
    //Set custom quality cuts for systematic checks
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit),trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kKaon,-1.);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa),trkQualityCut,cutKaCandidate,AliPID::kKaon,nsigmaKaTPC,nsigmaKaTOF);
  }else{
    //use default quality cuts std 2010 with crossed rows TPC
    Bool_t useCrossedRows = 1;
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit),AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kKaon,-1.,aodFilterBit,useCrossedRows);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate,nsigmaKa),cutKaCandidate,AliPID::kKaon,nsigmaKa,aodFilterBit,useCrossedRows);
  }

  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutK=task->AddTrackCuts(cutSetK);

  if(enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput()), monitorOpt.Data();
  }  

  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,kFALSE);
  /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
  /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP,kFALSE);
  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP,kFALSE);
  /* cos(theta) J     */ Int_t ctjID  = task->CreateValue(AliRsnMiniValue::kCosThetaJackson,kFALSE);
  /* cos(theta) J (MC)*/ Int_t ctjmID  = task->CreateValue(AliRsnMiniValue::kCosThetaJackson,kTRUE);
  /* cos(theta) T     */ Int_t cttID  = task->CreateValue(AliRsnMiniValue::kCosThetaTransversity,kFALSE);
  /* cos(theta) T (MC)*/ Int_t cttmID  = task->CreateValue(AliRsnMiniValue::kCosThetaTransversity,kTRUE);

  Double_t multbins[200];
  int j,nmult=0;
  for(j=0;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
  for(j=1;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
  for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
  //for(j=1;j<=18;j++){multbins[nmult]=10+5*j; nmult++;}
  
  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use    [11]={!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly, isMC,isMC,isMC,isMC,isMC, useMixLS,useMixLS};
  Int_t   useIM  [11]={ 1      ,  1     , 1      ,  1     ,  1     ,  1        ,  2      , 2           ,0       , 1        , 1        };
  TString name   [11]={"Unlike","Mixing","LikePP","LikeMM","Trues" ,"TruesFine","TruesMM","TruesFineMM","Res"   ,"MixingPP","MixingMM"};
  TString comp   [11]={"PAIR"  , "MIX"  ,"PAIR"  ,"PAIR"  , "TRUE" , "TRUE"    ,"TRUE"   ,"TRUE"       ,"TRUE"  ,"MIX"     ,"MIX"     };
  TString output [11]={"HIST","HIST","HIST","HIST","HIST","HIST"   ,"HIST" ,"HIST"     ,"HIST","HIST"  ,"HIST"  };
  Int_t   pdgCode[11]={333     , 333    ,333     ,333     , 333    , 333       ,333      ,333          ,333     , 333      ,333       };
  Char_t  charge1[11]={'+'     , '+'    ,'+'     ,'-'     , '+'    , '+'       ,'+'      , '+'         ,'+'     ,'+'       ,'-'       };
  Char_t  charge2[11]={'-'     , '-'    ,'+'     ,'-'     , '-'    , '-'       ,'-'      , '-'         ,'-'     ,'+'       ,'-'       };

  for(Int_t i=0;i<9;i++){
    if(!use[i]) continue;
    AliRsnMiniOutput* out=task->CreateOutput(Form("phi_%s%s",name[i].Data(),suffix),output[i].Data(),comp[i].Data());
    out->SetCutID(0,iCutK);
    out->SetCutID(1,iCutK);
    out->SetDaughter(0,AliRsnDaughter::kKaon);
    out->SetDaughter(1,AliRsnDaughter::kKaon);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(pdgCode[i]);
    out->SetMotherMass(1.019461);
    out->SetPairCuts(cutsPair);

    //axis X: invmass (or resolution)
    if(useIM[i]==1) out->AddAxis(imID,215,0.985,1.2);
    else if(useIM[i]==2) out->AddAxis(mmID,75,0.985,1.06);
    else out->AddAxis(diffID,200,-0.02,0.02);

    //axis Y: transverse momentum of pair as default - else chosen value
    if(yaxisVar==AliRsnMiniValue::kFirstDaughterPt) out->AddAxis(fdpt,100,0.,10.);
    else if(yaxisVar==AliRsnMiniValue::kSecondDaughterPt) out->AddAxis(sdpt,100,0.,10.);
    else if(yaxisVar==AliRsnMiniValue::kFirstDaughterP) out->AddAxis(fdp,100,0.,10.);
    else if(yaxisVar==AliRsnMiniValue::kSecondDaughterP)  out->AddAxis(sdp,100,0.,10.);
    else if(isMC && (i==5 || i==7)) out->AddAxis(ptID,300,0.,3.);//fine binning for efficiency weighting
    else out->AddAxis(ptID,200,0.,20.);//default use mother pt

    // axis Z: centrality-multiplicity
    //if(!isPP || MultBins) out->AddAxis(centID,nmult,multbins);//out->AddAxis(centID,100,0.,100.);
    //else out->AddAxis(centID,nmult,multbins);//out->AddAxis(centID,161,-0.5,160.5);
    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    // out->AddAxis(yID, 10, -0.5, 0.5);

    if (polarizationOpt.Contains("J")) out->AddAxis(ctjID,21,-1.,1);
    if (polarizationOpt.Contains("T")) out->AddAxis(cttID,21,-1.,1);
  }

  if(isMC){   
    //get mothers for phi PDG = 333
    AliRsnMiniOutput* outm=task->CreateOutput(Form("phi_Mother%s", suffix),"SPARSE","MOTHER");
    outm->SetDaughter(0,AliRsnDaughter::kKaon);
    outm->SetDaughter(1,AliRsnDaughter::kKaon);
    outm->SetMotherPDG(333);
    outm->SetMotherMass(1.019461);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID,215,0.985,1.2);
    outm->AddAxis(ptID,200,0.,20.);
    if(!isPP || MultBins) outm->AddAxis(centID,100,0.,100.);
    else outm->AddAxis(centID,161,-0.5,160.5);
    if (polarizationOpt.Contains("J")) outm->AddAxis(ctjmID,21,-1.,1.);
    if (polarizationOpt.Contains("T")) outm->AddAxis(cttmID,21,-1.,1.);

    AliRsnMiniOutput* outmf=task->CreateOutput(Form("phi_MotherFine%s", suffix),"SPARSE","MOTHER");
    outmf->SetDaughter(0,AliRsnDaughter::kKaon);
    outmf->SetDaughter(1,AliRsnDaughter::kKaon);
    outmf->SetMotherPDG(333);
    outmf->SetMotherMass(1.019461);
    outmf->SetPairCuts(cutsPair);
    outmf->AddAxis(imID,215,0.985,1.2);
    outmf->AddAxis(ptID,300,0.,3.);//fine binning for efficiency weighting
    if(!isPP || MultBins) outmf->AddAxis(centID,100,0.,100.);
    else outmf->AddAxis(centID,161,-0.5,160.5);
    if (polarizationOpt.Contains("J")) outmf->AddAxis(ctjmID,21,-1.,1.);
    if (polarizationOpt.Contains("T")) outmf->AddAxis(cttmID,21,-1.,1.);

    //get phase space of the decay from mothers
    AliRsnMiniOutput* outps=task->CreateOutput(Form("phi_phaseSpace%s", suffix),"HIST","TRUE");
    outps->SetDaughter(0,AliRsnDaughter::kKaon);
    outps->SetDaughter(1,AliRsnDaughter::kKaon);
    outps->SetCutID(0,iCutK);
    outps->SetCutID(1,iCutK);
    outps->SetMotherPDG(333);
    outps->SetMotherMass(1.019461);
    outps->SetPairCuts(cutsPair);
    outps->AddAxis(fdpt,100,0.,10.);
    outps->AddAxis(sdpt,100,0.,10.);
    outps->AddAxis(ptID,200,0.,20.);

    AliRsnMiniOutput* outpsf=task->CreateOutput(Form("phi_phaseSpaceFine%s", suffix),"HIST","TRUE");
    outpsf->SetDaughter(0,AliRsnDaughter::kKaon);
    outpsf->SetDaughter(1,AliRsnDaughter::kKaon);
    outpsf->SetCutID(0,iCutK);
    outpsf->SetCutID(1,iCutK);
    outpsf->SetMotherPDG(333);
    outpsf->SetMotherMass(1.019461);
    outpsf->SetPairCuts(cutsPair);
    outpsf->AddAxis(fdpt,30,0.,3.);
    outpsf->AddAxis(sdpt,30,0.,3.);
    outpsf->AddAxis(ptID,300,0.,3.);

    //get reflections
    if(checkReflex){
      AliRsnMiniOutput* outreflex=task->CreateOutput(Form("phi_reflex%s", suffix),"SPARSE","TRUE");
      outreflex->SetDaughter(0,AliRsnDaughter::kKaon);
      outreflex->SetDaughter(1,AliRsnDaughter::kKaon);
      outreflex->SetCutID(0,iCutK);
      outreflex->SetCutID(1,iCutK);
      outreflex->SetMotherPDG(333);
      outreflex->SetMotherMass(1.019461);
      outreflex->SetPairCuts(cutsPair);
      outreflex->AddAxis(imID,215,0.985,1.2);
      outreflex->AddAxis(ptID,200,0.,20.);
      if(!isPP) outreflex->AddAxis(centID,100,0.,100.);
      else outreflex->AddAxis(centID,400,0.5,400.5);
      if (polarizationOpt.Contains("J")) outreflex->AddAxis(ctjID,21,-1.,1.);
      if (polarizationOpt.Contains("T")) outreflex->AddAxis(cttID,21,-1.,1.);
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

    if(!customFilterBit){//ESD
      if(customQualityCutsID==3){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.0150+0.0500/pt^1.1");}
      else if(customQualityCutsID==4){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexXYPtDep("0.006+0.0200/pt^1.1");}
      else if(customQualityCutsID==5){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(5.);}
      else if(customQualityCutsID==6){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(0.2);}
      else if(customQualityCutsID==7){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(5.);}
      else if(customQualityCutsID==8){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(2.3);}
      else if(customQualityCutsID==9){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(60);}
      else if(customQualityCutsID==10){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(100);}
      else if(customQualityCutsID==11){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
      else if(customQualityCutsID==12){trkQualityCut->GetESDtrackCuts()->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}
      else if(customQualityCutsID==13){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(49.);}
      else if(customQualityCutsID==14){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(4.);}
      else if(customQualityCutsID==15){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(49.);}
      else if(customQualityCutsID==16){trkQualityCut->GetESDtrackCuts()->SetMaxChi2TPCConstrainedGlobal(25.);}
      else if(customQualityCutsID==17){trkQualityCut->GetESDtrackCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}
      else if(customQualityCutsID==56){trkQualityCut->GetESDtrackCuts()->SetMaxDCAToVertexZ(1.);}
      else if(customQualityCutsID==58){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterTPC(3.);}
      else if(customQualityCutsID==60){trkQualityCut->GetESDtrackCuts()->SetMinNCrossedRowsTPC(80);}
      else if(customQualityCutsID==64){trkQualityCut->GetESDtrackCuts()->SetMaxChi2PerClusterITS(25.);}
    }else{//AOD
      trkQualityCut->SetCheckOnlyFilterBit(kFALSE);
      if(customQualityCutsID==4){trkQualityCut->SetDCARPtFormula("0.006+0.0200/pt^1.1");}
      else if(customQualityCutsID==6){trkQualityCut->SetDCAZmax(0.2);}
      else if(customQualityCutsID==8){trkQualityCut->SetTrackMaxChi2(2.3);}
      else if(customQualityCutsID==10){trkQualityCut->SetMinNCrossedRowsTPC(100,kTRUE);}
      else if(customQualityCutsID==12){trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.9,kTRUE);}
      else if(customQualityCutsID==56){trkQualityCut->SetDCAZmax(1.);}
      else if(customQualityCutsID==58){trkQualityCut->SetTrackMaxChi2(3.5);}
      else if(customQualityCutsID==60){trkQualityCut->SetMinNCrossedRowsTPC(80,kTRUE);}
    }

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
