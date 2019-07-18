/***************************************************************************
              Anders Knospe: anders.knospe@cern.ch
  Macro to configure the resonance package for searches for rare resonances.

****************************************************************************/

AliRsnMiniAnalysisTask* AddTaskRhof0f2(
  TString lname,
  Bool_t isMC,
  Int_t system,
  Int_t EventCuts=0,
  Int_t isAOD=0,
  Int_t Strcut=2011,
  Int_t customQualityCutsID=0,
  Int_t TrackCuts1=0,
  Int_t TrackCuts2=0
){
  // ----- INITIALIZATION -----

  // retrieve analysis manager
  AliAnalysisManager* mgr=AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskRho", "No analysis manager to connect to.");
    return NULL;
  }

  // create the task and configure
  AliRsnMiniAnalysisTask* task=new AliRsnMiniAnalysisTask(lname,isMC);

  // trigger
  int trigger=EventCuts%10;
  if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
  else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMultV0);

  // multiplicity
  bool isPP=false;
  if(!system) isPP=true;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  if(isPP){
    if(MultBins==1) task->UseMultiplicity("AliMultSelection_V0M");
    else if(MultBins==2) task->UseMultiplicity("AliMultSelection_RefMult08");
    else task->UseMultiplicity("QUALITY");
  }else if(system==1) task->UseMultiplicity("AliMultSelection_V0A");
  else if(system==2) task->UseMultiplicity("AliMultSelection_V0M");
  else task->UseCentrality("V0M");

  // set event mixing options
  int nmix=5;
  if((EventCuts%10000)/1000==1) nmix=0;
  float maxDiffVzMix=1;
  float maxDiffMultMix=5;
  task->UseContinuousMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  ::Info("AddTaskRho", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));

  // vertex cuts
  float vtxZcut=10;
  Bool_t rejectPileUp=kTRUE;
  AliRsnCutPrimaryVertex* cutVertex=0;
  if(!MultBins || fabs(vtxZcut-10.)>1.e-10){
    cutVertex=new AliRsnCutPrimaryVertex("cutVertex",vtxZcut,0,kFALSE);
    if(!MultBins){
      cutVertex->SetCheckZResolutionSPD();
      cutVertex->SetCheckDispersionSPD();
      cutVertex->SetCheckZDifferenceSPDTrack();
    }
    if(0) cutVertex->SetCheckGeneratedVertexZ();
  }

  // other event selection cuts
  AliRsnCutEventUtils* cutEventUtils=0;
  if(1){
    cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
    if(!MultBins){
      cutEventUtils->SetCheckIncompleteDAQ();
      cutEventUtils->SetCheckSPDClusterVsTrackletBG();
    }else{
      cutEventUtils->SetRemovePileUppA2013(kFALSE);
      cutEventUtils->SetCheckAcceptedMultSelection();
    }
  }

  // set the check for pileup
  if(isPP && (!isMC) && cutVertex){
    cutVertex->SetCheckPileUp(rejectPileUp);
    ::Info("AddTaskRho", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
  }

  // define and fill cut set for event cuts
  AliRsnCutSet* eventCuts=0;
  if(cutEventUtils || cutVertex){
    eventCuts=new AliRsnCutSet("eventCuts",AliRsnTarget::kEvent);

    if(cutEventUtils && cutVertex){
      eventCuts->AddCut(cutEventUtils);
      eventCuts->AddCut(cutVertex);
      eventCuts->SetCutScheme(Form("%s&%s",cutEventUtils->GetName(),cutVertex->GetName()));
    }else if(cutEventUtils && !cutVertex){
      eventCuts->AddCut(cutEventUtils);
      eventCuts->SetCutScheme(Form("%s",cutEventUtils->GetName()));
    }else if(!cutEventUtils && cutVertex){
      eventCuts->AddCut(cutVertex);
      eventCuts->SetCutScheme(Form("%s",cutVertex->GetName()));
    }

    task->SetEventCuts(eventCuts);
  }

  // ----- EVENT-ONLY COMPUTATIONS -----
    
  Double_t multbins[1000];
  int j,nmult=0;
  if(!MultBins){
    for(j=0;j<=401;j++){multbins[nmult]=j-0.5; nmult++;}
  }else if(!trigger){
    for(j=0;j<=100;j++){multbins[nmult]=j; nmult++;}
  }else{
    for(j=0;j<10;j++){multbins[nmult]=0.0001*j; nmult++;}
    for(j=1;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
    for(j=1;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
    for(j=1;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
    for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
  }
  nmult--;

  //vertex
  Int_t vtxID=task->CreateValue(AliRsnMiniValue::kVz,kFALSE);
  AliRsnMiniOutput* outVtx=task->CreateOutput("eventVtx","HIST","EVENT");
  outVtx->AddAxis(vtxID,240,-12.0,12.0);

  //multiplicity or centrality
  Int_t multID=task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
  AliRsnMiniOutput* outMult=task->CreateOutput("eventMult","HIST","EVENT");
  outMult->AddAxis(multID,nmult+1,multbins);

  TH1F* hEventsVsMulti=new TH1F("hAEventsVsMulti","",nmult,multbins);
  task->SetEventQAHist("EventsVsMulti",hEventsVsMulti);//custom binning for fHAEventsVsMulti
    
  double ybins[1000];
  for(j=0;j<=240;j++) ybins[j]=-12+0.1*j;

  TH2F* hvz=new TH2F("hVzVsCent","",nmult,multbins, 240,ybins);
  task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

  for(j=0;j<=401;j++) ybins[j]=j-0.5;

  TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
  hmc->GetYaxis()->SetTitle("QUALITY");
  task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

  // ----- CONFIGURE -----

  cerr<<"configuring"<<endl;
  Config_pipi(task,lname,isMC,system,EventCuts,isAOD,Strcut,customQualityCutsID,TrackCuts1,TrackCuts2);
  cerr<<"done configuring"<<endl;
  
  // ----- CONTAINERS -----

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("AddTaskRho - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()),
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,
							  outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
   
  return task;
}


//=============================

Bool_t Config_pipi(
  AliRsnMiniAnalysisTask *task,
  TString     lname="pipi",
  Bool_t      isMC=kFALSE,
  Int_t       system=1,
  Int_t       EventCuts=0,
  Int_t       isAOD=0,
  Int_t       Strcut =2011,
  Int_t       customQualityCutsID=0,
  Int_t       TrackCutsPi=0,
  Int_t       TrackCuts2=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  




  // retrieve mass from PDG database
  Int_t pdg=TrackCuts2;
  TDatabasePDG* db=TDatabasePDG::Instance();
  TParticlePDG* part=db->GetParticle(pdg);
  Double_t mass=part->Mass();


  // set daughter cuts
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only

  AliRsnCutSetDaughterParticle* cutSetPi=0;
  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  if(!trkQualityCut) return kFALSE;
  
  if(SetCustomQualityCut(trkQualityCut,isAOD,customQualityCutsID,Strcut)){
    cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",CutTypePi,nsigmaPiTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
   
  }
  else{
    printf("Daughter Track cuts has been selected =================\n");
    return kFALSE;
  }

  
  // Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutPi=task->AddTrackCuts(cutSetPi);

  // monitoring
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);
  AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());

  
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
  /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges

    
  Bool_t  use    [15]={ 1      ,  1     , 1      ,  1      , isMC   ,   isMC   ,   isMC   ,  1       ,  1       , isMC    , isMC    , isMC   , isMC   , isMC   , isMC   };
  Int_t   useIM  [15]={ 1      ,  1     , 1      ,  1      ,  0     ,    0     ,    0     ,  1       ,  1       , 1       ,  1      , 1      ,  1     ,  1     ,  1     };
  TString name   [15]={"Unlike","Mixing","LikePP","LikeMM" ,"resrh0","resf0"   ,"resf2"   ,"MixingPP","MixingMM","f0T"   ,"f0M"     ,"RhoT"  ,"RhoM"  ,"f2T"   ,"f2M"  };
  TString comp   [15]={"PAIR"  , "MIX"  ,"PAIR"  ,"PAIR"   ,"TRUE"  , "TRUE"   , "TRUE"   ,"MIX"     ,"MIX"     ,"TRUE"   , "MOTHER","TRUE"  ,"MOTHER","TRUE"  ,"MOTHER"};
  TString output [15]={"HIST"  ,"HIST"  ,"HIST"  ,"HIST"   ,"HIST"  , "HIST"   , "HIST"   ,"HIST"    ,"HIST"    ,"HIST"   ,"HIST"   ,"HIST"  ,"HIST"  ,"HIST"  ,"HIST"  };
  // TString output [15]={"SPARSE","SPARSE","SPARSE","SPARSE" ,"SPARSE","SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE","SPARSE","SPARSE"};
  Char_t  charge1[15]={'+'     , '+'    ,'+'     ,'-'      , '+'    ,'+'       ,'+'       ,'+'       ,'-'       ,'+'      , '+'     ,'+'     ,'+'     , '+'    , '+'    };
  Char_t  charge2[15]={'-'     , '-'    ,'+'     ,'-'      , '-'    , '-'      , '-'      , '+'      ,'-'       ,'-'      , '-'     ,'-'     ,'-'     , '-'    , '-'    };
    Int_t PDGCode[15]={ 9010221,9010221 ,9010221 ,9010221  ,113     ,  9010221 ,  225     ,  9010221 ,9010221   ,9010221  , 9010221 , 113    ,  113   , 225    , 225    };
    Float_t Mass[15]={ 0.990   , 0.990  , 0.990  , 0.990   ,0.77526 , 0.990    ,  1.2755  , 0.990    , 0.990    ,0.990    ,0.990    , 0.77526, 0.77526, 1.2755 ,1.2755  };
  for(Int_t i=0;i<15;i++){
    if(!use[i]) continue;
    AliRsnMiniOutput *out=task->CreateOutput(Form("pipi_%s%s",name[i].Data(),suffix),output[i].Data(),comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kPion);
    out->SetDaughter(1,AliRsnDaughter::kPion);
    out->SetCutID(0,iCutPi);
    out->SetCutID(1,iCutPi);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(Mass[i]);
    ////out->SetPairCuts(cutsPair);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,173,0.27,2.);
    else out->AddAxis(diffID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }
 
  if(isMC) {

    AliRsnMiniOutput *outrho = task->CreateOutput(Form("rhoT_MM%s", suffix), "HIST", "TRUE");
    outrho->SetDaughter(0, AliRsnDaughter::kPion);
    outrho->SetDaughter(1, AliRsnDaughter::kPion);
    outrho->SetCutID(0,iCutPi);
    outrho->SetCutID(1,iCutPi);
    outrho->SetCharge(0,'+');
    outrho->SetCharge(1,'-');
    outrho->SetMotherPDG(113);
    outrho->SetMotherMass(0.77526);
    outrho->AddAxis(mmID, 173, 0.27, 2.);
    outrho->AddAxis(ptID, 200, 0.0, 20.0);
    outrho->AddAxis(centID,nmult,multbins);


    AliRsnMiniOutput *outf0 = task->CreateOutput(Form("f0T_MM%s", suffix), "HIST", "TRUE");
    outf0->SetDaughter(0, AliRsnDaughter::kPion);
    outf0->SetDaughter(1, AliRsnDaughter::kPion);
    outf0->SetCutID(0,iCutPi);
    outf0->SetCutID(1,iCutPi);
    outf0->SetCharge(0,'+');
    outf0->SetCharge(1,'-');
    outf0->SetMotherPDG(9010221);
    outf0->SetMotherMass(0.990);
    outf0->AddAxis(mmID, 173, 0.27, 2.);
    outf0->AddAxis(ptID, 200, 0.0, 20.0);
    outf0->AddAxis(centID,nmult,multbins);

    AliRsnMiniOutput *outf2 = task->CreateOutput(Form("f2T_MM%s", suffix), "HIST", "TRUE");
    outf2->SetDaughter(0, AliRsnDaughter::kPion);
    outf2->SetDaughter(1, AliRsnDaughter::kPion);
    outf2->SetCutID(0,iCutPi);
    outf2->SetCutID(1,iCutPi);
    outf2->SetCharge(0,'+');
    outf2->SetCharge(1,'-');
    outf2->SetMotherPDG(225);
    outf2->SetMotherMass(1.2755);
    outf2->AddAxis(mmID, 173, 0.27, 2.);
    outf2->AddAxis(ptID, 200, 0.0, 20.0);
    outf2->AddAxis(centID,nmult,multbins);

  

    //get phase space of the decay from mothers
    AliRsnMiniOutput *outps = task->CreateOutput(Form("rho_PS%s", suffix), "HIST", "TRUE");
    outps->SetDaughter(0, AliRsnDaughter::kPion);
    outps->SetCutID(0, iCutPi);
    outps->SetCharge(0,'+');
    outps->SetDaughter(1, AliRsnDaughter::kPion);
    outps->SetCutID(1, iCutPi);
    outps->SetCharge(0,'-');
    outps->SetMotherPDG(113);
    outps->SetMotherMass(0.77546);
    outps->SetPairCuts(cutsPair);
    outps->AddAxis(fdpt, 100, 0.0, 10.0);
    outps->AddAxis(sdpt, 100, 0.0, 10.0);
    outps->AddAxis(ptID, 200, 0.0, 20.0);

    AliRsnMiniOutput *outp = task->CreateOutput(Form("f0_PS%s", suffix), "HIST", "TRUE");
    outp->SetDaughter(0, AliRsnDaughter::kPion);
    outp->SetCutID(0, iCutPi);
    outp->SetCharge(0,'+');
    outp->SetDaughter(1, AliRsnDaughter::kPion);
    outp->SetCutID(1, iCutPi);
    outp->SetCharge(1,'-');
    outp->SetMotherPDG(9010221);
    outp->SetMotherMass(0.990);
    outp->SetPairCuts(cutsPair);
    outp->AddAxis(fdpt, 100, 0.0, 10.0);
    outp->AddAxis(sdpt, 100, 0.0, 10.0);
    outp->AddAxis(ptID, 200, 0.0, 20.0);

    AliRsnMiniOutput *outpp = task->CreateOutput(Form("f2_PS%s", suffix), "HIST", "TRUE");
    outpp->SetDaughter(0, AliRsnDaughter::kPion);
    outpp->SetCutID(0, iCutPi);
    outpp->SetCharge(0,'+');
    outpp->SetDaughter(1, AliRsnDaughter::kPion);
    outpp->SetCutID(1, iCutPi);
    outpp->SetCharge(1,'-');
    outpp->SetMotherPDG(225);
    outpp->SetMotherMass(1.2755);
    outpp->SetPairCuts(cutsPair);
    outpp->AddAxis(fdpt, 100, 0.0, 10.0);
    outpp->AddAxis(sdpt, 100, 0.0, 10.0);
    outpp->AddAxis(ptID, 200, 0.0, 20.0);
  }//end MC
  
  return kTRUE;
}



//-------------------------------------------------------  

Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t isAOD=0, Int_t customQualityCutsID = 0, Int_t trCut = 2011)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.

    if ((!trkQualityCut)){
        Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
        return kFALSE;
    }
    
    if(trCut == 2011){
        trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
        Printf(Form("::::: SetCustomQualityCut:: using standard 2011 track quality cuts"));
    }else if(trCut == 2015){
        trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
        trkQualityCut->GetESDtrackCuts()->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
        Printf(Form("::::: SetCustomQualityCut:: using standard 2015 track quality cuts"));
    }
    
    if(!isAOD){//ESD
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
        else if(customQualityCutsID){
            Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
            return kFALSE;
        }
    }else{//AOD
        if(customQualityCutsID) trkQualityCut->SetCheckOnlyFilterBit(kFALSE);
        if(customQualityCutsID==4){trkQualityCut->SetDCARPtFormula("0.006+0.0200/pt^1.1");}
        else if(customQualityCutsID==6){trkQualityCut->SetDCAZmax(0.2);}
        else if(customQualityCutsID==8){trkQualityCut->SetTrackMaxChi2(2.3);}
        else if(customQualityCutsID==10){trkQualityCut->SetMinNCrossedRowsTPC(100,kTRUE);}
        else if(customQualityCutsID==12){trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.9,kTRUE);}
        else if(customQualityCutsID==56){trkQualityCut->SetDCAZmax(1.);}
        else if(customQualityCutsID==58){trkQualityCut->SetTrackMaxChi2(3.5);}
        else if(customQualityCutsID==60){trkQualityCut->SetMinNCrossedRowsTPC(80,kTRUE);}
        else if(customQualityCutsID){
            Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
            return kFALSE;
        }
    }
    
    trkQualityCut->SetPtRange(0.15, 10000.0);
    trkQualityCut->SetEtaRange(-0.8, 0.8);
    
    Printf(Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));
    trkQualityCut->Print();
    return kTRUE;
}
