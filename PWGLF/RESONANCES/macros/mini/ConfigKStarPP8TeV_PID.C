/***************************************************************************
              fbellini@cern.ch - last modified on 02/07/2014
//modified for K*(892)0 at 8 TeV
 *** Configuration script for K*, anti-K* analysis of 2010 pp 7TeV datasets ***
This analysis task is used to extend the pT reach of the K* spectra published in Eur.
Phys. J. C72(2012)2183. 
****************************************************************************/
Bool_t ConfigKStarPP8TeV_PID
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    Int_t                  customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kQualityStd2010,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kQualityStd2010,
    Float_t                nsigmaPi = 2.0,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    TString                monitorOpt = "",
    Bool_t                 useMixLS = 0,
    Bool_t                 checkReflex = 0,
    AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt
)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetPi;
  AliRsnCutSetDaughterParticle * cutSetK;


 Float_t nsigmaPiTPC=fmod(nsigmaPi,1000.);
 Float_t nsigmaPiTOF=(nsigmaPi-fmod(nsigmaPi,1000.))/1000.;
 if(nsigmaPiTOF<1.e-10) nsigmaPiTOF=-1.;
 
 Float_t nsigmaKaTPC=fmod(nsigmaKa,1000.);
 Float_t nsigmaKaTOF=(nsigmaKa-fmod(nsigmaKa,1000.))/1000.;
 if(nsigmaKaTOF<1.e-10) nsigmaKaTOF=-1.;


  
  AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("myQualityCut");
  if (SetCustomQualityCut(trkQualityCut, customQualityCutsID, aodFilterBit)) {
    //Set custom quality cuts for systematic checks
    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0);
    cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi), trkQualityCut, cutPiCandidate, AliPID::kPion, nsigmaPiTPC,nsigmaPiTOF);
    cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutPiCandidate, nsigmaKa), trkQualityCut, cutKaCandidate, AliPID::kKaon, nsigmaKaTPC,nsigmaKaTOF);
  } else {

    //use default quality cuts std 2010 with crossed rows TPC
    Bool_t useCrossedRows = 1;
    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0, aodFilterBit, useCrossedRows);
    cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi), cutPiCandidate, AliPID::kPion, nsigmaPi, aodFilterBit, useCrossedRows);
    cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutPiCandidate, nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit, useCrossedRows);
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
  
  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("kstar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdgCode[i]);
    out->SetMotherMass(0.89594);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 90, 0.6, 1.5);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum of pair as default - else chosen value
    // if (yaxisVar==AliRsnMiniValue::kFirstDaughterPt)
    //   out->AddAxis(fdpt, 100, 0.0, 10.0);
    // else
    //   if (yaxisVar==AliRsnMiniValue::kSecondDaughterPt)
    // 	out->AddAxis(sdpt, 100, 0.0, 10.0);
    //   else
    // 	if (yaxisVar==AliRsnMiniValue::kFirstDaughterP)
    // 	  out->AddAxis(fdp, 100, 0.0, 10.0);
    // 	else
    // 	  if (yaxisVar==AliRsnMiniValue::kSecondDaughterP)
    // 	    out->AddAxis(sdp, 100, 0.0, 10.0);
    // 	  else
    
    out->AddAxis(ptID, 500, 0.0, 50.0); //default use mother pt
    
    // axis Z: centrality-multiplicity
    if (!isPP)
      out->AddAxis(centID, 100, 0.0, 100.0);
    else 
      out->AddAxis(centID, 400, 0.0, 400.0);
    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    // out->AddAxis(yID, 10, -0.5, 0.5);
  }   
  
  if (isMC){   
    //get mothers for K* PDG = 313
    AliRsnMiniOutput *outm = task->CreateOutput(Form("Ks_Mother%s", suffix), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kKaon);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(313);
    outm->SetMotherMass(0.89594);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID, 90, 0.6, 1.5);
    outm->AddAxis(ptID, 500, 0.0, 50.0);
    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
    
    //get mothers for antiK* PDG = -313
    AliRsnMiniOutput *outam = task->CreateOutput(Form("antiKs_Mother%s", suffix), "SPARSE", "MOTHER");
    outam->SetDaughter(0, AliRsnDaughter::kKaon);
    outam->SetDaughter(1, AliRsnDaughter::kPion);
    outam->SetMotherPDG(-313);
    outam->SetMotherMass(0.89594);
    outam->SetPairCuts(cutsPair);
    outam->AddAxis(imID, 90, 0.6, 1.5);
    outam->AddAxis(ptID, 500, 0.0, 50.0);
    if (!isPP){
      outam->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outam->AddAxis(centID, 400, 0.0, 400.0);
    }
    
    //get phase space of the decay from mothers
    AliRsnMiniOutput *outps = task->CreateOutput(Form("Ks_phaseSpace%s", suffix), "HIST", "TRUE");
    outps->SetDaughter(0, AliRsnDaughter::kKaon);
    outps->SetDaughter(1, AliRsnDaughter::kPion);
    outps->SetCutID(0, iCutK);
    outps->SetCutID(1, iCutPi);
    outps->SetMotherPDG(313);
    outps->SetMotherMass(0.89594);
    outps->SetPairCuts(cutsPair);
    outps->AddAxis(fdpt, 50, 0.0, 5.0);
    outps->AddAxis(sdpt, 50, 0.0, 5.0);
    outps->AddAxis(ptID, 500, 0.0, 50.0);
    
    AliRsnMiniOutput *outaps = task->CreateOutput(Form("antiKs_phaseSpace%s", suffix), "HIST", "TRUE");
    outaps->SetDaughter(0, AliRsnDaughter::kKaon);
    outaps->SetDaughter(1, AliRsnDaughter::kPion);
    outaps->SetCutID(0, iCutK);
    outaps->SetCutID(1, iCutPi);
    outaps->SetMotherPDG(-313);
    outaps->SetMotherMass(0.89594);
    outaps->SetPairCuts(cutsPair);
    outaps->AddAxis(fdpt, 50, 0.0, 5.0);
    outaps->AddAxis(sdpt, 50, 0.0, 5.0);
    outaps->AddAxis(ptID, 500, 0.0, 50.0);
   
    //get reflections
    if (checkReflex) { 

      AliRsnMiniOutput *outreflex = task->CreateOutput(Form("Ks_reflex%s", suffix), "SPARSE", "TRUE");
      outreflex->SetDaughter(0, AliRsnDaughter::kKaon);
      outreflex->SetDaughter(1, AliRsnDaughter::kPion);
      outreflex->SetCutID(0, iCutPi);
      outreflex->SetCutID(1, iCutK);
      outreflex->SetMotherPDG(313);
      outreflex->SetMotherMass(0.89594);
      outreflex->SetPairCuts(cutsPair);
      outreflex->AddAxis(imID, 90, 0.6, 1.5);
      outreflex->AddAxis(ptID, 500, 0.0, 50.0);
      if (!isPP){
	outreflex->AddAxis(centID, 100, 0.0, 100.0);
      }   else    { 
	outreflex->AddAxis(centID, 400, 0.0, 400.0);
      }
      
      AliRsnMiniOutput *outareflex = task->CreateOutput(Form("antiKs_reflex%s", suffix), "SPARSE", "TRUE");
      outareflex->SetDaughter(0, AliRsnDaughter::kKaon);
      outareflex->SetDaughter(1, AliRsnDaughter::kPion);
      outareflex->SetCutID(0, iCutPi);
      outareflex->SetCutID(1, iCutK);
      outareflex->SetMotherPDG(-313);
      outareflex->SetMotherMass(0.89594);
      outareflex->SetPairCuts(cutsPair);
      outareflex->AddAxis(imID, 90, 0.6, 1.5);
      outareflex->AddAxis(ptID, 500, 0.0, 50.0);
      if (!isPP){
	outareflex->AddAxis(centID, 100, 0.0, 100.0);
      }   else    { 
	outareflex->AddAxis(centID, 400, 0.0, 400.0);
      }

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
  
  /* NOTES FROM PRODUCTION LHC13b pass3 - AOD filtered with v5-03-Rev-20
  //(http://svnweb.cern.ch/world/wsvn/AliRoot/tags/v5-03-Rev-20/ANALYSIS/macros/AddTaskESDFilter.C)

  //filter bit 0: Cuts on primary tracks
  // AliESDtrackCuts* esdTrackCutsL = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();

  //filter bit 4: std but looser dca cut
  // AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  // esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
  // esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
  // esdTrackCutsH->SetDCAToVertex2D(kTRUE);

  //filter bit 5:  AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

   //filter bit 10: standard cuts with tight DCA cut, using cluster cut instead of crossed rows (a la 2010 default)
   //AliESDtrackCuts* esdTrackCutsH2Cluster = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 0); */



 if ((!trkQualityCut)){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }

  if(customQualityCutsID>=1 && customQualityCutsID<100 && customQualityCutsID!=2){
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    trkQualityCut->SetPtRange(0.15, 2000.0);
    trkQualityCut->SetEtaRange(-0.8, 0.8);
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
