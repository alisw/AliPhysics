/***************************************************************************
    fbellini@cern.ch, graham.richard.lee@cern.ch - created on 20/03/2014

 *** Configuration script for rho analysis of 2010 pp 7TeV datasets ***

****************************************************************************/
Bool_t ConfigRhoPP7TeV
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 38,
    Int_t                  customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPi1Candidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPi2Candidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
    Float_t                nsigmaPi1 = 3.0,
    Float_t                nsigmaPi2 = 3.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    TString                monitorOpt = "",
    Bool_t                 useMixLS = 0,
    AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt
)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetPi1;
  AliRsnCutSetDaughterParticle * cutSetPi2;
  
  AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("myQualityCut");
  if (SetCustomQualityCut(trkQualityCut, customQualityCutsID, aodFilterBit)) {
    //Set custom quality cuts for systematic checks
    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0);
    cutSetPi1 = new AliRsnCutSetDaughterParticle(Form("cutPi1_%i_%2.1fsigma",cutPi1Candidate, nsigmaPi1), trkQualityCut, cutPi1Candidate, AliPID::kPion, nsigmaPi1);
    cutSetPi2  = new AliRsnCutSetDaughterParticle(Form("cutPi2_%i_%2.1fsigma",cutPi2Candidate, nsigmaPi2), trkQualityCut, cutPi2Candidate, AliPID::kPion, nsigmaPi2);
  } else {
    //use default quality cuts std 2010 with crossed rows TPC
    Bool_t useCrossedRows = 1;
    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0, aodFilterBit, useCrossedRows);
    cutSetPi1 = new AliRsnCutSetDaughterParticle(Form("cutPi1_%i_%2.1fsigma",cutPi1Candidate, nsigmaPi1), cutPi1Candidate, AliPID::kPion, nsigmaPi1, aodFilterBit, useCrossedRows);
    cutSetPi2  = new AliRsnCutSetDaughterParticle(Form("cutPi2_%i_%2.1fsigma",cutPi2Candidate, nsigmaPi2), cutPi2Candidate, AliPID::kPion, nsigmaPi2, aodFilterBit, useCrossedRows);
  }
  
  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutPi1 = task->AddTrackCuts(cutSetPi1);
  Int_t iCutPi2 = task->AddTrackCuts(cutSetPi2);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetPi1->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetPi2->GetMonitorOutput()), monitorOpt.Data();
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
  Int_t   pdgCode [12] = {113       , 113       , 113      , 113      , 113     , 113     , 113     , -113    ,  113    , -113   ,  113      , 113         };
  Char_t  charge1 [12] = {'+'       , '-'       , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'    , '+'     , '-'    , '+'       , '-'         };
  Char_t  charge2 [12] = {'-'       , '+'       , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'    , '-'     , '+'    , '+'       , '-'         };
  Int_t   cutID1  [12] = { iCutPi1    ,  iCutPi1    ,  iCutPi1   ,  iCutPi1   ,  iCutPi1  ,  iCutPi1  ,  iCutPi1  ,   iCutPi1 ,  iCutPi1  , iCutPi1  , iCutPi1     , iCutPi1       };
  Int_t   cutID2  [12] = { iCutPi2   ,  iCutPi2   ,  iCutPi2  ,  iCutPi2  ,  iCutPi2 ,  iCutPi2 ,  iCutPi2 ,   iCutPi2,  iCutPi2 , iCutPi2 , iCutPi2    , iCutPi2      };
  
  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("rho_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kPion);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdgCode[i]);
    out->SetMotherMass(0.77549);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 180, 0.2, 2.0);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum of pair as default - else chosen value
    if (yaxisVar==AliRsnMiniValue::kFirstDaughterPt)
      out->AddAxis(fdpt, 100, 0.0, 10.0);
    else
      if (yaxisVar==AliRsnMiniValue::kSecondDaughterPt)
	out->AddAxis(sdpt, 100, 0.0, 10.0);
      else
	if (yaxisVar==AliRsnMiniValue::kFirstDaughterP)
	  out->AddAxis(fdp, 100, 0.0, 10.0);
	else
	  if (yaxisVar==AliRsnMiniValue::kSecondDaughterP)
	    out->AddAxis(sdp, 100, 0.0, 10.0);
	  else 
	    out->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
    
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
    //get mothers for rho PDG = 113
    AliRsnMiniOutput *outm = task->CreateOutput(Form("rho_Mother%s", suffix), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kPion);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(113);
    outm->SetMotherMass(0.77549);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID, 180, 0.2, 2.0);
    outm->AddAxis(ptID, 200, 0.0, 20.0);
    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
        
    //get phase space of the decay from mothers
    AliRsnMiniOutput *outps = task->CreateOutput(Form("rho_phaseSpace%s", suffix), "HIST", "TRUE");
    outps->SetDaughter(0, AliRsnDaughter::kPion);
    outps->SetDaughter(1, AliRsnDaughter::kPion);
    outps->SetCutID(0, iCutPi1);
    outps->SetCutID(1, iCutPi2);
    outps->SetMotherPDG(113);
    outps->SetMotherMass(0.77549);
    outps->SetPairCuts(cutsPair);
    outps->AddAxis(fdpt, 50, 0.0, 5.0);
    outps->AddAxis(sdpt, 50, 0.0, 5.0);
    outps->AddAxis(ptID, 200, 0.0, 20.0);

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
   //AliESDtrackCuts* esdTrackCutsH2Cluster = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 0);
   */

  if ((!trkQualityCut) || (customQualityCutsID<=0) || (customQualityCutsID>=AliRsnCutSetDaughterParticle::kNcustomQualityCuts)){
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
