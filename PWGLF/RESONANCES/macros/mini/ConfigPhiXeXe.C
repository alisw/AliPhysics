/***************************************************************************
              fbellini@cern.ch - last modified on 17/02/2014

// *** Configuration script for K*, anti-K* analysis with 2013 pPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/
Bool_t ConfigPhiXeXe(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPid = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, // pid cut set
    Float_t                nsigma = 3.0,          //nsigma of TPC PID cut
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    TString                monitorOpt = "NoSIGN",
    Bool_t                 useMixLS = 0,
    Bool_t                 checkReflex = 0,
    AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt
)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  //use default quality cuts std 2010 with crossed rows TPC
  Bool_t useCrossedRows = 1; 
  AliRsnCutSetDaughterParticle * cutSetKa = new AliRsnCutSetDaughterParticle("cutKa", cutPid, AliPID::kKaon, nsigma, aodFilterBit, useCrossedRows);
  cutSetKa->SetUse2011StdQualityCuts(kTRUE);
  Int_t iCutKa = task->AddTrackCuts(cutSetKa);

  //set daughter cuts
  Int_t iCut1 = iCutKa;
  Int_t iCut2 = iCutKa;

  //monitor single-track selection based on track quality cuts only
  AliRsnCutSetDaughterParticle * cutSetQuality = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, 10.0, aodFilterBit, useCrossedRows);
  Int_t iCutQuality = task->AddTrackCuts(cutSetQuality);
  
  //QA plots 
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQuality->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetKa->GetMonitorOutput(), monitorOpt.Data());
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
  /* pair pt res      */ Int_t resPt  = task->CreateValue(AliRsnMiniValue::kPairPtRes, kTRUE);
  /* pair y res       */ Int_t resY   = task->CreateValue(AliRsnMiniValue::kPairYRes, kTRUE);
  
  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0][1] = unlike +-, -+ (Minv, pT, multiplicity)
  // [2][3] = mixing +-, -+ (Minv, pT, multiplicity)
  // [4][5] = like ++, -- (Minv, pT, multiplicity)
  // [6][7] = MC true +-, -+ (Minv, pT, multiplicity)
  // [8][9] = MC true +-, -+ (Minv, pT, y)
  // [10][11] = mixing ++, -- (Minv, pT, multiplicity)
  Bool_t  use     [8] = {   !isMC,    !isMC,    !isMC,    !isMC,   isMC,    isMC,   useMixLS, useMixLS};
  TString name    [8] = {"Unlike", "Mixing", "LikePP", "LikeMM", "True", "TrueY", "MixingPP", "MixingMM"};
  TString comp    [8] = {"PAIR"  , "MIX"   , "PAIR"  , "PAIR"  , "TRUE", "TRUE" , "MIX"     , "MIX"};
  Int_t   pdgCode [8] = {333     , 333     , 333     , 333     , 333   , 333    , 333       , 333  };
  Char_t  charge1 [8] = {'+'     , '+'     , '+'     , '-'     , '+'   , '+'    , '+'       ,  '-' };
  Char_t  charge2 [8] = {'-'     , '-'     , '+'     , '-'     , '-'   , '-'    , '+'       ,  '-' };
  TString output  = "HIST";
  
  /*********************
     Data and MC true
    *******************/
  for (Int_t i = 0; i < 8; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("%s%s", name[i].Data(), suffix), output.Data(), comp[i].Data());
    out->SetCutID(0, iCut1);
    out->SetCutID(1, iCut2);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdgCode[i]);
    out->SetMotherMass(1.01995);
    out->SetPairCuts(cutsPair);
    
    // axis X: invmass
    out->AddAxis(imID, 350, 0.85, 1.2);
    // axis Y: transverse momentum of pair as default - else chosen value
    if (yaxisVar==AliRsnMiniValue::kFirstDaughterPt)
      out->AddAxis(fdpt, 100, 0.0, 5.0);
    else
      if (yaxisVar==AliRsnMiniValue::kSecondDaughterPt)
	out->AddAxis(sdpt, 100, 0.0, 5.0);
      else
	if (yaxisVar==AliRsnMiniValue::kFirstDaughterP)
	  out->AddAxis(fdp, 100, 0.0, 5.0);
	else
	  if (yaxisVar==AliRsnMiniValue::kSecondDaughterP)
	    out->AddAxis(sdp, 100, 0.0, 5.0);
	  else 
	    out->AddAxis(ptID, 100, 0.0, 5.0); //default use mother pt
    
    // axis Z: centrality or multiplicity or rapidity
    if (i==5) {
      out->AddAxis(yID, 400, -2.0, 2.0);
    } else {
      out->AddAxis(centID, 100, 0.0, 100.0);
    }
  }   
  
  /****************
     MONTECARLO
  ****************/
  if (isMC){   

    //Minv resolution, pair's pT and pair's y resolution vs pT vs y
    //Computed as (Xrec-Xgen)/Xgen, with a MC-true like computation
    TString nameR[3]    = {"Res", "ResPt", "ResY"};
    TString compR[3]    = {"TRUE" , "TRUE", "TRUE"};
    TString outputR[3]  = {"HIST","HIST", "HIST"};
    Int_t   pdgCodeR[3] = {333,   333, 333};
    Char_t  charge1R[3] = {'+', '+', '+'};
    Char_t  charge2R[3] = {'-', '-', '-'};

    for (Int_t j = 0; j < 3; j++) {
      AliRsnMiniOutput *outR = task->CreateOutput(Form("%s%s", nameR[j].Data(), suffix), outputR[j].Data(), compR[j].Data());
      outR->SetCutID(0, iCut1);
      outR->SetCutID(1, iCut2);
      outR->SetDaughter(0, AliRsnDaughter::kKaon);
      outR->SetDaughter(1, AliRsnDaughter::kKaon);
      outR->SetCharge(0, charge1R[j]);
      outR->SetCharge(1, charge2R[j]);
      outR->SetMotherPDG(pdgCodeR[j]);
      outR->SetMotherMass(1.01995);
      outR->SetPairCuts(cutsPair);
      // axis X: invmass resolution, pt resolution, y resolution
      if (j<2) outR->AddAxis(resID, 100, -0.01, 0.01);
      else if (j<4) outR->AddAxis(resPt, 100, -0.05, 0.05);
      else outR->AddAxis(resY, 100, -0.5, 0.5);
      //axis Y: mother pt
      outR->AddAxis(ptID, 100, 0.0, 5.0); 
      //axis Z: rapidity
      outR->AddAxis(yID, 400, -2.0, 2.0);
    }

    //get mothers for PDG = 333
    AliRsnMiniOutput *outm = task->CreateOutput(Form("Mother%s", suffix), "HIST", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kKaon);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(333);
    outm->SetMotherMass(1.01995);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID, 350, 0.85, 1.2);
    outm->AddAxis(ptID, 100, 0.0, 5.0);
    outm->AddAxis(centID, 100, 0.0, 100.0);
      
    //get mothers for PDG = 333
    AliRsnMiniOutput *outm2 = task->CreateOutput(Form("MotherY%s", suffix), "HIST", "MOTHER");
    outm2->SetDaughter(0, AliRsnDaughter::kKaon);
    outm2->SetDaughter(1, AliRsnDaughter::kKaon);
    outm2->SetMotherPDG(333);
    outm2->SetMotherMass(1.01995);
    outm2->SetPairCuts(cutsPair);
    outm2->AddAxis(imID, 350, 0.85, 1.2);
    outm2->AddAxis(ptID, 100, 0.0, 5.0);
    outm2->AddAxis(yID, 400, -2.0, 2.0);
    
    //get phase space of the decay from mothers
    AliRsnMiniOutput *outps = task->CreateOutput(Form("PhaseSpace%s", suffix), "HIST", "TRUE");
    outps->SetDaughter(0, AliRsnDaughter::kKaon);
    outps->SetDaughter(1, AliRsnDaughter::kKaon);
    outps->SetCutID(0, iCut1);
    outps->SetCutID(1, iCut2);
    outps->SetMotherPDG(333);
    outps->SetMotherMass(1.01995);
    outps->SetPairCuts(cutsPair);
    outps->AddAxis(fdpt, 50, 0.0, 5.0);
    outps->AddAxis(sdpt, 50, 0.0, 5.0);
    outps->AddAxis(ptID, 100, 0.0, 10.0);
    
    //get reflections
    //defined as MC-true like computation but checking what happens when 
    //pions are mis-identified as K and K are mis-identified as pions
    if (checkReflex) { 

      AliRsnMiniOutput *outreflex = task->CreateOutput(Form("Reflex%s", suffix), "HIST", "TRUE");
      outreflex->SetDaughter(0, AliRsnDaughter::kKaon);
      outreflex->SetDaughter(1, AliRsnDaughter::kKaon);
      outreflex->SetCutID(0, iCut1);
      outreflex->SetCutID(1, iCut2);
      outreflex->SetMotherPDG(333);
      outreflex->SetMotherMass(1.01995);
      outreflex->SetPairCuts(cutsPair);
      outreflex->AddAxis(imID, 350, 0.85, 1.2);
      outreflex->AddAxis(ptID, 100, 0.0, 5.0);
      outreflex->AddAxis(centID, 100, 0.0, 100.0);
      
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
