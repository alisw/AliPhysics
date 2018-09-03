/***************************************************************************
fbellini@cern.ch - last modified on 23/03/2018
Configuration script forphi analysis with 2017 XeXe runs
****************************************************************************/
#if !defined (__CINT__) || defined (__CLING__)
#include "AddMonitorOutput.C"
#endif

Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut = NULL, Int_t customQualityCutsID = 0, Int_t customFilterBit = 0);

Bool_t ConfigPhiVn(AliRsnMiniTaskPhiVn *task = 0x0, 
		     Bool_t                 isMC = kFALSE, 
		     TString                suffix = "",
		     AliRsnCutSet           *cutsPair = 0x0,
		     Int_t                  aodFilterBit = 5,
		     AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPid = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, // pid cut set
		     Float_t                nsigmaTPC = 2.0,          //nsigma of TPC PID cut
		     Float_t                nsigmaTOF = 3.0,          //nsigma of TOF PID cut
		     Float_t                ptMax = 15.0,
		     Bool_t                 enableMonitor = kTRUE,
		     Bool_t                 useMixLS = kFALSE,
		     Bool_t                 checkReflex = kFALSE)
{

  //use default quality cuts std 2010 with crossed rows TPC
  Bool_t useCrossedRows = 1; 
  AliRsnCutSetDaughterParticle * cutSetKa = new AliRsnCutSetDaughterParticle("cutKa", cutPid, AliPID::kKaon, nsigmaTPC, nsigmaTOF, aodFilterBit, useCrossedRows);
  cutSetKa->SetUse2011StdQualityCuts(kTRUE);
  Int_t icutKa = task->AddTrackCuts(cutSetKa);
  
  AliRsnCutSetDaughterParticle * cutSetKaNEW = new AliRsnCutSetDaughterParticle("cutNEW", cutPid, AliPID::kKaon, nsigmaTPC, nsigmaTOF, aodFilterBit, useCrossedRows);
  cutSetKaNEW->SetUse2011StdQualityCuts(kTRUE);
  Int_t icutKaNEW = task->AddTrackCuts(cutSetKaNEW);


  //set daughter cuts
  Int_t iCut1 = icutKa;
  Int_t iCut2 = icutKa;
  
  //monitor single-track selection based on track quality cuts only
  AliRsnCutSetDaughterParticle * cutSetQuality = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, 10.0, 10.0, aodFilterBit, useCrossedRows);
  Int_t icutQuality = task->AddTrackCuts(cutSetQuality);
  
#if defined (__CINT__) || !defined (__CLING__)
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
  //gROOT->LoadMacro("./AddMonitorOutput.C");
#endif
    
  //QA plots
  TString monitorOpt = "NoSIGN";
  if (enableMonitor){
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

  const Int_t nptbins = ptMax/0.1;
  
  // -- Create all needed outputs -----------------------------------------------------------------
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
    AliRsnMiniOutput *out = task->CreateOutput(Form("%s%s", name[i].Data(), suffix.Data()), output.Data(), comp[i].Data());
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
    out->AddAxis(imID, 420, 0.98, 1.4);

    // axis Y: transverse momentum of pair as default - else chosen value
    if (i==5) out->AddAxis(yID, 100, -0.5, 0.5);
    else out->AddAxis(ptID, nptbins, 0.0, ptMax); //default use mother pt
    
    // axis Z: centrality or multiplicity or rapidity
    out->AddAxis(centID, 100, 0.0, 100.0);
  }   
  
  if (!isMC) return kTRUE;   
  
  /****************
     MONTECARLO
  ****************/
  //Minv resolution, pair's pT and pair's y resolution vs pT vs y
  //Computed as (Xrec-Xgen)/Xgen, with a MC-true like computation
  TString nameR[4]    = {"Res", "ResCent", "ResPt", "ResY"};
  TString compR[4]    = {"TRUE" , "TRUE", "TRUE", "TRUE"};
  TString outputR[4]  = {"HIST", "HIST","HIST", "HIST"};
  Int_t   pdgCodeR[4] = {333, 333, 333, 333};
  Char_t  charge1R[4] = {'+', '+', '+', '+'};
  Char_t  charge2R[4] = {'-', '-', '-', '-'};

  for (Int_t j = 0; j < 4; j++) {
    AliRsnMiniOutput *outR = task->CreateOutput(Form("%s%s", nameR[j].Data(), suffix.Data()), outputR[j].Data(), compR[j].Data());
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
    switch (j) {
    case 0:
      outR->AddAxis(resID, 100, -0.01, 0.01);
      outR->AddAxis(ptID, nptbins, 0.0, ptMax); 
      outR->AddAxis(yID, 100, -0.5, 0.5);
      break;
    case 1:
      outR->AddAxis(resID, 100, -0.01, 0.01);
      outR->AddAxis(ptID, nptbins, 0.0, ptMax); 
      outR->AddAxis(centID, 100, 0.0, 100.0);
      break;
    case 2:
      outR->AddAxis(resPt, 100, -0.05, 0.05);
      outR->AddAxis(ptID, nptbins, 0.0, ptMax); 
      outR->AddAxis(yID, 100, -0.5, 0.5);
      break;
    case 3:
      outR->AddAxis(resY, 100, -0.05, 0.05);
      outR->AddAxis(ptID, nptbins, 0.0, ptMax); 
      outR->AddAxis(yID, 100, -0.5, 0.5);
      break;
    default:
      break;
    }
  }
  //get mothers for PDG = 333
  AliRsnMiniOutput *outm = task->CreateOutput(Form("Mother%s", suffix.Data()), "HIST", "MOTHER");
  outm->SetDaughter(0, AliRsnDaughter::kKaon);
  outm->SetDaughter(1, AliRsnDaughter::kKaon);
  outm->SetMotherPDG(333);
  outm->SetMotherMass(1.01995);
  outm->SetPairCuts(cutsPair);
  outm->AddAxis(imID, 420, 0.98, 1.4);
  outm->AddAxis(ptID, nptbins, 0.0, ptMax);
  outm->AddAxis(centID, 100, 0.0, 100.0);
      
  //get mothers for PDG = 333
  AliRsnMiniOutput *outm2 = task->CreateOutput(Form("MotherY%s", suffix.Data()), "HIST", "MOTHER");
  outm2->SetDaughter(0, AliRsnDaughter::kKaon);
  outm2->SetDaughter(1, AliRsnDaughter::kKaon);
  outm2->SetMotherPDG(333);
  outm2->SetMotherMass(1.01995);
  outm2->SetPairCuts(cutsPair);
  outm2->AddAxis(imID, 420, 0.98, 1.4);
  outm2->AddAxis(yID, 100, -0.5, 0.5);
  outm2->AddAxis(centID, 100, 0.0, 100.0);
	
  //get phase space of the decay from mothers
  AliRsnMiniOutput *outps = task->CreateOutput(Form("PhaseSpace%s", suffix.Data()), "HIST", "TRUE");
  outps->SetDaughter(0, AliRsnDaughter::kKaon);
  outps->SetDaughter(1, AliRsnDaughter::kKaon);
  outps->SetCutID(0, iCut1);
  outps->SetCutID(1, iCut2);
  outps->SetMotherPDG(333);
  outps->SetMotherMass(1.01995);
  outps->SetPairCuts(cutsPair);
  outps->AddAxis(fdpt, 50, 0.0, 5.0);
  outps->AddAxis(sdpt, 50, 0.0, 5.0);
  outps->AddAxis(ptID, nptbins, 0.0, ptMax);
    
  //get reflections
  //defined as MC-true like computation but checking what happens when 
  //pions are mis-identified as K and K are mis-identified as pions
  if (checkReflex) { 
    AliRsnMiniOutput *outreflex = task->CreateOutput(Form("Reflex2pi_%s", suffix.Data()), "HIST", "TRUE");
    outreflex->SetDaughter(0, AliRsnDaughter::kKaon);
    outreflex->SetDaughter(1, AliRsnDaughter::kKaon);
    outreflex->SetCutID(0, iCut1);
    outreflex->SetCutID(1, iCut2);
    outreflex->SetMotherPDG(333);
    outreflex->SetMotherMass(1.01995);
    outreflex->SetPairCuts(cutsPair);
    outreflex->AddAxis(imID, 420, 0.98, 1.4);
    outreflex->AddAxis(ptID, nptbins, 0.0, ptMax);
    outreflex->AddAxis(centID, 100, 0.0, 100.0);     
  }//end reflections

  
  return kTRUE;
}

//-------------------------------------------------------  
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID, Int_t customFilterBit)
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
  
  Printf("::::: SetCustomQualityCut:: using custom track quality cuts %i",customQualityCutsID);
  trkQualityCut->Print();
  return kTRUE;
}
