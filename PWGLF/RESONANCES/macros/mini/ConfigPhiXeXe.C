/***************************************************************************
fbellini@cern.ch - last modified on 16/04/2020
Configuration script for Phi analysis with 2017 Xe-Xe runs
****************************************************************************/
#if !defined (__CINT__) || defined (__CLING__)
#include "AddMonitorOutput.C"
#endif

Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut = NULL, Int_t customQualityCutsID = -1, Int_t customFilterBit = 0);
Bool_t UseCustomQualityCutsESD(AliRsnCutTrackQuality * trkQualityCut = NULL, Bool_t useGeoCut = 0, 
                              Bool_t useLowBcuts = 0, Bool_t useLooseDCAXYcut = 0, Bool_t useLooseDCAZcut = 0, 
                              Bool_t useLooseXRcut = 0, Bool_t useTightXEcut = 0, Bool_t useLooseXR2Cls =0, 
                              Bool_t useTightXR2Cls = 0, Bool_t useClsCut = 0, Bool_t useLooseChi2cut = 0);

Bool_t ConfigPhiXeXe(AliRsnMiniAnalysisTask *task = 0x0, 
		     Bool_t                 isMC = kFALSE, 
		     TString                suffix = "",
		     AliRsnCutSet           *cutsPair = 0x0,
		     Int_t                  aodFilterBit = 5,
		     Int_t                  customQualityCutsID = -1, //AliRsnCutSetDaughterParticle::kDisableCustom,
		     AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPid = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, // pid cut set
		     Float_t                nsigmaTPC = 2.0,          //nsigma of TPC PID cut
		     Float_t                nsigmaTOF = 3.0,          //nsigma of TOF PID cut
		     Float_t                ptMax = 12.0,
		     Bool_t                 enableMonitor = kTRUE,
		     Bool_t                 useMixLS = kFALSE,
		     Bool_t                 checkReflex = kFALSE)
{

  AliRsnCutSetDaughterParticle * cutSetQuality = NULL;
  AliRsnCutSetDaughterParticle * cutSetKa = NULL;

  AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality(Form("quality%i", customQualityCutsID));

  //Set custom quality cuts for systematic checks
  //use default quality cuts (std 2011)
  
  if (customQualityCutsID<0) {
    Printf("::::: ConfigPhiXeXe ::: Custom quality cuts disabled - using default");
    cutSetQuality  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, aodFilterBit, kTRUE);
    cutSetKa  = new AliRsnCutSetDaughterParticle("cutKa", cutPid, AliPID::kKaon, nsigmaTPC, nsigmaTOF, aodFilterBit, kTRUE);
    cutSetQuality->SetUse2011StdQualityCuts(kTRUE);
    cutSetKa->SetUse2011StdQualityCuts(kTRUE);
  }  else {
    if (customQualityCutsID>=100) {
      //Use ESD custom cuts -- since analysis releoad for final results
      Printf(Form("::::: ConfigPhiXeXe ::: Custom quality cuts enabled - using %i", customQualityCutsID));
      UseCustomQualityCutsESD(trkQualityCut, 
                              (customQualityCutsID==102 || customQualityCutsID == 103), //useGeoCut
                              (customQualityCutsID==101 || customQualityCutsID == 103), //useLowBcuts
                              !(customQualityCutsID%104), //useLooseDCAXYcut
                              !(customQualityCutsID%105), //useLooseDCAZcut
                              !(customQualityCutsID%106), //useLooseXRcut
                              !(customQualityCutsID%107), //useTightXRcut
                              !(customQualityCutsID%108), //useLooseXR2Cls
                              !(customQualityCutsID%109), //useTightXR2Cls
                              !(customQualityCutsID%110), //useClsCut
                              !(customQualityCutsID%111)); //useLooseChi2cut                           
    } else {
      //Use 2018 way of setting custom cuts -- used for preliminary QM2018, kept for bookeeping
      Printf(Form("::::: ConfigPhiXeXe ::: Custom quality cuts enabled (the OLD way) - using %i", customQualityCutsID));
      SetCustomQualityCut(trkQualityCut, customQualityCutsID, aodFilterBit);
    }
    cutSetQuality = new AliRsnCutSetDaughterParticle(trkQualityCut->GetName(), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, -1.0);
    cutSetKa = new AliRsnCutSetDaughterParticle("cutKa", trkQualityCut, cutPid, AliPID::kKaon, nsigmaTPC, nsigmaTOF);
  }

  if (!cutSetQuality || !cutSetKa) {
    Printf("ERROR: INVALID CUT SET in Configuration!");
    return kFALSE;
  }
  
  Int_t icutKa = task->AddTrackCuts(cutSetKa);
  Int_t icutQuality = task->AddTrackCuts(cutSetQuality);
  
  //set daughter cuts
  Int_t iCut1 = icutKa;
  Int_t iCut2 = icutKa;
  
#if defined (__CINT__) || !defined (__CLING__)
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
  //gROOT->LoadMacro("./AddMonitorOutputPhiXeXe.C");
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
  Bool_t  use     [8] = {   !isMC,    !isMC,    !isMC,    !isMC,   isMC,    isMC,   (!isMC & useMixLS), (!isMC & useMixLS)};
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
  Bool_t  usemc[4]     = {false, true, false, false};
  Int_t   pdgCodeR[4] = {333, 333, 333, 333};
  Char_t  charge1R[4] = {'+', '+', '+', '+'};
  Char_t  charge2R[4] = {'-', '-', '-', '-'};

  for (Int_t j = 0; j < 4; j++) {
    if (!usemc[j]) continue;
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

//-------------------------------------
Bool_t UseCustomQualityCutsESD(AliRsnCutTrackQuality * trkQualityCut, 
                              Bool_t useGeoCut, 
                              Bool_t useLowBcuts, 
                              Bool_t useLooseDCAXYcut,
                              Bool_t useLooseDCAZcut, 
                              Bool_t useLooseXRcut,
                              Bool_t useTightXRcut,
                              Bool_t useLooseXR2Cls,
                              Bool_t useTightXR2Cls,
                              Bool_t useClsCut,
                              Bool_t useLooseChi2cut)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.
  if (!trkQualityCut){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }
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
  trkQualityCut->SetDCARPtFormula("0.0105+0.0350/pt^1.1");
  trkQualityCut->SetDCAZmax(2.0); 
  if (useLowBcuts) trkQualityCut->SetDCARPtFormula("0.0119+0.049/pt"); // As in XeXe All Charged AN
  if (useLowBcuts & useGeoCut) ((AliESDtrackCuts  *)trkQualityCut->GetESDtrackCuts())->SetCutGeoNcrNcl(3, 130, 0.7, 0.85, 0.7); // As in XeXe All Charged AN
  if (useGeoCut) ((AliESDtrackCuts  *)trkQualityCut->GetESDtrackCuts())->SetCutGeoNcrNcl(3, 130, 0.7, 0.85, 0.7);  
  if (useLooseDCAXYcut) trkQualityCut->SetDCARmax(2.4);
  if (useLooseDCAZcut) trkQualityCut->SetDCAZmax(3.2);
  if (useLooseXRcut) trkQualityCut->SetMinNCrossedRowsTPC(60, kTRUE);
  if (useTightXRcut) trkQualityCut->SetMinNCrossedRowsTPC(80, kTRUE); 
  if (useLooseXR2Cls) trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.75, kTRUE);
  if (useTightXR2Cls) trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.85, kTRUE);  
  if (useClsCut) trkQualityCut->SetTPCminNClusters(70);
  if (useLooseChi2cut) trkQualityCut->SetTPCmaxChi2(3.5);
  trkQualityCut->SetPtRange(0.15, 20.0);
  trkQualityCut->SetEtaRange(-0.8, 0.8);  
  Printf("::::: SetCustomQualityCut:: using custom track quality cuts:");
  trkQualityCut->Print();
  return kTRUE;
}

//-------------------------------------------------------  
Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID, Int_t customFilterBit)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.
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

  if ((!trkQualityCut) || (customQualityCutsID<=0) || (customQualityCutsID>=AliRsnCutSetDaughterParticle::kNcustomQualityCuts)){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }

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
