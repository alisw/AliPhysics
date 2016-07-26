/***************************************************************************
              
               fbellini@cern.ch - inayat.rasool.bhat@cern.ch Created on 07/03/2016 
// *** Configuration script for Rho,  analysis with 2013 pPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/
Bool_t ConfigRhoPPb
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    Int_t                  customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,//kTOFpidKstarPbPb2010,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,//kTOFpidKstarPbPb2010,
    Float_t                nsigmaPi = 3.0,
    Float_t                nsigmaKa = 3.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    TString                monitorOpt = "",
    // Bool_t                 useMixLS = 0,
    Bool_t                 checkReflex = 1,
    AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt
)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetPi;
  AliRsnCutSetDaughterParticle * cutSetKa;
  
  AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("myQualityCut");
  if (SetCustomQualityCut(trkQualityCut, customQualityCutsID, aodFilterBit)) {
    //Set custom quality cuts for systematic checks
    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0);
    cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi), trkQualityCut, cutPiCandidate, AliPID::kPion, nsigmaPi);
    cutSetKa  = new AliRsnCutSetDaughterParticle(Form("cutKa%i_%2.1fsigma",cutKaCandidate, nsigmaKa), trkQualityCut, cutKaCandidate, AliPID::kPion, nsigmaKa);
  } else {
    //use defult quality cuts (std 2010 or 2011)
    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, aodFilterBit, kTRUE);
    cutSetQ->SetUse2011StdQualityCuts(kTRUE);
    cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi), cutPiCandidate, AliPID::kPion, nsigmaPi, aodFilterBit,kTRUE);
    cutSetPi->SetUse2011StdQualityCuts(kTRUE);
    cutSetKa  = new AliRsnCutSetDaughterParticle(Form("cutKa%i_%2.1fsigma",cutKaCandidate, nsigmaKa), cutKaCandidate, AliPID::kPion, nsigmaKa, aodFilterBit,kTRUE);
    cutSetKa->SetUse2011StdQualityCuts(kTRUE);
   }
  
  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutKa = task->AddTrackCuts(cutSetKa);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetKa->GetMonitorOutput()), monitorOpt.Data();
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


  Bool_t  use     [12] = {!IsMcTrueOnly ,!IsMcTrueOnly,!IsMcTrueOnly ,  isMC    ,  isMC    ,    isMC   ,  isMC    ,  isMC      ,  isMC      ,  isMC      ,  isMC     ,  isMC    };
  Bool_t  useIM   [12] = {      1       ,     1       ,        1     ,    1     ,    1     ,    1      ,    1     ,    1       ,    1       ,    1       ,    1      ,   0      };
  TString name    [12] = {    "Unlike"  ,   "LikePP"  ,  "LikeMM"    , "Rho"    , "Omega"  , "KShort"  , "KStar"  ,"AntiKStar" , "Eta"      ,"Etaprime"  ,  "Phi"    ,  "Res"   };
  TString comp    [12] = {    "PAIR"    ,    "PAIR"   ,  "PAIR"      , "TRUE"   , "TRUE"   , "TRUE"    , "TRUE"   , "TRUE"     , "TRUE"     , "TRUE"     , "TRUE"    , "TRUE"   };
  TString output  [12] = {    "SPARSE"  ,   "SPARSE"  ,  "SPARSE"    , "SPARSE" , "SPARSE" , "SPARSE"  , "SPARSE" , "SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"  , "SPARSE" };
  Int_t   pdgCode [12] = {    113       ,    113      ,  113         , 113      ,  223     ,  310      ,  313     ,  -313      ,  221       ,  331       ,   333     ,   113    };
  Char_t  charge1 [12] = {    '+'       ,    '+'      ,  '-'         , '+'      ,   '+'    ,   '+'     ,   '+'    ,   '+'      ,   '+'      ,   '+'      ,   '+'     ,   '+'    };
  Char_t  charge2 [12] = {    '-'       ,    '+'      ,  '-'         , '-'      ,   '-'    ,   '-'     ,   '-'    ,   '-'      ,   '-'      ,    '-'     ,    '-'    ,   '-'    };
  Int_t   cutID1  [12] = {   iCutKa    ,  iCutKa     , iCutKa       , iCutKa   , iCutKa   , iCutKa    , iCutKa   , iCutKa     , iCutKa     , iCutKa     , iCutKa    ,  iCutKa   };
  Int_t   cutID2  [12] = {   iCutPi   ,  iCutPi     , iCutPi       , iCutPi   , iCutPi   , iCutPi    , iCutPi   , iCutPi     , iCutPi     , iCutPi     , iCutPi    ,  iCutPi   };
   
  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("rho_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID2[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kPion);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdgCode[i]);
    out->SetMotherMass(0.77546);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 180, 0.2, 2.0);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum of pair as default - else chosen value
    if (yaxisVar==AliRsnMiniValue::kFirstDaughterPt)
      out->AddAxis(fdpt, 70, 0.0, 7.0);
    else
      if (yaxisVar==AliRsnMiniValue::kSecondDaughterPt)
	out->AddAxis(sdpt, 70, 0.0, 7.0);
      else
	if (yaxisVar==AliRsnMiniValue::kFirstDaughterP)
	  out->AddAxis(fdp, 70, 0.0, 7.0);
	else
	  if (yaxisVar==AliRsnMiniValue::kSecondDaughterP)
	    out->AddAxis(sdp, 70, 0.0, 7.0);
	  else 
	    out->AddAxis(ptID, 140, 0.0, 7.0); //default use mother pt
    
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
    //get mothers for particle Rho PDG = 113
    AliRsnMiniOutput *outm = task->CreateOutput(Form("Rho_Mother%s", suffix), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kPion);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(113);
    outm->SetMotherMass(0.77546);
    outm->SetPairCuts(cutsPair);
    outm->AddAxis(imID, 180, 0.2, 2.0);
    outm->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
    
    //get mothers for particle Omega and Omega Reflexion PDG = 223
   AliRsnMiniOutput *outomegam = task->CreateOutput(Form("Omega_Mother%s", suffix), "SPARSE", "MOTHER");
    outomegam->SetDaughter(0, AliRsnDaughter::kPion);
    outomegam->SetDaughter(1, AliRsnDaughter::kPion);
    outomegam->SetMotherPDG(223);
    outomegam->SetMotherMass(0.78265);
    outomegam->SetPairCuts(cutsPair);
    outomegam->AddAxis(imID, 180, 0.2, 2.0);
    outomegam->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outomegam->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outomegam->AddAxis(centID, 400, 0.0, 400.0);
    }

   //get mothers for particle KShort PDG = 310
    AliRsnMiniOutput *outkshortm = task->CreateOutput(Form("Kshort_Mother%s", suffix), "SPARSE", "MOTHER");
    outkshortm->SetDaughter(0, AliRsnDaughter::kPion);
    outkshortm->SetDaughter(1, AliRsnDaughter::kPion);
    outkshortm->SetMotherPDG(310);
    outkshortm->SetMotherMass(0.497614);
    outkshortm->SetPairCuts(cutsPair);
    outkshortm->AddAxis(imID, 180, 0.2, 2.0);
    outkshortm->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outkshortm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outkshortm->AddAxis(centID, 400, 0.0, 400.0);
    }
     //get mothers for particle KStar PDG = 313
      AliRsnMiniOutput *outkstarm = task->CreateOutput(Form("Kstartokapi_Mother%s", suffix), "SPARSE", "MOTHER");
    outkstarm->SetDaughter(0, AliRsnDaughter::kKaon);
    outkstarm->SetDaughter(1, AliRsnDaughter::kPion);
    outkstarm->SetMotherPDG(313);
    outkstarm->SetMotherMass(0.89166);
    outkstarm->SetPairCuts(cutsPair);
    outkstarm->AddAxis(imID, 180, 0.2, 2.0);
    outkstarm->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outkstarm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outkstarm->AddAxis(centID, 400, 0.0, 400.0);
    }
    //get mothers for particle Anti-KStar PDG = -313
     AliRsnMiniOutput *outantkstarm = task->CreateOutput(Form("Antikstartokapi_Mother%s", suffix), "SPARSE", "MOTHER");
    outantkstarm->SetDaughter(0, AliRsnDaughter::kKaon);
    outantkstarm->SetDaughter(1, AliRsnDaughter::kPion);
    outantkstarm->SetMotherPDG(-313);
    outantkstarm->SetMotherMass(0.89166);
    outantkstarm->SetPairCuts(cutsPair);
    outantkstarm->AddAxis(imID, 180, 0.2, 2.0);
    outantkstarm->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outantkstarm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outantkstarm->AddAxis(centID, 400, 0.0, 400.0);
    }
      //get mothers for particle KStar Template  PDG = 313
       AliRsnMiniOutput *outkstarmm = task->CreateOutput(Form("Kstartopipii_Mother%s", suffix), "SPARSE", "MOTHER");
    outkstarmm->SetDaughter(0, AliRsnDaughter::kPion);
    outkstarmm->SetDaughter(1, AliRsnDaughter::kPion);
    outkstarmm->SetMotherPDG(313);
    outkstarmm->SetMotherMass(0.89166);
    outkstarmm->SetPairCuts(cutsPair);
    outkstarmm->AddAxis(imID, 180, 0.2, 2.0);
    outkstarmm->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outkstarmm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outkstarmm->AddAxis(centID, 400, 0.0, 400.0);
    }
     //get mothers for particle Anti-KStar Template PDG = -313
     AliRsnMiniOutput *outantkstarmm = task->CreateOutput(Form("Antikstartopipi_Mother%s", suffix), "SPARSE", "MOTHER");
    outantkstarmm->SetDaughter(0, AliRsnDaughter::kPion);
    outantkstarmm->SetDaughter(1, AliRsnDaughter::kPion);
    outantkstarmm->SetMotherPDG(-313);
    outantkstarmm->SetMotherMass(0.89166);
    outantkstarmm->SetPairCuts(cutsPair);
    outantkstarmm->AddAxis(imID, 180, 0.2, 2.0);
    outantkstarmm->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outantkstarmm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outantkstarmm->AddAxis(centID, 400, 0.0, 400.0);
    }
    //get mothers for particle Eta PDG = 221
      AliRsnMiniOutput *outantetam = task->CreateOutput(Form("eta_Mother%s", suffix), "SPARSE", "MOTHER");
    outantetam->SetDaughter(0, AliRsnDaughter::kPion);
    outantetam->SetDaughter(1, AliRsnDaughter::kPion);
    outantetam->SetMotherPDG(221);
    outantetam->SetMotherMass(0.547862);
    outantetam->SetPairCuts(cutsPair);
    outantetam->AddAxis(imID, 180, 0.2, 2.0);
    outantetam->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outantetam->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outantetam->AddAxis(centID, 400, 0.0, 400.0);
    }
     //get mothers for particle Eta-prime PDG = 331
      AliRsnMiniOutput *outantetaprimem = task->CreateOutput(Form("etaprime_Mother%s", suffix), "SPARSE", "MOTHER");
    outantetaprimem->SetDaughter(0, AliRsnDaughter::kPion);
    outantetaprimem->SetDaughter(1, AliRsnDaughter::kPion);
    outantetaprimem->SetMotherPDG(331);
    outantetaprimem->SetMotherMass(0.95778);
    outantetaprimem->SetPairCuts(cutsPair);
    outantetaprimem->AddAxis(imID, 180, 0.2, 2.0);
    outantetaprimem->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outantetaprimem->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outantetaprimem->AddAxis(centID, 400, 0.0, 400.0);
    }
     //get mothers for particle Phi Template PDG = 333
     AliRsnMiniOutput *outantphim = task->CreateOutput(Form("phitopipi_Mother%s", suffix), "SPARSE", "MOTHER");
    outantphim->SetDaughter(0, AliRsnDaughter::kPion);
    outantphim->SetDaughter(1, AliRsnDaughter::kPion);
    outantphim->SetMotherPDG(333);
   // outantphim->SetMotherMass(1.019461);
    outantphim->SetPairCuts(cutsPair);
    outantphim->AddAxis(imID, 180, 0.2, 2.0);
    outantphim->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outantphim->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outantphim->AddAxis(centID, 400, 0.0, 400.0);
    }
    //get mothers for particle Phi PDG = 333
     AliRsnMiniOutput *outantphimm = task->CreateOutput(Form("phitoKK_Mother%s", suffix), "SPARSE", "MOTHER");
    outantphimm->SetDaughter(0, AliRsnDaughter::kKaon);
    outantphimm->SetDaughter(1, AliRsnDaughter::kKaon);
    outantphimm->SetMotherPDG(333);
    outantphimm->SetMotherMass(1.019461);
    outantphimm->SetPairCuts(cutsPair);
    outantphimm->AddAxis(imID, 180, 0.2, 2.0);
    outantphimm->AddAxis(ptID, 140, 0.0, 7.0);
    if (!isPP){
      outantphimm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outantphimm->AddAxis(centID, 400, 0.0, 400.0);
    }
   
    
    //get reflections 
    //Reflections for KStar , Anti-KStar , and Phi 
    if (checkReflex) { 

      AliRsnMiniOutput *outreflex = task->CreateOutput(Form("Ks_reflex%s", suffix), "SPARSE", "TRUE");
      outreflex->SetDaughter(0, AliRsnDaughter::kKaon);
      outreflex->SetDaughter(1, AliRsnDaughter::kPion);
      outreflex->SetCutID(0, iCutPi);
      outreflex->SetCutID(1, iCutKa);
      outreflex->SetMotherPDG(313);
      outreflex->SetMotherMass(0.77546);//0.89594);
      outreflex->SetPairCuts(cutsPair);
      outreflex->AddAxis(imID, 180, 0.2, 2.0);
      outreflex->AddAxis(ptID, 140, 0.0, 7.0);
      if (!isPP){
	outreflex->AddAxis(centID, 100, 0.0, 100.0);
      }   else    { 
	outreflex->AddAxis(centID, 400, 0.0, 400.0);
      }
      
      AliRsnMiniOutput *outareflex = task->CreateOutput(Form("antiKs_reflex%s", suffix), "SPARSE", "TRUE");
      outareflex->SetDaughter(0, AliRsnDaughter::kKaon);
      outareflex->SetDaughter(1, AliRsnDaughter::kPion);
      outareflex->SetCutID(0, iCutPi);
      outareflex->SetCutID(1, iCutKa);
      outareflex->SetMotherPDG(-313);
      outareflex->SetMotherMass(0.77546);//0.89594);
      outareflex->SetPairCuts(cutsPair);
      outareflex->AddAxis(imID, 180, 0.2, 2.0);
      outareflex->AddAxis(ptID, 140, 0.0, 7.0);
      if (!isPP){
	outareflex->AddAxis(centID, 100, 0.0, 100.0);
      }   else    { 
	outareflex->AddAxis(centID, 400, 0.0, 400.0);
      }

       AliRsnMiniOutput *outphiflex = task->CreateOutput(Form("phi_reflex%s", suffix), "SPARSE", "TRUE");
      outphiflex->SetDaughter(0, AliRsnDaughter::kKaon);
      outphiflex->SetDaughter(1, AliRsnDaughter::kKaon);
      outphiflex->SetCutID(0, iCutKa);
      outphiflex->SetCutID(1, iCutKa);
      outphiflex->SetMotherPDG(333);
      outphiflex->SetMotherMass(0.77546);//1.01946);
      outphiflex->SetPairCuts(cutsPair);
      outphiflex->AddAxis(imID, 180, 0.2, 2.0);
      outphiflex->AddAxis(ptID, 140, 0.0, 7.0);
      if (!isPP){
	outphiflex->AddAxis(centID, 100, 0.0, 100.0);
      }   else    { 
	outphiflex->AddAxis(centID, 400, 0.0, 400.0);
      }

    }//end reflections

      //get phase space of the decay from mothers
    AliRsnMiniOutput *outps = task->CreateOutput(Form("Rho_phaseSpace%s", suffix), "HIST", "TRUE");
    outps->SetDaughter(0, AliRsnDaughter::kPion);
    outps->SetDaughter(1, AliRsnDaughter::kPion);
    outps->SetCutID(0, iCutPi);
    outps->SetCutID(1, iCutPi);
    outps->SetMotherPDG(113);
    outps->SetMotherMass(0.77546);
    outps->SetPairCuts(cutsPair);
    outps->AddAxis(fdpt, 50, 0.0, 5.0);
    outps->AddAxis(sdpt, 50, 0.0, 5.0);
    outps->AddAxis(ptID, 100, 0.0, 7.0);
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
