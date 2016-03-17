/***************************************************************************
  rama.chandra.baral@cern.ch & sarita.sahoo@cern.ch - last modified on 20/08/2014

// *** Configuration script for L*, anti-L*, syst. analysis for pp and p-Pb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigLStar
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffixy,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,
    Float_t                nsigmaPr = 2.0,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableTrkSyst = kFALSE,
    Char_t                 DCAxyFormula[100] = "0.0182+0.035/pt^1.01",
    Double_t               dcazmax = 3.2,
    Double_t               minNcls = 70,
    Double_t               maxX2cls = 4.0,
    Float_t                TPCGlobalchi2 = 36,
    Double_t               ITSchi2 = 36,
    Double_t               minCrossedRows = 50.0,
    Double_t               maxClsCrossedRows = 0.8,

    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    Int_t                  signedPdg = 3124,

    TString                monitorOpt = "",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
    Bool_t                 useCrossedRows = kFALSE,
    AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt,
    Bool_t                 useMixLS = 0
)
{
  // manage suffix
  TString opt(suffixy);
  if (strlen(suffixy) > 0) opt = Form("_%s", suffixy);
  
  //cout<<"opt = "<<opt<<endl; return;

  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetP;
  AliRsnCutSetDaughterParticle * cutSetK;


  //vary track quality cuts for systematic checks
  if(enableTrkSyst) {
    AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("QualityCut");
    
    //trkQualityCut->DisableAll();//disable all cuts, filter bit, pT, eta, and DCAxy cuts will be reset later
    trkQualityCut->SetAODTestFilterBit(aodFilterBit);//reset the filter bit cut 
    trkQualityCut->SetCheckOnlyFilterBit(kFALSE);//tells the cut object to check all other cuts individually,

    if(useCrossedRows) {
      trkQualityCut->SetMinNCrossedRowsTPC(minCrossedRows, kTRUE);
      trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(maxClsCrossedRows, kTRUE); 
    } else {
      trkQualityCut->SetTPCminNClusters(minNcls);
    }

    trkQualityCut->SetDCARPtFormula(DCAxyFormula);
    trkQualityCut->SetDCAZmax(dcazmax);
    trkQualityCut->SetTPCmaxChi2(maxX2cls);
    trkQualityCut->SetRejectKinkDaughters(kTRUE);
    //trkQualityCut->SetSPDminNClusters(AliESDtrackCuts::kAny);
    trkQualityCut->SetSPDminNClusters(1);
    trkQualityCut->SetITSminNClusters(0);
    trkQualityCut->SetITSmaxChi2(ITSchi2);
    trkQualityCut->SetMaxChi2TPCConstrainedGlobal(TPCGlobalchi2);
    trkQualityCut->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);//already in defaults 2011
    trkQualityCut->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);//already in defaults 2011
    trkQualityCut->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);//already in defaults 2011
    trkQualityCut->SetPtRange(0.15, 20.0);
    trkQualityCut->SetEtaRange(-0.8, 0.8);
    trkQualityCut->Print();
    
    if(isPP) 
      {
	cutSetQ = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion,-1.0);
	cutSetP = new AliRsnCutSetDaughterParticle(Form("cutP%i_%2.1fsigma",cutPrCandidate, nsigmaPr), trkQualityCut, cutPrCandidate, AliPID::kProton, nsigmaPr);
	cutSetK = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa), trkQualityCut, cutKaCandidate, AliPID::kKaon, nsigmaKa);
      } 
   else 
     {    
       cutSetQ = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0);
       cutSetP = new AliRsnCutSetDaughterParticle(Form("cutP%i_%2.1fsigma",cutPrCandidate, nsigmaPr), trkQualityCut, cutPrCandidate, AliPID::kProton, nsigmaPr);
       cutSetP->SetUse2011StdQualityCuts(kTRUE);
       cutSetK = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa), trkQualityCut, cutKaCandidate, AliPID::kKaon, nsigmaKa);
       cutSetK->SetUse2011StdQualityCuts(kTRUE);
     }
  } 
  else {
    //default cuts 2010 for pp and 2011 for p-Pb
    if(isPP) {
      cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0, aodFilterBit,kFALSE);
      cutSetP  = new AliRsnCutSetDaughterParticle(Form("cutProton_%2.1fsigma",nsigmaPr), cutPrCandidate, AliPID::kProton, nsigmaPr, aodFilterBit,kFALSE);
      cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaon_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit,kFALSE);
    } else {
      cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, aodFilterBit,kFALSE);
      cutSetQ->SetUse2011StdQualityCuts(kTRUE);
      cutSetP = new AliRsnCutSetDaughterParticle(Form("cutProton2011_%2.1fsigma",nsigmaPr), cutPrCandidate, AliPID::kProton, nsigmaPr, aodFilterBit,kFALSE);
      cutSetP->SetUse2011StdQualityCuts(kTRUE);
      cutSetK = new AliRsnCutSetDaughterParticle(Form("cutKaon2011_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit,kFALSE);
      cutSetK->SetUse2011StdQualityCuts(kTRUE);
    }
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

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [12] = { !IsMcTrueOnly, !IsMcTrueOnly, !IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly, !IsMcTrueOnly , isMC  ,  isMC   , isMC   ,  isMC  , useMixLS , useMixLS};
  Bool_t  useIM   [12] = { 1            ,  1           ,  1           ,  1          ,  1          , 1           ,  1      ,   1     , 0      ,   0    , 1        , 1       };
  TString name    [12] = {"UnlikePM"    , "UnlikeMP"   , "MixingPM"   , "MixingMP"  , "LikePP"    , "LikeMM"    ,"TruesPM","TruesMP","ResPM" ,"ResMP" ,"MixingPP","MixingMM"};
  TString comp    [12] = {"PAIR"        , "PAIR"       , "MIX"        , "MIX"       , "PAIR"      , "PAIR"      , "TRUE"  , "TRUE"  , "TRUE" , "TRUE" , "MIX"    ,"MIX"   };
  TString output  [12] = {"SPARSE"      , "SPARSE"     , "SPARSE"     , "SPARSE"    , "SPARSE"    , "SPARSE"    , "SPARSE","SPARSE" ,"SPARSE","SPARSE","SPARSE"  ,"SPARSE"};
  Char_t  charge1 [12] = {'+'           , '-'          , '+'          , '-'         , '+'         , '-'         , '+'     ,  '-'    , '+'    ,  '-'   , '+'      , '-'    };
  Char_t  charge2 [12] = {'-'           , '+'          , '-'          , '+'         , '+'         , '-'         , '-'     ,  '+'    , '-'    ,  '+'   ,'+'       , '-'    };
  Int_t   cutID1  [12] = { iCutP        ,  iCutP       ,  iCutP       ,  iCutP      ,  iCutP      ,  iCutP      ,  iCutP  ,  iCutP  ,  iCutP ,  iCutP , iCutP    , iCutP };
  Int_t   cutID2  [12] = { iCutK        ,  iCutK       ,  iCutK       ,  iCutK      ,  iCutK      ,  iCutK      ,  iCutK  ,  iCutK  ,  iCutK ,  iCutK , iCutK    , iCutK };
  
  //TString output  [10] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };

  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("Lstar_%s%s", name[i].Data(), opt.Data()), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(signedPdg);
    out->SetMotherMass(1.51953);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 800, 1.4, 2.2);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum of pair as default - else chosen value
    if (yaxisVar==AliRsnMiniValue::kFirstDaughterPt) {
      out->AddAxis(fdpt, 100, 0.0, 10.0);
      out->AddAxis(sdpt, 100, 0.0, 10.0);  }
    else
      if (yaxisVar==AliRsnMiniValue::kFirstDaughterP) {
	out->AddAxis(fdp, 100, 0.0, 10.0);
	out->AddAxis(sdp, 100, 0.0, 10.0);  }
      else 
	out->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt

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
    // create output
    AliRsnMiniOutput *outm = task->CreateOutput(Form("Lstar_Mother%s", opt.Data()), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kProton);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(signedPdg);
    outm->SetMotherMass(1.51953);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 800, 1.4, 2.2);

    // axis Y: transverse momentum of pair as default - else chosen value
    if (yaxisVar==AliRsnMiniValue::kFirstDaughterPt) {
      outm->AddAxis(fdpt, 100, 0.0, 10.0);
      outm->AddAxis(sdpt, 100, 0.0, 10.0);  }
    else
      if (yaxisVar==AliRsnMiniValue::kFirstDaughterP) {
	outm->AddAxis(fdp, 100, 0.0, 10.0);
	outm->AddAxis(sdp, 100, 0.0, 10.0);  }
      else 
	outm->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt


    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
  }
  return kTRUE;
}
