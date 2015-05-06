/***************************************************************************
  sarita.sahoo@cern.ch - last modified on 12/02/2015

// *** Configuration script for Phi*, anti-Phi*, for service analysis task for pp and p-Pb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigPhiMassStudy
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableTrkSyst = kFALSE,
    Char_t                 DCAxyFormula[100] = "0.0105+0.035/pt^1.01",
    Double_t               dcazmax = 3.2,
    Double_t               minNcls = 70,
    Double_t               maxX2cls = 4.0,
    Double_t               globalX2cls = 36.0,
    Double_t               minCrossedRows = 50.0,
    Double_t               maxClsCrossedRows = 0.8,

    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    Int_t                  signedPdg = 333,

    TString                monitorOpt = "",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
    Bool_t                 useCrossedRows = kFALSE,
    const char*            yaxisVar = "PtDaughter_PDaughter_cent",  //yaxisVar = "PtDaughter_PDaughter_cent"
    Bool_t                 useMixLS = 0
)
{
  TString opt(yaxisVar);
  opt.ToUpper();

  Bool_t isPtDaughter  = opt.Contains("PTDAUGHTER") ;
  Bool_t isPDaughter   = opt.Contains("PDAUGHTER") ;
  Bool_t iscent        = opt.Contains("CENT") ;
  Bool_t iseta         = opt.Contains("ETA") ;
  Bool_t israpidity    = opt.Contains("RAPIDITY") ;



  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetK;


  //vary track quality cuts for systematic checks
  if(enableTrkSyst){
    AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("QualityCut");
    
    trkQualityCut->SetAODTestFilterBit(aodFilterBit);//reset the filter bit cut 
    trkQualityCut->SetCheckOnlyFilterBit(kFALSE);//tells the cut object to check all other cuts individually,
    trkQualityCut->SetDCARPtFormula(DCAxyFormula);
    trkQualityCut->SetDCAZmax(dcazmax);

    if(useCrossedRows) {       
      trkQualityCut->SetMinNCrossedRowsTPC(minCrossedRows, kTRUE);
      trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(maxClsCrossedRows, kTRUE); }
    else trkQualityCut->SetTPCminNClusters(minNcls);

    trkQualityCut->SetTPCmaxChi2(maxX2cls);
    trkQualityCut->SetITSmaxChi2(36);    // In Filter bit 
    trkQualityCut->SetMaxChi2TPCConstrainedGlobal(globalX2cls);  //  In the Filterbit

    trkQualityCut->SetPtRange(0.15, 20.0);
    trkQualityCut->SetEtaRange(-0.8, 0.8);
    
    trkQualityCut->Print();

    if(isPP) {
      cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0);
      cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa), trkQualityCut, cutKaCandidate, AliPID::kKaon, nsigmaKa); 
    } else {
      cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0);
      cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa), trkQualityCut, cutKaCandidate, AliPID::kKaon, nsigmaKa);
      cutSetK->SetUse2011StdQualityCuts(kTRUE);
    }
    
    
  }

  else
    {
      //default cuts 2010 for pp and 2011 for p-Pb
      if(isPP) {
	cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0, aodFilterBit, kFALSE);
	cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaon_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit, kFALSE);
      }
      else {
	cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, aodFilterBit, kFALSE);
	cutSetQ->SetUse2011StdQualityCuts(kTRUE);
	cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaon2011_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit, kFALSE);
	cutSetK->SetUse2011StdQualityCuts(kTRUE);
      }
    }
  
  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutK = task->AddTrackCuts(cutSetK);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput(), monitorOpt.Data());
  }

  Int_t         pdg  = signedPdg;
  TDatabasePDG *db   = TDatabasePDG::Instance();
  TParticlePDG *part = db->GetParticle(pdg);
  Double_t      mass = part->Mass();
  
  
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

  /* transv. momentum */ Int_t ptIDmc  = task->CreateValue(AliRsnMiniValue::kPt, kTRUE);
  /* 1st daughter pt  */ Int_t fdptmc  = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kTRUE);
  /* 2nd daughter pt  */ Int_t sdptmc  = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kTRUE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [12] = { !IsMcTrueOnly, !IsMcTrueOnly, !IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly, isMC    ,  isMC   , isMC   ,  isMC  , useMixLS , useMixLS};
  Bool_t  useIM   [12] = { 1            ,  1           ,  1           ,  1          ,  1          , 1           ,  1      ,   1     , 0      ,   0    , 1        , 1       };
  TString name    [12] = {"UnlikePM"    , "UnlikeMP"   , "MixingPM"   , "MixingMP"  , "LikePP"    , "LikeMM"    ,"TruesPM","TruesMP","ResPM" ,"ResMP" ,"MixingPP","MixingMM"};
  TString comp    [12] = {"PAIR"        , "PAIR"       , "MIX"        , "MIX"       , "PAIR"      , "PAIR"      , "TRUE"  , "TRUE"  , "TRUE" , "TRUE" , "MIX"    ,"MIX"   };
  TString output  [12] = {"SPARSE"      , "SPARSE"     , "SPARSE"     , "SPARSE"    , "SPARSE"    , "SPARSE"    , "SPARSE","SPARSE" ,"SPARSE","SPARSE","SPARSE"  ,"SPARSE"};
  Char_t  charge1 [12] = {'+'           , '-'          , '+'          , '-'         , '+'         , '-'         , '+'     ,  '-'    , '+'    ,  '-'   , '+'      , '-'    };
  Char_t  charge2 [12] = {'-'           , '+'          , '-'          , '+'         , '+'         , '-'         , '-'     ,  '+'    , '-'    ,  '+'   ,'+'       , '-'    };
  Int_t   cutID1  [12] = { iCutK        ,  iCutK       ,  iCutK       ,  iCutK      ,  iCutK      ,  iCutK      ,  iCutK  ,  iCutK  ,  iCutK ,  iCutK , iCutK    , iCutK };
  Int_t   cutID2  [12] = { iCutK        ,  iCutK       ,  iCutK       ,  iCutK      ,  iCutK      ,  iCutK      ,  iCutK  ,  iCutK  ,  iCutK ,  iCutK , iCutK    , iCutK };
  
  //TString output  [10] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };
  
  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("Phi_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(pdg);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 800, 0.8, 1.6);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum of pair as default - else chosen value
    out->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt

    if(isPtDaughter) {
      out->AddAxis(fdpt, 100, 0.0, 10.0);
      out->AddAxis(sdpt, 100, 0.0, 10.0);  
    }
    if(isPDaughter) {
      out->AddAxis(fdp, 100, 0.0, 10.0);
      out->AddAxis(sdp, 100, 0.0, 10.0);  
    }
      
    // axis Z: centrality-multiplicity
    if(iscent) {
      if (!isPP)	out->AddAxis(centID, 100, 0.0, 100.0);
      else       	out->AddAxis(centID, 400, 0.0, 400.0);      
    }

    // axis W: pseudorapidity
    if(iseta)       out->AddAxis(etaID, 20, -1.0, 1.0);
    
    // axis J: rapidity
    if(israpidity)  out->AddAxis(yID, 12, -0.6, 0.6);
    
  }   
 
  if (isMC){   
    // create output
    AliRsnMiniOutput *outm = task->CreateOutput(Form("Phi_Mother%s", suffix), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kKaon);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(pdg);
    outm->SetMotherMass(mass);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 800, 0.8, 1.6);
    // axis Y: transverse momentum of pair as default - else chosen value
    outm->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt

    if(isPtDaughter) {
      outm->AddAxis(fdpt, 100, 0.0, 10.0);
      outm->AddAxis(sdpt, 100, 0.0, 10.0);
    }
    
    if(isPDaughter) {
	outm->AddAxis(fdp, 100, 0.0, 10.0);
	outm->AddAxis(sdp, 100, 0.0, 10.0);
    }

    if(iscent) {
      if (!isPP)	outm->AddAxis(centID, 100, 0.0, 100.0);
      else       	outm->AddAxis(centID, 400, 0.0, 400.0);
    }
    // axis W: pseudorapidity

    if(iseta)       outm->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    if(israpidity)  outm->AddAxis(yID, 12, -0.6, 0.6);

    //create plot for generated Pt of mother vs generated Pt of daughters
    AliRsnMiniOutput *outPtGen = task->CreateOutput(Form("Phi_Mother_GenPt_%s", suffix), "SPARSE", "MOTHER");
    outPtGen->SetDaughter(0, AliRsnDaughter::kKaon);
    outPtGen->SetDaughter(1, AliRsnDaughter::kKaon);
    outPtGen->SetMotherPDG(pdg);
    outPtGen->SetMotherMass(mass);
    // pair cuts
    outPtGen->SetPairCuts(cutsPair);
    // binnings
    outPtGen->AddAxis(ptIDmc, 100, 0.0, 10.0); //mother pt - generated
    outPtGen->AddAxis(fdptmc, 100, 0.0, 10.0); //first daughter pt - generated
    outPtGen->AddAxis(sdptmc, 100, 0.0, 10.0); //second daughter pt - generated

    //create plot for reconstructed Pt of true mother vs generated Pt of daughters
    AliRsnMiniOutput *outPtTrueGen = task->CreateOutput(Form("Phi_True_GenPt_%s", suffix), "SPARSE", "TRUE");
    outPtTrueGen->SetDaughter(0, AliRsnDaughter::kKaon);
    outPtTrueGen->SetDaughter(1, AliRsnDaughter::kKaon);
    outPtTrueGen->SetCutID(0, iCutK);
    outPtTrueGen->SetCutID(1, iCutK);
    outPtTrueGen->SetCharge(0, '+');
    outPtTrueGen->SetCharge(1, '-');    
    outPtTrueGen->SetMotherPDG(pdg);
    outPtTrueGen->SetMotherMass(mass);
    // pair cuts
    outPtTrueGen->SetPairCuts(cutsPair);
    // binnings
    outPtTrueGen->AddAxis(ptID, 100, 0.0, 10.0); //true pt - generated
    outPtTrueGen->AddAxis(fdptmc, 100, 0.0, 10.0); //first daughter pt - generated
    outPtTrueGen->AddAxis(sdptmc, 100, 0.0, 10.0); //second daughter pt - generated

   //create plot for reconstructed Pt of true mother vs reconstructed Pt of daughters
    AliRsnMiniOutput *outPtTrueRec = task->CreateOutput(Form("Phi_True_RecPt_%s", suffix), "SPARSE", "TRUE");
    outPtTrueRec->SetDaughter(0, AliRsnDaughter::kKaon);
    outPtTrueRec->SetDaughter(1, AliRsnDaughter::kKaon);
    outPtTrueRec->SetCutID(0, iCutK);
    outPtTrueRec->SetCutID(1, iCutK);
    outPtTrueRec->SetCharge(0, '+');
    outPtTrueRec->SetCharge(1, '-');    
    outPtTrueRec->SetMotherPDG(pdg);
    outPtTrueRec->SetMotherMass(mass);
    // pair cuts
    outPtTrueRec->SetPairCuts(cutsPair);
    // binnings
    outPtTrueRec->AddAxis(ptID, 100, 0.0, 10.0); //mother pt - reconstructed
    outPtTrueRec->AddAxis(fdpt, 100, 0.0, 10.0); //first daughter pt - reconstructed
    outPtTrueRec->AddAxis(sdpt, 100, 0.0, 10.0); //second daughter pt - reconstructed
  }
  return kTRUE;
}
