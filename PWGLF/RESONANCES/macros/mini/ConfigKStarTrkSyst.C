/***************************************************************************
              fbellini@cern.ch - created on 19/09/2013

// *** Configuration script for K*, anti-K* syst. analysis with p-Pb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigKStarTrkSyst
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010,
    Float_t                nsigmaPi = 2.0,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    Bool_t                 useMixLS = 0,
    Int_t                  signedPdg = 313,
    TString                monitorOpt = "",
    AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt,
    Bool_t                 enableTrkSyst = kFALSE,
    Double_t               dcaxymax = 2.4,
    Double_t               dcazmax = 3.2,
    Double_t               minNcls = 70,
    Double_t               maxX2cls = 4.0
    // Double_t               minCrossedRows = 50.0,
    // Double_t               maxClsCrossedRows = 0.8
)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi), cutPiCandidate, AliPID::kPion, nsigmaPi, aodFilterBit);
  AliRsnCutSetDaughterParticle * cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate, nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit);
  AliRsnCutSetDaughterParticle * cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cut_quality"), AliRsnCutSetDaughterParticle::kQualityStd2011, 1e20, 1e20, aodFilterBit);
  
  if (enableTrkSyst) {
    ((AliRsnCutTrackQuality *)cutSetQ->GetQualityCut())->SetDCARmax(dcaxymax);
    ((AliRsnCutTrackQuality *)cutSetQ->GetQualityCut())->SetDCAZmax(dcazmax);
    ((AliRsnCutTrackQuality *)cutSetQ->GetQualityCut())->SetTPCminNClusters(minNcls);
    ((AliRsnCutTrackQuality *)cutSetQ->GetQualityCut())->SetTPCmaxChi2(maxX2cls);
    //cutSetQ->SetTPCminNcrossedRows(minCrossedRows); //TO BE IMPLEMENTED
    //cutSetQ->SetTPCminNclsPerCrossedRow(maxClsCrossedRows); //TO BE IMPLEMENTED
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

  Bool_t  use     [12] = { !IsMcTrueOnly,  !IsMcTrueOnly,  !IsMcTrueOnly,  !IsMcTrueOnly ,  !IsMcTrueOnly, !IsMcTrueOnly,  isMC   ,   isMC   ,  isMC   ,   isMC , useMixLS, useMixLS  };
  Bool_t  useIM   [12] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  1      ,   1      ,  0      ,   0 , 1    , 1     };
  TString name    [12] = {"UnlikePM", "UnlikeMP", "MixingPM", "MixingMP", "LikePP", "LikeMM", "TruesPM",  "TruesMP", "ResPM"  ,  "ResMP",  "MixingPP",  "MixingMM"  };
  TString comp    [12] = {"PAIR"   , "PAIR"   , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  ,  "TRUE"  , "TRUE"  ,  "TRUE", "MIX","MIX"};
  //TString output  [10] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };
  TString output  [12] = {"SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"  , "SPARSE"  , "SPARSE"  ,  "SPARSE"  , "SPARSE"  ,  "SPARSE", "SPARSE"  ,  "SPARSE"};
  Char_t  charge1 [12] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-'  , '+' , '-'};
  Char_t  charge2 [12] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+'  ,'+' , '-'   };
  Int_t   cutID1  [12] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK , iCutK, iCutK };
  Int_t   cutID2  [12] = { iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi ,  iCutPi ,  iCutPi ,   iCutPi ,  iCutPi ,   iCutPi, iCutPi, iCutPi };
  
  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("kstar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(signedPdg);
    out->SetMotherMass(0.89594);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 90, 0.6, 1.5);
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
    AliRsnMiniOutput *outm = task->CreateOutput(Form("kstar_Mother%s", suffix), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kKaon);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(signedPdg);
    outm->SetMotherMass(0.89594);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 90, 0.6, 1.5);
    outm->AddAxis(ptID, 100, 0.0, 10.0);
    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
  }
  return kTRUE;
}
