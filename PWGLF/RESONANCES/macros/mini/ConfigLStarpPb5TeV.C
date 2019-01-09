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

Bool_t ConfigLStarpPb5TeV
(  
 
 AliRsnMiniAnalysisTask *task, 
 Bool_t                 isMC, 
 Bool_t                 isPP,
 const char             *suffixy,
 AliRsnCutSet           *cutsPair,
 Int_t                  aodFilterBit = 5,
 Int_t                  customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar,
 Float_t                nsigmaPr = 3.0,
 Float_t                nsigmaKa = 3.0,
 Bool_t                 enableMonitor = kTRUE,
 Bool_t                 IsMcTrueOnly = kTRUE,
 Int_t                  signedPdg = 3124,
 
 TString                monitorOpt = "",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
 Bool_t                 useCrossedRows = kTRUE,
 const char*            yaxisVar = "PtDaughter_PDaughter_cent",  //yaxisVar = "PtDaughter_PDaughter_cent"
 //AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt,
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


  
  if(customQualityCutsID){
    AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("QualityCut");
    trkQualityCut->Print();
    if(isPP) {
      Bool_t useCrossedRows = 1;
      cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0, aodFilterBit, useCrossedRows);
      cutSetP  = new AliRsnCutSetDaughterParticle(Form("cutProton_%2.1fsigma",nsigmaPr), cutPrCandidate, AliPID::kProton, nsigmaPr, aodFilterBit, useCrossedRows);
      cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaon_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit, useCrossedRows);
    }
    else {
      Bool_t useCrossedRows = 1;
      cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, aodFilterBit, useCrossedRows);
      cutSetQ->SetUse2011StdQualityCuts(kTRUE);
      cutSetP = new AliRsnCutSetDaughterParticle(Form("cutProton2011_%2.1fsigma",nsigmaPr), cutPrCandidate, AliPID::kProton, nsigmaPr, aodFilterBit, useCrossedRows);
      cutSetP->SetUse2011StdQualityCuts(kTRUE);
      cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaon2011_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit, useCrossedRows);
      cutSetK->SetUse2011StdQualityCuts(kTRUE);
    }
  }
  else{
    
    Bool_t useCrossedRows = 1;
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit),AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.,aodFilterBit,useCrossedRows);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate,nsigmaKa),cutKaCandidate,AliPID::kKaon,nsigmaKa,aodFilterBit,useCrossedRows);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutPr%i_%2.1fsigma",cutPrCandidate,nsigmaPr),cutKaCandidate,AliPID::kProton,nsigmaPr,aodFilterBit,useCrossedRows);
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
  // return kTRUE;
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID    = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID   = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID    = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID  = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID   = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID     = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt  */ Int_t fdpt    = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt    = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p   */ Int_t fdp     = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p   */ Int_t sdp     = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);
  /* IM diff          */ Int_t IMdif   = task->CreateValue(AliRsnMiniValue::kInvMassDiff, kTRUE);
  //inverse cosine of the angle between daughter vector momenta

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
      out->AddAxis(imID, 1000, 1.4, 1.8);
    else
      out->AddAxis(resID, 500, -0.02, 0.02);
    
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
      out->AddAxis(centID, 20, 0.0, 100.0);
    else 
      out->AddAxis(centID, 400, 0.0, 400.0);
    



    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    // out->AddAxis(yID, 10, -0.5, 0.5);
    
  }   


  
  if (isMC){   


    AliRsnMiniOutput *outnew2 = task->CreateOutput(Form("Lstar_resolution2%s", opt.Data()), "SPARSE", "TRUE");
  
    outnew2->SetCutID(0, cutID1[0]);
    outnew2->SetCutID(1, cutID2[0]);
    outnew2->SetDaughter(0, AliRsnDaughter::kProton);
    outnew2->SetDaughter(1, AliRsnDaughter::kKaon);
    outnew2->SetCharge(0, charge1[0]);
    outnew2->SetCharge(1, charge2[0]);
    outnew2->SetMotherPDG(signedPdg);
    outnew2->SetMotherMass(1.51953);
    outnew2->SetPairCuts(cutsPair);

    outnew2->AddAxis(imID, 1000, 1.4, 1.8); // mass axis
    outnew2->AddAxis(ptID, 100, 0.0, 10.0); //pt axis
    //outnew2->AddAxis(resmass, 400, -0.2, 0.2);// mass resolution
    outnew2->AddAxis(IMdif, 1000, -0.5, 0.5);// mass diff (true -reco)
    /*
    AliRsnMiniOutput *outnew1 = task->CreateOutput(Form("Lstar_resolution%s%s", opt.Data()), "SPARSE", "TRUE");

    outnew1->SetCutID(0, cutID1[0]);
    outnew1->SetCutID(1, cutID2[0]);
    outnew1->SetDaughter(0, AliRsnDaughter::kProton);
    outnew1->SetDaughter(1, AliRsnDaughter::kKaon);

    outnew1->SetMotherPDG(signedPdg);
    outnew1->SetMotherMass(1.51953);
    outnew1->SetPairCuts(cutsPair);

    outnew1->AddAxis(imID, 800, 1.4, 2.2); // mass axis                                                                                 
    outnew1->AddAxis(ptID, 100, 0.0, 10.0); //pt axis                                                                                   
    outnew1->AddAxis(IMdif, 400, -0.5, 0.5);// mass diff (true -reco)
    */
    // create output
    AliRsnMiniOutput *outm = task->CreateOutput(Form("Lstar_Mother%s", opt.Data()), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kProton);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(signedPdg);
    outm->SetMotherMass(1.51953);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 1000, 1.4, 1.8);

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
      outm->AddAxis(centID, 20, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }


  }
  return kTRUE;
}
