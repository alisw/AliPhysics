/***************************************************************************
              subhash.singha@cern.ch - last modified on 20/01/2014

// *** Configuration script for K*, anti-K* analysis with 2010 PbPb runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
****************************************************************************/

Bool_t ConfigTPCanalysisKStarSyst
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  aodFilterBit = 5,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
    Float_t                nsigmaPi = 2.0,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableSyst = kFALSE,
    Char_t                 DCAxyFormula[100] = "0.0182+0.035/pt^1.01",
    Double_t               dcazmax = 3.2,
    Double_t               minNcls = 70,
    Double_t               maxX2cls = 5.0,
    Double_t               minCrossedRows = 50.0,
    Double_t               maxClsCrossedRows = 0.8,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    Int_t                  Pdg = 313,
    Int_t                  aodN = 0
)
{
  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetPi;
  AliRsnCutSetDaughterParticle * cutSetK;

  //vary track quality cuts for systematic checks
  if(enableSyst){
    AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("QualityCut");
    
    trkQualityCut->DisableAll();//disable all cuts, filter bit, pT, eta, and DCAxy cuts will be reset later
    trkQualityCut->SetAODTestFilterBit(aodFilterBit);//reset the filter bit cut 
    trkQualityCut->SetCheckOnlyFilterBit(kFALSE);//tells the cut object to check all other cuts individually,
    trkQualityCut->SetDCARPtFormula(DCAxyFormula);
    trkQualityCut->SetDCAZmax(dcazmax);
    trkQualityCut->SetMinNCrossedRowsTPC(minCrossedRows, kTRUE);
    trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(maxClsCrossedRows, kTRUE);
    trkQualityCut->SetTPCmaxChi2(maxX2cls);
    trkQualityCut->SetRejectKinkDaughters(kTRUE);
    trkQualityCut->SetSPDminNClusters(AliESDtrackCuts::kAny);
    trkQualityCut->SetITSmaxChi2(36);
    trkQualityCut->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);//already in defaults 2011
    trkQualityCut->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);//already in defaults 2011
    trkQualityCut->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);//already in defaults 2011
    //trkQualityCut->SetTPCminNClusters(70);
    trkQualityCut->SetPtRange(0.15, 20.0);
    trkQualityCut->SetEtaRange(-0.8, 0.8);
    
    trkQualityCut->Print();

    cutSetQ  = new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit), trkQualityCut, AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0);
    cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",cutPiCandidate, nsigmaPi), trkQualityCut, cutPiCandidate, AliPID::kPion, nsigmaPi);
    cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutPiCandidate, nsigmaKa), trkQualityCut, cutKaCandidate, AliPID::kKaon, nsigmaKa);

    
  }
  else
    {
      //default cuts 2011
      cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, aodFilterBit);
      cutSetQ->SetUse2011StdQualityCuts(kTRUE);
      cutSetPi = new AliRsnCutSetDaughterParticle(Form("cutPionTPCPbPb2011_%2.1fsigma",nsigmaPi), cutPiCandidate, AliPID::kPion, nsigmaPi, aodFilterBit);
      cutSetPi->SetUse2011StdQualityCuts(kTRUE);
      cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaonTPCPbPb2011_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit);
      cutSetK->SetUse2011StdQualityCuts(kTRUE);
    }
  
  
  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutK = task->AddTrackCuts(cutSetK);
  
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput());
  }  
  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --
  Bool_t  use     [12] = { !IsMcTrueOnly,  !IsMcTrueOnly,  !IsMcTrueOnly,  !IsMcTrueOnly ,  !IsMcTrueOnly, !IsMcTrueOnly,  isMC   ,   isMC   ,  isMC   ,   isMC,  !IsMcTrueOnly,  !IsMcTrueOnly   };
  Bool_t  useIM   [12] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  1      ,   1      ,  0      ,   0,    1       ,  1       };
  TString name    [12] = {"UnlikePM", "UnlikeMP", "MixingPM", "MixingMP", "LikePP", "LikeMM", "TruesPM",  "TruesMP", "ResPM"  ,  "ResMP",  "RotatePM",   "RotateMP"  };
  TString comp    [12] = {"PAIR"   , "PAIR"   , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  ,  "TRUE"  , "TRUE"  ,  "TRUE", "ROTATE2",   "ROTATE2"  };
  //TString output  [12] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };
  TString output  [12] = {"SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"  , "SPARSE"  , "SPARSE"  ,  "SPARSE"  , "SPARSE"  ,  "SPARSE" , "SPARSE"   , "SPARSE" };
  Char_t  charge1 [12] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-' ,   '+'      , '-'    };
  Char_t  charge2 [12] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+' , '-'      , '+'      };
  Int_t   cutID1  [12] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK, iCutK   ,  iCutK   };
  Int_t   cutID2  [12] = { iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi ,  iCutPi ,  iCutPi ,   iCutPi ,  iCutPi ,   iCutPi, iCutPi  ,  iCutPi };
  
  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    if(Pdg > 0) AliRsnMiniOutput *out = task->CreateOutput(Form("kstar1_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    if(Pdg < 0) AliRsnMiniOutput *out = task->CreateOutput(Form("kstar2_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(Pdg);//313
    out->SetMotherMass(0.89594);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 90, 0.6, 1.5);
    //else
    //out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID, 300, 0.0, 30.0);
    
    // axis Z: centrality-multiplicity
    if (!isPP)
      out->AddAxis(centID, 100, 0.0, 100.0);
    else 
      out->AddAxis(centID, 400, 0.0, 400.0);
    
    // axis W: pseudorapidity
    //    out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    //out->AddAxis(yID, 32, -0.8, 0.8);
    
  }
  
  if (isMC){
    // create output
    
    if(Pdg > 0) {AliRsnMiniOutput *outm = task->CreateOutput(Form("kstar_Mother1%s", suffix), "SPARSE", "MOTHER");}
    if(Pdg < 0) {AliRsnMiniOutput *outm = task->CreateOutput(Form("kstar_Mother2%s", suffix), "SPARSE", "MOTHER");}
    outm->SetDaughter(0, AliRsnDaughter::kKaon);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(Pdg);//313
    outm->SetMotherMass(0.89594);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 90, 0.6, 1.5);
    outm->AddAxis(ptID, 300, 0.0, 30.0);
    if (!isPP){
    outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else { 
     outm->AddAxis(centID, 400, 0.0, 400.0);
    }
    //outm->AddAxis(yID, 32, -0.8, 0.8);
  }
  return kTRUE;
}

