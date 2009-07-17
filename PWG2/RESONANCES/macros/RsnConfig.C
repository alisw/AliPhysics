AliRsnPairManager *RsnConfig_PHI(const char *name="PHI")
{
  return   CreatePairs(name, 333, AliAODTrack::kKaon, AliAODTrack::kKaon, 0, 10000);
}

AliRsnPairManager *RsnConfig_KSTAR(const char *name="KSTAR")
{
  return CreatePairs(name, 313, AliAODTrack::kPion, AliAODTrack::kKaon, 0, 10000);
}


AliRsnPairManager* CreatePairs
(
  const char       *pairMgrName,    // name for the pair manager
  Int_t             resonancePDG,   // PDG code of resonance (for true pairs)
  AliAODTrack::AODTrkPID_t type1,          // PID of one member of decay (+)
  AliAODTrack::AODTrkPID_t  type2,          // PID of other member of decay (-)
  Int_t             multMin,        // lower edge of multiplicity cut
  Int_t             multMax         // upper edge of multiplicity cut
)
{
//
// Creates an AliRsnPairMgr for a specified resonance, which contains:
// - signal (inv. mass)
// - event mixing (inv. mass)
// - like-sign (inv. mass)
// - true pairs (inv. mass, resolution)
//
// For all pairs, a binning in Pt and Eta is provided, and a cut in multiplicity
// which defines a multiplicity bin where the analysis is computed.
//
// Arguments define how the pair manager must be built, and are explained above
//

  AliRsnPairManager  *pairMgr  = new AliRsnPairManager(pairMgrName);
  //cout << "Creating " << pairMgrName << endl;

//   return pairMgr;

  // === PAIR DEFINITIONS =========================================================================

  // if particle #1 and #2 are different, two histograms must be built
  // for each scheme (signal, true, mixed, like-sign) exchanging both particles and signs
  Int_t i, j, nArray = 1;
  if (type1 != type2) nArray = 2;

  AliRsnPairDef *defUnlike[2] = {0, 0};
  AliRsnPairDef *defLikePP[2] = {0, 0};
  AliRsnPairDef *defLikeMM[2] = {0, 0};

  defUnlike[0] = new AliRsnPairDef(type1, '+', type2, '-', resonancePDG);
  defLikePP[0] = new AliRsnPairDef(type1, '+', type2, '+', resonancePDG);
  defLikeMM[0] = new AliRsnPairDef(type1, '-', type2, '-', resonancePDG);

  defUnlike[1] = new AliRsnPairDef(type2, '+', type1, '-', resonancePDG);
  defLikePP[1] = new AliRsnPairDef(type2, '+', type1, '+', resonancePDG);
  defLikeMM[1] = new AliRsnPairDef(type2, '-', type1, '-', resonancePDG);

  // === PAIR ANALYSIS ENGINES ====================================================================

  // define null (dummy) objects and initialize only the ones which are needed,
  // depending again on particle types;
  // array is organized as follows:
  // [0] - true pairs
  // [1] - signal
  // [2] - mixing
  // [3] - like PP
  // [4] - like MM
  AliRsnPair *noPIDnoCut[2][5]   = {0,0,0,0,0,0,0,0,0,0};
  AliRsnPair *noPIDwithCut[2][5] = {0,0,0,0,0,0,0,0,0,0};
  AliRsnPair *realisticPID[2][5] = {0,0,0,0,0,0,0,0,0,0};
  AliRsnPair *perfectPID[2][5]   = {0,0,0,0,0,0,0,0,0,0};

  for (i = 0; i < nArray; i++) {
    noPIDnoCut[i][0] = new AliRsnPair(AliRsnPair::kNoPID, defUnlike[i]);
    noPIDnoCut[i][1] = new AliRsnPair(AliRsnPair::kNoPID, defUnlike[i]);
    noPIDnoCut[i][2] = new AliRsnPair(AliRsnPair::kNoPIDMix, defUnlike[i]);
    noPIDnoCut[i][3] = new AliRsnPair(AliRsnPair::kNoPID, defLikePP[i]);
    noPIDnoCut[i][4] = new AliRsnPair(AliRsnPair::kNoPID, defLikeMM[i]);

    noPIDwithCut[i][0] = new AliRsnPair(AliRsnPair::kNoPID, defUnlike[i]);
    noPIDwithCut[i][1] = new AliRsnPair(AliRsnPair::kNoPID, defUnlike[i]);
    noPIDwithCut[i][2] = new AliRsnPair(AliRsnPair::kNoPIDMix, defUnlike[i]);
    noPIDwithCut[i][3] = new AliRsnPair(AliRsnPair::kNoPID, defLikePP[i]);
    noPIDwithCut[i][4] = new AliRsnPair(AliRsnPair::kNoPID, defLikeMM[i]);

    realisticPID[i][0] = new AliRsnPair(AliRsnPair::kRealisticPID, defUnlike[i]);
    realisticPID[i][1] = new AliRsnPair(AliRsnPair::kRealisticPID, defUnlike[i]);
    realisticPID[i][2] = new AliRsnPair(AliRsnPair::kRealisticPIDMix, defUnlike[i]);
    realisticPID[i][3] = new AliRsnPair(AliRsnPair::kRealisticPID, defLikePP[i]);
    realisticPID[i][4] = new AliRsnPair(AliRsnPair::kRealisticPID, defLikeMM[i]);

    perfectPID[i][0] = new AliRsnPair(AliRsnPair::kPerfectPID, defUnlike[i]);
    perfectPID[i][1] = new AliRsnPair(AliRsnPair::kPerfectPID, defUnlike[i]);
    perfectPID[i][2] = new AliRsnPair(AliRsnPair::kPerfectPIDMix, defUnlike[i]);
    perfectPID[i][3] = new AliRsnPair(AliRsnPair::kPerfectPID, defLikePP[i]);
    perfectPID[i][4] = new AliRsnPair(AliRsnPair::kPerfectPID, defLikeMM[i]);
  }

  // === CUTS =====================================================================================

  /*
  // cuts for tracks:
  // - probability to be a kaon (only for kaons)
  AliRsnCutStd *cutAssignedKaon = new AliRsnCut("cutAssignedKaon", "", AliRsnCut::kAssignedPID, AliAODTrack::kKaon);
  AliRsnCutStd *cutProbKaon = new AliRsnCut("cutProbKaon", "", AliRsnCut::kPIDProbForSpecies, (Int_t)AliAODTrack::kKaon);
  cutProbKaon->SetCutValues(AliRsnCut::kPIDProbForSpecies, 0.15, 1.0);
  */
  AliRsnCutSet *cutSetTrack = new AliRsnCutSet("tracks");
  //cutSetTrack->AddCut(cutAssignedKaon);
  //cutSetTrack->AddCut(cutProbKaon);
  //cutSetTrack->SetCutScheme("(!cutAssignedKaon)|(cutAssignedKaon&cutProbKaon)");

  // cuts on pairs:
  // - true daughters of the defined resonance (only for true pairs histogram)
  AliRsnCutStd *cutPairTrue    = new AliRsnCutStd("cutTrue", "", AliRsnCut::kTruePair, resonancePDG);
  AliRsnCutSet *cutSetPairTrue = new AliRsnCutSet("truePairs");
  cutSetPairTrue->AddCut(cutPairTrue);
  cutSetPairTrue->SetCutScheme("cutTrue");

  // cuts on events:
  // - multiplicity bin
  AliRsnCutStd *cutEventMult = new AliRsnCutStd("cutMult", "", AliRsnCut::kMultiplicity, multMin, multMax);
  AliRsnCutSet *cutSetEvent  = new AliRsnCutSet("multiplicity");
  cutSetEvent->AddCut(cutEventMult);
  cutSetEvent->SetCutScheme("cutMult");

  // define cut manager for NOT true pairs
  AliRsnCutMgr *cutMgr = new AliRsnCutMgr("default", "All pairs");
  cutMgr->SetCutSet(AliRsnCut::kEvent, cutSetEvent);

  // define cut manager for true pairs
  AliRsnCutMgr *cutMgrTrue = new AliRsnCutMgr("true", "True pairs");
  cutMgrTrue->SetCutSet(AliRsnCut::kEvent, cutSetEvent);
  cutMgrTrue->SetCutSet(AliRsnCut::kPair, cutSetPairTrue);

  // define cut manager for NOPID with kaon prob cut
  // define cut manager for NOT true pairs
  AliRsnCutMgr *cutMgrProb = new AliRsnCutMgr("probK", "All pairs with kaon probability cut");
  cutMgrProb->SetCutSet(AliRsnCut::kEvent, cutSetEvent);
  cutMgrProb->SetCutSet(AliRsnCut::kParticle, cutSetTrack);

  // define cut manager for true pairs
  AliRsnCutMgr *cutMgrTrueProb = new AliRsnCutMgr("true+probK", "True pairs with kaon probability cut");
  cutMgrTrueProb->SetCutSet(AliRsnCut::kEvent, cutSetEvent);
  cutMgrTrueProb->SetCutSet(AliRsnCut::kPair, cutSetPairTrue);
  cutMgrTrueProb->SetCutSet(AliRsnCut::kParticle, cutSetTrack);

  // add cuts to pair analysis
  for (i = 0; i < nArray; i++) {
    noPIDnoCut[i][0]->SetCutMgr(cutMgrTrue);
    noPIDwithCut[i][0]->SetCutMgr(cutMgrTrueProb);
    realisticPID[i][0]->SetCutMgr(cutMgrTrue);
    perfectPID[i][0]->SetCutMgr(cutMgrTrue);
    for (j = 1; j < 5; j++) {
      noPIDnoCut[i][j]->SetCutMgr(cutMgr);
      noPIDwithCut[i][j]->SetCutMgr(cutMgrProb);
      realisticPID[i][j]->SetCutMgr(cutMgr);
      perfectPID[i][j]->SetCutMgr(cutMgr);
    }
  }

  // === FUNCTIONS ================================================================================

  // define histogram templates
  AliRsnFunctionAxis *axisIM   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairInvMass,    1000,  0.0,   2.0);
  AliRsnFunctionAxis *axisPt   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairPt,           10,  0.0,  10.0);
  AliRsnFunctionAxis *axisEta  = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairEta,          20, -1.0,   1.0);
  AliRsnFunctionAxis *axisMult = new AliRsnFunctionAxis(AliRsnFunctionAxis::kEventMult,         8,  0.0, 200.0);

  // define functions axes
  AliRsnFunction *fcnIM = new AliRsnFunction;
  fcnIM->AddAxis(axisIM);
  fcnIM->AddAxis(axisPt);
  fcnIM->AddAxis(axisEta);
  fcnIM->AddAxis(axisMult);

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 5; j++) {
      noPIDnoCut[i][j]->AddFunction(fcnIM);
      noPIDwithCut[i][j]->AddFunction(fcnIM);
      realisticPID[i][j]->AddFunction(fcnIM);
      perfectPID[i][j]->AddFunction(fcnIM);
    }
  }

  // === ADD TO PAIR MANAGER ======================================================================

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 5; j++) {
      pairMgr->AddPair(noPIDnoCut[i][j]);
      pairMgr->AddPair(noPIDwithCut[i][j]);
      pairMgr->AddPair(realisticPID[i][j]);
      pairMgr->AddPair(perfectPID[i][j]);
    }
  }

  return pairMgr;
}
