AliRsnPairManager* RsnConfig
(
  AliRsnPair::EPairType  pidType,        // PID type (NoPID, RealisticPID or PerfectPID)
  const char            *pairMgrName,    // name for the pair manager
  Int_t                  resonancePDG,   // PDG code of resonance (for true pairs)
  AliPID::EParticleType  type1,          // PID of one member of decay (+)
  AliPID::EParticleType  type2,          // PID of other member of decay (-)
  Double_t               bbCut = 10000.0 // Bethe-Bloch TPC cut value (large --> excluded)
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

  // === NAME DEFINITIONS =========================================================================

  AliRsnPairManager  *pairMgr  = new AliRsnPairManager(pairMgrName);

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
  // [2] - like PP
  // [3] - like MM
  AliRsnPair *pair[2][4] = {0,0,0,0,0,0,0,0};

  for (i = 0; i < nArray; i++) {
    pair[i][0] = new AliRsnPair(pidType, defUnlike[i]);
    pair[i][1] = new AliRsnPair(pidType, defUnlike[i]);
    pair[i][2] = new AliRsnPair(pidType, defLikePP[i]);
    pair[i][3] = new AliRsnPair(pidType, defLikeMM[i]);
  }

  // === CUTS =====================================================================================

  // cuts for tracks:
  // -- primary track quality
  AliRsnCutESDPrimary *cutESDPrimary = new AliRsnCutESDPrimary("cutESDPrimary");
  cutESDPrimary->GetCuts()->SetMaxCovDiagonalElements(2.0, 2.0, 0.5, 0.5, 2.0);
  cutESDPrimary->GetCuts()->SetRequireSigmaToVertex(kTRUE);
  cutESDPrimary->GetCuts()->SetMaxNsigmaToVertex(4.0);
  cutESDPrimary->GetCuts()->SetRequireTPCRefit(kTRUE);
  cutESDPrimary->GetCuts()->SetAcceptKinkDaughters(kFALSE);
  cutESDPrimary->GetCuts()->SetMinNClustersTPC(50);
  cutESDPrimary->GetCuts()->SetMaxChi2PerClusterTPC(3.5);
  // -- Bethe-Bloch with kaon mass hypothesis
  AliRsnCutBetheBloch *cutKaonBB = new AliRsnCutBetheBloch("cutKaonBB", bbCut, AliPID::kKaon);
  cutKaonBB->SetCalibConstant(0, 0.76176e-1);
  cutKaonBB->SetCalibConstant(1, 10.632);
  cutKaonBB->SetCalibConstant(2, 0.13279e-4);
  cutKaonBB->SetCalibConstant(3, 1.8631);
  cutKaonBB->SetCalibConstant(4, 1.9479);

  // cuts on pairs:
  // -- true daughters of a phi resonance (only for true pairs histogram)cutSetPairTrue->AddCut(cutTrue);
  AliRsnCutStd *cutTruePair = new AliRsnCutStd("cutTrue", AliRsnCutStd::kTruePair, resonancePDG);

  // cuts on event:
  // -- none (specified for whole task)

  // cut set definition for all pairs
  AliRsnCutSet *cutSetParticle = new AliRsnCutSet("trackCuts");
  cutSetParticle->AddCut(cutESDPrimary);
  if (pidType == AliRsnPair::kNoPID && bbCut < 10.0)
  {
    Info("RsnConfig", "PID TYPE = NOPID AND BB CUT = %f -- Adding this cut", bbCut);
    cutSetParticle->AddCut(cutKaonBB);
    cutSetParticle->SetCutScheme("cutKaonBB&cutESDPrimary");
  }
  else {
    Info("RsnConfig", "PID TYPE = %d AND BB CUT = %f -- Excluding cut", pidType, bbCut);
    cutSetParticle->SetCutScheme("cutESDPrimary");
  }

  // cut set definition for true pairs
  AliRsnCutSet *cutSetPairTrue = new AliRsnCutSet("truePairs");
  cutSetPairTrue->AddCut(cutTruePair);
  cutSetPairTrue->SetCutScheme("cutTrue");

  // cut manager for all pairs
  AliRsnCutMgr *cutMgrAll = new AliRsnCutMgr("std", "All");
  cutMgrAll->SetCutSet(AliRsnCut::kParticle, cutSetParticle);

  // cut manager for all pairs
  AliRsnCutMgr *cutMgrTrue = new AliRsnCutMgr("true", "True");
  cutMgrTrue->SetCutSet(AliRsnCut::kParticle, cutSetParticle);
  cutMgrTrue->SetCutSet(AliRsnCut::kPair, cutSetPairTrue);

  for (i = 0; i < nArray; i++) {
    pair[i][0]->SetCutMgr(cutMgrTrue);
    pair[i][0]->SetOnlyTrue();
    for (j = 1; j < 4; j++) {
      pair[i][j]->SetCutMgr(cutMgrAll);
    }
  }

  // === FUNCTIONS ================================================================================

  // define all binnings
  AliRsnFunctionAxis *axisIM   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairInvMass,    1000,  0.0,   2.0);
  AliRsnFunctionAxis *axisPt   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairPt,          400,  0.0,  10.0);
  AliRsnFunctionAxis *axisEta  = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairEta,          10, -1.0,   1.0);
  AliRsnFunctionAxis *axisMult = new AliRsnFunctionAxis(AliRsnFunctionAxis::kEventMult,         8,  0.0, 200.0);

  // function #1: pt, mult
  AliRsnFunction *fcnPt = new AliRsnFunction;
  fcnPt->AddAxis(axisIM);
  fcnPt->AddAxis(axisPt);
  fcnPt->AddAxis(axisMult);
  // function #2: eta, mult
  AliRsnFunction *fcnEta = new AliRsnFunction;
  fcnEta->AddAxis(axisIM);
  fcnEta->AddAxis(axisEta);
  fcnEta->AddAxis(axisMult);

  // add functions to pairs
  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 4; j++) {
      pair[i][j]->AddFunction(fcnPt);
      pair[i][j]->AddFunction(fcnEta);
    }
  }

  // === ADD TO PAIR MANAGER ======================================================================

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 4; j++) {
      pairMgr->AddPair(pair[i][j]);
    }
  }

  return pairMgr;
}
