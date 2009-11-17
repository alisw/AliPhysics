//
// This is an example macro for creation of a pair manager
// for a resonance analysis with the PWG2/RESONANCES package.
// Its output is an AliRsnPairManager which will be added to the task.
// It will contain:
// - signal (inv. mass)
// - true pairs (inv. mass)
// - like-sign (inv. mass)
// - true pairs (inv. mass, resolution)
// All histogram are done w.r. to Pt or Eta and w.r. to multiplicity.
// The binnings are hard-coded. Change them according to your preferences.
//
// Arguments:
// -- pairMgrName : a name for the pair manager which must contain
//                  some keywords to choose amond the available cuts
//                  and PID selection available here (this is a personal customization)
// -- resonancePDG: PDG code of resonance (for true pairs)
// -- type1, type2: particle species for track 1 and 2 in the pair (using AliPID enumeration)
//
// NOTE:
// the keyword available here are the following:
// -- "NOPID": completely no PID analysis (only primary track cuts are applied)
// -- "BB"   : all tracks are used, but the TPC Bethe-Bloch cut is applied (cut value = 0.2)
// -- "PID"  : realistic PID is used
//
AliRsnPairManager* RsnConfig
(
  const char            *pairMgrName,    // name for the pair manager
  Int_t                  resonancePDG,   // PDG code of resonance (for true pairs)
  AliPID::EParticleType  type1,          // PID of one member of decay (+)
  AliPID::EParticleType  type2           // PID of other member of decay (-)
)
{
  // === NAME DEFINITIONS =========================================================================

  AliRsnPairManager  *pairMgr  = new AliRsnPairManager(pairMgrName);

  // examines the given name to define details about track selection and cuts
  TString               str(pairMgrName);
  AliRsnPair::EPairType pidType;
  Bool_t                useBBCut;
  if (str.Contains("NOPID"))
  {
    pidType  = AliRsnPair::kNoPID;
    useBBCut = kFALSE;
    Info("RsnConfig", "PID TYPE = No PID        -- BB CUT: not used");
  }
  else if (str.Contains("BB"))
  {
    pidType  = AliRsnPair::kNoPID;
    useBBCut = kTRUE;
    Info("RsnConfig", "PID TYPE = No PID        -- BB CUT: used");
  }
  else if (str.Contains("PID"))
  {
    pidType  = AliRsnPair::kRealisticPID;
    useBBCut = kFALSE;
    Info("RsnConfig", "PID TYPE = Realistic PID -- BB CUT: not used");
  }
  else
  {
    Error("RsnConfig", "Unrecognized keywords in the name. Can't continue");
    return 0x0;
  }

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
  cutESDPrimary->GetCuts()->SetMaxNsigmaToVertex(3.0);
  cutESDPrimary->GetCuts()->SetRequireTPCRefit(kTRUE);
  cutESDPrimary->GetCuts()->SetAcceptKinkDaughters(kFALSE);
  cutESDPrimary->GetCuts()->SetMinNClustersTPC(50);
  cutESDPrimary->GetCuts()->SetMaxChi2PerClusterTPC(3.5);
  // -- Bethe-Bloch with kaon mass hypothesis
  Double_t sigmaTPC = 0.065;
  AliRsnCutBetheBloch *cutKaonBB = new AliRsnCutBetheBloch("cutKaonBB", 3.0 * sigmaTPC, AliPID::kKaon);
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
  if (useBBCut)
  {
    cutSetParticle->AddCut(cutKaonBB);
    cutSetParticle->SetCutScheme("cutKaonBB&cutESDPrimary");
  }
  else {
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
  AliRsnFunctionAxis *axisIM   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairInvMass,    4000,  0.0,   2.0);
  AliRsnFunctionAxis *axisPt   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairPt,           50,  0.0,  10.0);
  AliRsnFunctionAxis *axisEta  = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairEta,          15, -1.5,   1.5);
  AliRsnFunctionAxis *axisMult = new AliRsnFunctionAxis(AliRsnFunctionAxis::kEventMult,         8,  0.0, 200.0);

  AliRsnFunction *fcn = new AliRsnFunction;
  fcn->AddAxis(axisIM);
  fcn->AddAxis(axisPt);
  fcn->AddAxis(axisEta);
  fcn->AddAxis(axisMult);

  // add functions to pairs
  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 4; j++) {
      pair[i][j]->AddFunction(fcn);
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
