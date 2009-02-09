//
// Creates an AliRsnPairMgr containing all invmass spectra
// for the PHI --> K+ K- resonance:
//
// - signal
// - event mixing
// - like-sign
// - true pairs
//
// When required, PDG code of phi is 333
//
// In order to allow analysis customization, some user-defined parameters
// are listed at the beginning of this macro, to define some details like
// the number of event mixing to do, and what PID to be used.
// Since they are many, they are hard-coded in the macro and the user should
// take care of them when launching it for an analysis.
// Moreover, here some cuts are defined, which are very general: if a user
// wants to add more specific cuts, he should take a look to the "CUTS" section
// of this macro.
//
// The idea of this macro is to be loaded and launched from another
// upper-level macro which prepares the AnalysisTask object to run it,
// which is defined in "CreateAnalysisManager.C" macro in this directory
//

AliRsnPairMgr* CreatePairsPhi(const char *name = "PHI")
{
  AliRsnPairMgr  *pairMgr  = new AliRsnPairMgr(name);

  // ========== USER CUSTOMIZATION VARIABLES ==========

  Int_t   iResPDG             = 333;
  Int_t   nMixEvents          = 10;
  Bool_t  boolUseNoPID        = kTRUE;
  Bool_t  boolUseRealisticPID = kTRUE;
  Bool_t  boolUsePerfectPID   = kTRUE;

  // ======= END USER CUSTOMIZATION VARIABLES =========

  // =========== DEFINE PAIRS ==============

  // decay tree definitions
  // for a PHI resonance (PDG = 333) decaying into K+ K-
  // and for related like-sign pairs
  AliRsnPairDef *defUnlike = new AliRsnPairDef(AliRsnPID::kKaon, '+', AliRsnPID::kKaon, '-', iResPDG);
  AliRsnPairDef *defLikePP = new AliRsnPairDef(AliRsnPID::kKaon, '+', AliRsnPID::kKaon, '+', iResPDG);
  AliRsnPairDef *defLikeMM = new AliRsnPairDef(AliRsnPID::kKaon, '-', AliRsnPID::kKaon, '-', iResPDG);

  // No PID
  AliRsnPair *pairUnlike_NoPID_Signal = new AliRsnPair(AliRsnPair::kNoPID, defUnlike);
  AliRsnPair *pairUnlike_NoPID_True   = new AliRsnPair(AliRsnPair::kNoPID, defUnlike);
  AliRsnPair *pairUnlike_NoPID_Mix    = new AliRsnPair(AliRsnPair::kNoPIDMix, defUnlike);
  AliRsnPair *pairLikePP_NoPID        = new AliRsnPair(AliRsnPair::kNoPID, defLikePP);
  AliRsnPair *pairLikeMM_NoPID        = new AliRsnPair(AliRsnPair::kNoPID, defLikeMM);
  // end No PID

  // Perfect PID
  AliRsnPair *pairUnlike_PerfectPID_Signal = new AliRsnPair(AliRsnPair::kPerfectPID, defUnlike);
  AliRsnPair *pairUnlike_PerfectPID_True   = new AliRsnPair(AliRsnPair::kPerfectPID, defUnlike);
  AliRsnPair *pairUnlike_PerfectPID_Mix    = new AliRsnPair(AliRsnPair::kPerfectPIDMix, defUnlike);
  AliRsnPair *pairLikePP_PerfectPID        = new AliRsnPair(AliRsnPair::kPerfectPID, defLikePP);
  AliRsnPair *pairLikeMM_PerfectPID        = new AliRsnPair(AliRsnPair::kPerfectPID, defLikeMM);
  // end Perfect PID

  // Perfect PID
  AliRsnPair *pairUnlike_RealisticPID_Signal = new AliRsnPair(AliRsnPair::kRealisticPID, defUnlike);
  AliRsnPair *pairUnlike_RealisticPID_True   = new AliRsnPair(AliRsnPair::kRealisticPID, defUnlike);
  AliRsnPair *pairUnlike_RealisticPID_Mix    = new AliRsnPair(AliRsnPair::kRealisticPIDMix, defUnlike);
  AliRsnPair *pairLikePP_RealisticPID        = new AliRsnPair(AliRsnPair::kRealisticPID, defLikePP);
  AliRsnPair *pairLikeMM_RealisticPID        = new AliRsnPair(AliRsnPair::kRealisticPID, defLikeMM);
  // end Realistic PID

  // =========== END DEFINE PAIRS ==============

  // =========== CUTS ==============

  // cuts on tracks:
  // - defined in 'AddRsnAnalysisTask.C' for single-step analysis

  // cuts on pairs:
  // - true daughters of a phi resonance (only for true pairs histogram)
  AliRsnCut *cutTruePair = new AliRsnCut("cutTrue", "cutTrue", AliRsnCut::kIsTruePair, iResPDG);

  // cut set definition for true pairs
  AliRsnCutSet *cutSetPairTrue = new AliRsnCutSet("truePairs");
  cutSetPairTrue->AddCut(cutTruePair);
  cutSetPairTrue->SetCutScheme("cutTrue");

  // define cut manager for true pairs
  AliRsnCutMgr *cutMgrTrue = new AliRsnCutMgr("true", "True pairs");
  cutMgrTrue->SetCutSet(AliRsnCut::kPair, cutSetPairTrue);

  // add cuts to pair analysis
  pairUnlike_NoPID_True->SetCutMgr(cutMgrTrue);
  pairUnlike_PerfectPID_True->SetCutMgr(cutMgrTrue);
  pairUnlike_RealisticPID_True->SetCutMgr(cutMgrTrue);

  // =========== END CUTS ==============

  // =========== FUNCTIONS ==============

  // define histogram templates
  AliRsnHistoDef *hdIM  = new AliRsnHistoDef(800, 0.0, 2.0);     // invmass
  AliRsnHistoDef *hdRES = new AliRsnHistoDef(200, -10.0, 10.0);  // resolution

  // functions
  AliRsnFunction *fcnIM  = new AliRsnFunction(AliRsnFunction::kInvMass, hdIM);      // invmass
  AliRsnFunction *fcnRES = new AliRsnFunction(AliRsnFunction::kResolution, hdRES);    // IM resolution

  // uncomment these lines when doing analysis in momentum bins
  // in this case, take care of the dimension and values in the template array
  Double_t mom[7] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0};
  fcnIM->SetBinningCut(AliRsnCut::kTransMomentum, 7, mom);
  fcnRES->SetBinningCut(AliRsnCut::kTransMomentum, 7, mom);

  pairUnlike_NoPID_Signal->AddFunction(fcnIM);
  pairUnlike_NoPID_True->AddFunction(fcnIM);
  pairUnlike_NoPID_Mix->AddFunction(fcnIM);
  pairLikePP_NoPID->AddFunction(fcnIM);
  pairLikeMM_NoPID->AddFunction(fcnIM);

  pairUnlike_PerfectPID_Signal->AddFunction(fcnIM);
  pairUnlike_PerfectPID_True->AddFunction(fcnIM);
  pairUnlike_PerfectPID_Mix->AddFunction(fcnIM);
  pairLikePP_PerfectPID->AddFunction(fcnIM);
  pairLikeMM_PerfectPID->AddFunction(fcnIM);

  pairUnlike_RealisticPID_Signal->AddFunction(fcnIM);
  pairUnlike_RealisticPID_True->AddFunction(fcnIM);
  pairUnlike_RealisticPID_Mix->AddFunction(fcnIM);
  pairLikePP_RealisticPID->AddFunction(fcnIM);
  pairLikeMM_RealisticPID->AddFunction(fcnIM);

  pairUnlike_NoPID_Signal->AddFunction(fcnRES);
  pairUnlike_PerfectPID_Signal->AddFunction(fcnRES);
  pairUnlike_RealisticPID_Signal->AddFunction(fcnRES);

  // =========== END FUNCTIONS =============

  if (boolUseNoPID) {
    pairMgr->AddPair(pairUnlike_NoPID_Signal);
    pairMgr->AddPair(pairUnlike_NoPID_True);
    pairMgr->AddPair(pairUnlike_NoPID_Mix);
    pairMgr->AddPair(pairLikePP_NoPID);
    pairMgr->AddPair(pairLikeMM_NoPID);
  }

  if (boolUsePerfectPID) {
    pairMgr->AddPair(pairUnlike_PerfectPID_Signal);
    pairMgr->AddPair(pairUnlike_PerfectPID_True);
    pairMgr->AddPair(pairUnlike_PerfectPID_Mix);
    pairMgr->AddPair(pairLikePP_PerfectPID);
    pairMgr->AddPair(pairLikeMM_PerfectPID);
  }

  if (boolUseRealisticPID) {
    pairMgr->AddPair(pairUnlike_RealisticPID_Signal);
    pairMgr->AddPair(pairUnlike_RealisticPID_True);
    pairMgr->AddPair(pairUnlike_RealisticPID_Mix);
    pairMgr->AddPair(pairLikePP_RealisticPID);
    pairMgr->AddPair(pairLikeMM_RealisticPID);
  }
  return pairMgr;
}
