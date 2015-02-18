//--------------------------------------------------------------------------
// Base macro for submitting muon Resolution analysis.
//
// It needs to load magnetic field, mapping, geometry (+alignment), and reconstruction parameters from the OCDB
// It is done by the task and OCDB storage locations can be parameterized in MuonResolution.C (default is raw://)
//
// The task reads ESDs
// Intermediate results are stored in a file chamberResolution_step<#step>.root
// Final results are stored in the file results.root
//
// Author: Philippe Pillot - SUBATECH Nantes
//--------------------------------------------------------------------------

//______________________________________________________________________________
void RunMuonResolution(TString smode = "local", TString inputFileName = "AliESDs.root",
		       TString rootVersion = "v5-34-08", TString aliphysicsVersion = "vAN-20150215", Int_t nSteps = 5,
		       Bool_t selectPhysics = kTRUE, Bool_t selectTrigger = kTRUE, Bool_t matchTrig = kTRUE,
		       Bool_t applyAccCut = kTRUE, Bool_t applyPDCACut = kTRUE, Double_t minMomentum = 0., Double_t minPt = 0.,
                       Bool_t isMC = kFALSE, Bool_t correctForSystematics = kTRUE, Int_t extrapMode = 1,
                       Bool_t shiftHalfCh = kFALSE, Bool_t shiftDE = kFALSE, Int_t nevents = 1234567890)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// - smode = "local" or "caf"
  /// - inputFileName = an ESD root file or a list of ESDs or a collection of ESDs or a dataset in proof mode
  /// - rootVersion = version of root package to enable on AAF (only used in proof mode)
  /// - alirootVersion = version of aliroot package to enable on AAF (only used in proof mode)
  /// - aliphysicsVersion = version of aliphysics package to enable on AAF (only used in proof mode)
  /// - nSteps = number of times to task is run (at each step it starts with the chamber resolution obtained in the previous one)
  /// - selectPhysics : apply or not the physics selection
  /// - selectTrigger : select only muon trigger events or not (the type of trigger can eventually be changed)
  /// - matchTrigger : select only the tracks matching the trigger or not
  /// - applyAccCut : select only the tracks passing the Rabs and the eta cut or not
  /// - minMomentum : select only the tracks with a total momentum above this value
  /// - if correctForSystematics == kTRUE: the systematic shifts of the residuals is included in the resolution
  /// - if extrapMode == 0: extrapolate from the closest cluster
  /// - if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  /// - nevents = maximum number of processed events
  
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  
  // compile analysis macro locally
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/MUON/dep/MuonResolution.C++g");
  MuonResolution(smode, inputFileName, rootVersion, aliphysicsVersion, nSteps, selectPhysics, selectTrigger, matchTrig, applyAccCut,
                 applyPDCACut, minMomentum, minPt, isMC, correctForSystematics, extrapMode, shiftHalfCh, shiftDE, nevents);
  
}

