//--------------------------------------------------------------------------
// Base macro for submitting muon Resolution analysis.
//
// In case it is not run with full aliroot, it needs the following libraries:
//  - libSTEERBase.so
//  - libESD.so
//  - libAOD.so
//  - libANALYSIS.so
//  - libANALYSISalice.so
//  - libGui.so
//  - libMinuit.so
//  - libProofPlayer.so
//  - libXMLParser.so
//  - libRAWDatabase.so
//  - libCDB.so
//  - libSTEER.so
//  - libMUONcore.so
//  - libMUONmapping.so
//  - libMUONcalib.so
//  - libMUONgeometry.so
//  - libMUONtrigger.so
//  - libMUONraw.so
//  - libMUONbase.so
//  - libMUONrec.so
//  - libCORRFW.so
//  - libPWGHFbase.so
//  - libPWGmuondep.so
//
// It also needs to load magnetic field, mapping, geometry (+alignment), and reconstruction parameters from the OCDB
//
// The task reads ESDs
// Intermediate results are stored in a file chamberResolution_step<#step>.root
// Final results are stored in the file results.root
//
// Author: Philippe Pillot - SUBATECH Nantes
//--------------------------------------------------------------------------

void LoadAlirootLocally(TString& extraLibs);

//______________________________________________________________________________
void RunMuonResolution(TString smode = "local", TString inputFileName = "AliESDs.root",
		       TString alirootVersion = "VO_ALICE@AliRoot::v4-20-12-AN", Int_t nSteps = 5,
		       Bool_t selectPhysics = kFALSE, Bool_t selectTrigger = kFALSE, Bool_t matchTrig = kTRUE,
		       Bool_t applyAccCut = kTRUE, Double_t minMomentum = 0., Bool_t correctForSystematics = kTRUE,
		       Int_t extrapMode = 1, Int_t nevents = 1234567890)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// - smode = "local" or "proof"
  /// - inputFileName = an ESD root file or a list of ESDs or a collection of ESDs or a dataset in proof mode
  /// - alirootVersion = version of aliroot package to enable on AAF (only used in proof mode)
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
  
  // Load libraries locally
  TString extraLibs = "RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec:CORRFW:PWGPPMUONdep";
  LoadAlirootLocally(extraLibs);
  
  // compile analysis macro locally
  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/MUON/dep/MuonResolution.C++g");
  MuonResolution(smode, inputFileName, alirootVersion, nSteps, selectPhysics, selectTrigger, matchTrig,
		 applyAccCut, minMomentum, correctForSystematics, extrapMode, nevents, extraLibs);
  
}

//______________________________________________________________________________
void LoadAlirootLocally(TString& extraLibs)
{
  /// Load libraries locally
  
  // Load common libraries
  gSystem->Load("libVMC");
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libGui.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  // Load additional libraries
  gSystem->Load("libProofPlayer");
  TObjArray* libs = extraLibs.Tokenize(":");
  for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
    gSystem->Load(Form("lib%s",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName()));
  delete libs;
  
  // Load lib for final mchview display
  gSystem->Load("libMUONgraphics");
  
  // Add AliRoot includes to compile the macro
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON/mapping");
  gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS/macros");
  
}

