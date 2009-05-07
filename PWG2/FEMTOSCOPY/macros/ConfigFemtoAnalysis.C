// ConfigFemtoAnalysis.C - macro to create the splitting/merging
// test with the pion correlation function.
// As a default the anti-splitting and anti-merging cuts are open
// and the two correlation functions:
// AliFemtoShareQualityCorrFctn and AliFemtoTPCInnerCorrFctn 
// can be used to study the splitting (former) and merging (latter) effect
// If ones needs to produce a "clean" sample with both effects removed, 
// one needs to change the cut values to the "reasonable" ones, or perform
// the full systematic analysis with the above-mentioned functions 

// Author: Adam Kisiel. Adam.Kisiel@cern.ch

// parameters:
//
// returns:
// a pointer to the created AliFemtoManager

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#endif

AliFemtoManager *ConfigFemtoAnalysis()
{
  double PionMass = 0.13956995;
  int chargePi = 1;
  
  // Set-up the reader for ALICE ESD
  AliFemtoEventReaderESDChain* Reader=new AliFemtoEventReaderESDChain();
  // Read only constrained momenta - primordial particles
  Reader->SetConstrained(true);
  Reader->SetReadTPCInner(true);
  
  // Setup the manager 
  AliFemtoManager* Manager=new AliFemtoManager();
  // Point to the data source - the reader
  Manager->SetEventReader(Reader);
  
  // Setup the analysis
  AliFemtoSimpleAnalysis* an =new AliFemtoSimpleAnalysis();
  // Number of events to construct the background
  an->SetNumEventsToMix(3);

  // The event selector
  AliFemtoBasicEventCut* mec = new AliFemtoBasicEventCut();
  // Accept events with the given multiplicity
  mec->SetEventMult(0,100000);
  // and z-vertex distance to the center of the TPC
  mec->SetVertZPos(-1000,1000);
	
  // The track selector
  AliFemtoESDTrackCut* dtc = new AliFemtoESDTrackCut();
  // We want positive pions
  dtc->SetPidProbPion(0.2,1.001);
  dtc->SetPidProbMuon(0.0,0.8);
  dtc->SetPidProbKaon(0.0,0.1);
  dtc->SetPidProbProton(0.0,0.1);
  dtc->SetMostProbablePion();
  dtc->SetCharge(chargePi);
  // so we set the correct mass
  dtc->SetMass(PionMass);
  // we select low pt
  dtc->SetPt(0.1,0.7);
  dtc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  dtc->SetminTPCncls(95);
  dtc->SetRemoveKinks(kTRUE);
  dtc->SetLabel(kFALSE);
  dtc->SetMaxITSChiNdof(3.0);
  dtc->SetMaxTPCChiNdof(2.0);
  dtc->SetMaxSigmaToVertex(3.0);

  AliFemtoCutMonitorParticleYPt *cutPass = new AliFemtoCutMonitorParticleYPt("cutPass", 0.13957);
  AliFemtoCutMonitorParticleYPt *cutFail = new AliFemtoCutMonitorParticleYPt("cutFail", 0.13957);
  dtc->AddCutMonitor(cutPass, cutFail);

  // Pair selector
  AliFemtoShareQualityTPCEntranceSepPairCut *sqpc = new AliFemtoShareQualityTPCEntranceSepPairCut();
  // remove split track pairs and pairs that share hits
  
  // Set maximim allowed "quality" for the pair
  //  1.0 - accept all pairs
  // -0.5 - reject all pairs
  // a reasonable value should lie between 0.0 and 0.5
  sqpc->SetShareQualityMax(1.0);

  // Set maximum allowed shared hits fraction per pair
  //  1.0 - accept all pairs
  //  0.0 - reject all pairs
  // a reasonable value is small but nno-zero (0.05)
  sqpc->SetShareFractionMax(1.0);

  // Set minimum allowed separation between nominal TPC entrance points
  // of the two tracks in the pair
  // 0.0 - accept all pairs
  // a reasonable value is 3.0 [cm]
  sqpc->SetTPCEntranceSepMinimum(0.0);
  sqpc->SetRemoveSameLabel(kFALSE);

  // Add the cuts to the analysis
  an->SetEventCut(mec);
  an->SetFirstParticleCut(dtc);
  an->SetSecondParticleCut(dtc);
  an->SetPairCut(sqpc);
  
  // Setup correlation functions
  // A simple qinv correlation function
  AliFemtoQinvCorrFctn *cqinv= new AliFemtoQinvCorrFctn("qinvcf",75,0.0,0.4);
  
  // A correlation function to monitor the splitting and cluster sharing in TPC
  AliFemtoShareQualityCorrFctn *csqqinv= new AliFemtoShareQualityCorrFctn("sqqinvcf",75,0.0,0.4);
  
  // A correlation function to monitor the distance at the entrance to the TPC
  AliFemtoTPCInnerCorrFctn *tpcin = new AliFemtoTPCInnerCorrFctn("tpcin",80, 0.0, 0.4);
  
  // add the correlation functions to the analysis
  an->AddCorrFctn(cqinv);
  an->AddCorrFctn(csqqinv);
  an->AddCorrFctn(tpcin);

  // Add the analysis to the manager
  Manager->AddAnalysis(an);	

  return Manager;
}
