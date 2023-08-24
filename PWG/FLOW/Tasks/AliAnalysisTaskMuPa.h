/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************** 
*   TBI add description eventually    * 
**************************************/ 

#ifndef ALIANALYSISTASKMUPA_H
#define ALIANALYSISTASKMUPA_H

#include <AliAnalysisTaskSE.h>
#include <AliVEvent.h>
#include <AliAODEvent.h>
#include <AliAODHeader.h>
#include <AliMCEvent.h>
#include <AliVTrack.h>
#include <AliVParticle.h>
#include <AliAODTrack.h>
#include <AliAODMCParticle.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TExMap.h>
#include <TComplex.h>
#include <TArrayD.h>
#include <TArrayI.h>
#include <TSystem.h>
#include <TF1.h>
#include <TF3.h>
#include <TStopwatch.h>
#include <TFormula.h>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

// Global variables:
const Int_t gCentralityEstimators = 4; // set here number of supported centrality estimators
const Int_t gKinematicVariables = 5; // number of supported kinematic variables: [phi,pt,eta,e,charge]
const Int_t gFilterBits = 17; // number of filterbits to scan
const Int_t gEventHistograms = 2; // total number of non-classified event histograms
const Int_t gParticleHistograms = 12; // total number of non-classified particle histograms, keep in sync. with eParticle 
const Int_t gCentralMultiplicity = 1; // multiplicities defined centrally, e.g. ref. mult.
const Int_t gWeights = 3; // phi, pt, eta
const Int_t gQAAnomalousEvents = 1; // |vertex| = 0; 
const Int_t gQASelfCorrelations = 3; // phi, pt, eta
const Int_t gQAEventCutCounter = 23; // see TString secc[gQAEventCutCounter] in .cxx
const Int_t gQAParticleCutCounter = 40; // see TString spcc[gQAParticleCutCounter] in .cxx . Used also for SequentialParticleCutCounter, via cloning
const Int_t gGenericCorrelations = 7; // correlations between various quantities (see .cxx for documentation)
const Int_t gMaxCorrelator = 12; // 
const Int_t gMaxHarmonic = 9; // 
const Int_t gMaxIndex = 300; // per order
const Int_t gMaxNoBinsKine = 1000; 
const Int_t gKineDependenceVariables = 2; // pt,eta => special treatment, since I need to calculate separately differential q-vectors here

// enums:
enum eBins {nBins,min,max};
enum eCentralityEstimator { V0M, SPDTracklets, CL0, CL1 }; 
enum eXYZ { X = 0, Y = 1, Z = 2 };
enum eBeforeAfter { BEFORE = 0, AFTER = 1 };
enum eRecoSim { RECO = 0, SIM = 1 };
enum eKinematics { PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4 };
enum eDefaultColors { COLOR = kBlack, FILLCOLOR = kGray };
enum eParticle { TPCNcls, TPCnclsS, TPCnclsFractionShared, TPCNCrossedRows, TPCChi2perNDF, TPCFoundFraction, Chi2TPCConstrainedVsGlobal, ITSNcls, ITSChi2perNDF, TPCNclsF, HasPointOnITSLayer, IsGlobalConstrained };
enum eEvent { MagneticField, PrimaryVertex };
enum eCentralMultiplicity { RefMultComb08 };
enum eqvectorKine { PTq = 0, ETAq = 1 };
enum eAsFunctionOf { AFO_INTEGRATED = 0, AFO_MULTIPLICITY = 1, AFO_CENTRALITY = 2, AFO_PT = 3, AFO_ETA = 4, eAsFunctionOf_N }; // prefix is needed, to avoid conflict with enum eKinematics

//================================================================================================================

class AliAnalysisTaskMuPa : public AliAnalysisTaskSE{
 public: 
  AliAnalysisTaskMuPa();
  AliAnalysisTaskMuPa(const char *name);
  virtual ~AliAnalysisTaskMuPa(); 
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  
  // 0) Methods called in the constructor:
  virtual void InitializeNonBuiltInTypes();
  virtual void InitializeArrays();
   virtual void InitializeArraysForQAHistograms();
   virtual void InitializeArraysForControlEventHistograms();
   virtual void InitializeArraysForControlParticleHistograms();
   virtual void InitializeArraysForQvectors();
   virtual void InitializeArraysForWeights();
   virtual void InitializeArraysForCentralityWeights();
   virtual void InitializeArraysForCorrelationsHistograms();
   virtual void InitializeArraysForNestedLoopsHistograms();
   virtual void InitializeArraysForToyNUA();
   virtual void InitializeArraysForInternalValidation();
   virtual void InitializeArraysForTest0();
   virtual void InitializeArraysForCommonLabels();
  virtual void DefaultConfiguration();
  virtual void DefaultBinning();
  virtual void DefaultCuts();

  // 1) Methods called in UserCreateOutputObjects():
  virtual void InsanityChecks(); 
  virtual void BookBaseProfile();
  virtual void BookAndNestAllLists();  
  virtual void BookQAHistograms();
  virtual void BookControlEventHistograms();
  virtual void BookControlParticleHistograms();
  virtual void BookQvectorHistograms();
  virtual void BookWeightsHistograms();
  virtual void BookCentralityWeightsHistograms();
  virtual void BookCorrelationsHistograms();
  virtual void BookNestedLoopsHistograms();
  virtual void BookToyNUAHistograms();
  virtual void BookInternalValidationHistograms();
  virtual void BookTest0Histograms();
  virtual void BookFinalResultsHistograms();
  virtual void StoreLabelsInPlaceholder(const char *source); 
  virtual Bool_t RetrieveCorrelationsLabels();

  // 2) Methods called in UserExec(Option_t *):
  virtual void InternalValidation();
  virtual void PrintEventInfo(AliVEvent *ave); // print event medatata (for AOD: fRun, fBunchCross, fOrbit, fPeriod). Enable via task->SetPrintEventInfo()  
  virtual void FillQAHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs);
  virtual void FillQAHistograms(AliAODEvent *aAOD, AliMCEvent *aMC);
  virtual void FilterEvent(AliVEvent *ave);
  virtual void UpdateHistogramBookingsWithRunInfo(AliVEvent *ave); // e.g. in user exec, I can use aAOD->GetRunNumber() to get run number. Clearly, I have to do this only once, therefore f
  virtual void FillControlEventHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs); // before or after event cuts, reco or sim
  virtual void FillControlParticleHistograms(AliVParticle *vParticle, const Int_t ba, const Int_t rs); // before or after particle cuts, reco or sim
  virtual void GlobalTracksAOD(AliAODEvent *aAOD);
  Bool_t SurvivesEventCuts(AliVEvent *ave);
  void EventCutCounter(AliVEvent *ave); // only for QA
  void SequentialEventCutCounter(AliVEvent *ave); // only for QA
  Bool_t SurvivesParticleCuts(AliVParticle *vParticle); // applied e.g. on TPC-only
  void ParticleCutCounter(AliVParticle *vParticle); // only for QA
  void SequentialParticleCutCounter(AliVParticle *vParticle); // only for QA
  virtual Double_t Weight(const Double_t &value, const char *variable);
  virtual Double_t CentralityWeight(const Double_t &value);
  virtual void CalculateCorrelations();
  virtual void CalculateNestedLoops(); // calculate all standard isotropic correlations  
  virtual Double_t CalculateCustomNestedLoop(TArrayI *harmonics); // calculate nested loop for the specified harmonics
  virtual void ResetEventByEventQuantities();
  virtual void ResetQ(); // reset the components of generic Q-vectors
  virtual void OnlineMonitoring();
  Bool_t SpecifiedEvent(AliVEvent *ave);
  virtual void MakeLookUpTable(AliAODEvent *aAOD, AliMCEvent *aMC);
  virtual void RandomIndices(AliVEvent *ave);
  virtual void CalculateTest0();
  virtual void CalculateKineCorrelations(const char* kc); 
  virtual void CalculateKineTest0(const char* kc);
  virtual Double_t CalculateKineCustomNestedLoop(TArrayI *harmonics, const char* kc, Int_t bin); // calculate custom nested loop for the specified harmonics, kine. variable, and bin

  // 3) Methods called in Terminate(Option_t *):
  //    a) Get pointers:
  virtual void GetPointers(TList *baseList);
   virtual void GetPointersForControlEventHistograms();
   virtual void GetPointersForControlParticleHistograms();
  virtual void ComparisonNestedLoopsVsCorrelations();

  // 4) Q-vectors:
  virtual TComplex Q(Int_t n, Int_t p);
  virtual TComplex One(Int_t n1);
  virtual TComplex Two(Int_t n1, Int_t n2);
  virtual TComplex Three(Int_t n1, Int_t n2, Int_t n3);
  virtual TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  virtual TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  virtual TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6);
  virtual TComplex Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7);
  virtual TComplex Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8);
  virtual TComplex Nine(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9);
  virtual TComplex Ten(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10);
  virtual TComplex Eleven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11);
  virtual TComplex Twelve(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8, Int_t n9, Int_t n10, Int_t n11, Int_t n12);
  virtual TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0); // Credits: Kristjan Gulbrandsen (gulbrand@nbi.dk) 
  virtual TComplex TheoreticalValue(TArrayI *harmonics, TArrayD *amplitudes, TArrayD *planes); // for the specified amplitudes and symmetry planes, return the theoretical value of correlator

  // 5) Setters and getters:
  void SetRealData(Bool_t rd){this->fRealData = rd;};
  void SetUseFisherYates(Bool_t ufy){this->fUseFisherYates = ufy;};
  void SetFixedNumberOfRandomlySelectedParticles(Int_t fnorsp)
  {
   this->fFixedNumberOfRandomlySelectedParticles = fnorsp;
   this->fUseFixedNumberOfRandomlySelectedParticles = kTRUE;
  };
  void SetDataTakingPeriod(const char *dtp) {this->fDataTakingPeriod = dtp;};
  void SetAODNumber(const char *an) {this->fAODNumber = an;};
  void SetRunNumber(const char *rn) {this->fRunNumber = rn;};
  void SetVerbose(Bool_t v) {this->fVerbose = v;};
  void SetRandomSeed(UInt_t rs) {this->fRandomSeed = rs;};
  void SetTrigger(const char *t) {this->fTrigger = t; this->fUseTrigger = kTRUE;};

  void SetControlEventHistogramsList(TList* const cehl) {this->fControlEventHistogramsList = cehl;};
  TList* GetControlEventHistogramsList() const {return this->fControlEventHistogramsList;} 

  // QA:
  void SetFillQAHistograms(Bool_t fqah) {this->fFillQAHistograms = fqah;};
  Bool_t GetFillQAHistograms() const {return this->fFillQAHistograms;};
  void SetFillQAHistogramsAll(Bool_t fqaha) {this->fFillQAHistogramsAll = fqaha;};
  Bool_t GetFillQAHistogramsAll() const {return this->fFillQAHistogramsAll;};
  void SetTerminateAfterQA(Bool_t taqa) {this->fTerminateAfterQA = taqa;};
  Bool_t GetTerminateAfterQA() const {return this->fTerminateAfterQA;};
  void SetQAFilterBits(TArrayI *fb){this->fQAFilterBits = fb;};
  void SetQACheckSelfCorrelations(Bool_t qacsc) {this->fQACheckSelfCorrelations = qacsc;};

  // Multiplicity:
  void SetMultiplicityBins(Int_t const nbins, Double_t min, Double_t max)
  {
   this->fMultiplicityBins[0] = nbins;
   this->fMultiplicityBins[1] = min;
   this->fMultiplicityBins[2] = max;
  };
  void SetSelectedTracksCuts(const Int_t min, const Int_t max)
  {
   this->fSelectedTracksCuts[0] = min;
   this->fSelectedTracksCuts[1] = max;
   this->fUseSelectedTracksCuts = kTRUE;
  }

  // Centrality:
  void SetCentralityBins(Int_t const nbins, Double_t min, Double_t max)
  {
   this->fCentralityBins[0] = nbins;
   this->fCentralityBins[1] = min;
   this->fCentralityBins[2] = max;
  };
  void SetCentralityCuts(const Double_t min, const Double_t max)
  {
   this->fCentralityCuts[0] = min;
   this->fCentralityCuts[1] = max;
   this->fUseCentralityCuts = kTRUE;
  }
  void SetCentralityEstimator(const char *ce) {this->fCentralityEstimator = ce;};
  TString GetCentralityEstimator() const {return this->fCentralityEstimator;};
  void SetTaskName(const char *tn) {this->fTaskName = tn;};

  void SetCentralityCorrelationsCuts(const char* firstEstimator, const char* secondEstimator, const Double_t cut, const int cutVersion)
  {
   // This is the cutVersion:
   // cutVersion = 0 (relative) : |(firstEstimator-secondEstimator)/(firstEstimator+secondEstimator)| > cut => reject the event
   // cutVersion = 1 (absolute) : |(firstEstimator-secondEstimator)| > cut => reject the event

   Int_t ce1 = -44;
   Int_t ce2 = -44;

   if(TString(firstEstimator).EqualTo("V0M")){ce1 = 0;} 
   else if (TString(firstEstimator).EqualTo("SPDTracklets")){ce1 = 1;} 
   else if (TString(firstEstimator).EqualTo("CL0")){ce1 = 2;} 
   else if (TString(firstEstimator).EqualTo("CL1")){ce1 = 3;} 
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }

   if(TString(secondEstimator).EqualTo("V0M")){ce2 = 0;} 
   else if (TString(secondEstimator).EqualTo("SPDTracklets")){ce2 = 1;} 
   else if (TString(secondEstimator).EqualTo("CL0")){ce2 = 2;} 
   else if (TString(secondEstimator).EqualTo("CL1")){ce2 = 3;} 
   else{ Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(2); }

   if(!(cutVersion==0||cutVersion==1)){ Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(3); }

   this->fUseCentralityCorrelationsCuts[ce1][ce2] = kTRUE;
   this->fCentralityCorrelationsCuts[ce1][ce2] = cut;
   this->fCentralityCorrelationCutVersion = cutVersion;

  } // void SetCentralityCorrelationsCuts(const char* firstEstimator, const char* secondEstimator, const Double_t cut)

  // Vertex:
  void SetVertexBins(Int_t const nbins, Double_t min, Double_t max)
  {
   this->fVertexBins[0] = nbins;
   this->fVertexBins[1] = min;
   this->fVertexBins[2] = max;
  };
  void SetVertexCuts(const char* sxyz, const Double_t min, const Double_t max)
  {
   Int_t xyz = -44;
   if(TString(sxyz).EqualTo("x")){xyz = 0;} 
   else if (TString(sxyz).EqualTo("y")){xyz = 1;} 
   else if (TString(sxyz).EqualTo("z")){xyz = 2;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fVertexCuts[xyz][0] = min;
   this->fVertexCuts[xyz][1] = max;
   this->fUseVertexCuts[xyz] = kTRUE;
  }
  void SetNContributorsBins(Int_t const nbins, Double_t min, Double_t max)
  {
   this->fNContributorsBins[0] = nbins;
   this->fNContributorsBins[1] = min;
   this->fNContributorsBins[2] = max;
  };
  void SetNContributorsCuts(const Int_t min, const Int_t max)
  {
   this->fNContributorsCuts[0] = min;
   this->fNContributorsCuts[1] = max;
   this->fUseNContributorsCuts = kTRUE;
  }
  void SetMinVertexDistance(const Double_t vd)
  {
   this->fMinVertexDistance = vd;
   this->fUseMinVertexDistanceCut = kTRUE;
  }

  // Remaining event distributions:
  void SetEventBins(const char* type, const Double_t nbins, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("MagneticField")){var = MagneticField;} 
   else if (TString(type).EqualTo("PrimaryVertex")){var = PrimaryVertex;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fEventBins[var][0] = nbins;
   this->fEventBins[var][1] = min;
   this->fEventBins[var][2] = max;
  }
  void SetEventCuts(const char* type, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("MagneticField")){var = MagneticField;} 
   else if (TString(type).EqualTo("PrimaryVertex")){var = PrimaryVertex;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); } 
   this->fEventCuts[var][0] = min;
   this->fEventCuts[var][1] = max;
   this->fUseEventCuts[var] = kTRUE;
  }

  // Generic correlations:
  void SetSelContrTreshold(const Float_t treshold)
  {
   this->fSelContrTreshold = treshold;
   this->fUseGenericCorrelationsCuts[0] = kTRUE;
  }

  void SetControlParticleHistogramsList(TList* const cphl) {this->fControlParticleHistogramsList = cphl;};
  TList* GetControlParticleHistogramsList() const {return this->fControlParticleHistogramsList;} 
  void SetFillControlParticleHistograms(Bool_t fcph) {this->fFillControlParticleHistograms = fcph;};
  void SetUseFakeTracks(const Bool_t uft) {this->fUseFakeTracks = uft;};

  void SetFilterBit(Int_t fb) {this->fFilterBit = fb;};
  Int_t GetFilterBit() const {return this->fFilterBit;};
  void SetUseDefaultFilterBitCuts(Int_t fb); // see .cxx
  void SetUseOnlyPrimaries(Bool_t uop) {this->fUseOnlyPrimaries = uop;};
  Int_t GetUseOnlyPrimaries() const {return this->fUseOnlyPrimaries;};
  void SetPrimaryDefinitionInMonteCarlo(const char *pdimc) {this->fPrimaryDefinitionInMonteCarlo = pdimc;};

  // Needed only for PID studies, setting FilterBit avoids double-counting:
  void SetFilterGlobalTracksAOD(const Bool_t fgta){this->fFilterGlobalTracksAOD = fgta;}; 

  // Kinematics:
  void SetKinematicsBins(const char* kv, const Double_t nbins, const Double_t min, const Double_t max) // used for control histograms. See also SetKineDependenceBins(...)
  {
   Int_t var = -44;
   if(TString(kv).EqualTo("phi")){var = PHI;} 
   else if (TString(kv).EqualTo("pt")){var = PT;} 
   else if (TString(kv).EqualTo("eta")){var = ETA;}
   else if (TString(kv).EqualTo("e")){var = E;}
   else if (TString(kv).EqualTo("charge")){var = CHARGE;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fKinematicsBins[var][0] = nbins;
   this->fKinematicsBins[var][1] = min;
   this->fKinematicsBins[var][2] = max;
  }

  void SetKineDependenceBins(const char* kv, TArrayD *iva) // used for custom binning of all histos related to results. See also SetKinematicsBins(...)
  {
   Int_t var = -44;
   if (TString(kv).EqualTo("pt")){var = PTq;} 
   else if (TString(kv).EqualTo("eta")){var = ETAq;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fUseCustomKineDependenceBins[var] = kTRUE;
   this->fKineDependenceBins[var] = iva;
  }

  void SetKinematicsCuts(const char* kc, const Double_t min, const Double_t max) // used for control histograms
  {
   Int_t var = -44;
   if(TString(kc).EqualTo("phi")){var = PHI;} 
   else if (TString(kc).EqualTo("pt")){var = PT;} 
   else if (TString(kc).EqualTo("eta")){var = ETA;}
   else if (TString(kc).EqualTo("e")){var = E;}
   else if (TString(kc).EqualTo("charge")){var = CHARGE;} // it is hardwired in SurvivesParticleCuts(...) that neutral particles are rejected, whenever this setting is called
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fKinematicsCuts[var][0] = min;
   this->fKinematicsCuts[var][1] = max;
   this->fUseKinematicsCuts[var] = kTRUE;
  }

  // DCA:
  void SetDCABins(const char* xyTz, const Double_t nbins, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(xyTz).EqualTo("xy")){var = 0;} 
   else if (TString(xyTz).EqualTo("z")){var = 1;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fDCABins[var][0] = nbins;
   this->fDCABins[var][1] = min;
   this->fDCABins[var][2] = max;
  }
  void SetDCACuts(const char* dc, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(dc).EqualTo("xy")){var = 0;} 
   else if (TString(dc).EqualTo("z")){var = 1;} 
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fDCACuts[var][0] = min;
   this->fDCACuts[var][1] = max;
   this->fUseDCACuts[var] = kTRUE;
  } 
  void SetPtDependentDCAxyCut(const char* parameterization)
  {
   // Supported pre-defined formulas for DCAxy vs. pT dependence: "definition-2010", "definition-2011", "definition-2015"
   // If none of them is used, then content of 'parameterization' is interpreted as user-defined formula for DCAxy vs. pT
   // The chosen parameterization is stored in 5th bin of fControlParticleHistogramsPro, so it can be retrieved afterward from the output of analysis.
   if(TString(parameterization).EqualTo("definition-2010")){this->fPtDependentDCAxyParameterization = "0.0182+0.0350/x^1.01";} 
   else if(TString(parameterization).EqualTo("definition-2011")){this->fPtDependentDCAxyParameterization = "0.0105+0.0350/x^1.1";} 
   else if(TString(parameterization).EqualTo("definition-2015")){this->fPtDependentDCAxyParameterization = "0.0105+0.0350/x^1.1";} 
   // else if ... 
   else // use user-supplied formula for parameterization. Use "x" to denote variable in the formula string, e.g. "1.44 + 2.*x"
   {
    this->fPtDependentDCAxyParameterization = parameterization;
   } 
    
   // Only if this setter is called, for whichever parameterization, pT dependent DCAxy cut will be used:
   this->fUsePtDependentDCAxyParameterization = kTRUE;

  } // void SetPtDependentDCAxyCut(const char* parameterization)

  // Remaining particle distributions:
  void SetParticleBins(const char* type, const Double_t nbins, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("TPCNcls")){var = TPCNcls;} 
   else if (TString(type).EqualTo("TPCnclsS")){var = TPCnclsS;}
   else if (TString(type).EqualTo("TPCnclsFractionShared")){var = TPCnclsFractionShared;}
   else if (TString(type).EqualTo("TPCNCrossedRows")){var = TPCNCrossedRows;}
   else if (TString(type).EqualTo("TPCChi2perNDF")){var = TPCChi2perNDF;}
   else if (TString(type).EqualTo("TPCFoundFraction")){var = TPCFoundFraction;}
   else if (TString(type).EqualTo("Chi2TPCConstrainedVsGlobal")){var = Chi2TPCConstrainedVsGlobal;}
   else if (TString(type).EqualTo("ITSNcls")){var = ITSNcls;}
   else if (TString(type).EqualTo("ITSChi2perNDF")){var = ITSChi2perNDF;}
   else if (TString(type).EqualTo("TPCNclsF")){var = TPCNclsF;}
   this->fParticleBins[var][0] = nbins;
   this->fParticleBins[var][1] = min;
   this->fParticleBins[var][2] = max;
  }
  void SetParticleCuts(const char* type, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("TPCNcls")){var = TPCNcls;} 
   else if (TString(type).EqualTo("TPCnclsS")){var = TPCnclsS;}
   else if (TString(type).EqualTo("TPCnclsFractionShared")){var = TPCnclsFractionShared;}
   else if (TString(type).EqualTo("TPCNCrossedRows")){var = TPCNCrossedRows;}
   else if (TString(type).EqualTo("TPCChi2perNDF")){var = TPCChi2perNDF;}
   else if (TString(type).EqualTo("TPCFoundFraction")){var = TPCFoundFraction;}
   else if (TString(type).EqualTo("Chi2TPCConstrainedVsGlobal")){var = Chi2TPCConstrainedVsGlobal;}
   else if (TString(type).EqualTo("ITSNcls")){var = ITSNcls;}
   else if (TString(type).EqualTo("ITSChi2perNDF")){var = ITSChi2perNDF;}
   else if (TString(type).EqualTo("TPCNclsF")){var = TPCNclsF;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fParticleCuts[var][0] = min;
   this->fParticleCuts[var][1] = max;
   this->fUseParticleCuts[var] = kTRUE;
  }  
  void SetUseParticleCuts(const char* type, const Bool_t b) // use this setter just to switch on/off sum cut, with was set either by default, or via SetParticleCuts(...)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("TPCNcls")){var = TPCNcls;} 
   else if (TString(type).EqualTo("TPCnclsS")){var = TPCnclsS;}
   else if (TString(type).EqualTo("TPCnclsFractionShared")){var = TPCnclsFractionShared;}
   else if (TString(type).EqualTo("TPCNCrossedRows")){var = TPCNCrossedRows;}
   else if (TString(type).EqualTo("TPCChi2perNDF")){var = TPCChi2perNDF;}
   else if (TString(type).EqualTo("TPCFoundFraction")){var = TPCFoundFraction;}
   else if (TString(type).EqualTo("Chi2TPCConstrainedVsGlobal")){var = Chi2TPCConstrainedVsGlobal;}
   else if (TString(type).EqualTo("ITSNcls")){var = ITSNcls;}
   else if (TString(type).EqualTo("ITSChi2perNDF")){var = ITSChi2perNDF;}
   else if (TString(type).EqualTo("TPCNclsF")){var = TPCNclsF;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   this->fUseParticleCuts[var] = b;
  }  
  void SetAtLeastOnePointInTheSPD(Bool_t alopits) {this->fAtLeastOnePointInTheSPD = alopits;};
  void SetIgnoreGlobalConstrained(Bool_t igc) {this->fIgnoreGlobalConstrained = igc;};

  void SetCalculateQvector(Bool_t cqv) {this->fCalculateQvector = cqv;};
  Bool_t GetCalculateQvector() const {return this->fCalculateQvector;};

  void SetCalculateCorrelations(Bool_t cc) {this->fCalculateCorrelations = cc;};
  Bool_t GetCalculateCorrelations() const {return this->fCalculateCorrelations;};

  void SetDoNotCalculateCorrelationsAsFunctionOf(const char* kc)
  {
   if(TString(kc).EqualTo("integrated")){fDoNotCalculateCorrelationsAsFunctionOf[AFO_INTEGRATED] = kTRUE;}
   else if(TString(kc).EqualTo("multiplicity")){fDoNotCalculateCorrelationsAsFunctionOf[AFO_MULTIPLICITY] = kTRUE;} 
   else if(TString(kc).EqualTo("centrality")){fDoNotCalculateCorrelationsAsFunctionOf[AFO_CENTRALITY] = kTRUE;} 
   else if(TString(kc).EqualTo("pt")){fDoNotCalculateCorrelationsAsFunctionOf[AFO_PT] = kTRUE;} 
   else if(TString(kc).EqualTo("eta")){fDoNotCalculateCorrelationsAsFunctionOf[AFO_ETA] = kTRUE;} 
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
  }

  void SetCalculateTest0(Bool_t c) {this->fCalculateTest0 = c;};
  Bool_t GetCalculateTest0() const {return this->fCalculateTest0;};
  void SetFileWithLabels(const char *externalFile) 
  {
   this->fFileWithLabels = new TString(externalFile); 
   if(gSystem->AccessPathName(fFileWithLabels->Data(),kFileExists)) // convention is opposite to Bash
   {
    cout<<Form("The file %s doesn't exist",externalFile)<<endl;
    cout<<__LINE__<<endl; exit(1);
   }
   this->StoreLabelsInPlaceholder("external");
  }

  void SetTest0List(TList* const list) {this->fTest0List = list;};
  TList* GetTest0List() const {return this->fTest0List;} 
  TProfile* GetTest0Pro(const Int_t order, const Int_t index, const Int_t var) {return this->fTest0Pro[order][index][var];}

  void SetCalculateNestedLoops(Bool_t cnl) {this->fCalculateNestedLoops = cnl;};
  Bool_t GetCalculateNestedLoops() const {return this->fCalculateNestedLoops;};

  void SetCalculateCustomNestedLoop(Bool_t ccnl) {this->fCalculateCustomNestedLoop = ccnl;};
  Bool_t GetCalculateCustomNestedLoop() const {return this->fCalculateCustomNestedLoop;};

  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  // Particle weights:
  void SetWeightsHist(TH1D* const hist, const char *variable); // .cxx
  TH1D* GetWeightsHist(const char *variable); // . cxx
  TH1D* GetHistogramWithWeights(const char *filePath, const char *variable); // .cxx
 
  // Centrality weights:
  void SetCentralityWeightsHist(TH1D* const hist); // .cxx
  TH1D* GetCentralityWeightsHist() const {return this->fCentralityWeightsHist;} 
  TH1D* GetHistogramWithCentralityWeights(const char *filePath, const char *estimator); // .cxx

  // Toy NUA:
  void SetToyNUACuts(const char* kc, const Double_t probability, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(kc).EqualTo("phi")){var = PHI;} 
   else if (TString(kc).EqualTo("pt")){var = PT;} 
   else if (TString(kc).EqualTo("eta")){var = ETA;}
   else if (TString(kc).EqualTo("e")){var = E;}
   else if (TString(kc).EqualTo("charge")){var = CHARGE;}
   else { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); }
   if(probability<0.||probability>1.){exit(1);}
   if(fToyNUACuts[var][1]>fToyNUACuts[var][2]){exit(1);}
   this->fToyNUACuts[var][0] = probability;
   this->fToyNUACuts[var][1] = min;
   this->fToyNUACuts[var][2] = max;
   this->fUseToyNUA[var] = kTRUE;
  }

  // Internal validation:
  void SetUseInternalValidation(Bool_t uiv, Int_t nEvts, Bool_t rescale) 
  {
   this->fUseInternalValidation = uiv; 
   this->fnEventsInternalValidation = nEvts;
   this->fRescaleWithTheoreticalInput = rescale;
  };
  void SetMultRangeInternalValidation(Int_t min, Int_t max) 
  {
   this->fMultRangeInternalValidation[0] = min; 
   this->fMultRangeInternalValidation[1] = max;
  };

  void SetHarmonicsOptionInternalValidation(const char* kc)
  {
   // If the setter is not called, fHarmonicsOptionInternalValidation is defaulted to "constant".
   // Supported options:
   // a) "constant": vn and psin are constant event-by-event, their values are set via SetInternalValidationAmplitudes and SetInternalValidationPlanes;
   // b) "correlated": vn's and psin's are sampled correlated event-by-event, see .cxx how pdf is implemented
   *(this->fHarmonicsOptionInternalValidation) = kc; // I can do it this way, since  DefaultConfiguration() is called in the constructor, and there fHarmonicsOptionInternalValidation is defined 
                                                     // as new TString("constant") already 
   if(! (fHarmonicsOptionInternalValidation->EqualTo("constant") || fHarmonicsOptionInternalValidation->EqualTo("correlated")) )
   { Red(__PRETTY_FUNCTION__); cout<<__LINE__<<endl; exit(1); } // bail out, if non-supported options is provided
  } // void SetHarmonicsOptionInternalValidation(const char* kc)
  void SetInternalValidationAmplitudes(TArrayD *iva){this->fInternalValidationAmplitudes = iva;};
  void SetInternalValidationPlanes(TArrayD *ivp){this->fInternalValidationPlanes = ivp;};

  // Utility:
  void Red(const char* text);
  void Green(const char* text);
  void Yellow(const char* text);
  void Blue(const char* text);
  TObject* GetObjectFromList(TList *list, Char_t *objectName); // see .cxx
  Int_t NumberOfNonEmptyLines(const char *externalFile);  

  // *.) Online monitoring:
  void SetUpdateOutputFile(const Int_t uf, const char *uqof)
  {
   // Example usage: task->SetUpdateOutputFile(44,"AnalysisResultsProgress.root");
   this->fOnlineMonitoring = kTRUE;
   this->fUpdateOutputFile = kTRUE;
   this->fUpdateFrequency = uf;
   this->fUpdateFile = new TString(uqof);
  };
  void SetMaxNumberOfEvents(const Int_t mnof, const char *uqof)
  {
   // Example usage: task->SetMaxNumberOfEvents(44,"AnalysisResultsBailOut.root");
   this->fOnlineMonitoring = kTRUE;
   this->fMaxNumberOfEvents = mnof;
   this->fBailOutFile = new TString(uqof);
  };

  // *.) Debugging:
  void SetProcessOnlySpecifiedEvent(Int_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period) 
  {
   this->fProcessOnlySpecifiedEvent = kTRUE; 
   this->fRun = run;
   this->fBunchCross = bunchCross;
   this->fOrbit = orbit;  
   this->fPeriod = period;
  }; // void SetProcessOnlySpecifiedEvent(UInt_t run, UShort_t bunchCross, UInt_t orbit, UInt_t period)
  virtual void SetPrintEventInfo(){this->fPrintEventInfo = kTRUE;}; // print event medatata (for AOD: fRun, fBunchCross, fOrbit, fPeriod)

 private:
  AliAnalysisTaskMuPa(const AliAnalysisTaskMuPa& aatmpf);
  AliAnalysisTaskMuPa& operator=(const AliAnalysisTaskMuPa& aatmpf);
  
  // 0) Base list:
  TList *fBaseList; // base list to hold all output object (a.k.a. grandmother of all lists)
  TProfile *fBasePro; // keeps flags relevant for the whole analysis
  Bool_t fRealData; // use SetRealData(...); in the steering macro
  Bool_t fUseFisherYates; // use SetUseFisherYates(kTRUE); in the steering macro to randomize particle indices
  TArrayI *fRandomIndices; // array to store random indices obtained from Fisher-Yates algorithm 
  TString fTaskName; // e.g. Form("Task=>%.1f-%.1f",centrMin,centrMax)
  TString fDataTakingPeriod; // the data taking period, use e.g. task->SetDataTakingPeriod("LHC10h")
  TString fAODNumber; // the AOD number, use e.g. task->SetAODNumber("AOD160")
  TString fRunNumber; // the run number, use e.g. task->SetRunNumber("000123456"). For sim, strip off 000.
  Bool_t fFillQAHistograms; // fill QA histograms (this shall be done only in one task, since these histos are heavy 2D objects). Additional loops over particles is performed.
  Bool_t fFillQAHistogramsAll; // if kFALSE, only most important QA histograms a filled
  Bool_t fTerminateAfterQA; // in UserExec(), bail out immediately after QA histograms are filled 
  Bool_t fVerbose; // print all additional info like Green(__PRETTY_FUNCTION__); etc.
  Int_t fEventCounter; // counter of all events, i.e. number of times UserExec() has been called
  UInt_t fRandomSeed; // argument to TRandom3 constructor. By default is 0, use SetRandomSeed(...) to change it
  TString fTrigger; // offline trigger, use e.g. task->SetTrigger("kMB")
  Bool_t fUseTrigger; // kFALSE by default. Set automatically when task->SetTrigger(...) is called
  Bool_t fUseFixedNumberOfRandomlySelectedParticles; // use or not fixed number of randomly selected particles in each event. Use always in combination with SetUseFisherYates(kTRUE)
  Int_t fFixedNumberOfRandomlySelectedParticles; // set here a fixed number of randomly selected particles in each event. Use always in combination with SetUseFisherYates(kTRUE)
  Bool_t fHistogramBookingsWithRunInfoWereUpdated; // makes sure that UpdateHistogramBookingsWithRunInfo() is called only once in UserExec()

  // 1) QA:
  TList *fQAList; // base list to hold all QA output object
  TProfile *fQAPro; // keeps flags relevant for the QA analysis
  TH1D *fQACentralityHist[gCentralityEstimators][2]; // centrality distribution [all supported estimators][before, after event cuts]
  TH2D *fQACentralityCorrHist[gCentralityEstimators][gCentralityEstimators][2]; // correlations of centrality distributions from all supported estimators [before, after event cuts]
  TH2D *fQAMultiplicityCorrHist[2]; // correlations between reference multiplicity and nSelected tracks [before, after event cuts]
  TH2D *fQAGenericCorrHist[gGenericCorrelations][2]; // correlations between various quantities (see .cxx for documentation) [before, after event cuts]
  TH1I *fQAFilterBitScan; // for each track in AOD, dump its filterbits
  TH2I *fQAIDvsFilterBit; // filterbit vs. atrack->ID()
  TH1D *fQAKinematicsFilterBits[gFilterBits][gKinematicVariables]; // kinematics [specified filter bit][phi,pt,eta,energy,charge] Use in combination with SetQAFilterBits(...)
  TArrayI *fQAFilterBits; // for these filterbits the kinematics in the previous line will be filled, use in combination with SetQAFilterBits(...)
  TH1I *fQAAnomalousEvents; // counter for anomalous events
  TH1D *fQASelfCorrelations[gQASelfCorrelations]; // check for self-correlations
  TH1D *fQASimRecoSelfCorrelations[gQASelfCorrelations]; // check for self-correlations between simulated and reconstructed particles
  Bool_t fQACheckSelfCorrelations; // kFALSE by default. If kTRUE, both fQASelfCorrelations[] and fQASimRecoSelfCorrelations[] are filled
  TH1I *fQAEventCutCounter; // counter for each event cut
  TH1I *fQASequentialEventCutCounter; // sequential counter for event cuts
  TH1I *fQAParticleCutCounter[2]; // counter for each particle cut. [0] = reco, [1] = sim
  TH1I *fQASequentialParticleCutCounter[2]; // sequential counter for each particle cut. [0] = reco, [1] = sim
  TH1I *fQATrigger[2]; // counter for triggers [0] = before event cuts, [1] = after event cuts

  // 2) Control event histograms:  
  TList *fControlEventHistogramsList; // list to hold all control event histograms
  TProfile *fControlEventHistogramsPro; // keeps flags relevant for the control event histograms

  //    Multiplicities:
  Double_t fMultiplicity;        // defined as a sum of track weights used to calculate Q-vectors (see below also fSelectedTracks) 
  TH1D *fMultiplicityHist;       // this is distribution of my multiplicity
  Double_t fMultiplicityBins[3]; // [nBins,minMultiplicity,maxMultiplicity]
  Int_t fSelectedTracks;         // this is an integer counter of tracks used to calculate Q-vectors, after all particle cuts have been applied. Calculated in FilterEvent(AliVEvent *ave)
  TH1I *fSelectedTracksHist;     // this is distribution fSelectedTracks. Also the counter of selected events for analysis (use e.g. in OnlineMonitoring() or InternalValidation())
  Int_t fSelectedTracksCuts[2];  // [min,max]
  Bool_t fUseSelectedTracksCuts; // kFALSE by default, set to kTRUE if task->SelectedTracksCuts(...) has been called
  TH1D *fCentralMultiplicityHist[gCentralMultiplicity][2]; // this is distribution of reference multiplicities, maintained centrally [before,after event cuts]
  Bool_t fUseGenericCorrelationsCuts[gGenericCorrelations]; // each of them is automatically set to kTRUE when setter is called, e.g. [0] with SetSelContrTreshold(...), etc.
  Double_t fSelContrTreshold;    // use via setter task->SetSelContrTreshold(...), applied as: fSelectedTracks > fSelContrTreshold * (Int_t)avtx->GetNContributors()

  //    Centrality:
  Double_t fCentrality;         // this is ebe centrality from default estimator
  TH1D *fCentralityHist[2];     //! centrality distribution from default estimator [before,after event cuts]
  Double_t fCentralityBins[3];  // [nBins,minCentrality,maxCentrality]
  TString fCentralityEstimator; // the default centrality estimator in this analysis, use e.g. task->SetCentralityEstimator("V0M")
  Double_t fCentralityCuts[2];  // [min,max]
  Bool_t fUseCentralityCuts;    // kFALSE by default, set to kTRUE if task->SetCentralityCuts(...) is called
  Bool_t fUseCentralityCorrelationsCuts[gCentralityEstimators][gCentralityEstimators]; // use specific cut below 
  Double_t fCentralityCorrelationsCuts[gCentralityEstimators][gCentralityEstimators]; // [firstEstimator][secondEstimator] |(firstEstimator-secondEstimator)/(firstEstimator+secondEstimator)| > cut => reject the event
  Int_t fCentralityCorrelationCutVersion; // see .cxx switch(fCentralityCorrelationCutVersion) for explanatio
  //    Vertex:
  TH1D *fVertexHist[2][2][3];      //! distribution of vertex components [before,after event cuts][reco,sim][x,y,z]
  Double_t fVertexBins[3];         // [nBins,minVertex,maxVertex]
  Double_t fVertexCuts[3][2];      // [x,y,z][min,max] vertex components
  Bool_t fUseVertexCuts[3];        // [x,y,z] use or not vertex cuts. Set the kTRUE only if task->SetVertexCuts(...) is called explicitly in steering macros
  TH1I *fNContributorsHist[2][2];  //! distribution of vertex components [before,after event cuts][reco,sim][x,y,z]
  Int_t fNContributorsBins[3];     // [nBins,min,max]
  Int_t fNContributorsCuts[2];     // [min,max]
  Bool_t fUseNContributorsCuts;    // by default kFALSE, set to kTRUE if task->SetNContributorsCuts(...) is called in steering macros
  Double_t fMinVertexDistance;     // if sqrt(vx^2+vy^2+vz^2) < fMinVertexDistance, the event is reject. This way, I remove suspicious events with |vertex| = 0.
  Bool_t fUseMinVertexDistanceCut; // Set the kTRUE only if task->SetMinVertexDistance(...) is called explicitly in steering macros

  //    All remaining event histograms: 
  TH1D *fEventHistograms[2][gEventHistograms]; //! [before,after event cuts][ type - see enum ]
  Double_t fEventBins[gEventHistograms][3]; // [nBins,min,max]
  Double_t fEventCuts[gEventHistograms][2]; // [type - see enum][min,max]
  Bool_t fUseEventCuts[gEventHistograms];   // if not set via setter, corresponding cut is kFALSE. Therefore, correspondig cut is open (default values are NOT used)

  // 3) Control particle histograms:  
  TList *fControlParticleHistogramsList; // list to hold all control histograms for particle distributions 
  TProfile *fControlParticleHistogramsPro; // keeps flags relevant for the control particle histograms
  Bool_t fFillControlParticleHistograms; // kTRUE by default, but for instance kFALSE when I do InternalValidation()
  TExMap *fSimReco; //! look up table between kine and reco particles (key = kine (Monte Carlo label), value = reco (track index in AOD))
  Bool_t fUseFakeTracks; // if kTRUE, the Monte Carlo particle is obtained as TMath:Abs(aRecoTrack->GetLabel())
  TExMap *fGlobalTracksAOD; //! global tracks in AOD
  Bool_t fFilterGlobalTracksAOD; // by default kFALSE, set to kTRUE when task->SetFilterGlobalTracksAOD(); is used. Neded only for PID studies, setting FilterBit avoids double-counting
  Int_t fFilterBit; // filter bit (its meaning can change from one production to another)
  Bool_t fUseOnlyPrimaries; // cut e.g. on AliAODTrack::kPrimary or aodmcParticle->IsPhysicalPrimary()
  TString fPrimaryDefinitionInMonteCarlo; // supported: "IsPhysicalPrimary" (default), "IsPrimary", ... Set via task->SetPrimaryDefinitionInMonteCarlo("...")

  //    Kinematics:
  TH1D *fKinematicsHist[2][2][gKinematicVariables]; // kinematics [before,after track cuts][reco,sim][phi,pt,eta,energy,charge]
  Double_t fKinematicsBins[gKinematicVariables][3]; // [phi,pt,eta,energy,charge][nBins,min,max]
  Double_t fKinematicsCuts[gKinematicVariables][2]; // [phi,pt,eta,energy,charge][min,max]
  Bool_t fUseKinematicsCuts[gKinematicVariables];   // if not set via setter, corresponding cut is kFALSE. Therefore, correspondig cut is open (default values are NOT used)
  //    DCA:
  TH1D *fDCAHist[2][2][2]; // distance of closest approach (DCA) [before,after][reco,sim][xy,z] // "xy" is transverse direction
  Double_t fDCABins[2][3]; // [xy,z][nBins,min,max]
  Double_t fDCACuts[2][2]; // [xy,z][min,max]
  Bool_t fUseDCACuts[2];   // [xy,z] if not set via setter, corresponding cut is kFALSE. Therefore, correspondig cut is open (default values are NOT used)
  Bool_t fUsePtDependentDCAxyParameterization; // use instead the pT dependent cut for DCAxy, see the setter SetPtDependentDCAxyCut(const char* parameterization)
  TString fPtDependentDCAxyParameterization; // math. parameterization for DCAxy vs. pT, see the setter SetPtDependentDCAxyCut(const char* parameterization)
  TFormula *fPtDependentDCAxyFormula; // the actual formula, used to evaluate for a given pT, the corresponding DCAxy, where the parameterization is given by fPtDependentDCAxyParameterization
  //    All remaining particle histograms: 
  TH1D *fParticleHist[2][2][gParticleHistograms]; //! distributions [before,after particle cuts][reco,sim][type - see enum]
  Double_t fParticleBins[gParticleHistograms][3]; // [nBins,min,max]
  Double_t fParticleCuts[gParticleHistograms][2]; // [type - see enum][min,max]
  Bool_t fUseParticleCuts[gParticleHistograms];   // set to kTRUE, only if the correspondign SetParticleCuts(...) was used in the steering macros
  Bool_t fAtLeastOnePointInTheSPD; // set to kTRUE via SetAtLeastOnePointInTheSPD( ... ), only tracks with one or two points in the SPD are taken
  Bool_t fIgnoreGlobalConstrained; // set to kTRUE by default, to avouid double counting in some cases.

  // 4) Q-vectors:
  TList *fQvectorList;        // list to hold all Q-vector objects       
  TProfile *fQvectorFlagsPro; // profile to hold all flags for Q-vector
  Bool_t fCalculateQvector;   // to calculate or not to calculate Q-vectors, that's a Boolean...
  Int_t fMaxHarmonic;         // 6 (not going beyond v6, if you change this value, change also fQvector[49][9]) 
  Int_t fMaxCorrelator;       // 8 (not going beyond 8-p correlations, if you change this value, change also fQvector[49][9]) 
  //TComplex fQvector[49][9];   //! Q-vector components [fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*8+1][8+1]  
  TComplex fQ[gMaxHarmonic*gMaxCorrelator+1][gMaxCorrelator+1]; //! generic Q-vector [fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*12+1][12+1]  
  TComplex fQvector[gMaxHarmonic*gMaxCorrelator+1][gMaxCorrelator+1]; //! "integrated" Q-vector [fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*12+1][12+1]  
  TComplex fqvector[2][gMaxNoBinsKine][gMaxHarmonic*gMaxCorrelator+1][gMaxCorrelator+1]; //! "differenttial" q-vector [kine var.][binNo][fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*12+1][12+1]  
  Int_t fqVectorEntries[2][gMaxNoBinsKine]; // count number of entries in each differential q-vector

  // 5) Particle weights: 
  TList *fWeightsList;          // list to hold all Q-vector objects       
  TProfile *fWeightsFlagsPro;   // profile to hold all flags for weights
  Bool_t fUseWeights[gWeights]; // use weights [phi,pt,eta]
  TH1D *fWeightsHist[gWeights]; // histograms holding weights [phi,pt,eta]

  // 6) Centrality weights: 
  TList *fCentralityWeightsList;        // list to hold all Q-vector objects       
  TProfile *fCentralityWeightsFlagsPro; // profile to hold all flags for CentralityWeights
  Bool_t fUseCentralityWeights;         // use centrality weights [V0M, SPDTracklets, CL0, CL1]
  TH1D *fCentralityWeightsHist;         // histograms holding centrality weights [V0M, SPDTracklets, CL0, CL1]

  // 7) Correlations:
  TList *fCorrelationsList;                                        // list to hold all correlations objects
  TProfile *fCorrelationsFlagsPro;                                 // profile to hold all flags for correlations
  Bool_t fCalculateCorrelations;                                   // calculate and store all correlations (see enum eAsFunctionOf for what is supported at the moment)
  Bool_t fDoNotCalculateCorrelationsAsFunctionOf[eAsFunctionOf_N]; // decide differentially which correlations shall NOT be calculated (by default, all are calculated)
  TProfile *fCorrelationsPro[4][gMaxHarmonic][5];                  //! multiparticle correlations [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...][0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta]
  Bool_t fUseCustomKineDependenceBins[gKineDependenceVariables];   // use or not custom binning for kine dependence => use setter SetKineDependenceBins(...)
  TArrayD *fKineDependenceBins[gKineDependenceVariables];          // custom binning kine dependence, SetKineDependenceBins(...). By default the same binning is used as in the corresponding control histograms. 

  // 8) Nested loops:
  TList *fNestedLoopsList;                       // list to hold all nested loops objects 
  TProfile *fNestedLoopsFlagsPro;                // profile to hold all flags for nested loops
  Bool_t fCalculateNestedLoops;                  // calculate and store correlations with nested loops, as a cross-check
  Bool_t fCalculateCustomNestedLoop;             // validate e-b-e all correlations with custom nested loop
  TProfile *fNestedLoopsPro[4][gMaxHarmonic][5]; //! multiparticle correlations from nested loops [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...][0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta]
  //TProfile *fNestedLoopsPerDemandPro[3];       // which correlator needs to be cross-checked with nested loops (no setter => recompile). [0=integrated,1=vs. multiplicity,2=vs. centrality]
  TArrayD *ftaNestedLoops[2];                    //! e-b-e container for nested loops [0=angles;1=product of all weights]   
  TArrayD *ftaNestedLoopsKine[gKineDependenceVariables][gMaxNoBinsKine][2]; //! e-b-e container for nested loops [0=pT,1=eta][kine.bin][0=angles;1=product of all weights]   

  // 9) Toy NUA:
  TList *fToyNUAList;                           // list to hold all Toy NUA objects
  TProfile *fToyNUAFlagsPro;                    // profile to hold all flags for Toy NUA
  Bool_t fUseToyNUA[gKinematicVariables];       // use toy NUA for particular kinematic variable
  Double_t fToyNUACuts[gKinematicVariables][3]; // stores probability [0] and NUA sector range min [1] and max [2]. Use task->SetToyNUACuts("variable",probability,min,max)

  //10) Internal validation:
  TList *fInternalValidationList;              // list to hold all objects for internal validation
  TProfile *fInternalValidationFlagsPro;       // profile to hold all flags for internal validation
  Bool_t fUseInternalValidation;               // use internal validation
  Bool_t fRescaleWithTheoreticalInput;         // if kTRUE, all measured correlators are rescaled with theoretical input, so that in profiles everything is at 1
  Int_t fnEventsInternalValidation;            // how many events will be sampled on-the-fly for internal validation
  TString *fHarmonicsOptionInternalValidation; // see .cxx for full documentation 
  TArrayD *fInternalValidationAmplitudes;      // 0 = v1, 1 = v2, etc.
  TArrayD *fInternalValidationPlanes;          // 0 = Psi1, 1 = Psi2, etc.
  Int_t fMultRangeInternalValidation[2];       // min and max values for uniform multiplicity distribution in on-the-fly analysis

  //11) Test0:  
  TList *fTest0List; // list to hold all objects for Test0
  TProfile *fTest0FlagsPro; // store all flags for Test0
  Bool_t fCalculateTest0; // calculate or not Test0 in general. Which one specifically, that's governed with fCalculateSpecificTest0[gMaxCorrelator][gMaxIndex]
  Bool_t fTest0LabelsWereStoredInPlaceholder; // Test0 labels were stored successfully in fTest0LabelsPlaceholder. From there, they will be extracted at run-time when booking
  TProfile *fTest0Pro[gMaxCorrelator][gMaxIndex][5]; //! [gMaxCorrelator][gMaxIndex][3] [order][index][0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta]
  TString *fFileWithLabels; // external file which specifies all labels of interest
  TString *fTest0Labels[gMaxCorrelator][gMaxIndex]; // all labels: k-p'th order is stored in k-1'th index. So yes, I also store 1-p
  TH1I *fTest0LabelsPlaceholder; // temporary workaround: store all Test0 labels in this histogram, until I find a better implementation

  // * Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results

  // * Common style:
  Int_t fBeforeAfterColor[2]; //! [0 = kRed,1 = kGreen] TBI 20210511 shall I use better enum here?

  // *.) Online monitoring:
  Bool_t fOnlineMonitoring; // enable online monitoring (not set excplicitly!), the flags below just refine it
  Bool_t fUpdateOutputFile; // update the output file after certain number of analysed events
  Int_t fUpdateFrequency;   // after how many events the output file will be updated
  TString *fUpdateFile;     // which file will be regularly updated with frequency fUpdateFrequency of events
  Int_t fMaxNumberOfEvents; // if this number of events is reached, write to external file fBailOutFile and bail out
  TString *fBailOutFile;    // in which file results will be bailed out after fMaxNumberOfEvents

  // *.) Debugging:
  Bool_t fProcessOnlySpecifiedEvent; // process only for the specified event
  Int_t fRun;                        // run number
  UShort_t fBunchCross;              // bunch crossing
  UInt_t fOrbit;                     // orbit
  UInt_t fPeriod;                    // period
  Bool_t fPrintEventInfo;            // print event medatata (for AOD: fRun, fBunchCross, fOrbit, fPeriod). Enabled indirectly via task->PrintEventInfo() 
 
  // Increase this counter in each new version:
  ClassDef(AliAnalysisTaskMuPa,36);

};

//================================================================================================================

#endif




