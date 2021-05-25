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
#include <AliAODTrack.h>
#include <AliAODEvent.h>
#include <AliVEvent.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TExMap.h>
#include <TComplex.h>
#include <TArrayD.h>
#include <TArrayI.h>
#include <TSystem.h>

// Global variables:
const Int_t gCentralityEstimators = 4; // set here number of supported centrality estimators
const Int_t gKinematicVariables = 5; // number of supported kinematic variables: [phi,pt,eta,e,charge]
const Int_t gFilterBits = 17; // number of filterbits to scan
const Int_t gEventHistograms = 2; // total number of non-classified event histograms
const Int_t gParticleHistograms = 9; // total number of non-classified particle histograms
const Int_t gCentralMultiplicity = 1; // multiplicities defined centrally, e.g. ref. mult.
const Int_t gWeights = 3; // phi, pt, eta

// enums:
enum eBins {nBins,min,max};
enum eCentralityEstimator { V0M, SPDTracklets, CL0, CL1 }; 
enum eXYZ { X = 0, Y = 1, Z = 2 };
enum eBeforeAfter { BEFORE = 0, AFTER = 1 };
enum eRecoSim { RECO = 0, SIM = 1 };
enum eKinematics { PHI = 0, PT = 1, ETA = 2, E = 3, CHARGE = 4 };
enum eDefaultColors { COLOR = kBlack, FILLCOLOR = kGray };
enum eParticle { TPCNcls, TPCnclsS, TPCnclsFractionShared, TPCNCrossedRows, TPCChi2perNDF, TPCFoundFraction, Chi2TPCConstrainedVsGlobal, ITSNcls, ITSChi2perNDF };
enum eEvent { MagneticField, PrimaryVertex };

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
   virtual void InitializeArraysForCorrelationsHistograms();
   virtual void InitializeArraysForNestedLoopsHistograms();
   virtual void InitializeArraysForCommonLabels();

  // 1) Methods called in UserCreateOutputObjects():
  virtual void InsanityChecks(); 
  virtual void BookBaseProfile();
  virtual void BookAndNestAllLists();  
  virtual void BookQAHistograms();
  virtual void BookControlEventHistograms();
  virtual void BookControlParticleHistograms();
  virtual void BookQvectorHistograms();
  virtual void BookWeightsHistograms();
  virtual void BookCorrelationsHistograms();
  virtual void BookNestedLoopsHistograms();
  virtual void BookFinalResultsHistograms();

  // 2) Methods called in UserExec(Option_t *):
  virtual void FillQAHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs);
  virtual void FilterEvent(AliVEvent *ave);
  virtual void FillControlEventHistograms(AliVEvent *ave, const Int_t ba, const Int_t rs); // before or after event cuts, reco or sim
  virtual void FillControlParticleHistograms(AliAODTrack *aTrack, const Int_t ba, const Int_t rs); // before or after particle cuts, reco or sim
  virtual void GlobalTracksAOD(AliAODEvent *aAOD);
  Bool_t SurvivesEventCuts(AliVEvent *ave);
  Bool_t SurvivesParticleCuts(AliAODTrack *aTrack); // applied e.g. on TPC-only
  virtual Double_t Weight(const Double_t &value, const char *variable);
  virtual void CalculateCorrelations();
  virtual void CalculateNestedLoops();
  virtual void ResetEventByEventQuantities();
  virtual void OnlineMonitoring();

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

  // 5) Setters and getters:
  void SetPeriod(const char *p) {this->fPeriod = p;};
  void SetAODNumber(const char *an) {this->fAODNumber = an;};

  void SetControlEventHistogramsList(TList* const cehl) {this->fControlEventHistogramsList = cehl;};
  TList* GetControlEventHistogramsList() const {return this->fControlEventHistogramsList;} 
  // QA:
  void SetFillQAHistograms(Bool_t fqah) {this->fFillQAHistograms = fqah;};
  Bool_t GetFillQAHistograms() const {return this->fFillQAHistograms;};
  void SetTerminateAfterQA(Bool_t taqa) {this->fTerminateAfterQA = taqa;};
  Bool_t GetTerminateAfterQA() const {return this->fTerminateAfterQA;};
  void SetQAFilterBits(TArrayI *fb){this->fQAFilterBits = fb;};

  // Multiplicity:
  void SetMultiplicityBins(Int_t const nbins, Double_t min, Double_t max)
  {
   this->fMultiplicityBins[0] = nbins;
   this->fMultiplicityBins[1] = min;
   this->fMultiplicityBins[2] = max;
  };
  void SetMultiplicityCuts(const Double_t min, const Double_t max)
  {
   this->fMultiplicityCuts[0] = min;
   this->fMultiplicityCuts[1] = max;
  }

  void SetSelectedTracksCuts(const Int_t min, const Int_t max)
  {
   this->fSelectedTracksCuts[0] = min;
   this->fSelectedTracksCuts[1] = max;
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
  }
  void SetCentralityEstimator(const char *ce) {this->fCentralityEstimator = ce;};

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
   this->fVertexCuts[xyz][0] = min;
   this->fVertexCuts[xyz][1] = max;
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
  }

  // Remaining event distributions:
  void SetEventBins(const char* type, const Double_t nbins, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("MagneticField")){var = MagneticField;} 
   else if (TString(type).EqualTo("PrimaryVertex")){var = PrimaryVertex;}
   this->fEventBins[var][0] = nbins;
   this->fEventBins[var][1] = min;
   this->fEventBins[var][2] = max;
  }
  void SetEventCuts(const char* type, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("MagneticField")){var = MagneticField;} 
   else if (TString(type).EqualTo("PrimaryVertex")){var = PrimaryVertex;}
   this->fEventCuts[var][0] = min;
   this->fEventCuts[var][1] = max;
  }
  void SetControlParticleHistogramsList(TList* const cphl) {this->fControlParticleHistogramsList = cphl;};
  TList* GetControlParticleHistogramsList() const {return this->fControlParticleHistogramsList;} 
  void SetFilterBit(Int_t fb) {this->fFilterBit = fb;};
  Int_t GetFilterBit() const {return this->fFilterBit;};
 
  // Kinematics:
  void SetKinematicsBins(const char* kv, const Double_t nbins, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(kv).EqualTo("phi")){var = PHI;} 
   else if (TString(kv).EqualTo("pt")){var = PT;} 
   else if (TString(kv).EqualTo("eta")){var = ETA;}
   else if (TString(kv).EqualTo("e")){var = E;}
   else if (TString(kv).EqualTo("charge")){var = CHARGE;}
   this->fKinematicsBins[var][0] = nbins;
   this->fKinematicsBins[var][1] = min;
   this->fKinematicsBins[var][2] = max;
  }

  void SetKinematicsCuts(const char* kc, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(kc).EqualTo("phi")){var = PHI;} 
   else if (TString(kc).EqualTo("pt")){var = PT;} 
   else if (TString(kc).EqualTo("eta")){var = ETA;}
   else if (TString(kc).EqualTo("e")){var = E;}
   else if (TString(kc).EqualTo("charge")){var = CHARGE;}
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
   this->fDCABins[var][0] = nbins;
   this->fDCABins[var][1] = min;
   this->fDCABins[var][2] = max;
  }
  void SetDCACuts(const char* dc, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(dc).EqualTo("xy")){var = 0;} 
   else if (TString(dc).EqualTo("z")){var = 1;} 
   this->fDCACuts[var][0] = min;
   this->fDCACuts[var][1] = max;
   //this->fUseDCACuts[var] = kTRUE; TBI 20210517 most likely obsolete
  }

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
   this->fParticleCuts[var][0] = min;
   this->fParticleCuts[var][1] = max;
  }
  void SetCalculateQvector(Bool_t cqv) {this->fCalculateQvector = cqv;};
  Bool_t GetCalculateQvector() const {return this->fCalculateQvector;};

  void SetCalculateCorrelations(Bool_t cc) {this->fCalculateCorrelations = cc;};
  Bool_t GetCalculateCorrelations() const {return this->fCalculateCorrelations;};

  void SetCalculateNestedLoops(Bool_t cnl) {this->fCalculateNestedLoops = cnl;};
  Bool_t GetCalculateNestedLoops() const {return this->fCalculateNestedLoops;};

  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  // Particle weights:
  void SetWeightsHist(TH1D* const hist, const char *variable); // .cxx
  TH1D* GetHistogramWithWeights(const char *filePath, const char *variable); // .cxx

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

 private:
  AliAnalysisTaskMuPa(const AliAnalysisTaskMuPa& aatmpf);
  AliAnalysisTaskMuPa& operator=(const AliAnalysisTaskMuPa& aatmpf);
  
  // 0) Base list:
  TList *fBaseList; // base list to hold all output object (a.k.a. grandmother of all lists)
  TProfile *fBasePro; // keeps flags relevant for the whole analysis
  TString fPeriod; // the data taking period, use e.g. task->SetPeriod("LHC10h")
  TString fAODNumber; // the AOD number, use e.g. task->SetPeriod("AOD160") (for "LHC10h")
  Bool_t fFillQAHistograms; // fill all QA histograms (this shall be done only in one task, since these histos are heavy 2D objects). Additional loop over particles is performed.
  Bool_t fTerminateAfterQA; // in UserExec(), bail out immediately after QA histograms are filled 

  // 1) QA:
  TList *fQAList; // base list to hold all QA output object
  TProfile *fQAPro; // keeps flags relevant for the QA analysis
  TH1D *fQACentralityHist[gCentralityEstimators][2]; // centrality distribution [all supported estimators][before, after event cuts]
  TH2D *fQACentralityCorrHist[gCentralityEstimators][gCentralityEstimators][2]; // correlations of centrality distributions from all supported estimators [before, after event cuts]
  TH1I *fQAFilterBitScan; // for each track in AOD, dump its filterbits
  TH2I *fQAIDvsFilterBit; // filterbit vs. atrack->ID()
  TH1D *fQAKinematicsFilterBits[gFilterBits][2][gKinematicVariables]; // kinematics [specified filter bit][reco,sim][phi,pt,eta,energy,charge] Use in combination with SetQAFilterBits(...)
  TArrayI *fQAFilterBits; // for these filterbits the kinematics in the previous line will be filled, use in combination with SetQAFilterBits(...)

  // 2) Control event histograms:  
  TList *fControlEventHistogramsList; // list to hold all control event histograms
  TProfile *fControlEventHistogramsPro; // keeps flags relevant for the control event histograms
  //    Multiplicity:
  Double_t fMultiplicity; // this is my ebe multiplicity, i.e. sum of track weights used to calculate Q-vectors (see below also fSelectedTracks) 
  TH1D *fMultiplicityHist; // this is distribution of my multiplicity
  TH1D *fCentralMultiplicityHist[gCentralMultiplicity]; // this is distribution of my multiplicity
  Double_t fMultiplicityBins[3]; // [nBins,minCentrality,maxCentrality]
  Double_t fMultiplicityCuts[2]; // [min,max]
  Int_t fSelectedTracks;         // this is counter of tracks used to calculate Q-vectors (see above also fMultiplicity) 
  TH1I *fSelectedTracksHist;     // this is distribution fSelectedTracks
  Int_t fSelectedTracksCuts[2];  // [min,max]
  //    Centrality:
  Double_t fCentrality;         // this is ebe centrality from default estimator
  TH1D *fCentralityHist[2];     //! centrality distribution from default estimator [before,after event cuts]
  Double_t fCentralityBins[3];  // [nBins,minCentrality,maxCentrality]
  TString fCentralityEstimator; // the default centrality estimator in this analysis, use e.g. task->SetCentralityEstimator("V0M")
  Double_t fCentralityCuts[2];  // [min,max]
  //    Vertex:
  TH1D *fVertexHist[2][2][3];     //! distribution of vertex components [before,after event cuts][reco,sim][x,y,z]
  Double_t fVertexBins[3];        // [nBins,minVertex,maxVertex]
  Double_t fVertexCuts[3][2];     // [x,y,z][min,max] vertex components
  TH1I *fNContributorsHist[2][2]; //! distribution of vertex components [before,after event cuts][reco,sim][x,y,z]
  Int_t fNContributorsBins[3];    // [nBins,min,max]
  Int_t fNContributorsCuts[2];    // [min,max]

  //    All remaining event histograms: 
  TH1D *fEventHistograms[2][gEventHistograms]; //! [before,after event cuts][ type - see enum ]
  Double_t fEventBins[gEventHistograms][3]; // [nBins,min,max]
  Double_t fEventCuts[gEventHistograms][2]; // [type - see enum][min,max]

  // 3) Control particle histograms:  
  TList *fControlParticleHistogramsList; // list to hold all control histograms for particle distributions
  TProfile *fControlParticleHistogramsPro; // keeps flags relevant for the control particle histograms
  TExMap *fGlobalTracksAOD; //! global tracks in AOD
  Int_t fFilterBit; // filter bit (its meaning can change from one production to another)
  //    Kinematics:
  TH1D *fKinematicsHist[2][2][gKinematicVariables]; // kinematics [before,after track cuts][reco,sim][phi,pt,eta,energy,charge]
  Double_t fKinematicsBins[gKinematicVariables][3]; // [phi,pt,eta,energy,charge][nBins,min,max]
  Double_t fKinematicsCuts[gKinematicVariables][2]; // [phi,pt,eta,energy,charge][min,max]
  Bool_t fUseKinematicsCuts[gKinematicVariables];   // if not set via setter, corresponding cut is kFALSE. Therefore, correspondig cut is open (default values are NOT used)
  //    DCA:
  TH1D *fDCAHist[2][2][2]; // distance of closest approach (DCA) [before,after][reco,sim][xy,z] // "xy" is transverse direction
  Double_t fDCABins[2][3]; // [xy,z][nBins,min,max]
  Double_t fDCACuts[2][2]; // [xy,z][min,max]
  //Bool_t fUseDCACuts[3]; // if not set via setter, corresponding cut is kFALSE. Therefore, correspondig cut is open (default values are NOT used)
  //    All remaining particle histograms: 
  TH1D *fParticleHist[2][2][gParticleHistograms]; //! distributions [before,after particle cuts][reco,sim][type - see enum]
  Double_t fParticleBins[gParticleHistograms][3]; // [nBins,min,max]
  Double_t fParticleCuts[gParticleHistograms][2]; // [type - see enum][min,max]

  // 4) Q-vectors:
  TList *fQvectorList;        // list to hold all Q-vector objects       
  TProfile *fQvectorFlagsPro; // profile to hold all flags for Q-vector
  Bool_t fCalculateQvector;   // to calculate or not to calculate Q-vector components, that's a Boolean...
  Int_t fMaxHarmonic;         // 6 (not going beyond v6, if you change this value, change also fQvector[49][9]) 
  Int_t fMaxCorrelator;       // 8 (not going beyond 8-p correlations, if you change this value, change also fQvector[49][9]) 
  TComplex fQvector[49][9];   //! Q-vector components [fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*8+1][8+1]  

  // 5) Particle weights: 
  TList *fWeightsList;                  // list to hold all Q-vector objects       
  TProfile *fWeightsFlagsPro;           // profile to hold all flags for weights
  Bool_t fUseWeights[gWeights]; // use weights [phi,pt,eta]
  TH1D *fWeightsHist[gWeights]; // histograms holding weights [phi,pt,eta]

  // 6) Correlations:
  TList *fCorrelationsList;            // list to hold all correlations objects
  TProfile *fCorrelationsFlagsPro;     // profile to hold all flags for correlations
  Bool_t fCalculateCorrelations;       // calculate and store correlations
  TProfile *fCorrelationsPro[4][6][3]; //! multiparticle correlations [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=6][0=integrated,1=vs. multiplicity,2=vs. centrality]

  // *) Nested loops:
  TList *fNestedLoopsList;               // list to hold all nested loops objects
  TProfile *fNestedLoopsFlagsPro;        // profile to hold all flags for nested loops
  Bool_t fCalculateNestedLoops;          // calculate and store correlations with nested loops, as a cross-check
  TProfile *fNestedLoopsPro[4][6][3];    //! multiparticle correlations from nested loops [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=6][0=integrated,1=vs. multiplicity,2=vs. centrality]
  //TProfile *fNestedLoopsPerDemandPro[3]; // which correlator needs to be cross-checked with nested loops (no setter => recompile). [0=integrated,1=vs. multiplicity,2=vs. centrality]
  TArrayD *ftaNestedLoops[2];            //! e-b-e container for nested loops [0=angles;1=product of all weights]   

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
 
  // Increase this counter in each new version:
  ClassDef(AliAnalysisTaskMuPa,3);

};

//================================================================================================================

#endif




