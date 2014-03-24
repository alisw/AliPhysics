/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

 /************************************ 
 * flow analysis with multi-particle *
 *           correlations            * 
 *                                   * 
 * author: Ante Bilandzic            * 
 *        (abilandzic@gmail.com)     *
 ************************************/ 

#ifndef ALIFLOWANALYSISWITHMULTIPARTICLECORRELATIONS_H
#define ALIFLOWANALYSISWITHMULTIPARTICLECORRELATIONS_H

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TComplex.h"
#include "TDirectoryFile.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TStopwatch.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"

class AliFlowAnalysisWithMultiparticleCorrelations{
 public:
  AliFlowAnalysisWithMultiparticleCorrelations();
  virtual ~AliFlowAnalysisWithMultiparticleCorrelations(); 
  // Member functions are grouped as:
  // 0.) Methods called in the constructor;
  // 1.) Method Init() and methods called in it (!);
  // 2.) Method Make() and methods called in it;
  // 3.) Method Finish() and methods called in it;
  // 4.) Method GetOutputHistograms() and methods called in it;
  // 5.) Setters and getters;
  // 6.) The rest.

  // 0.) Methods called in the constructor:
  virtual void InitializeArraysForControlHistograms(); 
  virtual void InitializeArraysForQvector(); 
  virtual void InitializeArraysForCorrelations(); 
  virtual void InitializeArraysForEbECumulants();
  virtual void InitializeArraysForWeights();
  virtual void InitializeArraysForQcumulants();
  virtual void InitializeArraysForDiffCorrelations(); 
  virtual void InitializeArraysForNestedLoops(); 

  // 1.) Method Init() and methods called in it (!):
  virtual void Init();
   virtual void CrossCheckSettings();
   virtual void BookAndNestAllLists(); 
   virtual void BookEverythingForBase();
   virtual void BookEverythingForControlHistograms();
   virtual void BookEverythingForQvector();
   virtual void BookEverythingForWeights();
   virtual void BookEverythingForCorrelations();
   virtual void BookEverythingForEbECumulants();
   virtual void BookEverythingForNestedLoops();
   virtual void BookEverythingForStandardCandles();
   virtual void BookEverythingForQcumulants();
   virtual void BookEverythingForDiffCorrelations();
   
  // 2.) Method Make() and methods called in it:
  virtual void Make(AliFlowEventSimple *anEvent);
   virtual Bool_t CrossCheckInternalFlags(AliFlowEventSimple *anEvent);
   virtual void CrossCheckPointersUsedInMake(); 
   virtual void FillControlHistograms(AliFlowEventSimple *anEvent);
   virtual void FillQvector(AliFlowEventSimple *anEvent);
   virtual void CalculateCorrelations(AliFlowEventSimple *anEvent);
   virtual void CalculateDiffCorrelations(AliFlowEventSimple *anEvent);
   virtual void CalculateEbECumulants(AliFlowEventSimple *anEvent);
   virtual void ResetQvector();
   virtual void CrossCheckWithNestedLoops(AliFlowEventSimple *anEvent);
   virtual void CrossCheckDiffWithNestedLoops(AliFlowEventSimple *anEvent);

  // 3.) Method Finish() and methods called in it:
  virtual void Finish();
   virtual void CrossCheckPointersUsedInFinish(); 
   virtual void CalculateStandardCandles();
   virtual void CalculateQcumulants();
   virtual void CalculateReferenceFlow();

  // 4.) Method GetOutputHistograms() and methods called in it: 
  virtual void GetOutputHistograms(TList *histList);
   virtual void GetPointersForControlHistograms(); 
   virtual void GetPointersForQvector(); 
   virtual void GetPointersForCorrelations(); 
   virtual void GetPointersForStandardCandles(); 
   virtual void GetPointersForQcumulants(); 
   virtual void GetPointersForDiffCorrelations(); 

  // 5.) Setters and getters:
  //  5.0.) Base list and internal flags:
  void SetHistList(TList* const hlist) {this->fHistList = hlist;} 
  TList* GetHistList() const {return this->fHistList;} 
  void SetInternalFlagsPro(TProfile* const ifp) {this->fInternalFlagsPro = ifp;};
  TProfile* GetInternalFlagsPro() const {return this->fInternalFlagsPro;}; 
  void SetMinNoRPs(Int_t const min) {fUseInternalFlags = kTRUE; this->fMinNoRPs = min;};
  Int_t GetMinNoRPs() const {return this->fMinNoRPs;};
  void SetMaxNoRPs(Int_t const max) {fUseInternalFlags = kTRUE; this->fMaxNoRPs = max;};
  Int_t GetMaxNoRPs() const {return this->fMaxNoRPs;};
  void SetExactNoRPs(Int_t const exact) {fUseInternalFlags = kTRUE; this->fExactNoRPs = exact;};
  Int_t GetExactNoRPs() const {return this->fExactNoRPs;};
  void SetAnalysisTag(const char *at) {this->fAnalysisTag = TString(at);};
  TString GetAnalysisTag() const {return this->fAnalysisTag;};
  void SetDumpThePoints(Bool_t const dtp, Int_t const max) {this->fDumpThePoints = dtp; this->fMaxNoEventsPerFile = max;};

  //  5.1.) Control histograms:  
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetControlHistogramsFlagsPro(TProfile* const chfp) {this->fControlHistogramsFlagsPro = chfp;};
  TProfile* GetControlHistogramsFlagsPro() const {return this->fControlHistogramsFlagsPro;}; 
  void SetFillControlHistograms(Bool_t const fch) {this->fFillControlHistograms = fch;};
  Bool_t GetFillControlHistograms() const {return this->fFillControlHistograms;};
  void SetFillKinematicsHist(Bool_t const fkh) {this->fFillKinematicsHist = fkh;};
  Bool_t GetFillKinematicsHist() const {return this->fFillKinematicsHist;};
  void SetFillMultDistributionsHist(Bool_t const mdh) {this->fFillMultDistributionsHist = mdh;};
  Bool_t GetFillMultDistributionsHist() const {return this->fFillMultDistributionsHist;};
  void SetFillMultCorrelationsHist(Bool_t const mch) {this->fFillMultCorrelationsHist = mch;};
  Bool_t GetFillMultCorrelationsHist() const {return this->fFillMultCorrelationsHist;};
  void SetnBins(const char *type, const char *variable, const Int_t nBins); // .cxx
  void SetMin(const char *type, const char *variable, const Double_t min); // .cxx
  void SetMax(const char *type, const char *variable, const Double_t max); // .cxx
  void SetnBinsMult(const char *type, const Int_t nBinsMult); // .cxx
  void SetMinMult(const char *type, const Double_t minMult); // .cxx
  void SetMaxMult(const char *type, const Double_t maxMult); // .cxx

  //  5.2.) Q-vectors:
  void SetQvectorList(TList* const qvl) {this->fQvectorList = qvl;};
  TList* GetQvectorList() const {return this->fQvectorList;} 
  void SetQvectorFlagsPro(TProfile* const qvfp) {this->fQvectorFlagsPro = qvfp;};
  TProfile* GetQvectorFlagsPro() const {return this->fQvectorFlagsPro;}; 
  void SetCalculateQvector(Bool_t const cqv) {this->fCalculateQvector = cqv;};
  Bool_t GetCalculateQvector() const {return this->fCalculateQvector;};
  void SetCalculateDiffQvectors(Bool_t const cdqv) {this->fCalculateDiffQvectors = cdqv;};
  Bool_t GetCalculateDiffQvectors() const {return this->fCalculateDiffQvectors;};

  //  5.3.) Correlations:
  void SetCorrelationsList(TList* const cl) {this->fCorrelationsList = cl;};
  TList* GetCorrelationsList() const {return this->fCorrelationsList;} 
  void SetCorrelationsFlagsPro(TProfile* const cfp) {this->fCorrelationsFlagsPro = cfp;};
  TProfile* GetCorrelationsFlagsPro() const {return this->fCorrelationsFlagsPro;}; 
  void SetCalculateCorrelations(Bool_t const cc) {this->fCalculateCorrelations = cc;};
  Bool_t GetCalculateCorrelations() const {return this->fCalculateCorrelations;};
  void SetCalculateIsotropic(Bool_t const ci) {this->fCalculateIsotropic = ci;};
  Bool_t GetCalculateIsotropic() const {return this->fCalculateIsotropic;};
  void SetCalculateSame(Bool_t const cs) {this->fCalculateSame = cs;};
  Bool_t GetCalculateSame() const {return this->fCalculateSame;};
  void SetSkipZeroHarmonics(Bool_t const szh) {this->fSkipZeroHarmonics = szh;};
  Bool_t GetSkipZeroHarmonics() const {return this->fSkipZeroHarmonics;};
  void SetCalculateSameIsotropic(Bool_t const csi) {this->fCalculateSameIsotropic = csi;};
  Bool_t GetCalculateSameIsotropic() const {return this->fCalculateSameIsotropic;};
  void SetCalculateAll(Bool_t const ca) {this->fCalculateAll = ca;};
  Bool_t GetCalculateAll() const {return this->fCalculateAll;};
  void SetDontGoBeyond(Int_t const dgb) {this->fDontGoBeyond = dgb;};
  Int_t GetDontGoBeyond() const {return this->fDontGoBeyond;};
  void SetCalculateOnlyForHarmonicQC(Bool_t const cofhqc) {this->fCalculateOnlyForHarmonicQC = cofhqc;};
  Bool_t GetCalculateOnlyForHarmonicQC() const {return this->fCalculateOnlyForHarmonicQC;};
  void SetCalculateOnlyForSC(Bool_t const cofsc) {this->fCalculateOnlyForSC = cofsc;};
  Bool_t GetCalculateOnlyForSC() const {return this->fCalculateOnlyForSC;};
  void SetCalculateOnlyCos(Bool_t const coc) {this->fCalculateOnlyCos = coc;};
  Bool_t GetCalculateOnlyCos() const {return this->fCalculateOnlyCos;};
  void SetCalculateOnlySin(Bool_t const cos) {this->fCalculateOnlySin = cos;};
  Bool_t GetCalculateOnlySin() const {return this->fCalculateOnlySin;};

  //  5.4.) Event-by-event cumulants:
  void SetEbECumulantsList(TList* const ebecl) {this->fEbECumulantsList = ebecl;};
  TList* GetEbECumulantsList() const {return this->fEbECumulantsList;} 
  void SetEbECumulantsFlagsPro(TProfile* const ebecfp) {this->fEbECumulantsFlagsPro = ebecfp;};
  TProfile* GetEbECumulantsFlagsPro() const {return this->fEbECumulantsFlagsPro;}; 
  void SetCalculateEbECumulants(Bool_t const cebec) {this->fCalculateEbECumulants = cebec;};
  Bool_t GetCalculateEbECumulants() const {return this->fCalculateEbECumulants;};

  //  5.5.) Weights: 
  void SetWeightsList(TList* const wl) {this->fWeightsList = (TList*)wl->Clone();};
  TList* GetWeightsList() const {return this->fWeightsList;} 
  void SetWeightsFlagsPro(TProfile* const wfp) {this->fWeightsFlagsPro = wfp;};
  TProfile* GetWeightsFlagsPro() const {return this->fWeightsFlagsPro;}; 
  void SetWeightsHist(TH1D* const hist, const char *type, const char *variable); // .cxx
  
  //  5.6.) Nested loops:
  void SetNestedLoopsList(TList* const nll) {this->fNestedLoopsList = nll;} 
  TList* GetNestedLoopsList() const {return this->fNestedLoopsList;} 
  void SetNestedLoopsFlagsPro(TProfile* const nlfp) {this->fNestedLoopsFlagsPro = nlfp;};
  TProfile* GetNestedLoopsFlagsPro() const {return this->fNestedLoopsFlagsPro;}; 
  void SetCrossCheckWithNestedLoops(Bool_t const ccwnl) {this->fCrossCheckWithNestedLoops = ccwnl;};
  Bool_t GetCrossCheckWithNestedLoops() const {return this->fCrossCheckWithNestedLoops;};
  void SetCrossCheckDiffWithNestedLoops(Bool_t const ccdwnl) {this->fCrossCheckDiffWithNestedLoops = ccdwnl;};
  Bool_t GetCrossCheckDiffWithNestedLoops() const {return this->fCrossCheckDiffWithNestedLoops;};
  void SetCrossCheckDiffCSCOBN(Int_t const cs, Int_t const co, Int_t const bn)  
  {
   this->fCrossCheckDiffCSCOBN[0] = cs; // cos/sin
   this->fCrossCheckDiffCSCOBN[1] = co; // correlator order [1p,2p,3p,4p]
   this->fCrossCheckDiffCSCOBN[2] = bn; // bin number
  };

  // 5.7.) 'Standard candles':
  void SetStandardCandlesList(TList* const scl) {this->fStandardCandlesList = scl;} 
  TList* GetStandardCandlesList() const {return this->fStandardCandlesList;} 
  void SetStandardCandlesFlagsPro(TProfile* const scfp) {this->fStandardCandlesFlagsPro = scfp;};
  TProfile* GetStandardCandlesFlagsPro() const {return this->fStandardCandlesFlagsPro;}; 
  void SetCalculateStandardCandles(Bool_t const csc) {this->fCalculateStandardCandles = csc;};
  Bool_t GetCalculateStandardCandles() const {return this->fCalculateStandardCandles;};
  void SetPropagateErrorSC(Bool_t const pesc) {this->fPropagateErrorSC = pesc;};
  Bool_t GetPropagateErrorSC() const {return this->fPropagateErrorSC;};
  void SetStandardCandlesHist(TH1D* const sch) {this->fStandardCandlesHist = sch;};
  TH1D* GetStandardCandlesHist() const {return this->fStandardCandlesHist;}; 
  void SetProductsSCPro(TProfile2D* const psc) {this->fProductsSCPro = psc;};
  TProfile2D* GetProductsSCPro() const {return this->fProductsSCPro;}; 

  //  5.8.) Q-cumulants:
  void SetQcumulantsList(TList* const qcl) {this->fQcumulantsList = qcl;};
  TList* GetQcumulantsList() const {return this->fQcumulantsList;} 
  void SetQcumulantsFlagsPro(TProfile* const qcfp) {this->fQcumulantsFlagsPro = qcfp;};
  TProfile* GetQcumulantsFlagsPro() const {return this->fQcumulantsFlagsPro;}; 
  void SetCalculateQcumulants(Bool_t const cqc) {this->fCalculateQcumulants = cqc;};
  Bool_t GetCalculateQcumulants() const {return this->fCalculateQcumulants;};
  void SetHarmonicQC(Int_t const hqc) {this->fHarmonicQC = hqc;};
  Int_t GetHarmonicQC() const {return this->fHarmonicQC;};
  void SetPropagateErrorQC(Bool_t const peqc) {this->fPropagateErrorQC = peqc;};
  Bool_t GetPropagateErrorQC() const {return this->fPropagateErrorQC;};
  void SetQcumulantsHist(TH1D* const qch) {this->fQcumulantsHist = qch;};
  TH1D* GetQcumulantsHist() const {return this->fQcumulantsHist;}; 
  void SetReferenceFlowHist(TH1D* const rfh) {this->fReferenceFlowHist = rfh;};
  TH1D* GetReferenceFlowHist() const {return this->fReferenceFlowHist;}; 
  void SetProductsQCPro(TProfile2D* const pqc) {this->fProductsQCPro = pqc;};
  TProfile2D* GetProductsQCPro() const {return this->fProductsQCPro;}; 

  //  5.9.) Differential correlations:
  void SetDiffCorrelationsList(TList* const dcl) {this->fDiffCorrelationsList = dcl;};
  TList* GetDiffCorrelationsList() const {return this->fDiffCorrelationsList;} 
  void SetDiffCorrelationsFlagsPro(TProfile* const cdfp) {this->fDiffCorrelationsFlagsPro = cdfp;};
  TProfile* GetDiffCorrelationsFlagsPro() const {return this->fDiffCorrelationsFlagsPro;}; 
  void SetCalculateDiffCorrelations(Bool_t const cdc) {this->fCalculateDiffCorrelations = cdc;};
  Bool_t GetCalculateDiffCorrelations() const {return this->fCalculateDiffCorrelations;};
  void SetDiffHarmonics(Int_t order, Int_t *harmonics); // see implementation in .cxx file 
  void SetCalculateDiffCos(Bool_t const cdc) {this->fCalculateDiffCos = cdc;};
  Bool_t GetCalculateDiffCos() const {return this->fCalculateDiffCos;};
  void SetCalculateDiffSin(Bool_t const cds) {this->fCalculateDiffSin = cds;};
  Bool_t GetCalculateDiffSin() const {return this->fCalculateDiffSin;};
  void SetCalculateDiffCorrelationsVsPt(Bool_t const cdcvspt) {this->fCalculateDiffCorrelationsVsPt = cdcvspt;};
  Bool_t GetCalculateDiffCorrelationsVsPt() const {return this->fCalculateDiffCorrelationsVsPt;};
  void SetUseDefaultBinning(Bool_t const udb) {this->fUseDefaultBinning = udb;};
  Bool_t GetUseDefaultBinning() const {return this->fUseDefaultBinning;};
  void SetnDiffBins(Int_t const ndb) {this->fnDiffBins = ndb;};
  Int_t GetnDiffBins() const {return this->fnDiffBins;};
  void SetRangesDiffBins(Double_t* const rdb) {this->fRangesDiffBins = rdb;};
  Double_t* GetRangesDiffBins() const {return this->fRangesDiffBins;};

  // 6.) The rest:
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);
  virtual TComplex Q(Int_t n, Int_t p);
  virtual TComplex p(Int_t n, Int_t p);
  virtual TComplex q(Int_t n, Int_t p);
  virtual TComplex One(Int_t n1);
  virtual TComplex Two(Int_t n1, Int_t n2);
  virtual TComplex Three(Int_t n1, Int_t n2, Int_t n3);
  virtual TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  virtual TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  virtual TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6);
  virtual TComplex Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7);
  virtual TComplex Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8);
  virtual TComplex OneDiff(Int_t n1);
  virtual TComplex TwoDiff(Int_t n1, Int_t n2);
  virtual TComplex ThreeDiff(Int_t n1, Int_t n2, Int_t n3);
  virtual TComplex FourDiff(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  virtual Double_t Weight(const Double_t &value, const char *type, const char *variable); // value, [RP,POI], [phi,pt,eta]
  virtual Double_t CastStringToCorrelation(const char *string, Bool_t numerator);
  virtual Double_t Covariance(const char *x, const char *y, TProfile2D *profile2D, Bool_t bUnbiasedEstimator = kFALSE);
  virtual TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0); // Credits: Kristjan Gulbrandsen (gulbrand@nbi.dk) 
  virtual void CalculateProductsOfCorrelations(AliFlowEventSimple *anEvent, TProfile2D *profile2D);
  static void DumpPointsForDurham(TGraphErrors *ge);
  static void DumpPointsForDurham(TH1D *h);
  static void DumpPointsForDurham(TH1F *h);
  virtual void DumpThePoints(AliFlowEventSimple *anEvent);
  TH1D* GetHistogramWithWeights(const char *filePath, const char *listName, const char *type, const char *variable);

 private:
  AliFlowAnalysisWithMultiparticleCorrelations(const AliFlowAnalysisWithMultiparticleCorrelations& afawQc);
  AliFlowAnalysisWithMultiparticleCorrelations& operator=(const AliFlowAnalysisWithMultiparticleCorrelations& afawQc); 
  // Data members are grouped as:
  // 0.) Base list and internal flags;
  // 1.) Control histograms;  
  // 2.) Q-vectors;
  // 3.) Correlations;
  // 4.) Event-by-event cumulants;
  // 5.) Weights;
  // 6.) Nested loops;
  // 7.) 'Standard candles';
  // 8.) Q-cumulants;
  // 9.) Differential correlations. 

  // 0.) Base list and internal flags:
  TList* fHistList;            // base list to hold all output object (a.k.a. grandmother of all lists)
  TProfile *fInternalFlagsPro; // profile to hold all internal flags and settings
  Bool_t fUseInternalFlags;    // use internal flags (automatically set if some internal flag is used)
  Int_t fMinNoRPs;             // minimum number of RPs required for the analysis 
  Int_t fMaxNoRPs;             // maximum number of RPs allowed for the analysis 
  Int_t fExactNoRPs;           // exact (randomly shuffled) number of RPs selected for the analysis 
  Bool_t fPropagateError;      // prevent error propagation if something strange happens during calculations 
  TString fAnalysisTag;        // tag internally this analysis
  Bool_t fDumpThePoints;       // dump the data points into the external file 
  Int_t fMaxNoEventsPerFile;   // maximum number of events to be dumped in a single file

  // 1.) Control histograms:  
  TList *fControlHistogramsList;        // list to hold all 'control histograms' objects
  TProfile *fControlHistogramsFlagsPro; // profile to hold all flags for control histograms
  Bool_t fFillControlHistograms;        // fill or not control histograms (by default they are all filled)
  Bool_t fFillKinematicsHist;           // fill or not fKinematicsHist[2][3]
  Bool_t fFillMultDistributionsHist;    // fill or not TH1D *fMultDistributionsHist[3]    
  Bool_t fFillMultCorrelationsHist;     // fill or not TH2D *fMultCorrelationsHist[3]  
  TH1D *fKinematicsHist[2][3];          // [RP,POI][phi,pt,eta] distributions
  TH1D *fMultDistributionsHist[3];      // multiplicity distribution [RP,POI,reference multiplicity]
  TH2I *fMultCorrelationsHist[3];       // [RP vs. POI, RP vs. refMult, POI vs. refMult]  
  Int_t fnBins[2][3];                   // [RP,POI][phi,pt,eta], corresponds to fKinematicsHist[2][3]
  Double_t fMin[2][3];                  // [RP,POI][phi,pt,eta], corresponds to fKinematicsHist[2][3]
  Double_t fMax[2][3];                  // [RP,POI][phi,pt,eta], corresponds to fKinematicsHist[2][3]
  Int_t fnBinsMult[3];                  // [RP,POI,REF], corresponds to fMultDistributionsHist[3]   
  Double_t fMinMult[3];                 // [RP,POI,REF], corresponds to fMultDistributionsHist[3]   
  Double_t fMaxMult[3];                 // [RP,POI,REF], corresponds to fMultDistributionsHist[3]   
  
  // 2.) Q-vectors:
  TList *fQvectorList;           // list to hold all Q-vector objects       
  TProfile *fQvectorFlagsPro;    // profile to hold all flags for Q-vector
  Bool_t fCalculateQvector;      // to calculate or not to calculate Q-vector components, that's a Boolean...
  TComplex fQvector[49][9];      // Q-vector components [fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*8+1][8+1]  
  Bool_t fCalculateDiffQvectors; // to calculate or not to calculate p- and q-vector components, that's a Boolean...  
  TComplex fpvector[100][49][9]; // p-vector components [bin][fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*8+1][8+1] TBI hardwired 100
  TComplex fqvector[100][49][9]; // q-vector components [bin][fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*8+1][8+1] TBI hardwired 100

  // 3.) Correlations:
  TList *fCorrelationsList;           // list to hold all correlations objects
  TProfile *fCorrelationsFlagsPro;    // profile to hold all flags for correlations
  TProfile *fCorrelationsPro[2][8];   // multi-particle correlations [0=cos,1=sin][1p,2p,...,8p]
  Bool_t fCalculateCorrelations;      // calculate and store correlations
  Int_t fMaxHarmonic;                 // 6 (not going beyond v6, if you change this value, change also fQvector[49][9]) 
  Int_t fMaxCorrelator;               // 8 (not going beyond 8-p correlations, if you change this value, change also fQvector[49][9]) 
  Bool_t fCalculateIsotropic;         // calculate only isotropic correlations
  Bool_t fCalculateSame;              // calculate only 'same abs harmonics' correlations TBI 
  Bool_t fSkipZeroHarmonics;          // skip correlations which have some of the harmonicc equal to zero
  Bool_t fCalculateSameIsotropic;     // calculate all isotropic correlations in 'same abs harmonic' TBI this can be implemented better
  Bool_t fCalculateAll;               // calculate all possible correlations 
  Int_t  fDontGoBeyond;               // do not go beyond fDontGoBeyond-p correlators
  Bool_t fCalculateOnlyForHarmonicQC; // calculate only isotropic correlations in |fHarmonicQC|
  Bool_t fCalculateOnlyForSC;         // calculate only correlations needed for 'standard candles'
  Bool_t fCalculateOnlyCos;           // calculate only 'cos' correlations
  Bool_t fCalculateOnlySin;           // calculate only 'sin' correlations

  // 4.) Event-by-event cumulants:
  TList *fEbECumulantsList;         // list to hold all e-b-e cumulants objects
  TProfile *fEbECumulantsFlagsPro;  // profile to hold all flags for e-b-e cumulants
  TProfile *fEbECumulantsPro[2][8]; // multi-particle e-b-e cumulants [0=cos,1=sin][1p,2p,...,8p]
  Bool_t fCalculateEbECumulants;    // calculate and store e-b-e cumulants
 
  // 5.) Weights:
  TList *fWeightsList;        // list to hold all weights objects
  TProfile *fWeightsFlagsPro; // profile to hold all flags for weights
  Bool_t fUseWeights[2][3];   // use weights [RP,POI][phi,pt,eta]
  TH1D *fWeightsHist[2][3];   // histograms holding weights [RP,POI][phi,pt,eta]

  // 6.) Nested loops:
  TList *fNestedLoopsList;               // list to hold all nested loops objects
  TProfile *fNestedLoopsFlagsPro;        // profile to hold all flags for nested loops
  Bool_t fCrossCheckWithNestedLoops;     // cross-check results with nested loops
  Bool_t fCrossCheckDiffWithNestedLoops; // cross-check differential correlators with nested loops
  Int_t fCrossCheckDiffCSCOBN[3];        // [0=cos,1=sin][1p,2p,...,4p][binNo]
  TProfile *fNestedLoopsResultsCosPro;   // profile to hold nested loops results (cosine)
  TProfile *fNestedLoopsResultsSinPro;   // profile to hold nested loops results (sinus)
  TProfile *fNestedLoopsDiffResultsPro;  // profile to hold differential nested loops results // TBI

  // 7.) 'Standard candles':
  TList *fStandardCandlesList;        // list to hold all 'standard candles' objects
  TProfile *fStandardCandlesFlagsPro; // profile to hold all flags fo 'standard candles'
  Bool_t fCalculateStandardCandles;   // calculate and store 'standard candles'
  Bool_t fPropagateErrorSC;           // propagate and store error for 'standard candles'
  TH1D *fStandardCandlesHist;         // histogram to hold results for 'standard candles'
  TProfile2D *fProductsSCPro;         // 2D profile to hold products of correlations needed for SC error propagation

  // 8.) Q-cumulants:
  TList *fQcumulantsList;        // list to hold all Q-cumulants objects
  TProfile *fQcumulantsFlagsPro; // profile to hold all flags for Q-cumulants
  Bool_t fCalculateQcumulants;   // calculate and store Q-cumulants
  Int_t fHarmonicQC;             // calculate Q-cumulants in this harmonic (default is 2) 
  Bool_t fPropagateErrorQC;      // propagate and store error for Q-cumulants
  TH1D *fQcumulantsHist;         // two- and multi-particle Q-cumulants
  TH1D *fReferenceFlowHist;      // reference flow from two- and multi-particle Q-cumulants
  TProfile2D *fProductsQCPro;    // 2D profile to hold products of correlations needed for QC error propagation

  // 9.) Differential correlations:
  TList *fDiffCorrelationsList;          // list to hold all correlations objects
  TProfile *fDiffCorrelationsFlagsPro;   // profile to hold all flags for correlations
  Bool_t fCalculateDiffCorrelations;     // calculate and store differential correlations
  Bool_t fCalculateDiffCos;              // calculate and store differential cosine correlations (kTRUE by default)
  Bool_t fCalculateDiffSin;              // calculate and store differential sinus correlations (kFALSE by default)
  Bool_t fCalculateDiffCorrelationsVsPt; // calculate differential correlations vs pt (default), or vs eta
  Bool_t fUseDefaultBinning;             // use default binning in pt or in eta
  Int_t fnDiffBins;                      // number of differential bins in pt or in eta (when non-default binning is used)
  Double_t *fRangesDiffBins;             // ranges for differential bins in pt or in eta (when non-default binning is used)
  Int_t fDiffHarmonics[4][4];            // harmonics for differential correlations [order][{n1},{n1,n2},...,{n1,n2,n3,n4}] 
  TProfile *fDiffCorrelationsPro[2][4];  // multi-particle correlations [0=cos,1=sin][1p,2p,3p,4p]
  UInt_t fDiffBinNo;                     // differential bin number

  ClassDef(AliFlowAnalysisWithMultiparticleCorrelations,1);

};

//================================================================================================================

#endif





