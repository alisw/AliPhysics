/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/****************************************
 * analysis task for flow analysis with *
 *     multi-particle correlations      * 
 *                                      * 
 * author: Ante Bilandzic               *
 *         (abilandzic@gmail.com)       * 
 ***************************************/

#ifndef ALIANALYSISTASKMULTIPARTICLECORRELATIONS_H
#define ALIANALYSISTASKMULTIPARTICLECORRELATIONS_H

#include "AliAnalysisTaskSE.h"
#include "AliFlowAnalysisWithMultiparticleCorrelations.h"
#include "AliFlowEventSimple.h"

//================================================================================================================

class AliAnalysisTaskMultiparticleCorrelations : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskMultiparticleCorrelations();
  AliAnalysisTaskMultiparticleCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskMultiparticleCorrelations(){}; 
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  // Internal flags:
  void SetMinNoRPs(Int_t const min) {fUseInternalFlags = kTRUE; this->fMinNoRPs = min;};
  Int_t GetMinNoRPs() const {return this->fMinNoRPs;};
  void SetMaxNoRPs(Int_t const max) {fUseInternalFlags = kTRUE; this->fMaxNoRPs = max;};
  Int_t GetMaxNoRPs() const {return this->fMaxNoRPs;};
  void SetExactNoRPs(Int_t const exact) {fUseInternalFlags = kTRUE; this->fExactNoRPs = exact;};
  Int_t GetExactNoRPs() const {return this->fExactNoRPs;};
  void SetAnalysisTag(const char *at) {this->fAnalysisTag = TString(at);};
  TString GetAnalysisTag() const {return this->fAnalysisTag;};
  void SetDumpThePoints(Bool_t const dtp, Int_t const max) {this->fDumpThePoints = dtp; this->fMaxNoEventsPerFile = max;};

  // Control histograms:
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

  // Q-vectors:
  void SetCalculateQvector(Bool_t const cqv) {this->fCalculateQvector = cqv;};
  Bool_t GetCalculateQvector() const {return this->fCalculateQvector;};
  void SetCalculateDiffQvectors(Bool_t const cdqv) {this->fCalculateDiffQvectors = cdqv;};
  Bool_t GetCalculateDiffQvectors() const {return this->fCalculateDiffQvectors;};

  // Weights:              
  void SetWeightsHist(TH1D* const hist, const char *type, const char *variable); // .cxx
  TH1D* GetHistogramWithWeights(const char *filePath, const char *listName, const char *type, const char *variable)
  {
   AliFlowAnalysisWithMultiparticleCorrelations *mpc;
   return mpc->GetHistogramWithWeights(filePath,listName,type,variable);
  };

  // Correlations:
  void SetCalculateCorrelations(Bool_t const cc) {this->fCalculateCorrelations = cc;};
  Bool_t GetCalculateCorrelations() const {return this->fCalculateCorrelations;};
  void SetCalculateIsotropic(Bool_t const ci) {this->fCalculateIsotropic = ci;};
  Bool_t GetCalculateIsotropic() const {return this->fCalculateIsotropic;};
  void SetCalculateSame(Bool_t const csh) {this->fCalculateSame = csh;};
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

  // Event-by-event cumulants:
  void SetCalculateEbECumulants(Bool_t const cebec) {this->fCalculateEbECumulants = cebec;};
  Bool_t GetCalculateEbECumulants() const {return this->fCalculateEbECumulants;};

  // Nested loops:
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

  // 'Standard candles':
  void SetCalculateStandardCandles(Bool_t const csc) {this->fCalculateStandardCandles = csc;};
  Bool_t GetCalculateStandardCandles() const {return this->fCalculateStandardCandles;};
  void SetPropagateErrorSC(Bool_t const pesc) {this->fPropagateErrorSC = pesc;};
  Bool_t GetPropagateErrorSC() const {return this->fPropagateErrorSC;};

  // Q-cumulants:
  void SetCalculateQcumulants(Bool_t const cqc) {this->fCalculateQcumulants = cqc;};
  Bool_t GetCalculateQcumulants() const {return this->fCalculateQcumulants;};
  void SetHarmonicQC(Int_t const hqc) {this->fHarmonicQC = hqc;};
  Int_t GetHarmonicQC() const {return this->fHarmonicQC;};
  void SetPropagateErrorQC(Bool_t const peqc) {this->fPropagateErrorQC = peqc;};
  Bool_t GetPropagateErrorQC() const {return this->fPropagateErrorQC;};

  // Differential correlations:
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

 private:
  AliAnalysisTaskMultiparticleCorrelations(const AliAnalysisTaskMultiparticleCorrelations& aatqc);
  AliAnalysisTaskMultiparticleCorrelations& operator=(const AliAnalysisTaskMultiparticleCorrelations& aatqc);
  
  AliFlowEventSimple *fEvent; // the input event
  AliFlowAnalysisWithMultiparticleCorrelations *fMPC; // "multi-particle correlations" object
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // Internal flags:
  Bool_t fUseInternalFlags;  // use internal flags (automatically set if some internal flag is used)
  Int_t fMinNoRPs;           // minimum number of RPs required for the analysis 
  Int_t fMaxNoRPs;           // maximum number of RPs allowed for the analysis 
  Int_t fExactNoRPs;         // exact (randomly shuffled) number of RPs selected for the analysis 
  TString fAnalysisTag;      // tag internally this analysis
  Bool_t fDumpThePoints;     // dump the data points into the external file
  Int_t fMaxNoEventsPerFile; // maximum number of events to be dumped in a single file

  // Control histograms:
  Bool_t fFillControlHistograms;     // fill or not control histograms (by default they are filled)
  Bool_t fFillKinematicsHist;        // fill or not fKinematicsHist[2][3]
  Bool_t fFillMultDistributionsHist; // fill or not TH1D *fMultDistributionsHist[3]    
  Bool_t fFillMultCorrelationsHist;  // fill or not TH2D *fMultCorrelationsHist[3] 
  Int_t fnBins[2][3];                // [RP,POI][phi,pt,eta], corresponds to fKinematicsHist[2][3]
  Double_t fMin[2][3];               // [RP,POI][phi,pt,eta], corresponds to fKinematicsHist[2][3]
  Double_t fMax[2][3];               // [RP,POI][phi,pt,eta], corresponds to fKinematicsHist[2][3]
  Int_t fnBinsMult[3];               // [RP,POI,REF], corresponds to fMultDistributionsHist[3]   
  Double_t fMinMult[3];              // [RP,POI,REF], corresponds to fMultDistributionsHist[3]   
  Double_t fMaxMult[3];              // [RP,POI,REF], corresponds to fMultDistributionsHist[3]  

  // Q-vectors:
  Bool_t fCalculateQvector;      // to calculate or not to calculate Q-vector components, that's a Boolean...
  Bool_t fCalculateDiffQvectors; // to calculate or not to calculate p- and q-vector components, that's a Boolean...

  // Weights:
  Bool_t fUseWeights[2][3]; // use weights [RP,POI][phi,pt,eta]
  TH1D *fWeightsHist[2][3]; // histograms holding weights [RP,POI][phi,pt,eta]

  // Correlations:
  Bool_t fCalculateCorrelations;      // calculate and store correlations, or perhaps not, if the weather is bad...
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

  // Event-by-event cumulants:
  Bool_t fCalculateEbECumulants; // calculate and store event-by-event cumulants

  // Nested loops:
  Bool_t fCrossCheckWithNestedLoops;     // cross-check results with nested loops
  Bool_t fCrossCheckDiffWithNestedLoops; // cross-check differential correlators with nested loops
  Int_t fCrossCheckDiffCSCOBN[3];        // [0=cos,1=sin][1p,2p,...,4p][binNo]

  // 'Standard candles':
  Bool_t fCalculateStandardCandles; // calculate and store 'standard candles'
  Bool_t fPropagateErrorSC;         // propagate and store error for 'standard candles'

  // Q-cumulants:
  Bool_t fCalculateQcumulants; // calculate and store Q-cumulants
  Int_t fHarmonicQC;           // calculate Q-cumulants in this harmonic (default is 2) 
  Bool_t fPropagateErrorQC;    // propagate and store error for Q-cumulants

  // Differential correlations:
  Bool_t fCalculateDiffCorrelations;     // calculate and store differential correlations
  Bool_t fCalculateDiffCos;              // calculate and store differential cosine correlations (kTRUE by default)
  Bool_t fCalculateDiffSin;              // calculate and store differential sinus correlations (kFALSE by default)
  Bool_t fCalculateDiffCorrelationsVsPt; // calculate differential correlations vs pt (default), or vs eta
  Bool_t fUseDefaultBinning;             // use default binning in pt or in eta
  Int_t fnDiffBins;                      // number of differential bins in pt or in eta (when non-default binning is used)
  Double_t *fRangesDiffBins;             // ranges for differential bins in pt or in eta (when non-default binning is used)

  ClassDef(AliAnalysisTaskMultiparticleCorrelations,1); 

};

//================================================================================================================

#endif











