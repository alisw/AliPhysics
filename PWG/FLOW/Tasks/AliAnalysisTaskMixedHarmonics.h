/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/**************************************
 * analysis task for mixed harmomics  * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/

#ifndef ALIANALYSISTASKMIXEDHARMONICS_H
#define ALIANALYSISTASKMIXEDHARMONICS_H


#include "AliAnalysisTaskSE.h"
#include "TF1.h"

class TString;
class TList;
class AliFlowEventSimple;
class AliFlowAnalysisWithMixedHarmonics;
class AliAODEvent;
class AliVVertex;
class AliMultSelection;
class AliAnalysisUtils;

//================================================================================================================

class AliAnalysisTaskMixedHarmonics : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskMixedHarmonics();
  AliAnalysisTaskMixedHarmonics(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskMixedHarmonics(){
	  if(fV0CutPU)      delete fV0CutPU;
	  if(fSPDCutPU)     delete fSPDCutPU;
	  if(fMultCutPU)    delete fMultCutPU;
	  if(fCenCutLowPU)  delete fCenCutLowPU; 
	  if(fCenCutHighPU) delete fCenCutHighPU;   
  }; 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  // common:
  void SetHarmonic(Int_t const h) {this->fHarmonic = h;};
  Int_t GetHarmonic() const {return this->fHarmonic;}; 
  void SetNoOfMultipicityBins(Int_t const nomb) {this->fNoOfMultipicityBins = nomb;};
  Int_t GetNoOfMultipicityBins() const {return this->fNoOfMultipicityBins;};   
  void SetMultipicityBinWidth(Double_t const mbw) {this->fMultipicityBinWidth = mbw;};
  Double_t GetMultipicityBinWidth() const {return this->fMultipicityBinWidth;};   
  void SetMinMultiplicity(Double_t const mm) {this->fMinMultiplicity = mm;};
  Double_t GetMinMultiplicity() const {return this->fMinMultiplicity;}; 
  void SetOppositeChargesPOI(Bool_t const ocp) {this->fOppositeChargesPOI = ocp;};
  Bool_t GetOppositeChargesPOI() const {return this->fOppositeChargesPOI;}; 
  void SetEvaluateDifferential3pCorrelator(Bool_t const ed3pc) {this->fEvaluateDifferential3pCorrelator = ed3pc;};
  Bool_t GetEvaluateDifferential3pCorrelator() const {return this->fEvaluateDifferential3pCorrelator;};       
  void SetCorrectForDetectorEffects(Bool_t const cfde) {this->fCorrectForDetectorEffects = cfde;};
  Bool_t GetCorrectForDetectorEffects() const {return this->fCorrectForDetectorEffects;}; 
  void SetPrintOnTheScreen(Bool_t const pots) {this->fPrintOnTheScreen = pots;};
  Bool_t GetPrintOnTheScreen() const {return this->fPrintOnTheScreen;};  
  void SetCalculateVsM(Bool_t const cvm) {this->fCalculateVsM = cvm;};
  Bool_t GetCalculateVsM() const {return this->fCalculateVsM;};  
  void SetShowBinLabelsVsM(Bool_t const sblvm) {this->fShowBinLabelsVsM = sblvm;};
  Bool_t GetShowBinLabelsVsM() const {return this->fShowBinLabelsVsM;};     
  void SetCalculateVsZDC(Bool_t const cvzdc) {this->fCalculateVsZDC = cvzdc;};
  Bool_t GetCalculateVsZDC() const {return this->fCalculateVsZDC;};
  
  // particle weights:
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
  void  SetRejectPileUp(Bool_t  pileup) {this->fRejectPileUp = pileup;}
  void  SetRejectPileUpTight(Bool_t pileupT) {this->fRejectPileUpTight = pileupT;}
  void  SetFillQAHistograms(Bool_t const fillQA) {this->fFillQAHistograms = fillQA;};
  
  // set fIs2018Data
  void SetIs2018Data(Bool_t const is2018) {this->fIs2018Data = is2018;};
  Bool_t GetIs2018Data() {return this->fIs2018Data;};

  // set pileup for 2018
  void SetupPileUpRemovalFunctions18qPass3();
  void SetupPileUpRemovalFunctions18rPass3();
  
 private:
  AliAnalysisTaskMixedHarmonics(const AliAnalysisTaskMixedHarmonics& aatmh);
  AliAnalysisTaskMixedHarmonics& operator=(const AliAnalysisTaskMixedHarmonics& aatmh);
  
  Bool_t CheckEventIsPileUp(AliAODEvent* faod);
  Bool_t PileUpMultiVertex(const AliAODEvent* faod);
  Bool_t CheckEventIsPileUp2018(AliAODEvent* faod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);

  AliMultSelection*   fMultSelection;    //! MultSelection (RUN2 centrality estimator)
  AliAnalysisUtils*    fAnalysisUtil;    //! Event selection


  AliFlowEventSimple *fEvent; // the input event
  AliFlowAnalysisWithMixedHarmonics *fMH; // mixed harmonics object
  TList *fListHistos; // collection of output 
  // common:  
  Int_t fHarmonic; // integer n in cos[n(phi1+phi2-2phi3)]
  Int_t fNoOfMultipicityBins; // number of multiplicity bins
  Double_t fMultipicityBinWidth; // width of multiplicity bin
  Double_t fMinMultiplicity; // minimal multiplicity
  Bool_t fOppositeChargesPOI; // two POIs, psi1 and psi2, in correlator <<cos[psi1+psi2-2phi3)]>> will be taken with opposite charges 
  Bool_t fEvaluateDifferential3pCorrelator; // evaluate <<cos[psi1+psi2-2phi3)]>>, where psi1 and psi2 are two POIs      
  Bool_t fCorrectForDetectorEffects; // correct 3-p correlator for detector effects  
  Bool_t fPrintOnTheScreen; // print or not the final results on the screen
  Bool_t fCalculateVsM; // calculate correlators vs multiplicity
  Bool_t fShowBinLabelsVsM; // in histograms holding results vs multiplicity show bin labels in the format M_lowEdge \leq M < M_upperEdge  
  Bool_t fCalculateVsZDC; // calculate correlators vs ZDC-P energy

  // particle weights:
  Bool_t fUseParticleWeights; // use any particle weights
  Bool_t fUsePhiWeights; // use phi weights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUseEtaWeights; // use eta weights  

  Bool_t fRejectPileUp; //
  Bool_t fRejectPileUpTight; //
  Bool_t fFillQAHistograms; //

  TH2F   *fTPCvsGlobalTrkBefore; //!  Global vs TPC tracks for QA
  TH2F   *fTPCvsGlobalTrkAfter; //!  Global vs TPC tracks for QA
  TH2F   *fTPCvsESDTrk; //!  Global vs TPC tracks for QA

  TList    *fWeightsList; // list with weights

  // bool variable for if the data set is 2018
  Bool_t fIs2018Data;
  // Functions for Pile Up Event Removal 2018 period:
  TF1                   *fV0CutPU;      //
  TF1                   *fSPDCutPU;     //
  TF1                   *fMultCutPU;    //
  TF1                   *fCenCutLowPU;  //
  TF1                   *fCenCutHighPU; //
  
  // QA histograms
  TH2F          *fHistTPConlyVsCL1Before; //!
  TH2F          *fHistTPConlyVsCL1After;  //!
  TH2F          *fHistTPConlyVsV0MBefore; //!
  TH2F          *fHistTPConlyVsV0MAfter;  //!
  TH2F          *fHistCentCL0VsV0MBefore; //!
  TH2F          *fHistCentCL0VsV0MAfter;  //!  
  TH2F          *fHistTPCVsESDTrkBefore;  //!
  TH2F          *fHistTPCVsESDTrkAfter;   //!




  ClassDef(AliAnalysisTaskMixedHarmonics, 1); 
};

//================================================================================================================

#endif










