/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/**************************************
 * analysis task for Q-cumulants      * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/

#ifndef ALIANALYSISTASKQCUMULANTS_H
#define ALIANALYSISTASKQCUMULANTS_H

#include "AliAnalysisTask.h"

class AliESDEvent;
class AliAODEvent;
class AliCFManager;
class AliFlowAnalysisWithQCumulants;
class AliFlowEventSimpleMaker;
class TFile;

//================================================================================================================

class AliAnalysisTaskQCumulants : public AliAnalysisTask{
 public:
  AliAnalysisTaskQCumulants();
  AliAnalysisTaskQCumulants(const char *name, Bool_t useWeights=kFALSE);
  virtual ~AliAnalysisTaskQCumulants(){}; 
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
 
 private:
  AliAnalysisTaskQCumulants(const AliAnalysisTaskQCumulants& aatqc);
  AliAnalysisTaskQCumulants& operator=(const AliAnalysisTaskQCumulants& aatqc);
  
  AliFlowEventSimple *fEvent;          // the input event
  AliFlowAnalysisWithQCumulants *fQCA; // Q-cumulant Analysis (QCA) object
  TList  *fListHistos;                 // collection of output 

  Bool_t fUseWeights;                  // use any weights
  Bool_t fUsePhiWeights;               // use phi weights
  Bool_t fUsePtWeights;                // use pt weights
  Bool_t fUseEtaWeights;               // use eta weights  
  TList *fListWeights;                 // list with weights
     
  ClassDef(AliAnalysisTaskQCumulants, 1); 
};

//================================================================================================================

#endif











