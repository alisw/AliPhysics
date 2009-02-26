/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/**************************************
 * analysis task for cumulant method  * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/

#ifndef ALIANALYSISTASKCUMULANTS_H
#define ALIANALYSISTASKCUMULANTS_H

#include "AliAnalysisTask.h"

class AliESDEvent;
class AliAODEvent;
class AliCFManager;
class AliFlowAnalysisWithCumulants;
class AliFlowEventSimpleMaker;
class TFile;

//================================================================================================================

class AliAnalysisTaskCumulants : public AliAnalysisTask{
 public:
  AliAnalysisTaskCumulants();
  AliAnalysisTaskCumulants(const char *name, Bool_t useWeights=kFALSE);
  virtual ~AliAnalysisTaskCumulants(){}; 
  
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
  AliAnalysisTaskCumulants(const AliAnalysisTaskCumulants& aatc);
  AliAnalysisTaskCumulants& operator=(const AliAnalysisTaskCumulants& aatc);

  AliFlowEventSimple *fEvent;          // the input event
  AliFlowAnalysisWithCumulants *fGFCA; // Generating Function Cumulant (GFC) analysis object
  TList *fListHistos;                  // collection of output 

  Bool_t fUseWeights;                  // use any weights
  Bool_t fUsePhiWeights;               // use phi weights
  Bool_t fUsePtWeights;                // use pt weights
  Bool_t fUseEtaWeights;               // use eta weights  
  TList  *fListWeights;                // list with weights
           
  ClassDef(AliAnalysisTaskCumulants, 1); 
};

//================================================================================================================

#endif











