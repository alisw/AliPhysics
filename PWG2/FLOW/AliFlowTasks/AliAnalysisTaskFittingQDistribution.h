/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/**************************************
 *    analysis task for fitting       * 
 *         q-distribution             *
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/

#ifndef ALIANALYSISTASKFITTINGQDISTRIBUTION_H
#define ALIANALYSISTASKFITTINGQDISTRIBUTION_H

#include "TString.h"
#include "AliAnalysisTask.h"

class AliESDEvent;
class AliAODEvent;
class AliCFManager;
class AliFittingQDistribution;
class AliFlowEventSimpleMaker;
class TFile;

//================================================================================================================

class AliAnalysisTaskFittingQDistribution : public AliAnalysisTask{
 public:
  AliAnalysisTaskFittingQDistribution();
  AliAnalysisTaskFittingQDistribution(const char *name, Bool_t QAon = kFALSE, Bool_t useWeights=kFALSE);
  virtual ~AliAnalysisTaskFittingQDistribution(){}; 
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void           SetAnalysisType(TString type) {this->fAnalysisType = type;}
  TString GetAnalysisType() const {return this->fAnalysisType;}
  
  void SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr;} 
  AliCFManager* GetCFManager1()           {return this->fCFManager1;}
  void SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr;} 
  AliCFManager* GetCFManager2()           {return this->fCFManager2;} 
  void          SetQAList1(TList* list)   {this->fQAInt = list; }
  TList*        GetQAList1()              {return this->fQAInt; }
  void          SetQAList2(TList* list)   {this->fQADiff = list; }
  TList*        GetQAList2()              {return this->fQADiff; }
  void          SetQAOn(Bool_t kt)        {this->fQA = kt; }
  Bool_t        GetQAOn()                 {return this->fQA; }
  
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
 
 private:
  AliAnalysisTaskFittingQDistribution(const AliAnalysisTaskFittingQDistribution& aatfqd);
  AliAnalysisTaskFittingQDistribution& operator=(const AliAnalysisTaskFittingQDistribution& aatfqd);

  AliESDEvent *fESD;                      // ESD object
  AliAODEvent* fAOD;                      // AOD object
  AliFittingQDistribution* fFQDA;         // Fitting q-distribution Analysis (FQDA) object
  AliFlowEventSimpleMaker* fEventMaker;   // FlowEventSimple maker object
  TString fAnalysisType;                  // string to select which kind of input to analyse (ESD, AOD or MC)
  AliCFManager* fCFManager1;              // correction framework manager
  AliCFManager* fCFManager2;              // correction framework manager
  TList  *fListHistos;                    // collection of output 
     
  TList*        fQAInt;                   // QA histogram list
  TList*        fQADiff;                  // QA histogram list

  Bool_t fQA;                             // flag to set the filling of the QA histograms  
  
  Bool_t       fUseWeights;               // use any weights
  Bool_t       fUsePhiWeights;            // phi weights
  TList*       fListWeights;              // list with weights                                               
                                                           
  ClassDef(AliAnalysisTaskFittingQDistribution, 1); 
};

//================================================================================================================

#endif











