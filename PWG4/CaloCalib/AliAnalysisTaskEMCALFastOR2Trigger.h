/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIANALYSISTASKEMCALFASTOR2TRIGGER_H
#define ALIANALYSISTASKEMCALFASTOR2TRIGGER_H

class TList;
class TTree;
class AliEMCALFastORPatch;
class AliCaloCalibPedestal;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALFastOR2Trigger : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEMCALFastOR2Trigger();
  AliAnalysisTaskEMCALFastOR2Trigger(const char *name);
  virtual ~AliAnalysisTaskEMCALFastOR2Trigger();
  
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  
  void SetCheckDeadClusters(Bool_t c) {fCheckDeadClusters = c;}
  Bool_t GetCheckDeadClusters() {return fCheckDeadClusters;}
  
  void SetPedestal(AliCaloCalibPedestal *pds) {fPedestal = pds;}
  AliCaloCalibPedestal* GetPedestal() {return fPedestal;}
  
  void SetnCol(Int_t n) {nCol = n;}
  Int_t GetnCol() {return nCol;}
  void SetnRow(Int_t n) {nRow = n;}
  Int_t GetnRow() {return nRow;}
  void SetshiftCol(Int_t n) {shiftCol = n;}
  Int_t GetshiftCol() {return shiftCol;}
  void SetshiftRow(Int_t n) {shiftRow = n;}
  Int_t GetshiftRow() {return shiftRow;}
  
  void SetTriggerClustersName(TString name) {fTriggerClustersName = name;}
  
  Int_t GetMinL0Time() {return fMinL0Time;}
  void SetMinL0Time(Int_t t) {fMinL0Time = t;}
  Int_t GetMaxL0Time() {return fMaxL0Time;}
  void SetMaxL0Time(Int_t t) {fMaxL0Time = t;}
  
  Bool_t GetTimeCutOn() {return fTimeCutOn;}
  
  void SetTimeCutOn(Bool_t yes) {fTimeCutOn = yes;}
  
private:
  Bool_t                fCheckDeadClusters;
  Int_t                 nCol;
  Int_t                 nRow;
  Int_t                 shiftCol;
  Int_t                 shiftRow;
  AliCaloCalibPedestal *fPedestal;
  TString               fTriggerClustersName;
  Int_t                 fMinL0Time;
  Int_t                 fMaxL0Time;
  Bool_t                fTimeCutOn;
  
  AliAnalysisTaskEMCALFastOR2Trigger (const AliAnalysisTaskEMCALFastOR2Trigger&); // not implemented
  AliAnalysisTaskEMCALFastOR2Trigger operator=(const AliAnalysisTaskEMCALFastOR2Trigger&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALFastOR2Trigger, 1); 
};

#endif

