#ifndef ALIANALYSISTASKAODCENTRALITYMAKER_H
#define ALIANALYSISTASKAODCENTRALITYMAKER_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskAODCentralityMaker
// AliAnalysisTaskSE to make AOD centrality
// Author: Alberica Toia, CERN, Alberica.Toia@cern.ch
//*************************************************************************

#include "AliAnalysisTaskSE.h"
class AliAODCentrality;
class AliAODHeader;


class AliAnalysisTaskAODCentralityMaker : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskAODCentralityMaker();
  AliAnalysisTaskAODCentralityMaker(const char *name);
  virtual ~AliAnalysisTaskAODCentralityMaker();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetDeltaAODFileName(const char* name) {fDeltaAODFileName=name;}
  const char* GetDeltaAODFileName() const {return fDeltaAODFileName.Data();}

  void SetMCInput() {fIsMCInput = kTRUE;}

 private:


  AliAnalysisTaskAODCentralityMaker(const AliAnalysisTaskAODCentralityMaker &source);
  AliAnalysisTaskAODCentralityMaker& operator=(const AliAnalysisTaskAODCentralityMaker& source); 
  AliAODCentrality *fAODCentrality;    // AOD centrality pointer   
  TString       fDeltaAODFileName;     // Name of output file
  AliAODHeader* fAODHeader;            // Header for replaction
  
  Bool_t   fIsMCInput;          // true when input is MC


  ClassDef(AliAnalysisTaskAODCentralityMaker,1); // AliAnalysisTaskSE to make AOD centrality
};

#endif

