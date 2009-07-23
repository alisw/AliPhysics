#ifndef ALIANALYSISTASKCALOFILTER_H
#define ALIANALYSISTASKCALOFILTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskCaloFilter.h  $ */

//////////////////////////////////////////////////////////
// Filter the ESDCaloClusters and ESDCaloCells of EMCAL,
// PHOS or both, creating the corresponing AODCaloClusters
// and AODCaloCells.
// Keep also the AODHeader information and the vertex.
// Needed for calorimeter calibration.
// Copy of AliAnalysisTaskESDfilter.
// Author: Gustavo Conesa Balbastre (INFN - Frascati)
//////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCaloFilter : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskCaloFilter();
  AliAnalysisTaskCaloFilter(const char* name);
  virtual ~AliAnalysisTaskCaloFilter() {;}
  
  // Implementation of interface methods
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() {Init();}
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  
  void SetCalorimeter(TString calo) {fCalorimeter = calo;}
  TString GetCalorimeter() const {return fCalorimeter;}
  
  void CreateAODFromESD();
  void CreateAODFromAOD();

 private:
  AliAnalysisTaskCaloFilter(const AliAnalysisTaskCaloFilter&);
  AliAnalysisTaskCaloFilter& operator=(const AliAnalysisTaskCaloFilter&);
  
  TString fCalorimeter; //Calorimeter to filter
  
  ClassDef(AliAnalysisTaskCaloFilter, 1); // Analysis task for standard ESD filtering
};

#endif
