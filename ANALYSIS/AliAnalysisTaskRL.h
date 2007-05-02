#ifndef ALIANALYSISTASKRL_H
#define ALIANALYSISTASKRL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Panos Christakoglou, 31/05/2006

//============================================================================
//   AliAnalysysTaskRL - Class representing a basic analysis task. 
//   Access to the run loader
//============================================================================

#include "AliAnalysisTask.h"

class TTree;

class AliRunLoader;
class AliHeader;
class AliStack;

class AliAnalysisTaskRL : public AliAnalysisTask {
 public:  
  AliAnalysisTaskRL();
  AliAnalysisTaskRL(const char *name, const char *title);
  virtual ~AliAnalysisTaskRL();
  
  virtual Bool_t Notify();

 protected:
  Bool_t        GetEntry(Long64_t ientry);

  AliRunLoader *GetRunLoader();
  AliHeader    *GetHeader();
  AliStack     *GetStack();

 private:
  void DeleteRunLoader();

  TTree        *fTree; //! pointer to the ESD tree
  AliRunLoader *fRunLoader; //! pointer to the RunLoader if galice.root was opened
  Bool_t        fKinematicsLoaded; // determines if the stack is properly loaded (AliRunLoader::LoadKinematics() succeeded), this is needed because the GetStack returnes a invalid stack object when the function failed
  Bool_t fHeaderLoaded; // determines if the header is properly loaded
  
  AliAnalysisTaskRL(const AliAnalysisTaskRL&);
  AliAnalysisTaskRL& operator=(const AliAnalysisTaskRL&);
 
  ClassDef(AliAnalysisTaskRL,2);  // Class describing an analysis task
};
#endif
