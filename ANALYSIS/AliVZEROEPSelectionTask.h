#ifndef ALIVZEROEPSELECTIONTASK_H
#define ALIVZEROEPSELECTIONTASK_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliVZEROEPSelectionTask
//   author: Cvetan Cheshkov
//   30/01/2012
//   This analysis task reads the OADB and
//   provides the parameters needed to flatten
//   the VZERO event plane in AliEventplane
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TProfile;

class AliOADBContainer;
class AliEventplane;

class AliVZEROEPSelectionTask : public AliAnalysisTaskSE {

 public:
  AliVZEROEPSelectionTask();
  AliVZEROEPSelectionTask(const char *name);
  virtual ~AliVZEROEPSelectionTask();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  void SetUserParams(const char* inFileName, const char* listName);
  void UseVZEROCentrality() {fUseVZEROCentrality = kTRUE;}
  void SetEventplaneParams(AliEventplane *esdEP,Float_t percentile);  
 private:

  void SetHistograms(TList *list);
  void SetParamsFromOADB();
   
  AliVZEROEPSelectionTask(const AliVZEROEPSelectionTask& ep);
  AliVZEROEPSelectionTask& operator= (const AliVZEROEPSelectionTask& ep); 

  Int_t    fRunNumber;			// runnumber
  Bool_t   fUserParams;		        // in case one wants to use custom flatenning params
  Bool_t   fUseVZEROCentrality;         // use VZERO centrality estimator instead of SPD
  AliOADBContainer* fVZEROEPContainer;	// VZERO event-plane OADB Container

  TProfile *fX2In[11];                   // Profile histogram for Q^2_x (read from input file)
  TProfile *fY2In[11];                   // Profile histogram for Q^2_y (read from input file)
  TProfile *fX2Y2In[11];                 // Profile histogram for Q^2_x*Q^2_y (read from input file)
  TProfile *fCos8PsiIn[11];              // Profile histogram for Cos(8*Psi) (read from input file)

  ClassDef(AliVZEROEPSelectionTask,2) 
};

#endif

