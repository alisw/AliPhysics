#ifndef ALIVZEROQATASK_H
#define ALIVZEROQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//___________________________________________________________________________
//
//   An analysis task to check the VZERO data in simulated data
//
//___________________________________________________________________________

#include <TTree.h> 
#include "AliAnalysisTask.h" 

class AliESD; 
class TH1I; 

class AliVZEROQATask : public AliAnalysisTask {

public:
  AliVZEROQATask(const char *name);
  virtual  ~AliVZEROQATask();
   
  virtual void Exec(Option_t * opt = "");
  virtual void ConnectInputData(Option_t *); 
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

public:
  
  TTree   * fChain;             //! pointer to the analyzed TTree or TChain
  AliESD  * fESD;               //! declaration of leave types

  TObjArray * fOutputContainer; //! output data container

// Histograms

  TH1I    * fhVZERONbPMA;       //! histo of V0A PMs
  TH1I    * fhVZERONbPMC;       //! histo of V0C PMs
  TH1I    * fhVZEROMultA;       //! histo of multiplicity in V0A 
  TH1I    * fhVZEROMultC;       //! histo of multiplicity in V0C
  
   
  ClassDef(AliVZEROQATask, 0); // a VZERO analysis task 
};
#endif // ALIVZEROQATASK_H
