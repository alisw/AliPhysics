#ifndef ALIT0QATASK_H
#define ALIT0QATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the T0 data in simulated data
//
//*-- Alla Maevskaya
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"  

class AliESD ; 
class TNtuple ;
class TH1D ; 
class TH1I ; 

class AliT0QATask : public AliAnalysisTask {

public:
  AliT0QATask(const char *name) ;
  virtual ~AliT0QATask() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  // Histograms
  TH1F * fhT01;
  TH1F * fhT02;
  TH1F * fhT03;   
  
  ClassDef(AliT0QATask, 0); // a T0 photon analysis task 
};
#endif // ALIT0QATASK_H
