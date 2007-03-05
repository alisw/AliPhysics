#ifndef ALIANACALOTRIGGER_H
#define ALIANACALOTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the PHOS photon data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"  

class AliESD ; 
class TNtuple ;
class TH1D ; 
class TH1I ; 

class AliAnaCaloTrigger : public AliAnalysisTask {

public:
  AliAnaCaloTrigger(const char *name) ;
  virtual ~AliAnaCaloTrigger() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

  TString GetCalorimeter()     const   {return fCalorimeter ; }
  void    SetCalorimeter(TString calo) {fCalorimeter = calo ; }

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  TString fCalorimeter ; // "PHOS" or "EMCAL"

  // Histograms
  TNtuple * fNtTrigger22 ;
  TNtuple * fNtTriggerNN ;

  ClassDef(AliAnaCaloTrigger, 0); // a PHOS photon analysis task 
};
#endif // ALIANACALOTRIGGER_H
