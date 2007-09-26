#ifndef ALIANACALOTRIGGER_H
#define ALIANACALOTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the PHOS photon data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////


#include "AliAnalysisTask.h"  

class AliESDEvent ; 
class TNtuple ;
class TH1D ; 
class TH1I ; 
class TChain;

class AliAnaCaloTrigger : public AliAnalysisTask {

public:
  AliAnaCaloTrigger() ;
  AliAnaCaloTrigger(const char *name) ;
  AliAnaCaloTrigger(const AliAnaCaloTrigger & trig) ;
  AliAnaCaloTrigger & operator=(const AliAnaCaloTrigger& source);
  virtual ~AliAnaCaloTrigger() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

  TString GetCalorimeter()     const   {return fCalorimeter ; }
  void    SetCalorimeter(TString calo) {fCalorimeter = calo ; }

private:
  TChain   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  TString fCalorimeter ; // "PHOS" or "EMCAL"

  // Histograms
  TNtuple * fNtTrigger22 ;
  TNtuple * fNtTriggerNN ;

  ClassDef(AliAnaCaloTrigger, 0); // a PHOS photon analysis task 
};
#endif // ALIANACALOTRIGGER_H
