#ifndef ALIQACHECKERBASE_H
#define ALIQACHECKERBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/*
  Base class for detectors quality assurance checkers 
  Compares Data made by QADataMakers with reference data
  Y. Schutz CERN August 2007
*/


// --- ROOT system ---
#include <TNamed.h>
#include "AliQA.h"
class TFile ; 
class TH1 ; 
class TList ; 

// --- Standard library ---

// --- AliRoot header files ---

class AliQACheckerBase: public TNamed {

public:
  AliQACheckerBase(const char * name = "", const char * title = "") ;          // ctor
  AliQACheckerBase(const AliQACheckerBase& qac) ;   
  AliQACheckerBase& operator = (const AliQACheckerBase& qac) ;
  virtual ~AliQACheckerBase() {;} // dtor

  void   Init(const AliQA::DETECTORINDEX det) ; 
  void   Run(AliQA::ALITASK tsk, TList * list=0x0); 
  void   SetRefandData(TDirectory * ref, TDirectory * data=NULL) { fRefSubDir = ref ;  fDataSubDir = data ; }

protected:
  virtual const Double_t Check() ;
  virtual const Double_t Check(TList * list) ;
  const Double_t DiffC(const TH1 * href, const TH1 * hin) const ;   
  const Double_t DiffK(const TH1 * href, const TH1 * hin) const ;   
  void           Finish() const ; 

  TDirectory * fDataSubDir ; //! directory for the current task directory in the current detector directory in the data file
  TDirectory * fRefSubDir  ; //! directory for the current task directory in the current detector directory in the reference file

  ClassDef(AliQACheckerBase,1)  // description 

};

#endif // AliQUALASSCHECKERBASE_H
