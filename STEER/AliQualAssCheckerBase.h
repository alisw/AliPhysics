#ifndef ALIQUALASSCHECKERBASE_H
#define ALIQUALASSCHECKERBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/*
  Base class for detectors quality assurance checkers 
  Compares Data made by QualAssDataMakers with reference data
  Y. Schutz CERN August 2007
*/


// --- ROOT system ---
#include <TNamed.h>
#include "AliQualAss.h"
class TFile ; 
class TH1 ; 

// --- Standard library ---

// --- AliRoot header files ---

class AliQualAssCheckerBase: public TNamed {

public:
  AliQualAssCheckerBase(const char * name = "", const char * title = "") ;          // ctor
  AliQualAssCheckerBase(const AliQualAssCheckerBase& qac) ;   
  AliQualAssCheckerBase& operator = (const AliQualAssCheckerBase& qac) ;
  virtual ~AliQualAssCheckerBase() {;} // dtor

  void   Run(AliQualAss::ALITASK tsk); 
  void   Init(const AliQualAss::DETECTORINDEX det) ; 
  void   SetRefandData(TDirectory * ref, TDirectory * data) { fRefSubDir = ref ;  fDataSubDir = data ; }

protected:
  virtual const Double_t Check() ;
  const Double_t DiffC(const TH1 * href, const TH1 * hin) const ;   
  const Double_t DiffK(const TH1 * href, const TH1 * hin) const ;   
  void           Finish() const ; 

  TDirectory * fDataSubDir ; //! directory for the current task directory in the current detector directory in the data file
  TDirectory * fRefSubDir  ; //! directory for the current task directory in the current detector directory in the reference file

  ClassDef(AliQualAssCheckerBase,1)  // description 

};

#endif // AliQUALASSCHECKERBASE_H
