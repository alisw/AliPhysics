#ifndef ALIQACHECKERBASE_H
#define ALIQACHECKERBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Base class for detectors quality assurance checkers 
//  Compares Data made by QADataMakers with reference data
//  Y. Schutz CERN August 2007
//


// --- ROOT system ---
#include <TNamed.h>
#include "AliQAv1.h"
class TFile ; 
class TH1 ; 
class TObjArray ; 
class TDirectory ; 
class TNtupleD ;

// --- Standard library ---

// --- AliRoot header files ---

class AliQACheckerBase: public TNamed {

public:
  AliQACheckerBase(const char * name = "", const char * title = "") ;          // ctor
  AliQACheckerBase(const AliQACheckerBase& qac) ;   
  AliQACheckerBase& operator = (const AliQACheckerBase& qac) ;
  virtual ~AliQACheckerBase() ; // dtor

  virtual void   Init(const AliQAv1::DETECTORINDEX_t det)   { AliQAv1::Instance(det) ; }
  void           Run(AliQAv1::ALITASK_t tsk, TObjArray ** list = NULL); 
  void           Run(AliQAv1::ALITASK_t /*tsk*/, TNtupleD ** /*nt*/) {;} 
  void           SetHiLo(Float_t * hiValue, Float_t * lowValue) ; 
  void           SetRefandData(TDirectory * ref, TObjArray ** refOCDB, TDirectory * data=NULL) { fRefSubDir = ref ;  fRefOCDBSubDir = refOCDB, fDataSubDir = data ; }

protected:
  virtual      Double_t * Check(AliQAv1::ALITASK_t index) ;
  virtual      Double_t * Check(AliQAv1::ALITASK_t, TObjArray **) ; 

  Double_t     DiffC(const TH1 * href, const TH1 * hin) const ;   
  Double_t     DiffK(const TH1 * href, const TH1 * hin) const ;   
  void         Finish() const ; 
  virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const ;	

  TDirectory  * fDataSubDir     ; //! directory for the current task directory in the current detector directory in the data file
  TDirectory  * fRefSubDir      ; //! directory for the current task directory in the current detector directory in the reference file
  TObjArray   ** fRefOCDBSubDir ; //! Entry in OCDB for the current detector 
  Float_t     * fLowTestValue   ; // array of lower bounds for INFO, WARNING, ERROR, FATAL   
  Float_t     * fUpTestValue    ; // array of upper bounds for INFO, WARNING, ERROR, FATAL   

  ClassDef(AliQACheckerBase,1)  // description 

};

#endif // AliQUALASSCHECKERBASE_H
