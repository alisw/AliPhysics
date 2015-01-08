#ifndef ALIPHOSSDigitizer_H
#define ALIPHOSSDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.29  2007/10/10 09:05:10  schutz
 * Changing name QualAss to QA
 *
 * Revision 1.28  2007/09/30 17:08:20  schutz
 * Introducing the notion of QA data acquisition cycle (needed by online)
 *
 * Revision 1.27  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.26  2006/08/28 10:01:56  kharlov
 * Effective C++ warnings fixed (Timur Pocheptsov)
 *
 * Revision 1.25  2005/11/30 18:56:26  schutz
 * Small corrections to fix compilation errors
 *
 * Revision 1.24  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Class for making SDigits in PHOS      
// A Summable Digits is the sum of all hits originating 
// from one primary in one active cell
//*--
//*-- Author: Dmitri Peressounko(SUBATECH & KI)


// --- ROOT system ---
#include "TNamed.h"
#include "AliConfig.h"
class TFile ;


// --- Standard library ---

// --- AliRoot header files ---
//class AliPHOSQADataMaker ; 

class AliPHOSSDigitizer: public TNamed {

public:
  AliPHOSSDigitizer() ;          // ctor
  AliPHOSSDigitizer(const char * alirunFileName, const char * eventFolderName = AliConfig::GetDefaultEventFolderName()) ; 
  AliPHOSSDigitizer(const AliPHOSSDigitizer& sd) ;   
  AliPHOSSDigitizer& operator = (const AliPHOSSDigitizer& sd) ;

  virtual ~AliPHOSSDigitizer(); // dtor

  virtual void   Digitize(Option_t *option); 
  Int_t          GetSDigitsInRun() const {return fSDigitsInRun ;}  
  virtual void   Print(const Option_t * = "") const ;
  void           SetEventFolderName(TString name) { fEventFolderName = name ; }
  void           SetEventRange(Int_t first=0, Int_t last=-1) {fFirstEvent=first; fLastEvent=last; }

  Bool_t operator == (const AliPHOSSDigitizer & sd) const ;

 
private:

  void     Init() ;
  void     InitParameters() ;
  void     PrintSDigits(Option_t * option) ;
  void     Unload() const ;


private:
  Float_t fPrimThreshold ;  // To store primari if Elos > threshold
  Bool_t  fDefaultInit;     //! Says if the task was created by defaut ctor (only parameters are initialized)
  TString fEventFolderName; // event folder name
  Bool_t  fInit ;           //! tells if initialisation wennt OK, will revent exec if not
  Int_t   fSDigitsInRun ;   //! Total number of sdigits in one run
  Int_t   fFirstEvent;      // first event to process
  Int_t   fLastEvent;       // last  event to process

  ClassDef(AliPHOSSDigitizer,6)  // description 

};

#endif // AliPHOSSDigitizer_H
