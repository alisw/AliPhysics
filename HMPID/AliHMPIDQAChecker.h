#ifndef ALIHMPIDQACHECKER_H
#define ALIHMPIDQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for HMPID
//


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQACheckerBase.h"

class AliHMPIDQAChecker: public AliQACheckerBase {

public:
  AliHMPIDQAChecker() ;          // ctor
  AliHMPIDQAChecker(const AliHMPIDQAChecker& qac) ; // cpy ctor   
  virtual ~AliHMPIDQAChecker() ; // dtor
  
  virtual void Check(Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * recoParam) ;
  
  Double_t CheckEntries(TObjArray * list) const ;
  Double_t CheckRec(TObjArray *listrec, TObjArray *listref) const ;
  Double_t CheckSim(TObjArray *listsim, TObjArray *listref) const ;

private:
  AliHMPIDQAChecker& operator= (const AliHMPIDQAChecker&); // Not implemented
  Bool_t        fNoReference ; //! flag telling if reference data hqve been found or not  
  TObjArray *   fQARefRec ;    //! Reference data from OCDB 
      
  ClassDef(AliHMPIDQAChecker,1)  // description 

};

#endif // AliHMPIDQAChecker_H
