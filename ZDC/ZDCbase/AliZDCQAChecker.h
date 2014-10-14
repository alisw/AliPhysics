#ifndef ALIZDCQACHECKER_H
#define ALIZDCQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
//  						    //
//  Checks the quality assurance.                   //
//  C. Oppedisano Chiara.Oppedisano@to.infn.it      //
//  						    //
//////////////////////////////////////////////////////

#include "AliQAv1.h"
#include "AliQACheckerBase.h"

class TObjArray;

class AliZDCQAChecker: public AliQACheckerBase {

public:
  AliZDCQAChecker();             // ctor
  virtual ~AliZDCQAChecker() {;} // dtor

 protected:

  virtual void Check(Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list,
      const AliDetectorRecoParam * /*recoParam*/); 
  void SetupHisto(const TObjArray& messages, TH1& histo, Float_t& code);
 
  void    GetThresholds();
  void    PrintThresholds();

 private:  
  AliZDCQAChecker(const AliZDCQAChecker& qac); // cpy ctor   
  AliZDCQAChecker& operator= (const AliZDCQAChecker & /*checker*/);

  TObjArray *fQAThresholds;    		//! Reference data from OCDB 
  Double_t    fZDCQAThr_ZNCTDCRefThr; 	// TDC reference value for QA checks
  Double_t    fZDCQAThr_ZPCTDCRefThr; 	// TDC reference value for QA checks
  Double_t    fZDCQAThr_ZNATDCRefThr; 	// TDC reference value for QA checks
  Double_t    fZDCQAThr_ZPATDCRefThr; 	// TDC reference value for QA checks
  Double_t    fZDCQAThr_ZEM1TDCRefThr; 	// TDC reference value for QA checks
  Double_t    fZDCQAThr_ZEM2TDCRefThr; 	// TDC reference value for QA checks
  
  ClassDef(AliZDCQAChecker,1)  // description 

};

#endif // AliZDCQAChecker_H
