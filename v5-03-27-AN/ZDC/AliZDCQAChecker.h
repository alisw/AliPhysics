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
  AliZDCQAChecker() : AliQACheckerBase("ZDC","ZDC Quality Assurance Data Maker") {;}          // ctor
  virtual ~AliZDCQAChecker() {;} // dtor

 protected:

  virtual void Check(Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list,
      const AliDetectorRecoParam * /*recoParam*/); 
  void SetupHisto(const TObjArray& messages, TH1& histo, Float_t& code);

 private:  
  AliZDCQAChecker(const AliZDCQAChecker& qac); // cpy ctor   
  AliZDCQAChecker& operator= (const AliZDCQAChecker & /*checker*/);
  ClassDef(AliZDCQAChecker,1)  // description 

};

#endif // AliZDCQAChecker_H
