#ifndef ALITOFQACHECKER_H
#define ALITOFQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Checks the quality assurance. 
  By analysis of the histograms & comparing with reference data
  Author S.Arcelli
*/

class TFile ; 
class TH1F ; 
class TH1I ; 
#include "AliQACheckerBase.h"

class AliTOFQAChecker: public AliQACheckerBase {

public:
  AliTOFQAChecker() : AliQACheckerBase("TOF","TOF Quality Assurance Data Maker") {;}          // ctor
  AliTOFQAChecker(const AliTOFQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()) {;} // cpy ctor   
  AliTOFQAChecker& operator = (const AliTOFQAChecker& qac) ;
  virtual ~AliTOFQAChecker() {;} // dtor

private:
  
  ClassDef(AliTOFQAChecker,1)  // description 

};

#endif // AliTOFQAChecker_H
