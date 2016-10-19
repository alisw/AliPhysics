#ifndef ALITOFQACHECKER_H
#define ALITOFQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 // 
//  Checks the quality assurance.                                  //
//  By analysis of the histograms & comparing with reference data  //
//  Author S.Arcelli                                               //
//                                                                 // 
/////////////////////////////////////////////////////////////////////

#include "AliQAv1.h"
#include "AliQACheckerBase.h"

//class TFile ; 
//class TH1F ; 
//class TH1I ; 

class AliTOFQAChecker: public AliQACheckerBase {

public:
 AliTOFQAChecker() : AliQACheckerBase("TOF","TOF Quality Assurance Data Maker"), fCheckGolden(kTRUE) {;}          // ctor
 AliTOFQAChecker(const AliTOFQAChecker& qac) : AliQACheckerBase(qac.GetName(), qac.GetTitle()), fCheckGolden(qac.fCheckGolden) {;} // cpy ctor

  AliTOFQAChecker& operator = (const AliTOFQAChecker& qac);
  virtual ~AliTOFQAChecker() {;} // dtor
  void SetCheckGolden(Bool_t doit) {fCheckGolden = doit; return;}

 protected:
  virtual void Check(Double_t * test, AliQAv1::ALITASK_t /*index*/, TObjArray ** list,
		   const AliDetectorRecoParam * recoParam=0);
  Int_t  CheckRaws(TH1* histo, Int_t specie);

 private:
  Bool_t fCheckGolden; //flag to enable automatic checks on the full list of histos

  ClassDef(AliTOFQAChecker, 3)  // description 

};

#endif // AliTOFQAChecker_H
