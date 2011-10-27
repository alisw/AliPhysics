#ifndef ALIQAChecker_H
#define ALIQAChecker_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// class for running the Quality Assurance Checker
// to run type:
//   AliQAChecker qac;
//   qac.Run();
//   qac.SelectDetectors("list of detectors") ; 
//   qac.SelectTargets("list of targets among Hits, Digits, ESD") ; 
//

#include <TNamed.h>
#include <TFile.h>  
class TNtupleD ;

#include "AliQAv1.h"
#include "AliRecoParam.h"
class AliCDBEntry ; 
class AliRunInfo ;
class AliQACheckerBase ; 
class AliDetectorRecoParam ; 

class AliQAChecker: public TNamed {
public:
  AliQAChecker(const char* name = "AliQAChecker", 
		    const char* title = "Quality Assurance checker for Raws, Hits, Digits and ESDs");
  AliQAChecker(const AliQAChecker& qac);
  AliQAChecker& operator = (const AliQAChecker& qac);
  virtual  ~AliQAChecker();

  static  AliQAChecker *     Instance() ;
  AliQACheckerBase *         GetDetQAChecker(Int_t det) ; 
  Bool_t Run(const char * fileName = NULL, AliDetectorRecoParam * recoParam = NULL) ;
  Bool_t Run(AliQAv1::DETECTORINDEX_t det, AliQAv1::TASKINDEX_t task, TObjArray ** list, AliDetectorRecoParam * recoParam = NULL);
  Bool_t Run(AliQAv1::DETECTORINDEX_t det, AliQAv1::TASKINDEX_t task, TNtupleD ** list, AliDetectorRecoParam * recoParam = NULL);
  void   SetRunInfo(AliRunInfo * ei) {fRunInfo = ei;}
  Int_t  GetRunNumber() const { return fRun ; } 
  void   SetRunNumber(Int_t run) { fRun = run ; } 

private:

  void LoadRunInfoFromGRP() ; 

  static AliQAChecker *       fgQAChecker ;             // pointer to the instance of the singleton
  TFile *                     fDataFile ;               //! Data file to check
  AliRunInfo *                fRunInfo ;                //! Event info object 
  Bool_t                      fRunInfoOwner;            //! owns fRunInfo or not
  TFile *                     fRefFile ;                //! Reference Data file 
  TString                     fFoundDetectors ;         //! detectors for which the Quality assurance could be done
  AliQACheckerBase *          fCheckers[AliQAv1::kNDET] ; //! list of detectors checkers
  AliRecoParam::EventSpecie_t fEventSpecie ;            //! event specie deduced from the GRP data
	Int_t                       fRun ;                    //! run number
  ClassDef(AliQAChecker, 1)  // class for running generation, simulation and digitization
};
#endif
