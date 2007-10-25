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

#include "AliQA.h"
class AliQACheckerBase ; 

class AliQAChecker: public TNamed {
public:
  AliQAChecker(const char* name = "AliQAChecker", 
		    const char* title = "Quality Assurance checker for Raws, Hits, Digits and ESDs");
  AliQAChecker(const AliQAChecker& qac);
  AliQAChecker& operator = (const AliQAChecker& qac);
  virtual ~AliQAChecker();

  static  AliQAChecker * Instance() ;
  AliQACheckerBase *     GetDetQAChecker(Int_t det) ; 
  TDirectory *           GetRefSubDir(const char * det, const char * task) ;
  static TFile *         GetQAResultFile() ;
  static const char *    GetQAResultFileName() { return fgQAResultFileName.Data() ; }
  void                   SetQAResultDirName(const char * name) ; 
  void                   SetRefDirName(const char * name) ; 

  virtual Bool_t Run(const char * fileName = NULL) ;
  virtual Bool_t Run(AliQA::DETECTORINDEX det, AliQA::TASKINDEX task, TList * list);

private:
  TFile *      GetDataFile(const char * fileName) ; 

  static AliQAChecker *fgQAChecker ; // pointer to the instance of the singleton
  TFile * fDataFile ;                     //! Data file to check
  static TFile * fgQAResultFile ;         //! File where to find the QA result
  static TString fgQAResultDirName ;      //! directory where to find the QA result
  static TString fgQAResultFileName ;     //! file name where to find the QA result
  TString fRefDirName ;                   //! name of directory where to find the reference data file
  TString fRefName ;                      //! file name where to find the reference data
  TString fFoundDetectors ;               //! detectors for which the Quality assurance could be done
  AliQACheckerBase * fCheckers[AliQA::kNDET] ; //! list of detectors checkers
  ClassDef(AliQAChecker, 1)  // class for running generation, simulation and digitization
};
#endif
