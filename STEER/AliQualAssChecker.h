#ifndef ALIQUALASSChecker_H
#define ALIQUALASSChecker_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// class for running the Quality Assurance Checker
// to run type:
//   AliQualAssChecker qac;
//   qac.Run();
//   qac.SelectDetectors("list of detectors") ; 
//   qac.SelectTargets("list of targets among Hits, Digits, ESD") ; 
//

#include <TNamed.h>
#include <TFile.h>  

#include "AliQualAss.h"
class AliQualAssCheckerBase ; 

class AliQualAssChecker: public TNamed {
public:
  AliQualAssChecker(const char* name = "AliQualAssChecker", 
		    const char* title = "Quality Assurance checker for Hits, Digits and ESDs");
  AliQualAssChecker(const AliQualAssChecker& qac);
  AliQualAssChecker& operator = (const AliQualAssChecker& qac);
  virtual ~AliQualAssChecker();

  AliQualAssCheckerBase * GetDetQualAssChecker(Int_t det) ; 
  TDirectory *            GetRefSubDir(const char * det, const char * task) ;
  static TFile *          GetQAResultFile() ;
  static const char *     GetQAResultFileName() { return fgQAResultFileName.Data() ; }
  void                    SetQAResultDirName(const char * name) ; 
  void                    SetRefDirName(const char * name) ; 

  virtual Bool_t Run();
  
private:
//   AliRunLoader*  LoadRun(const char* mode = "UPDATE") const;
  TFile *      GetDataFile() ; 

  TFile * fDataFile ;                //! Data file to check
  static TFile * fgQAResultFile ;    //! File where to find the QA result
  static TString fgQAResultDirName ; //! directory where to find the QA result
  static TString fgQAResultFileName ;//! file name where to find the QA result
  TString fRefDirName ;              //! name of directory where to find the reference data file
  TString fRefName ;                 //! file name where to find the reference data
  TString fFoundDetectors ;          //! detectors for which the Quality assurance could be done
  AliQualAssCheckerBase * fCheckers[AliQualAss::kNDET] ; //! list of detectors checkers
  ClassDef(AliQualAssChecker, 1)  // class for running generation, simulation and digitization
};

#endif
