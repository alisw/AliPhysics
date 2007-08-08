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
class AliRunLoader; 

class AliQualAssChecker: public TNamed {
public:
  AliQualAssChecker(const char* name = "AliQualAssChecker", 
		    const char* title = "Quality Assurance checker for Hits, Digits and ESDs");
  AliQualAssChecker(const AliQualAssChecker& qac);
  AliQualAssChecker& operator = (const AliQualAssChecker& qac);
  virtual ~AliQualAssChecker();

  static TFile *      GetDataFile()    { return TFile::Open(AliQualAss::GetOutputName()) ; }
  static TFile *      GetRefFile()     { return TFile::Open(GetRefFileName()) ; }
  static const char * GetRefFileName() { return fgRefName.Data() ; }
  static TFile * GetOutFile() ;
  static const char * GetOutFileName() { return fgOutName.Data() ; }
  void   SetOutDir(const char * outDir) ; 
  void   SetRefDir(const char * refDir) ; 
  void   SetGAliceFile(const char* fileName) ;

  virtual Bool_t Run();
  
private:
  AliRunLoader*  LoadRun(const char* mode = "UPDATE") const;
   
  TString        fGAliceFileName ;    // name of the galice file
  static TFile * fgOutFile ;          // File where to find the QA result
  static TString fgOutDir ;           // directory where to find the QA result
  static TString fgOutName ;          // file name where to find the QA result
  static TString fgRefDir ;           // directory where to find the reference data
  static TString fgRefName ;          // file name where to find the reference data
  Bool_t         fStopOnError;        // stop or continue on errors

  ClassDef(AliQualAssChecker, 1)  // class for running generation, simulation and digitization
};

#endif
