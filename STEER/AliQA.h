#ifndef ALIQA_H
#define ALIQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Quality Assurance Object
//

#include <TNamed.h> 
class TFile ; 

#include "AliLog.h"

class AliQA : public TNamed {
public:

  enum DETECTORINDEX {
    kNULLDET=-1, kITS, kTPC, kTRD, kTOF, kPHOS, kHMPID, kEMCAL, kMUON, kFMD,
    kZDC, kPMD, kT0, kVZERO, kACORDE, kHLT, kNDET
  };
  enum ALITASK {
    kNULLTASK=-1, kRAW, kSIM, kREC, kESD, kANA, kNTASK
  };
  enum QABIT {
    kNULLBit=-1, kINFO, kWARNING, kERROR, kFATAL, kNBIT
  };
  
  enum TASKINDEX {
    kRAWS, kHITS, kSDIGITS, kDIGITS, kRECPOINTS, kTRACKSEGMENTS, kRECPARTICLES, kESDS
  };
  
 // Creators - destructors
  AliQA(); // beware singleton, not to be used
  AliQA(const ALITASK tsk) ;
  AliQA(const DETECTORINDEX det) ;
  AliQA(const AliQA& qa) ;   
  AliQA& operator = (const AliQA& qa) ;
  virtual ~AliQA();
 
  static  AliQA *   Instance() ;
  static  AliQA *   Instance(const DETECTORINDEX det) ;
  static  AliQA *   Instance(const ALITASK tsk) ;
  const Bool_t           AddQAData2CDB(const char * defSto) const ;
  const Bool_t           CheckFatal() const ;
  static void            Close() ; 
  static const char *    GetAliTaskName(ALITASK tsk) ;
  static const TString   GetDetName(DETECTORINDEX det) { return fgDetNames[det] ; }
  static const TString   GetTaskName(TASKINDEX tsk) { return fgTaskNames[tsk] ; }
  static const char *    GetDetName(Int_t det) ;
  static const char *    GetQADataFileName() { return fgQADataFileName.Data() ; }
  static TFile *         GetQADataFile(const char * name, const Int_t run, const Int_t cycle) ; 
  static TFile *		 GetQADataFile(const char * fileName) ;
  static TFile *         GetQAResultFile() ; 
  static TFile *         GetQARefFile() ; 
  static const char  *   GetQAResultFileName() { return (fgQAResultDirName + fgQAResultFileName).Data() ; }
  static const char  *   GetQARefFileName() { return (fgQARefDirName + fgQARefFileName).Data() ; }
  const Bool_t           IsSet(DETECTORINDEX det, ALITASK tsk, QABIT bit) const ;
  void                   Set(QABIT bit) ;
  static void			 SetQAResultDirName(const char * name) ; 
  static void            SetQARefDir(const char * name) ; 
  void                   Show() const { ShowStatus(fDet) ; }
  void                   ShowAll() const ;

private:      

  const Bool_t         CheckRange(DETECTORINDEX det) const ;
  const Bool_t         CheckRange(ALITASK tsk) const ;
  const Bool_t         CheckRange(QABIT bit) const ;
  const char *         GetBitName(QABIT bit) const ;
  const ULong_t        GetStatus(DETECTORINDEX det) const  { return fQA[det] ;}
  void                 Finish() const ;  
  const ULong_t        Offset(ALITASK tsk) const ;
  virtual void         ShowStatus(DETECTORINDEX det) const ;
  void                 ResetStatus(DETECTORINDEX det) { fQA[det] = 0 ; }
  void                 Set(DETECTORINDEX det) { fDet = det ;}
  void                 Set(ALITASK tsk) { fTask = tsk ; AliInfo(Form("Ready to set QA status in %s", GetAliTaskName(tsk) )) ; }
  void                 SetStatus(DETECTORINDEX det, UShort_t status) { fQA[det] = status ; }
  void                 SetStatusBit(DETECTORINDEX det, ALITASK tsk, QABIT bit) ;

  static AliQA *fgQA					; // pointer to the instance of the singleton
  Int_t              fNdet				; // number of detectors
  ULong_t    *       fQA				; //[fNdet] the status word 4 bits for SIM, REC, ESD, ANA each
  DETECTORINDEX      fDet				; //!  the current detector (ITS, TPC, ....)
  ALITASK            fTask				; //!  the current environment (SIM, REC, ESD, ANA)
  static TString     fgDetNames[]		; //! list of detector names   
  static TFile *     fgQADataFile		; //! the output file where the quality assurance maker store their results
  static TString     fgQADataFileName	; //! the name of the file where the quality assurance maker store their results
  static TFile *     fgQARefFile		; //! the output file where the quality assurance maker store their results
  static TString     fgQARefDirName		; //! name of directory where to find the reference data file
  static TString     fgQARefFileName	; //! file name where to find the reference data
  static TFile *     fgQAResultFile     ; //! File where to find the QA result
  static TString     fgQAResultDirName  ; //! the location of the output file where the QA results are stored  
  static TString     fgQAResultFileName ; //! the output file where the QA results are stored  
  static TString     fgTaskNames[]		; //! list of tasks names   

 ClassDef(AliQA,1)  //ALICE Quality Assurance Object
};
#endif
