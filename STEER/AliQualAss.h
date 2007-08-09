#ifndef ALIQualAss_H
#define ALIQualAss_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Quality Assurance Object
//

#include <TNamed.h>
#include "AliLog.h" 
class TFile ; 

class AliQualAss : public TNamed {
public:

  enum DETECTORINDEX {
    kNULLDET=-1, kITS, kTPC, kTRD, kTOF, kPHOS, kHMPID, kEMCAL, kMUON, kFMD,
    kZDC, kPMD, kT0, kVZERO, kACORDE, kHLT, kNDET
  };
  enum ALITASK {
    kNULLTASK=-1, kSIM, kREC, kESD, kANA, kNTASK
  };
  enum QABIT {
    kNULLBit=-1, kINFO, kWARNING, kERROR, kFATAL, kNBIT
  };
  
  enum TASKINDEX {
    kHITS, kSDIGITS, kDIGITS, kRECPOINTS, kTRACKSEGMENTS, kRECPARTICLES, kESDS
  };
  
 // Creators - destructors
  AliQualAss(); // beware singleton, not to be used
  AliQualAss(ALITASK tsk) ;
  AliQualAss(const AliQualAss& qa) ;   
  AliQualAss& operator = (const AliQualAss& qa) ;
  virtual ~AliQualAss();
 
  static  AliQualAss *   Instance() ;
  static  AliQualAss *   Instance(DETECTORINDEX det) ;
  static  AliQualAss *   Instance(ALITASK tsk) ;
  const Bool_t           CheckFatal() const ;
  static const char *    GetOutputName() { return fgOutputName.Data() ; }
  static TFile *         GetQADMOutFile() ; 
  void                   Set(QABIT bit) ;
  void                   Show() const { ShowStatus(fDet) ; }
  void                   ShowAll() const ;

private:      

  const Bool_t         CheckRange(DETECTORINDEX det) const ;
  const Bool_t         CheckRange(ALITASK tsk) const ;
  const Bool_t         CheckRange(QABIT bit) const ;
  const char *         GetDetectorName(DETECTORINDEX det) const ;
  const char *         GetTaskName(ALITASK tsk) const ;
  const char *         GetBitName(QABIT bit) const ;
  const ULong_t        GetStatus(DETECTORINDEX det) const  { return fQA[det] ;}
  void                 Finish() const ;  
  const Bool_t         IsSet(DETECTORINDEX det, ALITASK tsk, QABIT bit) const ;
  const ULong_t        Offset(ALITASK tsk) const ;
  virtual void         ShowStatus(DETECTORINDEX det) const ;
  void                 ResetStatus(DETECTORINDEX det) { fQA[det] = 0 ; }
  void                 Set(DETECTORINDEX det) { fDet = det ;}
  void                 Set(ALITASK tsk) { fTask = tsk ; AliInfo(Form("Ready to set QA status in %s\n", GetTaskName(tsk) )) ; }
  void                 SetStatus(DETECTORINDEX det, UShort_t status) { fQA[det] = status ; }
  void                 SetStatusBit(DETECTORINDEX det, ALITASK tsk, QABIT bit) ;

  static AliQualAss *fgQA         ; // pointer to the instance of the singleton
  Int_t              fNdet        ; // # of detectors
  ULong_t *          fQA          ; //[kNDET] the status word 4 bits for SIM, REC, ESD, ANA each
  DETECTORINDEX      fDet         ; //!  the current detector (ITS, TPC, ....)
  ALITASK            fTask        ; //!  the current environment (SIM, REC, ESD, ANA)
  static TFile *     fgOutput     ; //! the output file where the quality assurance maker store their results
  static TString     fgOutputName ; //! the output name file where the quality assurance maker store their results
    
  ClassDef(AliQualAss,1)  //ALICE Quality Assurance Obbject
};
#endif
