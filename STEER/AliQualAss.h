#ifndef ALIQualAss_H
#define ALIQualAss_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Quality Assurance Object
//

#include <TNamed.h> 
class TFile ; 

#include "AliLog.h"

class AliQualAss : public TNamed {
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
  AliQualAss(); // beware singleton, not to be used
  AliQualAss(const ALITASK tsk) ;
  AliQualAss(const DETECTORINDEX det) ;
  AliQualAss(const AliQualAss& qa) ;   
  AliQualAss& operator = (const AliQualAss& qa) ;
  virtual ~AliQualAss();
 
  static  AliQualAss *   Instance() ;
  static  AliQualAss *   Instance(const DETECTORINDEX det) ;
  static  AliQualAss *   Instance(const ALITASK tsk) ;
  const Bool_t           CheckFatal() const ;
  static const char *    GetAliTaskName(ALITASK tsk) ;
  static const char *    GetDataName() { return fgDataName.Data() ; }
  static const TString   GetDetName(DETECTORINDEX det) { return fgDetNames[det] ; }
  static const TString   GetTaskName(TASKINDEX tsk) { return fgTaskNames[tsk] ; }
  static const char *    GetDetName(Int_t det) ;
  static TFile *         GetQADMOutFile(const char * name, const Int_t run, const Int_t cycle) ; 
  void                   Set(QABIT bit) ;
  void                   Show() const { ShowStatus(fDet) ; }
  void                   ShowAll() const ;
  void                   print() { printf("%d %x\n", kNDET, fQA) ; } 

private:      

  const Bool_t         CheckRange(DETECTORINDEX det) const ;
  const Bool_t         CheckRange(ALITASK tsk) const ;
  const Bool_t         CheckRange(QABIT bit) const ;
  const char *         GetBitName(QABIT bit) const ;
  const ULong_t        GetStatus(DETECTORINDEX det) const  { return fQA[det] ;}
  void                 Finish() const ;  
  const Bool_t         IsSet(DETECTORINDEX det, ALITASK tsk, QABIT bit) const ;
  const ULong_t        Offset(ALITASK tsk) const ;
  virtual void         ShowStatus(DETECTORINDEX det) const ;
  void                 ResetStatus(DETECTORINDEX det) { fQA[det] = 0 ; }
  void                 Set(DETECTORINDEX det) { fDet = det ;}
  void                 Set(ALITASK tsk) { fTask = tsk ; AliInfo(Form("Ready to set QA status in %s", GetAliTaskName(tsk) )) ; }
  void                 SetStatus(DETECTORINDEX det, UShort_t status) { fQA[det] = status ; }
  void                 SetStatusBit(DETECTORINDEX det, ALITASK tsk, QABIT bit) ;

  static AliQualAss *fgQA         ; // pointer to the instance of the singleton
  ULong_t    *       fQA          ; //[kNDET] the status word 4 bits for SIM, REC, ESD, ANA each
  DETECTORINDEX      fDet         ; //!  the current detector (ITS, TPC, ....)
  ALITASK            fTask        ; //!  the current environment (SIM, REC, ESD, ANA)
  static TFile *     fgDataFile   ; //! the output file where the quality assurance maker store their results
  static TString     fgDataName   ; //! the name of the file where the quality assurance maker store their results
  static TString     fgDetNames[] ; //! list of detector names   
  static TString     fgTaskNames[]; //! list of tasks names   
 ClassDef(AliQualAss,1)  //ALICE Quality Assurance Object
};
#endif
