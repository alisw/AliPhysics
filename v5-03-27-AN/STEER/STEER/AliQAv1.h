#ifndef ALIQAV1_H
#define ALIQAV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliQAv1.h 31503 2009-03-16 11:01:16Z schutz $ */

//
// Quality Assurance Object
//

#include <TNamed.h> 
#include <TMath.h> 
class TFile ; 
class TH1 ; 

#include "AliLog.h"
#include "AliRecoParam.h"

class AliQAv1 : public TNamed {
public:
  
  enum DETECTORINDEX_t {
    kNULLDET=-1, kITS, kTPC, kTRD, kTOF, kPHOS, kHMPID, kEMCAL, kMUON, kFMD,
    kZDC, kPMD, kT0, kVZERO, kACORDE, kHLT, kGLOBAL, kCORR, kNDET
    #ifdef MFT_UPGRADE
    , kMFT
    #endif 
    };
  enum ALITASK_t {
    kNULLTASK=-1, kRAW, kSIM, kREC, kESD, kANA, kNTASK };
  enum QABIT_t {
    kNULLBit=-1, kINFO, kWARNING, kERROR, kFATAL, kNBIT };
  enum TASKINDEX_t {
    kNULLTASKINDEX=-1, kRAWS, kHITS, kSDIGITS, kDIGITS, kDIGITSR, kRECPOINTS, kTRACKSEGMENTS, kRECPARTICLES, kESDS, kNTASKINDEX };
  enum MODE_t {kNULLMODE=0, kSIMMODE, kRECMODE, kQAMODE};
  
  // Creators - destructors
  AliQAv1(); // beware singleton, not to be used
  AliQAv1(const Int_t qalength, ULong_t * qa, const Int_t eslength, Bool_t * es) ;
  AliQAv1(const ALITASK_t tsk) ;
  AliQAv1(const DETECTORINDEX_t det) ;
  AliQAv1(const AliQAv1& qa) ;   
  AliQAv1& operator = (const AliQAv1& qa) ;
  virtual ~AliQAv1();
  
  static  AliQAv1 *      Instance() ;
  static  AliQAv1 *      Instance(const Int_t qalength, ULong_t * qa, const Int_t eslength, Bool_t * es) ;
  static  AliQAv1 *      Instance(const DETECTORINDEX_t det) ;
  static  AliQAv1 *      Instance(const ALITASK_t tsk) ;
  static  AliQAv1 *      Instance(const TASKINDEX_t tsk) ;
  Bool_t                 CheckFatal() const ;
  static void            Close() ; 
  static const char *    GetAliTaskName(ALITASK_t tsk) ;
  static  TH1 *          GetData(TObjArray** list, Int_t index, AliRecoParam::EventSpecie_t eventSpecie) ; 
  Bool_t *               GetEventSpecies() { return fEventSpecies ; }
  static const TString   GetExpert() { return fgkExpert ; }
  static       UInt_t    GetExpertBit() { return fgkExpertBit ; }
  static       UInt_t    GetImageBit()  { return fgkImageBit ; }
  static const char *    GetImageFileName() { return fgkImageFileName.Data() ; }
  static const char *    GetImageFileFormat() { return fImageFileFormat.Data() ; }
  static const TString   GetLabLocalFile() { return fgkLabLocalFile ; } 
  static const TString   GetLabLocalOCDB() { return fgkLabLocalOCDB ; } 
  static const TString   GetLabAliEnOCDB() { return fgkLabAliEnOCDB ; } 
  static DETECTORINDEX_t GetDetIndex(const char * name) ; 
  static const TString   GetDetName(DETECTORINDEX_t det) { return fgDetNames[det] ; }
  static const char *    GetDetName(Int_t det) ;
  static const TString   GetGRPPath() { return fgGRPPath ; }  
  ULong_t *              GetQA() { return fQA ; }
  static       UInt_t    GetQABit() { return fgkQABit ; }
  static TFile *         GetQADataFile(const char * name, Int_t run) ; 
  static TFile *	     GetQADataFile(const char * fileName) ;
  static const char *    GetQADataFileName(const char * name, Int_t run) 
  {return Form("%s.%s.%d.root", name, fgQADataFileName.Data(), run)  ; }
  static const char *    GetQADataFileName() { return fgQADataFileName.Data() ; }
  static Int_t           GetQADebugLevel() { return fgkQADebugLevel ; }
  static const char *    GetQAName() { return fgkQAName ; } 
  static const char *    GetQACorrName() { return fgkQACorrNtName ; }
  static  UShort_t       GetMaxQAObj() { return fgkMaxQAObjects ; }
  static TFile *         GetQAResultFile() ; 
  static const char  *   GetQAResultFileName() { return (fgQAResultDirName + fgQAResultFileName).Data() ; }
  static const char  *   GetQARefDefaultStorage() { return fgkQARefOCDBDefault.Data() ; }
  static const char  *   GetQARefFileName() { return fgQARefFileName ; }
  static const char  *   GetQARefStorage() { return fgQARefDirName.Data() ; }
  QABIT_t                GetQAStatusBit(AliRecoParam::EventSpecie_t es = AliRecoParam::kDefault, DETECTORINDEX_t det=kNULLDET, ALITASK_t tsk=kNULLTASK) const ; 
  static const char  *   GetRefOCDBDirName() { return fgkRefOCDBDirName.Data() ; }
  static const char  *   GetRefDataDirName() { return fgRefDataDirName.Data() ; }
  //
  static       UInt_t    GetClonedBit()               { return fgkQAClonedBit; }
  static       UInt_t    GetForbidCloningBit()        { return fgkForbidCloningBit; }
  static       UInt_t    GetOrigHistoKeptBit()        { return fgkOrigHistoKeptBit; }
  //
  static     TASKINDEX_t GetTaskIndex(const char * name) ; 
  static       TString   GetTaskName(UInt_t tsk) { return fgTaskNames[tsk] ; }
  static const char *    GetModeName(MODE_t mode) { return (mode == kSIMMODE || mode == kRECMODE || mode == kQAMODE) ? (fgModeNames[mode]).Data() : "" ; }
  Bool_t                 IsEventSpecieSet(AliRecoParam::EventSpecie_t es) const 
  {Int_t ibit=0; while(es!=1<<ibit) ++ibit; return fEventSpecies[ibit];}
  Bool_t                 IsEventSpecieSet(Int_t es) const { return fEventSpecies[es] ; }
  Bool_t                 IsSet(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es, QABIT_t bit) const ;
  Bool_t                 IsSet(DETECTORINDEX_t det, ALITASK_t tsk, Int_t es, QABIT_t bit) const ;
  Bool_t                 IsSetAny(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es) const ;
  Bool_t                 IsSetAny(DETECTORINDEX_t det, AliRecoParam::EventSpecie_t es) const ;
  Long64_t               Merge(const TCollection * list) ; 
  static MODE_t          Mode(TASKINDEX_t task) ;
  void                   Set(QABIT_t bit, AliRecoParam::EventSpecie_t es) ;
  void                   Set(QABIT_t bit, Int_t es) ;
  void                   SetEventSpecie(AliRecoParam::EventSpecie_t es) 
//    {Int_t ibit=0; while(es!=1<<ibit) ++ibit; fEventSpecies[ibit] = kTRUE ; }
  {if(es==AliRecoParam::kDefault) return; Int_t ibit=0; while(es!=1<<ibit) ++ibit; fEventSpecies[ibit] = kTRUE ; }
  static void            SetQADebug() { AliLog::SetGlobalDebugLevel(GetQADebugLevel()); }
  static void            SetQAResultDirName(const char * name) ; 
  static void            SetQARefStorage(const char * name) ; 
  static void            SetQARefDataDirName(AliRecoParam::EventSpecie_t es) { fgRefDataDirName = AliRecoParam::GetEventSpecieName(es) ; }
  static void            SetQARefDataDirName(Int_t es) { fgRefDataDirName = AliRecoParam::GetEventSpecieName(es) ; }
  void                   Show(DETECTORINDEX_t det = kNULLDET) const ;
  void                   ShowAll() const ;
  void                   ShowStatus(DETECTORINDEX_t det, ALITASK_t tsk=kNULLTASK, AliRecoParam::EventSpecie_t es=AliRecoParam::kDefault) const ;
  void                   UnSet(QABIT_t bit, AliRecoParam::EventSpecie_t es) ;
  void                   UnSet(QABIT_t bit, Int_t es) ;
  
private:      
  
  Bool_t                CheckRange(DETECTORINDEX_t det) const ;
  Bool_t                CheckRange(ALITASK_t tsk) const ;
  Bool_t                CheckRange(QABIT_t bit) const ;
  Bool_t                CheckRange(AliRecoParam::EventSpecie_t es) const ;
  const char *          GetBitName(QABIT_t bit) const ;
  ULong_t               GetStatus(DETECTORINDEX_t det, AliRecoParam::EventSpecie_t es) const  { return fQA[det*fNEventSpecies+AliRecoParam::AConvert(es)] ;}
  void                  Finish() const ;  
  ULong_t               Offset(ALITASK_t tsk) const ;
  void                  ShowASCIIStatus(AliRecoParam::EventSpecie_t es, DETECTORINDEX_t det, ALITASK_t tsk, ULong_t status) const ; 
  void                  Reset(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es) ;  

  void                  ResetStatus(DETECTORINDEX_t det) ; 
  void                  Set(DETECTORINDEX_t det) { fDet = det ;}
  void                  Set(ALITASK_t tsk) { fTask = tsk ; AliDebug(GetQADebugLevel(), Form("Ready to set QA status in %s", GetAliTaskName(tsk) )) ; }
  void                  SetStatus(DETECTORINDEX_t det, AliRecoParam::EventSpecie_t es, ULong_t status) { fQA[det*fNEventSpecies+AliRecoParam::AConvert(es)] = status ; }
  void                  SetStatusBit(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es, QABIT_t bit) ;
  void                  UnSetStatusBit(DETECTORINDEX_t det, ALITASK_t tsk, AliRecoParam::EventSpecie_t es, QABIT_t bit) ;
  
  static AliQAv1 *       fgQA		                  ; // pointer to the instance of the singleton
  Int_t                  fNdet     	              ; // number of detectors
  Int_t                  fNEventSpecies           ; // number of Event Species (see AliRecoParam)
  Int_t                  fLengthQA                ; // Auxiliary length of fQA
  ULong_t    *           fQA		                  ; //[fLengthQA]  the status word 4 bits for SIM, REC, ESD, ANA each
  DETECTORINDEX_t fDet		                        ; //! the current detector (ITS, TPC, ....)
  ALITASK_t       fTask	                          ; //! the current environment (SIM, REC, ESD, ANA)
  AliRecoParam::EventSpecie_t fEventSpecie        ; //! the current event specie
  static TString         fgDetNames[]	            ; //! list of detector names   
  static TString         fgGRPPath                ; //! path of the GRP object in OCDB
  static TFile *         fgQADataFile	            ; //! the output file where the quality assurance maker store their results
  static TString         fgQADataFileName         ; //! the name of the file where the quality assurance maker store their results
  static TFile *         fgQARefFile	            ; //! the output file where the quality assurance maker store their results
  static TString         fgQARefDirName	          ; //! name of directory where to find the reference data file
  static TString         fgQARefFileName          ; //! file name where to find the reference data
  static TFile *         fgQAResultFile           ; //! File where to find the QA result
  static TString         fgQAResultDirName        ; //! the location of the output file where the QA results are stored  
  static TString         fgQAResultFileName       ; //! the output file where the QA results are stored  
  static TString         fgRTNames[]	            ; //! list of Run Type names   
  static TString         fgTaskNames[]            ; //! list of tasks names   
  static TString         fgModeNames[]	          ; //! list of mode names   
  static const TString   fgkExpert                ; //! name for the expert directory
  static const UInt_t    fgkExpertBit             ; //! TObject bit identifing the object as "expert"
  static const UInt_t    fgkImageBit              ; //! TObject bit identifing the object to be plotted on the QA image
  static const TString   fgkLabLocalFile          ; //! label to identify a file as local 
  static const TString   fgkLabLocalOCDB          ; //! label to identify a file as local OCDB 
  static const TString   fgkLabAliEnOCDB          ; //! label to identify a file as AliEn OCDB 
  static const TString   fgkRefFileName           ; //! name of Reference File Name 
  static const UInt_t    fgkQABit                 ; //! bit in the QA data object which is set when Checker does not return 0
  static const Int_t     fgkQADebugLevel          ; //! debug level used for QA verbosity
  static const TString   fgkQAName                ; //! name of QA object 
  static const TString   fgkQACorrNtName          ; //! name of QA Correlation Ntuple
  static const TString   fgkRefOCDBDirName        ; //! name of Reference directory name in OCDB  	
  static       TString   fgRefDataDirName         ; //! name of Reference directory name in OCDB for data  	
  static const TString   fgkQARefOCDBDefault      ; //! default storage for QA in OCDB 
  Bool_t *               fEventSpecies            ; //[fNEventSpecies] list of event species encountered in a run
  static const TString   fgkImageFileName           ; //! name of the file that contains all the QA images
  static const TString   fImageFileFormat         ; //! format of the file that contains all the QA images
  static const UShort_t  fgkMaxQAObjects          ;//! maximum number of QA objects allowed dor each task (Raws, digits,....) 
  //
  static const UInt_t    fgkQAClonedBit           ;//! flag that the histrogram was cloned per trigger class
  static const UInt_t    fgkForbidCloningBit      ;//! flag that the histogram is forbiden to clone per trigger class
  static const UInt_t    fgkOrigHistoKeptBit      ;//! flag that the histogram was clonned but original is kept
 ClassDef(AliQAv1,3)  //ALICE Quality Assurance Object
};
#endif
