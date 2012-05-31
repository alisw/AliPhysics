#ifndef ALIQAMANAGER_H
#define ALIQAMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliQAManager.h 30796 2009-01-28 11:05:10Z schutz $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running the QA makers                                           //
//                                                                           //1
//   AliQAManager qas;                                                //
//   qas.Run(AliQAv1::kRAWS, rawROOTFileName);                                 //
//   qas.Run(AliQAv1::kHITS);                                                  //
//   qas.Run(AliQAv1::kSDIGITS);                                               //
//   qas.Run(AliQAv1::kDIGITS);                                                //
//   qas.Run(AliQAv1::kRECPOINTS);                                             //
//   qas.Run(AliQAv1::kESDS);                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include "AliQAv1.h"
#include "AliLoader.h"
#include "AliCDBManager.h"
#include "AliRecoParam.h"
#include "AliEventInfo.h"
 
class TCanvas ; 
class AliESDEvent ; 
class AliDetectorRecoParam ; 
class AliESDEvent ; 
class AliQADataMaker ;
class AliRawReader ;  
class AliRunLoader ; 
class AliCorrQADataMakerRec ;

class AliQAManager : public AliCDBManager {
public:
  static void      Destroy() ;
  void             EndOfCycle(TObjArray * detArray=0x0) ; 
  void             EndOfCycle(TString detectors) ; 
	UInt_t           GetCurrentEvent() const { return fCurrentEvent ; }
	TObjArray *      GetFromOCDB(AliQAv1::DETECTORINDEX_t det, AliQAv1::TASKINDEX_t task, const Char_t * year) const ; 
  const AliEventInfo *   GetEventInfo() const { return fEventInfo ; }
  AliRecoParam::EventSpecie_t GetEventSpecieFromESD() ;
  TCanvas **       GetImage(Char_t * detName) ;
  const Char_t *   GetMode(){ return fMode.Data() ; }
  AliQAv1     *    GetQA(UInt_t run, UInt_t evt) ; 
	AliQADataMaker * GetQADataMaker(const Int_t iDet) ; 
	void             Increment(const AliQAv1::TASKINDEX_t taskIndex = AliQAv1::kNULLTASKINDEX) ;
	void             InitQADataMaker(UInt_t run, TObjArray * detArray=0x0) ;
  Bool_t           IsSaveData() { return fSaveData ; } 
	Bool_t           IsSelected(const Char_t * detName)  ;
	Bool_t           Merge(Int_t runNumber = -1, const char *fileName = NULL) const ;  
  void             MergeCustom() const ;
  Bool_t           MergeXML(const Char_t * collection, const Char_t * subFile = 0, const Char_t * outFile = 0) ; 
  static           AliQAManager * QAManager(AliQAv1::MODE_t = AliQAv1::kNULLMODE, TMap *entryCache = NULL, Int_t run = -1) ;
  static           AliQAManager * QAManager(AliQAv1::TASKINDEX_t task) ;  
	void             Reset(const Bool_t sameCycle = kFALSE) ;  
  void             ResetDetectors(AliQAv1::TASKINDEX_t task, AliQAv1::DETECTORINDEX_t det=AliQAv1::kNULLDET) ; 
	TString          Run(const Char_t * detectors, const AliQAv1::TASKINDEX_t taskIndex=AliQAv1::kNULLTASKINDEX, Bool_t const sameCycle = kFALSE, const Char_t * fileName = NULL) ; 
	TString          Run(const Char_t * detectors, AliRawReader * rawReader, Bool_t const sameCycle = kFALSE) ; 
	TString          Run(const Char_t * detectors, const Char_t * filename, Bool_t const sameCycle = kFALSE) ;
	void             RunOneEvent(AliRawReader * rawReader) ; 
	void             RunOneEventInOneDetector(Int_t det, TTree * tree) ; 
	void             RunOneEvent(AliESDEvent *& esd, AliESDEvent *& hltesd)  ;
	Bool_t           Save2OCDB(const Int_t runNumber, AliRecoParam::EventSpecie_t es, const Char_t * year = "08", const Char_t * detectors = "ALL") const ; 
	void             SetActiveDetectors(TString aDet) { fDetectors = aDet ;  }
  void             SetCheckerExternParam(AliQAv1::DETECTORINDEX_t det, TList * parameterList) ;  
	void             SetCycleLength(const AliQAv1::DETECTORINDEX_t det, const Int_t cycle) { fQACycles[det] = cycle ; }
	void             SetWriteExpert(const AliQAv1::DETECTORINDEX_t det) { fQAWriteExpert[det] = kTRUE ; }
  void             SetEventInfo(AliEventInfo *info) { fEventInfo = info ;} 
	void             SetEventRange(UInt_t first, UInt_t last) { fFirstEvent = first ; fMaxEvents = last - first + 1 ; }    
  void             SetEventSpecie(AliRecoParam::EventSpecie_t es) ; 
	void             SetFirsEvent(UInt_t first) { fFirstEvent = first ; }      
	void             SetMaxEvents(UInt_t max) { fMaxEvents = max ; }      
	void             SetNewCycle() { fCycleSame = kTRUE ; }
  void             SetPrintImage(Bool_t opt = kTRUE) { fPrintImage = opt ; }
	void             SetRecoParam(const Int_t det, const AliDetectorRecoParam *par) ;
	void             SetRunLoader(AliRunLoader * rl) { fRunLoader = rl ; }
  void             SetSaveData(Bool_t opt = kTRUE ) { fSaveData = opt ; }
	void             SetTasks(TString tasks) { fTasks = tasks ; }
  void             SetWriteExpert() ; 
  void             ShowQA() ; 
  
private: 
  AliQAManager() ; 
	AliQAManager(AliQAv1::MODE_t mode, const Char_t * gAliceFilename = "galice.root") ; 
	AliQAManager(const AliQAManager & qas) ; 
	AliQAManager & operator = (const AliQAManager & qas) ; 
  ~AliQAManager() ; 
  
	Bool_t			DoIt(const AliQAv1::TASKINDEX_t taskIndex) ;
	AliLoader * GetLoader(Int_t iDet) ; 
	Int_t       GetQACycles(const Int_t iDet) const { return fQACycles[iDet] ; }
	Bool_t 		  InitQA(const AliQAv1::TASKINDEX_t taskIndex, const  Char_t * fileName = NULL) ;
	Bool_t      InitRunLoader() ; 
	Bool_t      Finish(const AliQAv1::TASKINDEX_t taskIndex) ;
	Bool_t      MergeData(const Int_t runNumber, const char *fileName = NULL) const ;  
	Bool_t      MergeResults(const Int_t runNumber) const ;  
	Bool_t      SaveIt2OCDB(const Int_t runNumber, TFile * inputFile, const Char_t * year, AliRecoParam::EventSpecie_t es) const ;  

 	static AliQAManager*        fgQAInstance;                   // AliQAManager instance
	UInt_t                      fCurrentEvent ;                 //! event counter
	Bool_t                      fCycleSame ;                    //! true if 2 consecutive data making for a same detector   
	TString                     fDetectors ;                    //! list of active detectors 
	TString                     fDetectorsW ;                   //! list of active detectors with QA implemented 
	AliESDEvent *               fESD ;                          //! current ESD
	TTree *                     fESDTree ;                      //! current ESD Tree
  AliEventInfo *              fEventInfo ;                    //! info on the current event  
	TString                     fGAliceFileName ;               //! name of the galice file
	UInt_t                      fFirstEvent ;                   //! first event to process
	Long64_t                    fMaxEvents ;                    //! number of events to process
	TString                     fMode ;                         //! sim or rec
	Long64_t                    fNumberOfEvents ;               //! number of events in the run 
  AliRecoParam                fRecoParam;                     //! container for the reco-param objects for detectors
	UInt_t                      fRunNumber ;                    //! current run number
	AliRawReader *              fRawReader ;                    //! current raw reader object 
	Bool_t                      fRawReaderDelete ;              //! tells if the rawReader has been created by this
	AliRunLoader *              fRunLoader ;                    //! current run loader object
	TString                     fTasks ;                        //! list of QA tasks to be performed
	static const UInt_t         fgkNDetectors = AliQAv1::kNDET ;  //! number of detectors    
	AliLoader      *            fLoader[fgkNDetectors];         //! array of detectors loader
	AliQADataMaker *            fQADataMaker[fgkNDetectors];    //! array of QA data maker objects
	Int_t                       fQACycles[fgkNDetectors];       //! array of QA cycle length
	Bool_t                      fQAWriteExpert[fgkNDetectors];  //! array of QA cycle length
  AliRecoParam::EventSpecie_t fEventSpecie ;                  //! type of event 
  Bool_t                      fPrintImage ;                   //! flag to print the images or not
  Bool_t                      fSaveData ;                     //! flag to sve the QA data or not   
    
  ClassDef(AliQAManager, 2)      // class for running the QA makers
};

#endif
