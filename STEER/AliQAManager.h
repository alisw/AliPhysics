#ifndef ALIQAMANAGER_H
#define ALIQAMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliQAManager.h 30796 2009-01-28 11:05:10Z schutz $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running the QA makers                                           //
//                                                                           //
//   AliQAManager qas;                                                //
//   qas.Run(AliQA::kRAWS, rawROOTFileName);                                 //
//   qas.Run(AliQA::kHITS);                                                  //
//   qas.Run(AliQA::kSDIGITS);                                               //
//   qas.Run(AliQA::kDIGITS);                                                //
//   qas.Run(AliQA::kRECPOINTS);                                             //
//   qas.Run(AliQA::kESDS);                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include "AliQA.h"
#include "AliLoader.h"
#include "AliCDBManager.h"
#include "AliRecoParam.h"
 
class AliESDEvent ; 
class AliDetectorRecoParam ; 
class AliESDEvent ; 
class AliQADataMaker ;
class AliRawReader ;  
class AliRunLoader ; 
class AliCorrQADataMakerRec ;

class AliQAManager : public AliCDBManager {
private:
  AliQAManager() ; 
	AliQAManager(const Char_t * mode, const Char_t * gAliceFilename = "galice.root") ; 
	AliQAManager(const AliQAManager & qas) ; 
	AliQAManager & operator = (const AliQAManager & qas) ; 

public:
	virtual ~AliQAManager() ; 
	void        EndOfCycle(TObjArray * detArray=0x0) ; 
  void        EndOfCycle(TString detectors) ; 
	UInt_t      GetCurrentEvent() const { return fCurrentEvent ; }
	TObjArray * GetFromOCDB(AliQA::DETECTORINDEX_t det, AliQA::TASKINDEX_t task, const char * year) const ; 
  AliQA     * GetQA(UInt_t run, UInt_t evt) ; 
	AliQADataMaker * GetQADataMaker(const Int_t iDet) ; 
	void        Increment() ;
	void        InitQADataMaker(UInt_t run, TObjArray * detArray=0x0) ;
	Bool_t      Merge(Int_t runNumber = -1) const ;  
  void        MergeCustom() const ;
  Bool_t      MergeXML(const char * collection, const char * subFile = 0, const char * outFile = 0) ; 
  static      AliQAManager * QAManager(const Char_t * mode="", TMap *entryCache = NULL, Int_t run = -1) ;
	void        Reset(const Bool_t sameCycle = kFALSE) ;  
	TString     Run(const char * detectors, const AliQA::TASKINDEX_t taskIndex=AliQA::kNULLTASKINDEX, Bool_t const sameCycle = kFALSE, const char * fileName = NULL) ; 
	TString     Run(const char * detectors, AliRawReader * rawReader, Bool_t const sameCycle = kFALSE) ; 
	TString     Run(const char * detectors, const char * filename, Bool_t const sameCycle = kFALSE) ;
	void        RunOneEvent(AliRawReader * rawReader) ; 
	void        RunOneEventInOneDetector(Int_t det, TTree * tree) ; 
	void        RunOneEvent(AliESDEvent *& esd)  ;
	Bool_t      Save2OCDB(const Int_t runNumber, AliRecoParam::EventSpecie_t es, const char * year = "08", const char * detectors = "ALL") const ; 
	void        SetActiveDetectors(TString aDet) { fDetectors = aDet ;  }
	void        SetCycleLength(const AliQA::DETECTORINDEX_t det, const Int_t cycle) { fQACycles[det] = cycle ; }
	void        SetWriteExpert(const AliQA::DETECTORINDEX_t det) { fQAWriteExpert[det] = kTRUE ; }
	void        SetEventRange(UInt_t first, UInt_t last) { fFirstEvent = first ; fMaxEvents = last - first + 1 ; }    
  void        SetEventSpecie(AliRecoParam::EventSpecie_t es) ; 
	void        SetFirsEvent(UInt_t first) { fFirstEvent = first ; }      
	void        SetMaxEvents(UInt_t max) { fMaxEvents = max ; }      
	void        SetNewCycle() { fCycleSame = kTRUE ; }
	void        SetRecoParam(const Int_t det, const AliDetectorRecoParam *par) ;
	void        SetRunLoader(AliRunLoader * rl) { fRunLoader = rl ; }
	void        SetTasks(TString tasks) { fTasks = tasks ; }
  void        SetWriteExpert() ; 
private: 
	Bool_t			  DoIt(const AliQA::TASKINDEX_t taskIndex) ;
	AliLoader   * GetLoader(Int_t iDet) ; 
	Int_t         GetQACycles(const Int_t iDet) const { return fQACycles[iDet] ; }
	Bool_t			  InitQA(const AliQA::TASKINDEX_t taskIndex, const  char * fileName = NULL) ;
	Bool_t        InitRunLoader() ; 
	Bool_t        IsSelected(const char * detName)  ;
	Bool_t        Finish(const AliQA::TASKINDEX_t taskIndex) ;
	Bool_t        MergeData(const Int_t runNumber) const ;  
	Bool_t        MergeResults(const Int_t runNumber) const ;  
	Bool_t        SaveIt2OCDB(const Int_t runNumber, TFile * inputFile, const char * year, AliRecoParam::EventSpecie_t es) const ;  

 	static AliQAManager* fgQAInstance; // AliQAManager instance
	UInt_t                  fCurrentEvent ;                 //! event counter
	Bool_t                  fCycleSame ;                    //! true if 2 consecutive data making for a same detector   
	TString                 fDetectors ;                    //! list of active detectors 
	TString                 fDetectorsW ;                   //! list of active detectors with QA implemented 
	AliESDEvent *           fESD ;                          //! current ESD
	TTree *                 fESDTree ;                      //! current ESD Tree
	TString                 fGAliceFileName ;               //! name of the galice file
	UInt_t                  fFirstEvent ;                   //! first event to process
	Long64_t                fMaxEvents ;                    //! number of events to process
	const Char_t *          fMode ;                         //! sim or rec
	Long64_t                fNumberOfEvents ;               //! number of events in the run 
  AliRecoParam            fRecoParam;                     //! container for the reco-param objects for detectors
	UInt_t                  fRunNumber ;                    //! current run number
	AliRawReader *          fRawReader ;                    //! current raw reader object 
	Bool_t                  fRawReaderDelete ;              //! tells if the rawReader has been created by this
	AliRunLoader *          fRunLoader ;                    //! current run loader object
	TString                 fTasks ;                        //! list of QA tasks to be performed
	static const UInt_t     fgkNDetectors = AliQA::kNDET ;  //! number of detectors    
	AliLoader      *        fLoader[fgkNDetectors];         //! array of detectors loader
	AliQADataMaker *        fQADataMaker[fgkNDetectors];    //! array of QA data maker objects
	Int_t                   fQACycles[fgkNDetectors];       //! array of QA cycle length
	Bool_t                  fQAWriteExpert[fgkNDetectors];  //! array of QA cycle length
	AliRecoParam::EventSpecie_t fEventSpecie ;              //! event specie, see AliRecoParam::EventSpecie_t 
  ClassDef(AliQAManager, 0)      // class for running the QA makers
};

#endif
