#ifndef ALIQADATAMAKERSTEER_H
#define ALIQADATAMAKERSTEER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running the QA makers                                           //
//                                                                           //
//   AliQADataMakerSteer qas;                                                //
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
#include "AliRecoParam.h"
 
class AliESDEvent ; 
class AliDetectorRecoParam ; 
class AliESDEvent ; 
class AliQADataMaker ;
class AliRawReader ;  
class AliRunLoader ; 
class AliCorrQADataMakerRec ;

class AliQADataMakerSteer: public TNamed {
public:
	AliQADataMakerSteer(char * mode, const char* gAliceFilename = "galice.root", 
						const char * name = "AliQADataMakerSteer", 
						const char * title = "QA makers") ; 
	AliQADataMakerSteer(const AliQADataMakerSteer & qas) ; 
	AliQADataMakerSteer & operator = (const AliQADataMakerSteer & qas) ; 
	virtual ~AliQADataMakerSteer() ; 
	void        EndOfCycle(TObjArray * detArray=0x0) ; 
	UInt_t      GetCurrentEvent() const { return fCurrentEvent ; }
	TObjArray * GetFromOCDB(AliQA::DETECTORINDEX_t det, AliQA::TASKINDEX_t task, const char * year) const ; 
	AliQADataMaker * GetQADataMaker(const Int_t iDet) ; 
	void        Increment() ;
	void        InitQADataMaker(UInt_t run, const AliRecoParam & par, TObjArray * detArray=0x0) ;
	Bool_t      Merge(const Int_t runNumber = -1 ) const ;  
	void        Reset(const Bool_t sameCycle = kFALSE) ;  
	TString     Run(const char * detectors, const AliQA::TASKINDEX_t taskIndex=AliQA::kNULLTASKINDEX, Bool_t const sameCycle = kFALSE, const char * fileName = NULL) ; 
	TString     Run(const char * detectors, AliRawReader * rawReader, Bool_t const sameCycle = kFALSE) ; 
	TString     Run(const char * detectors, const char * filename, Bool_t const sameCycle = kFALSE) ;
	void        RunOneEvent(AliRawReader * rawReader) ; 
	void        RunOneEventInOneDetector(Int_t det, TTree * tree) ; 
	void        RunOneEvent(AliESDEvent *& esd)  ;
	Bool_t      Save2OCDB(const Int_t runNumber, const char * year = "08", const Int_t cycleNumber=0, const char * detectors = "ALL") const ; 
	void        SetActiveDetectors(TString aDet) { fDetectors = aDet ;  }
	void        SetCycleLength(const AliQA::DETECTORINDEX_t det, const Int_t cycle) { fQACycles[det] = cycle ; }
	void        SetEventRange(UInt_t first, UInt_t last) { fFirstEvent = first ; fMaxEvents = last - first + 1 ; }      
	void        SetFirsEvent(UInt_t first) { fFirstEvent = first ; }      
	void        SetMaxEvents(UInt_t max) { fMaxEvents = max ; }      
	void        SetNewCycle() { fCycleSame = kTRUE ; }
	void        SetRecoParam(const char* detector, AliDetectorRecoParam *par) ;
	void        SetRunLoader(AliRunLoader * rl) { fRunLoader = rl ; }
	void        SetTasks(TString tasks) { fTasks = tasks ; }

private: 
	Bool_t			  DoIt(const AliQA::TASKINDEX_t taskIndex) ;
	AliLoader   * GetLoader(Int_t iDet) ; 
	const Int_t   GetQACycles(const Int_t iDet) { return fQACycles[iDet] ; }
	Bool_t			  Init(const AliQA::TASKINDEX_t taskIndex, const  char * fileName = NULL) ;
	Bool_t        InitRunLoader() ; 
	Bool_t        IsSelected(const char * detName)  ;
	Bool_t        Finish(const AliQA::TASKINDEX_t taskIndex) ;
	Bool_t        MergeData(const Int_t runNumber) const ;  
	Bool_t        MergeResults(const Int_t runNumber) const ;  
	Bool_t        SaveIt2OCDB(const Int_t runNumber, TFile * inputFile, const char * year) const ;  

 
	UInt_t                  fCurrentEvent ;                 //! event counter
	Bool_t                  fCycleSame ;                    //! true if 2 consecutive data making for a same detector   
	TString                 fDetectors ;                    //! list of active detectors 
	TString                 fDetectorsW ;                   //! list of active detectors with QA implemented 
	AliESDEvent *           fESD ;                          //! current ESD
	TTree *                 fESDTree ;                      //! current ESD Tree
	TString                 fGAliceFileName ;               //! name of the galice file
	UInt_t                  fFirstEvent ;                   //! first event to process
	Int_t                   fMaxEvents ;                    //! number of events to process
	char *                  fMode ;                         //! sim or rec
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
	
  ClassDef(AliQADataMakerSteer, 0)      // class for running the QA makers
};

#endif
