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
 
class AliQADataMaker ;
class AliRawReader ;  
class AliRunLoader ; 
class AliESDEvent ; 

class AliQADataMakerSteer: public TNamed {
public:
	AliQADataMakerSteer(const char* gAliceFilename = "galice.root", 
						const char * name = "AliQADataMakerSteer", 
						const char * title = "QA makers") ; 
	AliQADataMakerSteer(const AliQADataMakerSteer & qas) ; 
	AliQADataMakerSteer & operator = (const AliQADataMakerSteer & qas) ; 
	virtual ~AliQADataMakerSteer() ; 
	Bool_t Merge() ;  
    void   Reset() ;  
	Bool_t Run(const char * detectors, const AliQA::TASKINDEX taskIndex, const char * fileName = NULL) ; 
	Bool_t Run(const char * detectors, AliRawReader * rawReader) ; 
	void   SetCycleLength(const AliQA::DETECTORINDEX det, const Int_t cycle) { fQACycles[det] = cycle ; }
    void   SetRunLoader(AliRunLoader * rl) { fRunLoader = rl ; }
private: 
	Bool_t			 DoIt(const AliQA::TASKINDEX taskIndex) ;
	AliLoader      * GetLoader(Int_t iDet) ; 
	const Int_t      GetQACycles(const Int_t iDet) { return fQACycles[iDet] ; }
	AliQADataMaker * GetQADataMaker(Int_t iDet) ; 
	Bool_t			 Init(const AliQA::TASKINDEX taskIndex, const  char * fileName = NULL) ;
	Bool_t           InitRunLoader() ; 
	Bool_t           IsSelected(const char * detName)  ;
	Bool_t           Finish(const AliQA::TASKINDEX taskIndex) ;

 
	Bool_t			   fCycleSame ;                    //! true if 2 consecutive data making for a same detector   
    TString            fDetectors ;                    //! list of active detectors 
	AliESDEvent *      fESD ;                          //! current ESD
	TTree *            fESDTree ;                      //! current ESD Tree
	Bool_t             fFirst ;                        //! to search the detector QA data maker only once
	TString            fGAliceFileName ;               //! name of the galice file
	UInt_t             fRunNumber ;                    //! current run number
	Long64_t           fNumberOfEvents ;               //! number of events in the run 
	AliRawReader     * fRawReader ;                    //! current raw reader object 
	Bool_t             fRawReaderDelete ;              //! tells if the rawReader has been created by this
	AliRunLoader *     fRunLoader ;                    //! current run loader object
	static const UInt_t fgkNDetectors = AliQA::kNDET ; //! number of detectors    
	AliLoader      *   fLoader[fgkNDetectors];         //! array of detectors loader
	AliQADataMaker *   fQADataMaker[fgkNDetectors];    //! array of QA data maker objects
	Int_t              fQACycles[fgkNDetectors];       //! array of QA cycle length
	
  ClassDef(AliQADataMakerSteer, 0)      // class for running the QA makers
};

#endif
