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
class AliRawReaderRoot ;  
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
    void   Reset() ;  
	Bool_t Run(const AliQA::TASKINDEX taskIndex, const char * fileName = NULL) ; 
	void   SetCycleLength(const AliQA::DETECTORINDEX det, const Int_t cycle) { fQACycles[det] = cycle ; }

private: 
	AliLoader      * GetLoader(Int_t iDet) ; 
	const Int_t      GetQACycles(const Int_t iDet) { return fQACycles[iDet] ; }
	AliQADataMaker * GetQADataMaker(Int_t iDet) ; 
	Bool_t			 Init(const AliQA::TASKINDEX taskIndex, const  char * fileName = NULL) ;
	Bool_t           InitRunLoader() ; 
	Bool_t           Finish(const AliQA::TASKINDEX taskIndex) ;

 
	Bool_t			   fCycleSame ;                    //! true if 2 consecutive data making for a same detector   
	AliESDEvent *      fESD ;                          //! current ESD
	TTree *            fESDTree ;                      //! current ESD Tree
	Bool_t             fFirst ;                        //! to search the detector QA data maker only once
	TString            fGAliceFileName ;               //! name of the galice file
	UInt_t             fRunNumber ;                    //! current run number
	Long64_t           fNumberOfEvents ;               //! number of events in the run 
	AliRawReaderRoot * fRawReader ;                    //! current raw reader object 
	AliRunLoader *     fRunLoader ;                    //! current run loader object
	static const UInt_t fgkNDetectors = AliQA::kNDET ; //! number of detectors    
	AliLoader      *   fLoader[fgkNDetectors];         //! array of detectors loader
	AliQADataMaker *   fQADataMaker[fgkNDetectors];    //! array of QA data maker objects
	Int_t              fQACycles[fgkNDetectors];       //! array of QA cycle length
	
  ClassDef(AliQADataMakerSteer, 0)      // class for running the QA makers
};

#endif
