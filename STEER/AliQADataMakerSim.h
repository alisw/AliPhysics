#ifndef AliQADataMakerSim_H
#define AliQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//  Base Class:
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  Y. Schutz CERN July 2007


// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMaker.h"

class AliQADataMakerSim: public AliQADataMaker {
  
public:
	
	AliQADataMakerSim(const char * name="", const char * title="") ;          // ctor
	AliQADataMakerSim(const AliQADataMakerSim& qadm) ;   
	AliQADataMakerSim& operator = (const AliQADataMakerSim& qadm) ;
	virtual ~AliQADataMakerSim() ; // dtor
  
	virtual Int_t Add2DigitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE)    { return Add2List(hist, index, fDigitsQAList, expert) ; }
	virtual Int_t Add2ESDsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/)      { return -1 ; } 
	virtual Int_t Add2HitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE)      { return Add2List(hist, index, fHitsQAList, expert) ; }
	virtual Int_t Add2RecPointsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/) { return -1 ; } 
	virtual Int_t Add2RawsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/, const Bool_t /*saveForCorr = kFALSE*/)      { return -1 ; }  
	virtual Int_t Add2SDigitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE)   { return Add2List(hist, index, fSDigitsQAList, expert) ; }
	virtual void        Exec(AliQA::TASKINDEX_t task, TObject * data) ;
	virtual void        EndOfCycle() ;
	virtual void        EndOfCycle(AliQA::TASKINDEX_t task) ;
	virtual void        EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray ** ) {AliInfo("To be implemented by detectors");} 
	virtual TH1 *       GetDigitsData(const Int_t index)    { return dynamic_cast<TH1 *>(GetData(fDigitsQAList, index)) ; }
	virtual TH1 *       GetESDsData(const Int_t /*index*/)      { return NULL ; }
	virtual TH1 *       GetHitsData(const Int_t index)      { return dynamic_cast<TH1 *>(GetData(fHitsQAList, index)) ; }
	virtual TH1 *       GetRecPointsData(const Int_t /*index*/) { return NULL ; }
	virtual TH1 *       GetRawsData(const Int_t /*index*/)      { return NULL ; } 
	virtual TH1 *       GetSDigitsData(const Int_t index)   { return dynamic_cast<TH1 *>(GetData(fSDigitsQAList, index)) ; }
	virtual TObjArray** Init(AliQA::TASKINDEX_t task, Int_t cycles = -1) ;
	virtual void        Init(AliQA::TASKINDEX_t task, TObjArray ** list, Int_t run, Int_t cycles = -1) ;
	virtual void        StartOfCycle(Int_t run = -1) ;
	virtual void        StartOfCycle(AliQA::TASKINDEX_t task, Int_t run, const Bool_t sameCycle = kFALSE) ;

protected: 
	
	virtual void   InitDigits()                     {AliInfo("To be implemented by detectors");}
	virtual void   InitESDs()                       {AliWarning("Call not valid") ; } 
	virtual void   InitHits()                       {AliInfo("To be implemented by detectors");}
	virtual void   InitRecPoints()                  {AliWarning("Call not valid") ; } 
	virtual void   InitRaws()                       {AliWarning("Call not valid") ; }
	virtual void   InitSDigits()                    {AliInfo("To be implemented by detectors");}
	virtual void   MakeESDs(AliESDEvent * )         {AliWarning("Call not valid") ; }
	virtual void   MakeHits(TClonesArray * )        {AliInfo("To be implemented by detectors");} 
	virtual void   MakeHits(TTree * )               {AliInfo("To be implemented by detectors");} 
	virtual void   MakeDigits(TClonesArray * )      {AliInfo("To be implemented by detectors");} 
	virtual void   MakeDigits(TTree * )             {AliInfo("To be implemented by detectors");} 
	virtual void   MakeRaws(AliRawReader *)         {AliWarning("Call not valid") ; } 
	virtual void   MakeRecPoints(TTree * )          {AliWarning("Call not valid") ; } 
	virtual void   MakeSDigits(TClonesArray * )     {AliInfo("To be implemented by detectors");} 
	virtual void   MakeSDigits(TTree * )            {AliInfo("To be implemented by detectors");} 
	virtual void   StartOfDetectorCycle()           {AliInfo("To be implemented by detectors");} 

	TObjArray * *    fDigitsQAList ;    //! list of the digits QA data objects
	TObjArray * *    fHitsQAList ;      //! list of the hits QA data objects
	TObjArray * *   fSDigitsQAList ;   //! list of the sdigits QA data objects
  
 ClassDef(AliQADataMakerSim,1)  // description 

};

#endif // AliQADataMakerSim_H

