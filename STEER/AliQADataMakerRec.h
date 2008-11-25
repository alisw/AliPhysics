#ifndef ALIQADATAMAKERREC_H
#define ALIQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Base Class:
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  Y. Schutz CERN July 2007
//

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMaker.h"
#include "AliQA.h"

class AliQADataMakerRec: public AliQADataMaker {
  
public:
	
	AliQADataMakerRec(const char * name="", const char * title="") ;          // ctor
	AliQADataMakerRec(const AliQADataMakerRec& qadm) ;   
	AliQADataMakerRec& operator = (const AliQADataMakerRec& qadm) ;
	virtual ~AliQADataMakerRec() ; // dtor
  
 	virtual Int_t Add2DigitsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/)    { return -1 ; } 
	virtual Int_t Add2ESDsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE)  { return Add2List(hist, index, fESDsQAList, expert) ; }
	virtual Int_t Add2HitsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/)       { return -1 ; }  
	virtual Int_t Add2RecPointsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE)  { return Add2List(hist, index, fRecPointsQAList, expert) ; }
	virtual Int_t Add2RawsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t saveForCorr = kFALSE)  { return Add2List(hist, index, fRawsQAList, expert, saveForCorr) ; }
	virtual Int_t Add2SDigitsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/)   { return -1 ; } 
	virtual void        Exec(AliQA::TASKINDEX_t task, TObject * data) ;
	virtual void        EndOfCycle() ;
	virtual void        EndOfCycle(AliQA::TASKINDEX_t task) ;
	virtual void        EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray * ) {AliInfo("To be implemented by detectors");} 
	virtual TH1 *       GetDigitsData(const Int_t /*index*/)    { return NULL ; } 
	virtual TH1 *       GetESDsData(const Int_t index)      { return dynamic_cast<TH1 *>(GetData(fESDsQAList, index)) ; }
	virtual TH1 *       GetHitsData(const Int_t /*index*/)      { return NULL ; }
	virtual TH1 *       GetRecPointsData(const Int_t index) { return dynamic_cast<TH1 *>(GetData(fRecPointsQAList, index)) ; }
	virtual TH1 *       GetRawsData(const Int_t index)     { return fRawsQAList ? dynamic_cast<TH1 *>(GetData(fRawsQAList, index)) : NULL ; }
 	virtual TH1 *       GetSDigitsData(const Int_t /*index*/)   { return NULL ; }  
	virtual TObjArray * Init(AliQA::TASKINDEX_t task, Int_t cycles = -1) ;
	virtual void        Init(AliQA::TASKINDEX_t task, TObjArray * list, Int_t run, Int_t cycles = -1) ;
	virtual void        StartOfCycle(Int_t run = -1) ;
	virtual void        StartOfCycle(AliQA::TASKINDEX_t task, Int_t run, const Bool_t sameCycle = kFALSE) ;

	virtual void        SetRecoParam(const AliDetectorRecoParam *param) { fRecoParam = param; }

protected: 

	virtual void   InitDigits()                        {AliWarning("Call not valid") ; }
	virtual void   InitESDs()                          {AliInfo("To be implemented by detectors");}
	virtual void   InitHits()                          {AliWarning("Call not valid") ; }
	//virtual void   InitRecParticles()                {AliInfo("To be implemented by detectors");}
	virtual void   InitRecPoints()                     {AliInfo("To be implemented by detectors");}
	virtual void   InitRaws()                          {AliInfo("To be implemented by detectors");}
	virtual void   InitSDigits()                       {AliWarning("Call not valid") ; }
	//virtual void   InitTrackSegments()               {AliInfo("To ne implemented by detectors");}
	virtual void   MakeESDs(AliESDEvent * )            {AliInfo("To be implemented by detectors");} 
	virtual void   MakeHits(TClonesArray * )           {AliWarning("Call not valid") ; }
	virtual void   MakeHits(TTree * )                  {AliWarning("Call not valid") ; }  
	virtual void   MakeDigits(TClonesArray * )         {AliWarning("Call not valid") ; }    
	virtual void   MakeDigits(TTree * )                {AliWarning("Call not valid") ; }   
 	//virtual void   MakeRecParticles(TClonesArray * ) {AliInfo("To be implemented by detectors");} 
	virtual void   MakeRaws(AliRawReader *)            {AliInfo("To be implemented by detectors");} 
	virtual void   MakeRecPoints(TTree * )             {AliInfo("To be implemented by detectors");} 
	virtual void   MakeSDigits(TClonesArray * )        {AliWarning("Call not valid") ; }     
	virtual void   MakeSDigits(TTree * )               {AliWarning("Call not valid") ; }    
	virtual void   StartOfDetectorCycle()              {AliInfo("To be implemented by detectors");} 

	TObjArray *                 fESDsQAList ;      //! list of the ESDs QA data objects
	TObjArray *                 fRawsQAList ;      //! list of the raws QA data objects
	TObjArray *                 fRecPointsQAList ; //! list of the RecPoints QA data objects
  TObject   *                 fObject ;          //! This is used by Corr only to hold its Ntuple. Allows to write something else than TH1 to QA data file
  const AliDetectorRecoParam *fRecoParam;        //! const pointer to the reco parameters to be used in the reco QA
  
 ClassDef(AliQADataMakerRec,2)  // description 

};

#endif // ALIQADATAMAKERREC_H
