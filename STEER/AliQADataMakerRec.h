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

class AliQADataMakerRec: public AliQADataMaker {
  
public:
	
	AliQADataMakerRec(const char * name="", const char * title="") ;          // ctor
	AliQADataMakerRec(const AliQADataMakerRec& qadm) ;   
	AliQADataMakerRec& operator = (const AliQADataMakerRec& qadm) ;
	virtual ~AliQADataMakerRec() {;} // dtor
  
 	virtual Int_t Add2DigitsList(TH1 * /*hist*/, const Int_t /*index*/)    { return -1 ; } 
	virtual Int_t Add2ESDsList(TH1 * hist, const Int_t index)  { return Add2List(hist, index, fESDsQAList) ; }
	virtual Int_t Add2HitsList(TH1 * /*hist*/, const Int_t /*index*/)       { return -1 ; }  
	virtual Int_t Add2RecPointsList(TH1 * hist, const Int_t index)  { return Add2List(hist, index, fRecPointsQAList) ; }
	virtual Int_t Add2RawsList(TH1 * hist, const Int_t index)  { return Add2List(hist, index, fRawsQAList) ; }
	virtual Int_t Add2SDigitsList(TH1 * /*hist*/, const Int_t /*index*/)   { return -1 ; } 
	virtual void        Exec(AliQA::TASKINDEX task, TObject * data) ;
	virtual void        EndOfCycle(AliQA::TASKINDEX task) ;
	virtual TH1 *       GetDigitsData(const Int_t /*index*/)    { return NULL ; } 
	virtual TH1 *       GetESDsData(const Int_t index)      { return dynamic_cast<TH1 *>(GetData(fESDsQAList, index)) ; }
	virtual TH1 *       GetHitsData(const Int_t /*index*/)      { return NULL ; }
	virtual TH1 *       GetRecPointsData(const Int_t index) { return dynamic_cast<TH1 *>(GetData(fRecPointsQAList, index)) ; }
	virtual TH1 *       GetRawsData(const Int_t index)     { return dynamic_cast<TH1 *>(GetData(fRawsQAList, index)) ; }
 	virtual TH1 *       GetSDigitsData(const Int_t /*index*/)   { return NULL ; }  
	virtual TObjArray * Init(AliQA::TASKINDEX task, Int_t run, Int_t cycles = -1) ;
	virtual void        Init(AliQA::TASKINDEX task, TObjArray * list, Int_t run, Int_t cycles = -1) ;
	virtual void        StartOfCycle(AliQA::TASKINDEX task, const Bool_t sameCycle = kFALSE) ;

protected: 

	virtual void   EndOfDetectorCycle(AliQA::TASKINDEX, TObjArray * ) {AliInfo("To be implemented by detectors");} 
	virtual void   InitDigits()                        {AliFatal("Call not valid") ; }
	virtual void   InitESDs()                          {AliInfo("To be implemented by detectors");}
	virtual void   InitHits()                          {AliFatal("Call not valid") ; }
	//virtual void   InitRecParticles()                {AliInfo("To be implemented by detectors");}
	virtual void   InitRecPoints()                     {AliInfo("To be implemented by detectors");}
	virtual void   InitRaws()                          {AliInfo("To be implemented by detectors");}
	virtual void   InitSDigits()                       {AliFatal("Call not valid") ; }
	//virtual void   InitTrackSegments()               {AliInfo("To ne implemented by detectors");}
	virtual void   MakeESDs(AliESDEvent * )            {AliInfo("To be implemented by detectors");} 
	virtual void   MakeHits(TClonesArray * )           {AliFatal("Call not valid") ; }
	virtual void   MakeHits(TTree * )                  {AliFatal("Call not valid") ; }  
	virtual void   MakeDigits(TClonesArray * )         {AliFatal("Call not valid") ; }    
	virtual void   MakeDigits(TTree * )                {AliFatal("Call not valid") ; }   
 	//virtual void   MakeRecParticles(TClonesArray * ) {AliInfo("To be implemented by detectors");} 
	virtual void   MakeRaws(AliRawReader *)            {AliInfo("To be implemented by detectors");} 
	virtual void   MakeRecPoints(TTree * )             {AliInfo("To be implemented by detectors");} 
	virtual void   MakeSDigits(TClonesArray * )        {AliFatal("Call not valid") ; }     
	virtual void   MakeSDigits(TTree * )               {AliFatal("Call not valid") ; }    
	virtual void   StartOfDetectorCycle()              {AliInfo("To be implemented by detectors");} 

	TObjArray *    fESDsQAList ;      //! list of the ESDs QA data objects
	TObjArray *    fRawsQAList ;      //! list of the raws QA data objects
	TObjArray *    fRecPointsQAList ; //! list of the RecPoints QA data objects
  
 ClassDef(AliQADataMakerRec,1)  // description 

};

#endif // ALIQADATAMAKERREC_H
