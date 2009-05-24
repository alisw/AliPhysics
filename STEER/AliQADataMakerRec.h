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
class TNtupleD ; 
// --- Standard library ---

// --- AliRoot header files ---
class AliDetectorRecoParam ;
#include "AliQADataMaker.h"
#include "AliQAv1.h"

class AliQADataMakerRec: public AliQADataMaker {
  
public:
	
	AliQADataMakerRec(const char * name="", const char * title="") ;          // ctor
	AliQADataMakerRec(const AliQADataMakerRec& qadm) ;   
	AliQADataMakerRec& operator = (const AliQADataMakerRec& qadm) ;
	virtual ~AliQADataMakerRec() ; // dtor
  
 	virtual Int_t Add2DigitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)    
    { return Add2List(hist, index, fDigitsQAList, expert, image) ; }
	virtual Int_t Add2ESDsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)                      
    { return Add2List(hist, index, fESDsQAList, expert, image) ; }
	virtual Int_t Add2HitsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/, const Bool_t /*image = kFALSE*/)      
    { return -1 ; }  
	virtual Int_t Add2RecPointsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)                 
    { return Add2List(hist, index, fRecPointsQAList, expert, image) ; }
	virtual Int_t Add2RawsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE, const Bool_t saveForCorr = kFALSE)  { 
    return Add2List(hist, index, fRawsQAList, expert, image, saveForCorr) ; }
	virtual Int_t Add2SDigitsList(TH1 * /*hist*/, const Int_t /*index*/, const Bool_t /*expert = kFALSE*/, const Bool_t /*image = kFALSE*/)   { return -1 ; } 
	virtual void        Exec(AliQAv1::TASKINDEX_t task, TObject * data) ;
	virtual void        EndOfCycle() ;
	virtual void        EndOfCycle(AliQAv1::TASKINDEX_t task) ;
	virtual void        EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** ) {AliInfo("To be implemented by detectors");} 
	virtual TH1 *       GetDigitsData(const Int_t index   )  { return dynamic_cast<TH1 *>(GetData(fDigitsQAList, index)) ; } 
	virtual TH1 *       GetESDsData(const Int_t index)       { return dynamic_cast<TH1 *>(GetData(fESDsQAList, index)) ; }
	virtual TH1 *       GetHitsData(const Int_t /*index*/)   { return NULL ; }
  virtual const AliDetectorRecoParam * GetRecoParam() { return fRecoParam ; }

	virtual TH1 *       GetRecPointsData(const Int_t index)  { return dynamic_cast<TH1 *>(GetData(fRecPointsQAList, index)) ; }
	virtual TH1 *       GetRawsData(const Int_t index)       { return dynamic_cast<TH1 *>(GetData(fRawsQAList, index))  ; }
 	virtual TH1 *       GetSDigitsData(const Int_t /*index*/)   { return NULL ; }  
  virtual void        MakeImage(AliQAv1::TASKINDEX_t task) ; 
	virtual TObjArray** Init(AliQAv1::TASKINDEX_t task, Int_t cycles = -1) ;
	virtual void        Init(AliQAv1::TASKINDEX_t task, TObjArray ** list, Int_t run, Int_t cycles = -1) ;
	virtual void        StartOfCycle(Int_t run = -1) ;
	virtual void        StartOfCycle(AliQAv1::TASKINDEX_t task, Int_t run, const Bool_t sameCycle = kFALSE) ;
	virtual void        SetRecoParam(const AliDetectorRecoParam *param) { fRecoParam = param; }

protected: 

	virtual void   InitDigits()                        {AliInfo("To be implemented by detectors");}
	virtual void   InitESDs()                          {AliInfo("To be implemented by detectors");}
  virtual void   InitRecoParams() ; 
	virtual void   InitHits()                          {AliWarning("Call not valid") ; }
	//virtual void   InitRecParticles()                {AliInfo("To be implemented by detectors");}
	virtual void   InitRecPoints()                     {AliInfo("To be implemented by detectors");}
	virtual void   InitRaws()                          {AliInfo("To be implemented by detectors");}
	virtual void   InitSDigits()                       {AliWarning("Call not valid") ; }
	//virtual void   InitTrackSegments()               {AliInfo("To ne implemented by detectors");}
	virtual void   MakeESDs(AliESDEvent * )            {AliInfo("To be implemented by detectors");} 
	virtual void   MakeHits(TClonesArray * )           {AliWarning("Call not valid") ; }
	virtual void   MakeHits(TTree * )                  {AliWarning("Call not valid") ; }  
	virtual void   MakeDigits(TClonesArray * )         {AliInfo("To be implemented by detectors");}   
	virtual void   MakeDigits(TTree * )                {AliInfo("To be implemented by detectors");}   
 	//virtual void   MakeRecParticles(TClonesArray * ) {AliInfo("To be implemented by detectors");} 
	virtual void   MakeRaws(AliRawReader *)            {AliInfo("To be implemented by detectors");} 
	virtual void   MakeRecPoints(TTree * )             {AliInfo("To be implemented by detectors");} 
	virtual void   MakeSDigits(TClonesArray * )        {AliWarning("Call not valid") ; }     
	virtual void   MakeSDigits(TTree * )               {AliWarning("Call not valid") ; }    
	virtual void   StartOfDetectorCycle()              {AliInfo("To be implemented by detectors");} 

	TObjArray * *               fDigitsQAList ;    //! list of the digits QA data objects
	TObjArray * *               fESDsQAList ;      //! list of the ESDs QA data objects
	TObjArray * *               fRawsQAList ;      //! list of the raws QA data objects
	TObjArray * *               fRecPointsQAList ; //! list of the RecPoints QA data objects
  TNtupleD  * *               fCorrNt ;          //! This is used by Corr only to hold its Ntuple. 
  const AliDetectorRecoParam *fRecoParam;        //! const pointer to the reco parameters to be used in the reco QA
  
 ClassDef(AliQADataMakerRec,3)  // description 

};

#endif // ALIQADATAMAKERREC_H
