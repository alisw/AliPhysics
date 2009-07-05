#ifndef ALIQADATAMAKER_H
#define ALIQADATAMAKER_H
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
#include <TH1.h>
#include <TObjArray.h>
#include <TNamed.h>  
//class TCanvas ; 
class TClonesArray;
class TDirectory;
class TFile;  
class TObject; 
class TTree; 
class AliESDEvent;
class AliRawReader;
class AliDetectorRecoParam;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQAv1.h"
#include "AliRecoParam.h" 

class AliQADataMaker: public TNamed {
  
public:
	
	AliQADataMaker(const Char_t * name="", const Char_t * title="") ;          // ctor
	AliQADataMaker(const AliQADataMaker& qadm) ;   
	virtual ~AliQADataMaker() ; // dtor
  
	virtual Int_t Add2DigitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)          = 0 ; 
	virtual Int_t Add2ESDsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)            = 0 ; 
	virtual Int_t Add2HitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)            = 0 ; 
	virtual Int_t Add2RecPointsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)       = 0 ; 
	virtual Int_t Add2RawsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE, const Bool_t saveForCorr = kFALSE)            = 0 ; 
	virtual Int_t Add2SDigitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE)         = 0 ; 
	virtual void        Exec(AliQAv1::TASKINDEX_t, TObject * data)             = 0 ;
	virtual void        EndOfCycle()                                           = 0 ;
	virtual void        EndOfCycle(AliQAv1::TASKINDEX_t )                      = 0 ;
	void                Finish() const ; 
	virtual TH1 *       GetDigitsData(const Int_t index)                       = 0 ; 
	virtual TH1 *       GetESDsData(const Int_t index)                         = 0 ; 
	virtual TH1 *       GetHitsData(const Int_t index)                         = 0 ; 
	virtual TH1 *       GetRecPointsData(const Int_t index)                    = 0 ; 
	virtual TH1 *       GetRawsData(const Int_t index)                         = 0 ; 
	virtual TH1 *       GetSDigitsData(const Int_t index)                      = 0 ; 
	const Char_t *      GetDetectorDirName() const { return fDetectorDirName.Data() ; }
  TList *             GetParameterList() const { return fParameterList[AliRecoParam::AConvert(fEventSpecie)] ; }
  virtual const AliDetectorRecoParam * GetRecoParam() { return NULL ; }
	Int_t               Increment() { return ++fCycleCounter ; } 
	virtual TObjArray** Init(AliQAv1::TASKINDEX_t, Int_t cycles = -1)                                 = 0 ;
  TObjArray*          Init(AliQAv1::TASKINDEX_t, AliRecoParam::EventSpecie_t es, Int_t cycles = -1) ;
	virtual void        Init(AliQAv1::TASKINDEX_t, TObjArray ** list, Int_t run, Int_t cycles = -1)   = 0 ;
	virtual void        InitRaws()          = 0 ; 
  virtual void        InitRecPoints()     = 0 ; 
  Bool_t              IsCycleDone() const { return fCycleCounter > fCycle ? kTRUE : kFALSE ; }
  Bool_t              IsValidEventSpecie(Int_t eventSpecieIndex, TObjArray ** list) ; 
	void                Reset() { fCycleCounter = 0 ; }
	void                SetCycle(Int_t nevts) { fCycle = nevts ; } 
  void                SetWriteExpert() { fWriteExpert = kTRUE ; }
	virtual void        StartOfCycle(Int_t run = -1)                                                   = 0 ;
	virtual void        StartOfCycle(AliQAv1::TASKINDEX_t, Int_t run, const Bool_t sameCycle = kFALSE) = 0 ;
  void                UnSetWriteExpert() { fWriteExpert = kFALSE ; }
  Bool_t              WriteExpert() { return fWriteExpert ; }
  void                SetEventSpecie(AliRecoParam::EventSpecie_t es) { fEventSpecie = es ; }
  void                SetEventSpecie(Int_t es) { fEventSpecie = AliRecoParam::Convert(es) ; }
  virtual void        SetRecoParam(const AliDetectorRecoParam *) {;}

  virtual void        InitRecPointsForTracker() {;} // needed by AliGlobalQADataMaker

protected: 

	Int_t          Add2List(TH1 * hist, const Int_t index, TObjArray ** list, const Bool_t expert = kFALSE, const Bool_t image = kFALSE, const Bool_t saveForCorr = kFALSE) ;
  TH1 *          CloneMe(TH1 * hist, Int_t specie) const ; 
	virtual void   DefaultEndOfDetectorCycle(AliQAv1::TASKINDEX_t task ) ; 
	virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list ) = 0 ; 
	TObject *      GetData(TObjArray ** list, const Int_t index) ;
	virtual void   InitDigits()        = 0 ; 
	virtual void   InitESDs()          = 0 ; 
	virtual void   InitHits()          = 0 ; 
  //virtual void   InitRecParticles()  = 0 ; 
	virtual void   InitSDigits()       = 0 ; 
  //virtual void   InitTrackSegments()  = 0 ; 
	virtual void   MakeESDs(AliESDEvent * )          = 0 ; 
	virtual void   MakeHits()         = 0 ; 
	virtual void   MakeHits(TTree * )                = 0 ;  
	virtual void   MakeDigits()       = 0 ;  
	virtual void   MakeDigits(TTree * )              = 0 ; 
  //virtual void   MakeRecParticles( ) = 0 ; 
	virtual void   MakeRaws(AliRawReader *)          = 0 ; 
	virtual void   MakeRecPoints(TTree * )           = 0 ; 
	virtual void   MakeSDigits()      = 0 ;  
	virtual void   MakeSDigits(TTree * )             = 0 ;  
  //virtual void   MakeTrackSegments(TTree * )		 = 0 ;  
	void           ResetCycle() { fCurrentCycle++ ; fCycleCounter = 0 ; } 
	virtual void   StartOfDetectorCycle()            = 0 ;
	
	TFile *        fOutput ;          //! output root file
	TDirectory *   fDetectorDir ;     //! directory for the given detector in the file
	TString        fDetectorDirName ; //! detector directory name in the quality assurance data file
	Int_t          fCurrentCycle ;    //! current cycle number
	Int_t          fCycle ;           //! length (# events) of the QA data acquisition cycle  
	Int_t          fCycleCounter ;    //! cycle counter
  Bool_t         fWriteExpert ;     //! flag to write or not the expert QA data
  TList **       fParameterList ;   //! list of QA data parameters
	Int_t          fRun ;             //! run number
  AliRecoParam::EventSpecie_t fEventSpecie ; //! event specie, see AliRecoParam
  TClonesArray * fDigitsArray ;    //! array to hold the sdigits
  
private:
	AliQADataMaker& operator = (const AliQADataMaker& /*qadm*/); // Not implemented

  
 ClassDef(AliQADataMaker,3)  // description 

};

#endif // AliQADataMaker_H
