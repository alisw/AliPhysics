#ifndef ALIQUALASSDATAMAKER_H
#define ALIQUALASSDATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/*
  Base Class:
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
#include <TH1.h>
#include <TList.h>
#include <TNamed.h>  
class TFile;  
class TDirectory;
class TObject; 
class TTree; 
class AliESDEvent;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQualAss.h"

class AliQualAssDataMaker: public TNamed {
  
public:
  
  AliQualAssDataMaker(const char * name="", const char * title="") ;          // ctor
  AliQualAssDataMaker(const AliQualAssDataMaker& qadm) ;   
  AliQualAssDataMaker& operator = (const AliQualAssDataMaker& qadm) ;
  virtual ~AliQualAssDataMaker() {;} // dtor
  
  virtual void        Exec(AliQualAss::TASKINDEX, TObject * data) ;
  void                Finish(AliQualAss::TASKINDEX task) const ; 
  static const char * GetDetectorDirName() { return fDetectorDirName.Data() ; }
  TList *             Init(AliQualAss::TASKINDEX) ;
  Int_t               Add2DigitsList(TH1 * hist, Int_t index)    { return Add2List(hist, index, fDigitsQAList) ; }
  Int_t               Add2ESDsList(TH1 * hist, Int_t index)      { return Add2List(hist, index, fESDsQAList) ; }
  Int_t               Add2HitsList(TH1 * hist, Int_t index)      { return Add2List(hist, index, fHitsQAList) ; }
  Int_t               Add2RecPointsList(TH1 * hist, Int_t index) { return Add2List(hist, index, fRecPointsQAList) ; }
  Int_t               Add2RawsList(TH1 * hist, Int_t index)      { return Add2List(hist, index, fRawsQAList) ; }
  Int_t               Add2SDigitsList(TH1 * hist, Int_t index)   { return Add2List(hist, index, fSDigitsQAList) ; }
  TH1 *               GetDigitsData(Int_t index)    { return dynamic_cast<TH1 *>(GetData(fDigitsQAList, index)) ; }
  TH1 *               GetESDsData(Int_t index)      { return dynamic_cast<TH1 *>(GetData(fESDsQAList, index)) ; }
  TH1 *               GetHitsData(Int_t index)      { return dynamic_cast<TH1 *>(GetData(fHitsQAList, index)) ; }
  TH1 *               GetRecPointsData(Int_t index) { return dynamic_cast<TH1 *>(GetData(fRecPointsQAList, index)) ; }
  TH1 *               GetRawsData(Int_t index)      { return dynamic_cast<TH1 *>(GetData(fRawsQAList, index)) ; }
  TH1 *               GetSDigitsData(Int_t index)   { return dynamic_cast<TH1 *>(GetData(fSDigitsQAList, index)) ; }

protected: 

  Int_t          Add2List(TH1 * hist, Int_t index, TList * list) { list->AddAt(hist, index) ; return list->LastIndex() ; }
  TObject *      GetData(TList * list, Int_t index)  { return list->At(index) ; } 
  virtual void   InitDigits()        {AliInfo("To ne implemented by detectors");}
  virtual void   InitESDs()          {AliInfo("To ne implemented by detectors");}
  virtual void   InitHits()          {AliInfo("To ne implemented by detectors");}
  //virtual void   InitRecParticles()  {AliInfo("To ne implemented by detectors");}
  virtual void   InitRecPoints()     {AliInfo("To ne implemented by detectors");}
  virtual void   InitRaws()          {AliInfo("To ne implemented by detectors");}
  virtual void   InitSDigits()       {AliInfo("To ne implemented by detectors");}
  //virtual void   InitTrackSegments() {AliInfo("To ne implemented by detectors");}
  virtual void   MakeESDs(AliESDEvent * )          {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeHits(TObject * )              {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeDigits(TObject * )            {AliInfo("To ne implemented by detectors");} 
  //  virtual void   MakeRecParticles(TClonesArray * ) {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeRaws(TObject * )                {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeRecPoints(TTree * )           {AliInfo("To ne implemented by detectors");} 
  virtual void   MakeSDigits(TObject * )           {AliInfo("To ne implemented by detectors");} 
  //virtual void   MakeTrackSegments(TTree * )       {AliInfo("To ne implemented by detectors");} 

  TFile *        fOutput ;          //! output root file
  TDirectory *   fDetectorDir ;     //! directory for the given detector in the file
  static TString fDetectorDirName ; //! detector directory name in the quality assurance data file
  TList *        fDigitsQAList ;    //! list of the digits QA data objects
  TList *        fESDsQAList ;      //! list of the ESDs QA data objects
  TList *        fHitsQAList ;      //! list of the hits QA data objects
  TList *        fRawsQAList ;      //! list of the raws QA data objects
  TList *        fRecPointsQAList ; //! list of the recpoints QA data objects
  TList *        fSDigitsQAList ;   //! list of the sdigits QA data objects
 ClassDef(AliQualAssDataMaker,1)  // description 

};

#endif // AliQualAssDataMaker_H
