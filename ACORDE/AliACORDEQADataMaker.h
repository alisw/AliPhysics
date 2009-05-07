#ifndef ALIACORDEQADATAMAKER_H
#define ALIACORDEQADATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliACORDEQADataMaker.h 25659 2008-05-030 15:13:46Z ldg $ */


//  Produces the data needed to calculate the quality assurance. 
//  ACORDE QA for Hits, Digits, RAW and ESD's


//  Authors:
//
//  Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//  Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP)
//  Arturo Fernandez Tellez <afernan@mail.cern.ch (FCFM-BUAP)
//
//  Created: June 13th 2008
//---




#include "AliQADataMaker.h"

class AliACORDEQADataMaker: public AliQADataMaker {

public:
  AliACORDEQADataMaker() ;          // constructor
  AliACORDEQADataMaker(const AliACORDEQADataMaker& qadm) ;   
  AliACORDEQADataMaker& operator = (const AliACORDEQADataMaker& qadm) ;
  virtual ~AliACORDEQADataMaker() {;} // destructor

private:
  virtual void   InitHits() ;      //book hit QA histo 
  virtual void   InitDigits() ;    //book Digit QA histo
  virtual void   InitRaws() ;    //book Digit QA histo
  virtual void   InitRecPoints();  //book cluster QA histo
  virtual void   InitESDs() ;      //book ESD QA histo 
  virtual void   MakeHits(TTree * hits) ;       //Fill hit QA histo
  virtual void   MakeRaws(AliRawReader* rawReader) ;
  virtual void   MakeDigits(TTree* digitsTree) ;   //Fill Digit QA histo
  virtual void   MakeRecPoints(TTree * clusters)    ;  //Fill cluster QA histo
  virtual void   MakeESDs(AliESDEvent * esd) ;         //Fill hit QA histo

/*******/

  virtual Int_t  Add2DigitsList(TH1*, Int_t){return 0;};
  virtual Int_t  Add2HitsList(TH1*, Int_t){return 0;};
  virtual Int_t  Add2RecPointsList(TH1*, Int_t){return 0;};
  virtual Int_t  Add2RawsList(TH1*, Int_t){return 0;};
  virtual Int_t  Add2SDigitsList(TH1*, Int_t){return 0;};
  virtual void   Exec(AliQAv1::TASKINDEX_t, TObject*){};
  virtual void   EndOfCycle(AliQAv1::TASKINDEX_t){};
  virtual Int_t  Add2ESDsList(TH1*, Int_t){return 0;};
  virtual TH1*   GetDigitsData(Int_t){return 0;};
  virtual TH1*   GetESDsData(Int_t){return 0;};
  virtual TH1*   GetHitsData(Int_t){return 0;};
  virtual TH1*   GetRecPointsData(Int_t){return 0;};
  virtual TH1*   GetRawsData(Int_t){return 0;};
  virtual TH1*   GetSDigitsData(Int_t){return 0;};
  virtual TObjArray* Init(AliQAv1::TASKINDEX_t, Int_t, Int_t){return 0;};
  virtual void   Init(AliQAv1::TASKINDEX_t, TObjArray*, Int_t, Int_t){};
  virtual void   StartOfCycle(AliQAv1::TASKINDEX_t, Bool_t){};
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray*){};
  virtual void   InitSDigits(){};
  virtual void   MakeHits(TClonesArray*){};
  virtual void   MakeDigits(TClonesArray*){};
  virtual void   MakeSDigits(TClonesArray*){};
  virtual void   MakeSDigits(TTree*){};
  virtual void   StartOfDetectorCycle() ;
/******/
  ClassDef(AliACORDEQADataMaker,1)  // description 

};

#endif // AliACORDEQADataMaker_H
