#ifndef ALIT0QADATAMAKER_H
#define ALIT0QADATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  A. Mastroserio



// --- ROOT system ---


class TH1F ; 
class TH2F ;

// --- Standard library ---
#include <TString.h>
// --- AliRoot header files ---

#include "AliQADataMaker.h"

class AliT0QADataMaker: public AliQADataMaker {

public:
  AliT0QADataMaker() ;          // ctor
  AliT0QADataMaker(const AliT0QADataMaker& qadm) ;   
  AliT0QADataMaker& operator = (const AliT0QADataMaker& qadm) ;
  virtual ~AliT0QADataMaker() {;} // dtor

private:
  virtual void   InitHits() ;      //book hit QA histo 
  virtual void   InitDigits() ;    //book Digit QA histo
  virtual void   InitRecPoints();  //book cluster QA histo
  virtual void   InitESDs() ;      //book ESD QA histo 
  virtual void   MakeHits(TObject * hits) ;       //Fill hit QA histo
  virtual void   MakeDigits(TObject * digits) ;   //Fill Digit QA histo
  virtual void   MakeRecPoints(TTree * clusters)    ;  //Fill cluster QA histo
  virtual void   MakeESDs(AliESDEvent * esd) ;         //Fill hit QA histo

  TH1F *fhHitsTime[24];
  TH1F *fhHitsEff;

  TH1F *fhDigCFD[24];
  TH1F *fhDigLEDamp[24];
  TH1F *fhDigQTC[24];
  TH1F *fhDigMean;
  TH1F *fhDigEff;

  TH1F *fhRecCFD[24];
  TH1F *fhRecLEDamp[24];
  TH1F *fhRecQTC[24];
  TH1F *fhRecMean;
  TH1F *fhRecEff;

  TH1F *fhESDMean;
  TH1F *fhESDVertex;

  ClassDef(AliT0QADataMaker,1)  // description 

};

#endif // AliT0QADataMaker_H
