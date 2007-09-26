#ifndef ALIPHOSQUALASSDATAMAKER_H
#define ALIPHOSQUALASSDATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQualAssDataMaker.h"

class AliPHOSQualAssDataMaker: public AliQualAssDataMaker {

public:
  AliPHOSQualAssDataMaker() ;          // ctor
  AliPHOSQualAssDataMaker(const AliPHOSQualAssDataMaker& qadm) ;   
  AliPHOSQualAssDataMaker& operator = (const AliPHOSQualAssDataMaker& qadm) ;
  virtual ~AliPHOSQualAssDataMaker() {;} // dtor
  
private:
  virtual void   InitHits() ; 
  virtual void   InitESDs() ; 
  virtual void   InitDigits() ; 
  virtual void   InitRecParticles() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitTrackSegments() ; 
  virtual void   InitSDigits() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeHits(TObject * hits) ;
  virtual void   MakeDigits(TObject * digits) ; 
  // virtual void   MakeRecParticles(TTree * recpar) ; 
  virtual void   MakeRecPoints(TTree * recpo) ; 
  virtual void   MakeSDigits(TObject * sigits) ; 
  //virtual void   MakeTrackSegments(TTree *ts ) ; 
  
  TH1F  * fhHits ;            //! hits energy histogram
  TH1I  * fhHitsMul ;         //! hits multiplicity histogram
  TH1I  * fhDigits ;          //! digits energy histogram
  TH1I  * fhDigitsMul ;       //! digits multiplicity histogram
  TH1F  * fhSDigits ;         //! sdigits energy histogram
  TH1I  * fhSDigitsMul ;      //! sdigits multiplicity histogram
  TH1F  * fhEmcRecPoints ;    //! Emc recpoints energy histogram
  TH1I  * fhEmcRecPointsMul ; //! emc recpoints multiplicity histogram
  TH1F  * fhCpvRecPoints ;    //! cpv recpoints energy histogram
  TH1I  * fhCpvRecPointsMul ; //! cpv recpoints multiplicity histogram
  TH1F  * fhTrackSegments ;   //! tracksegments energy histogram
  TH1I  * fhTrackSegmentsMul ;//! tracksegments multiplicity histogram
  TH1F  * fhRecParticles ;    //! recparticles energy histogram
  TH1I  * fhRecParticlesMul ; //! recparticles multiplicity histogram
  TH1F  * fhESDs ;            //! ESDs energy histogram
  TH1I  * fhESDsMul ;         //! ESDs multiplicity histogram

  ClassDef(AliPHOSQualAssDataMaker,1)  // description 

};

#endif // AliPHOSQualAssDataMaker_H
