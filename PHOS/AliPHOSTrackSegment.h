#ifndef ALIPHOSSUBTRACK_H
#define ALIPHOSSUBTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////
//  Short description                          //
//  Version SUBATECH                           //
//  Author Dmitri Peressounko RRC KI           //
//      comment: contains pairs (triplets) of  //  
//               EMC+PPSD(+PPSD) clusters, and //
//               evaluates particle type,      // 
//               energy, etc                   //
/////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h"
#include "TVector3.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"

const static Int_t kGAMMA         = 0 ; 
const static Int_t kELECTRON      = 1 ;
const static Int_t kNEUTRAL       = 2 ;  
const static Int_t kCHARGEDHADRON = 3 ;  

class AliPHOSTrackSegment : public TObject  {

public:

  AliPHOSTrackSegment() {} ;       // ctor 
  AliPHOSTrackSegment(AliPHOSEmcRecPoint * EmcRecPoint , AliPHOSPpsdRecPoint * PpsdUp, 
		      AliPHOSPpsdRecPoint * PpsdLow  ) ; // ctor
  AliPHOSTrackSegment(const AliPHOSTrackSegment & ts) ;  // ctor                   
  virtual ~AliPHOSTrackSegment() ; // dtor 

  void Copy(TObject & obj) ;  
  virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
  virtual void   Draw(Option_t * option="") ;
  virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);
  Int_t GetPartType() ;                    // Returns 0 - gamma, 1 - e+, e- ;  2 - neutral hadron ; 3 - charged hadron
  Float_t GetEnergy(){ return fEmcRecPoint->GetTotalEnergy() ;}   // Returns energy in EMC
  
  Float_t GetDistanceInPHOSPlane(void) ;   // Computes in PHOS plane the relative position between EMC and PPSD clusters 
  virtual Int_t  GetPHOSMod(void) {return fEmcRecPoint->GetPHOSMod();  }
  TVector3 GetMomentumDirection() ;        // Returns the momentum direction
  void GetPosition( TVector3 & pos ) ;     // Returns positions of hit
  virtual  void  Paint(Option_t * option="");
  void Print() ;
  void SetDispersionCutOff(Float_t Dcut) {fCutOnDispersion = Dcut ; }    
  
  
private:
  
  AliPHOSEmcRecPoint  * fEmcRecPoint ;
  AliPHOSPpsdRecPoint * fPpsdLow ;
  AliPHOSPpsdRecPoint * fPpsdUp ;
  
  Float_t fCutOnDispersion ;   

  ClassDef(AliPHOSTrackSegment,1)  // description , version 1

};

#endif // AliPHOSSUBTRACK_H
