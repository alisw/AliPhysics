#ifndef ALIPHOSSUBTRACK_H
#define ALIPHOSSUBTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Track segment in PHOS
//  Can be : 1 EmcRecPoint
//           1 EmcRecPoint + 1 PPSD
//           1 EmcRecPoint + 1 PPSD + 1 PPSD     
//                  
//*-- Author:  Dmitri Peressounko (RRC KI & SUBATECH)

// --- ROOT system ---

#include "TObject.h"
#include "TVector3.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"

class AliPHOSTrackSegment : public TObject  {

public:

  AliPHOSTrackSegment() {}       // ctor 
  AliPHOSTrackSegment(AliPHOSEmcRecPoint * EmcRecPoint , AliPHOSPpsdRecPoint * PpsdUp, 
		      AliPHOSPpsdRecPoint * PpsdLow  ) ; // ctor
  AliPHOSTrackSegment(const AliPHOSTrackSegment & ts) ;  // ctor                   
  virtual ~AliPHOSTrackSegment() {} // dtor 

  void Copy(TObject & obj) ;  
  virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
  virtual void   Draw(Option_t * option="") ;
  virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);
  Float_t GetEnergy(){ return fEmcRecPoint->GetTotalEnergy() ;}   // Returns energy in EMC
  
  Float_t GetDistanceInPHOSPlane(void) ;   // Computes in PHOS plane the relative position between EMC and PPSD clusters 
  virtual Int_t  GetPHOSMod(void) {return fEmcRecPoint->GetPHOSMod();  }
  TVector3 GetMomentumDirection() ;        // Returns the momentum direction
  void GetPosition( TVector3 & pos ) ;     // Returns positions of hit
  Int_t * GetPrimariesEmc(Int_t & number) ;
  Int_t * GetPrimariesPpsdLow(Int_t & number) ;
  Int_t * GetPrimariesPpsdUp(Int_t & number) ;
  AliPHOSEmcRecPoint * GetEmcRecPoint() const { return fEmcRecPoint ; } 
  AliPHOSPpsdRecPoint * GetPpsdLow() const { return fPpsdLow ; } 
  AliPHOSPpsdRecPoint * GetPpsdUp() const { return fPpsdUp ; } 
  virtual  void  Paint(Option_t * option="");
  void Print() ;
  
  
private:
  
  AliPHOSEmcRecPoint  * fEmcRecPoint ; //! The EMC reconstructed point
  AliPHOSPpsdRecPoint * fPpsdLow ;     //! The PPSD reconstructed point from the lower layer
  AliPHOSPpsdRecPoint * fPpsdUp ;      //! The PPSD reconstructed point from the upper layer
  Int_t fEmcRecPointId ; // The EMC reconstructed point Id in the list
  Int_t fPpsdLowId ;     // The PPSD reconstructed point from the lower layer Id in the list
  Int_t fPpsdUpId ;      // The PPSD reconstructed point from the upper layer Id in the list

  ClassDef(AliPHOSTrackSegment,1)  // Track segment in PHOS

};

#endif // AliPHOSSUBTRACK_H
