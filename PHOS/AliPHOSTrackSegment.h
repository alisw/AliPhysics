#ifndef ALIPHOSTRACKSEGMENT_H
#define ALIPHOSTRACKSEGMENT_H
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

  AliPHOSTrackSegment() {
    // default ctor
  } 
  AliPHOSTrackSegment(AliPHOSEmcRecPoint * EmcRecPoint , 
		      AliPHOSPpsdRecPoint * PpsdUp, 
		      AliPHOSPpsdRecPoint * PpsdLow  ) ; // ctor
  AliPHOSTrackSegment(const AliPHOSTrackSegment & ts) ;  // ctor                   
  virtual ~AliPHOSTrackSegment() {
    // dtor
  } 

  void Copy(TObject & obj) ;  
  virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
  virtual void   Draw(Option_t * option="") ;
  virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);
  Float_t GetEnergy() ;   // Returns energy in EMC
  
  Float_t GetDistanceInPHOSPlane(void) ;   // Computes in PHOS plane the relative position between EMC and PPSD clusters 
  virtual Int_t  GetPHOSMod(void) ; 
  TVector3 GetMomentumDirection() ;        // Returns the momentum direction
  void    GetPosition( TVector3 & pos ) ;  // Returns positions of hit
  Int_t * GetPrimariesEmc(Int_t & number) ;
  Int_t * GetPrimariesPpsdLow(Int_t & number) ;
  Int_t * GetPrimariesPpsdUp(Int_t & number) ;
  AliPHOSEmcRecPoint *   GetEmcRecPoint() const ;  
  Int_t   GetIndexInList() const { 
    // returns the index of this in the list  
    return fIndexInList ; 
  } 
  Int_t                   GetEmcRecPointIndex() {
    // return the index of the EMC recpoint at the origine of this 
    return fEmcRecPoint ; 
  }
  AliPHOSPpsdRecPoint *  GetPpsdLowRecPoint() const ;
  Int_t                  GetPpsdLowRecPointIndex() {
    // returns the PPSD low recpoint at the origin of this 
    return fPpsdLowRecPoint ; 
  }
  AliPHOSPpsdRecPoint *  GetPpsdUpRecPoint() const ; 
  Int_t                  GetPpsdUpRecPointIndex() {
    // returns the PPSD up recpoint at the origin of this 
    return fPpsdUpRecPoint ; 
  }
  virtual void  Paint(Option_t * option="");
  void    Print() ;
  void    SetIndexInList(Int_t val) { 
    // sets the index value  of this in the list 
    fIndexInList = val ; 
  } 
 
 private:
  
  Int_t fEmcRecPoint ;     // The EMC reconstructed point index in array stored in TreeR/PHOSEmcRP
  Int_t fIndexInList ;     // the index of this TrackSegment in the list stored in TreeR (to be set by analysis)
  Int_t fPpsdLowRecPoint ; // The PPSD reconstructed point from the lower layer index in array stored in TreeR/PHOSPpsdRP
  Int_t fPpsdUpRecPoint ;  // The PPSD reconstructed point from the upper layer index in array stored in TreeR/PHOSPpsdRP

  typedef TClonesArray TrackSegmentsList ; 

  ClassDef(AliPHOSTrackSegment,1)  // Track segment in PHOS

};

#endif // ALIPHOSTRACKSEGMENT_H
