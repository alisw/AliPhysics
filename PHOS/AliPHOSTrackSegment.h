#ifndef ALIPHOSSUBTRACK_H
#define ALIPHOSSUBTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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



class AliPHOSTrackSegment : public TObject  {

public:

  AliPHOSTrackSegment() {} ; // ctor 
  AliPHOSTrackSegment(AliPHOSEmcRecPoint * EmcRecPoint , AliPHOSPpsdRecPoint * PpsdUp, 
                  AliPHOSPpsdRecPoint * PpsdLow  ) ;                    
  virtual ~AliPHOSTrackSegment() ; // dtor 

  virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
  virtual void   Draw(Option_t * option="") ;
  virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);
  Int_t GetPartType() ;          // Returns 0 - gamma, 1 - e+, e- ;  2 - neutral hadron ; 3 - charged hadron
  Float_t GetEnergy(){ return fEmcRecPoint->GetTotalEnergy() ;}   // Returs energy in EMC
  
  Float_t GetDistanceInPHOSPlane(void) ;    // computes in PHOS plane the relative position between EMC and PPSD clusters 
  virtual Int_t  GetPHOSMod(void) {return fEmcRecPoint->GetPHOSMod();  }
  Bool_t GetMomentumDirection( TVector3 & dir ) ;   // True if determined
  void GetPosition( TVector3 & pos ) ;              // Returns positions of hits
  virtual  void  Paint(Option_t * option="");
  void Print() ;
  void SetDispersionCutOff(Float_t Dcut) {fCutOnDispersion = Dcut ; }    
  
  
private:
  
  AliPHOSEmcRecPoint  * fEmcRecPoint ;
  AliPHOSPpsdRecPoint * fPpsdLow ;
  AliPHOSPpsdRecPoint * fPpsdUp ;
  
  Float_t fCutOnDispersion ;   
  
public:

  ClassDef(AliPHOSTrackSegment,1)  // description , version 1

};

#endif // AliPHOSSUBTRACK_H
