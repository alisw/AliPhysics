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

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSRecPoint.h"

class AliPHOSTrackSegment : public TObject  {

public:

  AliPHOSTrackSegment() {} 
  AliPHOSTrackSegment(AliPHOSEmcRecPoint * EmcRecPoint , 
		      AliPHOSRecPoint * PpsdUp, 
		      AliPHOSRecPoint * PpsdLow  ) ; // ctor
  AliPHOSTrackSegment(const AliPHOSTrackSegment & ts) ;  // ctor                   
  virtual ~AliPHOSTrackSegment() {  } 

  void Copy(TObject & obj) ;  

  Int_t   GetIndexInList() const {    return fIndexInList ;   } 
  Int_t   GetEmcIndex()const {  return fEmcRecPoint ;   }
  Int_t   GetPpsdIndex()const{  return fPpsdLowRecPoint;}
  Int_t   GetCpvIndex()const {  return fPpsdUpRecPoint; }

  virtual void  Print(Option_t * option) ;
  void    SetIndexInList(Int_t val){ fIndexInList = val ;     } 
  void    SetCpvRecPoint(AliPHOSRecPoint * PpsdUpRecPoint ); //sets PPSD up Rec Point

  typedef TClonesArray TrackSegmentsList ; 
 
 private:
  
  Int_t fEmcRecPoint ;     // The EMC reconstructed point index in array stored in TreeR/PHOSEmcRP
  Int_t fIndexInList ;     // the index of this TrackSegment in the list stored in TreeR (to be set by analysis)
  Int_t fPpsdLowRecPoint ; // The PPSD reconstructed point from the lower layer index in array stored in TreeR/PHOSPpsdRP
  Int_t fPpsdUpRecPoint ;  // The PPSD reconstructed point from the upper layer index in array stored in TreeR/PHOSPpsdRP
  
  ClassDef(AliPHOSTrackSegment,1)  // Track segment in PHOS

};

#endif // ALIPHOSTRACKSEGMENT_H
