#ifndef ALIEMCALTRACKSEGMENT_H
#define ALIEMCALTRACKSEGMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Track segment in EMCAL
//  Can be any combination of : 1 PRERecPoint, ECRecPoint and HCRecPoint     
//                  
//*-- Author:  Dmitri Peressounko (RRC KI & SUBATECH)
//             Adapted from PHOS by Y. Schutz (SUBATECH)
// --- ROOT system ---

#include "TObject.h"
class TClonesArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALTowerRecPoint.h" 
class AliEMCALRecPoint ; 

class AliEMCALTrackSegment : public TObject  {

public:

  AliEMCALTrackSegment() {} 
  AliEMCALTrackSegment(AliEMCALTowerRecPoint * ec, AliEMCALTowerRecPoint * pre, AliEMCALTowerRecPoint * hc) ; 
  AliEMCALTrackSegment(const AliEMCALTrackSegment & ts) ;  // ctor                   
  virtual ~AliEMCALTrackSegment() {  } 

  void Copy(TObject & obj) ;  

  Int_t   GetIndexInList() const {  return fIndexInList ; } 
  Int_t   GetPREIndex()    const {  return fPRERecPoint ; }
  Int_t   GetECAIndex()    const {  return fECARecPoint; }
  Int_t   GetHCAIndex()    const {  return fHCARecPoint; }

  virtual void  Print(Option_t * option) const;
  void SetIndexInList(Int_t val){ fIndexInList = val ;     } 
  void SetPRERecPoint(AliEMCALRecPoint * pre ) ; 
  void SetHCARecPoint(AliEMCALRecPoint * hc ) ; 

  typedef TClonesArray TrackSegmentsList ; 
 
 private:
  
  Int_t fPRERecPoint ; // The PRE reconstructed point index in array stored in TreeR/EMCALPRERP
  Int_t fECARecPoint ; // The EC reconstructed point index in array stored in TreeR/EMCALECRP
  Int_t fHCARecPoint ; // The HC reconstructed point index in array stored in TreeR/EMCALHCRP
  Int_t fIndexInList ; // The index of this TrackSegment in the list stored in TreeR (to be set by analysis)
  
  ClassDef(AliEMCALTrackSegment,1)  // Track segment in EMCAL

};

#endif // ALIEMCALTRACKSEGMENT_H
