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

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALTowerRecPoint.h"

class AliEMCALTrackSegment : public TObject  {

public:

  AliEMCALTrackSegment() {} 
  AliEMCALTrackSegment(AliEMCALTowerRecPoint * ec, AliEMCALTowerRecPoint * pre, AliEMCALTowerRecPoint * hc) ; 
  AliEMCALTrackSegment(const AliEMCALTrackSegment & ts) ;  // ctor                   
  virtual ~AliEMCALTrackSegment() {  } 

  void Copy(TObject & obj) ;  

  Int_t   GetIndexInList() const {  return fIndexInList ; } 
  Int_t   GetPREIndex()    const {  return fPRERecPoint ; }
  Int_t   GetECIndex()     const {  return fECRecPoint; }
  Int_t   GetHCIndex()     const {  return fHCRecPoint; }

  virtual void  Print(Option_t * option) const;
  void SetIndexInList(Int_t val){ fIndexInList = val ;     } 
  void SetPRERecPoint(AliEMCALRecPoint * pre ) ; 
  void SetHCRecPoint(AliEMCALRecPoint * hc ) ; 

  typedef TClonesArray TrackSegmentsList ; 
 
 private:
  
  Int_t fPRERecPoint ; // The PRE reconstructed point index in array stored in TreeR/EMCALPRERP
  Int_t fECRecPoint ;  // The EC reconstructed point index in array stored in TreeR/EMCALECRP
  Int_t fHCRecPoint ;  // The HC reconstructed point index in array stored in TreeR/EMCALHCRP
  Int_t fIndexInList ; // The index of this TrackSegment in the list stored in TreeR (to be set by analysis)
  
  ClassDef(AliEMCALTrackSegment,1)  // Track segment in EMCAL

};

#endif // ALIEMCALTRACKSEGMENT_H
