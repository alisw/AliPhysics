#ifndef ALIEMCALTRACKSEGMENT_H
#define ALIEMCALTRACKSEGMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Track segment in EMCAL    
//                  
//*-- Author:  Dmitri Peressounko (RRC KI & SUBATECH)
//             Adapted from PHOS by Y. Schutz (SUBATECH)
// --- ROOT system ---

#include "TObject.h"
class TClonesArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALRecPoint.h" 
class AliEMCALRecPoint ; 

class AliEMCALTrackSegment : public TObject  {

public:

  AliEMCALTrackSegment() {} 
  AliEMCALTrackSegment(AliEMCALRecPoint * ec) ; 
  AliEMCALTrackSegment(const AliEMCALTrackSegment & ts) ;  // ctor                   
  virtual ~AliEMCALTrackSegment() {  } 

  void Copy(TObject & obj) ;  

  Int_t   GetIndexInList() const {  return fIndexInList ; } 
  Int_t   GetECAIndex()    const {  return fECARecPoint; }

  virtual void  Print(Option_t * option) const;
  void SetIndexInList(Int_t val){ fIndexInList = val ;     } 

  typedef TClonesArray TrackSegmentsList ; 
 
 private:
  Int_t fECARecPoint ; // The EC reconstructed point index in array stored in TreeR/EMCALECRP
  Int_t fIndexInList ; // The index of this TrackSegment in the list stored in TreeR (to be set by analysis)
  
  ClassDef(AliEMCALTrackSegment,1)  // Track segment in EMCAL

};

#endif // ALIEMCALTRACKSEGMENT_H
