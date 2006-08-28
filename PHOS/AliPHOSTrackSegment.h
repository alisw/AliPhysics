#ifndef ALIPHOSTRACKSEGMENT_H
#define ALIPHOSTRACKSEGMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.22  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Track segment in PHOS
//  Can be : 1 EmcRecPoint
//           1 EmcRecPoint + 1 CPV
//                  
//*-- Author:  Dmitri Peressounko (RRC KI & SUBATECH)

// --- ROOT system ---

#include "TObject.h"
class TClonesArray ; 

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSRecPoint ; 
class AliPHOSEmcRecPoint ; 
class AliPHOSCpvRecPoint ; 
class AliESDtrack ; 

class AliPHOSTrackSegment : public TObject  {

public:

  AliPHOSTrackSegment() ;
  AliPHOSTrackSegment(AliPHOSEmcRecPoint * EmcRecPoint , 
		      AliPHOSRecPoint * Cpv) ;
  AliPHOSTrackSegment(AliPHOSEmcRecPoint * EmcRecPoint , 
		      AliPHOSRecPoint * Cpv, Int_t track) ;
  AliPHOSTrackSegment(const AliPHOSTrackSegment & ts) ;  // ctor                   
  virtual ~AliPHOSTrackSegment() {  } 

  void Copy(TObject & obj) const;  

  Int_t   GetIndexInList() const {  return fIndexInList ;   } 
  Int_t   GetEmcIndex()    const {  return fEmcRecPoint ;   }
  Int_t   GetCpvIndex()    const {  return fCpvRecPoint; }
  Int_t   GetTrackIndex()  const {  return fTrack; }

  virtual void  Print(const Option_t * = "") const;
  void    SetIndexInList(Int_t val){ fIndexInList = val ;     } 
  void    SetCpvRecPoint(AliPHOSRecPoint * CpvRecPoint ); //sets CPV Rec Point

  typedef TClonesArray TrackSegmentsList ; 
 
 private:
  
  Int_t fEmcRecPoint ;     // The EMC reconstructed point index in array stored in TreeR/PHOSEmcRP
  Int_t fIndexInList ;     // the index of this TrackSegment in the list stored in TreeR (to be set by analysis)
  Int_t fCpvRecPoint ;     // The CPV reconstructed point in array stored in TreeR/PHOSCpvRP
  Int_t fTrack ;           // The charged track index (from global tracking) in ESD file 

  ClassDef(AliPHOSTrackSegment,1)  // Track segment in PHOS

};

#endif // ALIPHOSTRACKSEGMENT_H
