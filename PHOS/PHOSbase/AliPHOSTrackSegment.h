#ifndef ALIPHOSTRACKSEGMENT_H
#define ALIPHOSTRACKSEGMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.23  2006/08/28 10:01:56  kharlov
 * Effective C++ warnings fixed (Timur Pocheptsov)
 *
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
  AliPHOSTrackSegment(AliPHOSEmcRecPoint * EmcRecPoint , 
		      AliPHOSRecPoint * Cpv, Int_t track, 
                      Float_t dx, Float_t dz) ;
  AliPHOSTrackSegment(const AliPHOSTrackSegment & ts) ;  // ctor                   
  virtual ~AliPHOSTrackSegment() {  } 

  void    SetCpvDistance(Float_t x,Float_t z){fDcpv[0]=x ; fDcpv[1]=z ; }

  void Copy(TObject & obj) const;  

  Int_t   GetIndexInList() const {  return fIndexInList ;   } 
  Int_t   GetEmcIndex()    const {  return fEmcRecPoint ;   }
  Int_t   GetCpvIndex()    const {  return fCpvRecPoint; }
  Int_t   GetTrackIndex()  const {  return fTrack; }
  Float_t GetCpvDistance(const Option_t* dr="r") const ;

  virtual void  Print(const Option_t * = "") const;
  void    SetIndexInList(Int_t val){ fIndexInList = val ;     } 
  void    SetCpvRecPoint(AliPHOSRecPoint * CpvRecPoint ); //sets CPV Rec Point

  typedef TClonesArray TrackSegmentsList ; 
 
private:
  AliPHOSTrackSegment & operator = (const AliPHOSTrackSegment & /*ts*/);
 private:
  
  Int_t fEmcRecPoint ;     // The EMC reconstructed point index in array stored in TreeR/PHOSEmcRP
  Int_t fIndexInList ;     // the index of this TrackSegment in the list stored in TreeR (to be set by analysis)
  Int_t fCpvRecPoint ;     // The CPV reconstructed point in array stored in TreeR/PHOSCpvRP
  Int_t fTrack ;           // The charged track index (from global tracking) in ESD file 
  Float_t fDcpv[2] ;       // Distance to projection of CPV cluster

 ClassDef(AliPHOSTrackSegment,1)  // Track segment in PHOS

};

#endif // ALIPHOSTRACKSEGMENT_H
