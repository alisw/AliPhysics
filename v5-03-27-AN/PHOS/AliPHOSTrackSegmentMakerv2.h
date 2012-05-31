#ifndef ALIPHOSTRACKSEGMENTMAKERV2_H
#define ALIPHOSTRACKSEGMENTMAKERV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* History of cvs commits:
 *
 * $Log$
 * Revision 1.3  2007/04/25 19:39:42  kharlov
 * Track extracpolation improved
 *
 * Revision 1.2  2007/04/01 19:16:52  kharlov
 * D.P.: Produce EMCTrackSegments using TPC/ITS tracks (no CPV)
 *
 * Revision 1.50  2007/03/06 06:54:48  kharlov
 * DP:Calculation of cluster properties dep. on vertex added
 *
 * Revision 1.49  2007/02/01 13:59:11  hristov
 * Forward declaration
 *
 * Revision 1.48  2006/08/28 10:01:56  kharlov
 * Effective C++ warnings fixed (Timur Pocheptsov)
 *
 * Revision 1.47  2005/11/17 12:35:27  hristov
 * Use references instead of objects. Avoid to create objects when they are not really needed
 *
 * Revision 1.46  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation version 1 of algorithm class to construct PHOS track segments
// Associates EMC and CPV lusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)
#include <vector>

// --- ROOT system ---
#include <TVector3.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSTrackSegmentMaker.h"

class AliPHOSEmcRecPoint ;
class AliPHOSCpvRecPoint ;
class AliESDtrack ;
class TClonesArray;

class  AliPHOSTrackSegmentMakerv2 : public AliPHOSTrackSegmentMaker {

public:

  AliPHOSTrackSegmentMakerv2() ;                     
  AliPHOSTrackSegmentMakerv2(AliPHOSGeometry *geom);
  AliPHOSTrackSegmentMakerv2(const AliPHOSTrackSegmentMakerv2 & tsm);

  virtual ~ AliPHOSTrackSegmentMakerv2() ; // dtor
  
  virtual void   Clusters2TrackSegments(Option_t *option); // Does the job
          void   FillOneModule() ;       // Finds range in which RecPoints belonging current PHOS module are

          void   MakeLinks() ;           //Evaluates distances(links) between EMC and CPV
          void   MakePairs() ;           //Finds pairs(triplets) with smallest link
  virtual void   Print(const Option_t * = "") const ;
  //Switch to "on flyght" mode, without writing to TreeR and file  
  void SetWriting(Bool_t toWrite = kFALSE){fWrite = toWrite;} 
  virtual void   SetMaxTPCDistance(Float_t r){ fRtpc = r ;} //Maximal distance 
                                                               //between EMCrp and extrapolation of TPC track
 virtual const char * Version() const { return "tsm-v2" ; }  

  AliPHOSTrackSegmentMakerv2 & operator = (const AliPHOSTrackSegmentMakerv2 & )  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }

  virtual TClonesArray * GetTrackSegments() const { return fTrackSegments; }

private:

  void  GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * EmcClu , AliESDtrack* track,
                               Float_t &dx, Float_t &dz ) const ; // see R0
  void    Init() ;
  void    InitParameters() ;
  void    PrintTrackSegments(Option_t *option) ;

protected: 
  struct TrackInPHOS_t {
     AliESDtrack* track ;
     Float_t      x ; //Raw X intersection coordinate
     Float_t      z ; //Raw Z intersection coordinate
  } ;

private:  

  Bool_t  fDefaultInit;               //! Says if the task was created by defaut ctor (only parameters are initialized)
  Bool_t  fWrite ;                   // Write Tracks to TreeT  
 
  Int_t fNTrackSegments ; // number of track segments found 

  Float_t fRtpc ;        // Maximum distance between a EMC RecPoint and extrapolation of a TPC track   
  
  TVector3 fVtx ;        //! Vertex in current position

  TClonesArray * fLinkUpArray  ;  //!
  std::vector<TrackInPHOS_t> fTPCtracks[5];  //!
  Int_t fNtpcTracks[5] ;  //!
  Int_t fEmcFirst;     //! Index of first EMC RecPoint belonging to currect PHOS module
  Int_t fEmcLast ;     //!
  Int_t fModule ;      //! number of module being processed

  TClonesArray * fTrackSegments; // Array with found track-segments

  ClassDef( AliPHOSTrackSegmentMakerv2,2)  // Implementation version 1 of algorithm class to make PHOS track segments 

 };

#endif // AliPHOSTRACKSEGMENTMAKERV2_H
