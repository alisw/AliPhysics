#ifndef ALIPHOSTRACKSEGMENTMAKERV1_H
#define ALIPHOSTRACKSEGMENTMAKERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Implementation version 1 of algorithm class to construct PHOS track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSTrackSegmentMaker.h"

class AliPHOSEmcRecPoint ;
class AliPHOSRecPoint ;


class  AliPHOSTrackSegmentMakerv1 : public AliPHOSTrackSegmentMaker {

public:

  AliPHOSTrackSegmentMakerv1() ;                     
  AliPHOSTrackSegmentMakerv1(const char* headerFile,const char* branchTitle = 0) ;                     
  AliPHOSTrackSegmentMakerv1(const AliPHOSTrackSegmentMakerv1 & tsm) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    abort() ; 
  }
   
  virtual ~ AliPHOSTrackSegmentMakerv1() ; // dtor
  
  virtual char*  GetRecPointsBranch    (void)const{return (char*)fRecPointsBranchTitle.Data() ;}
  virtual char*  GetTrackSegmentsBranch(void)const{return (char*)fTSBranchTitle.Data() ;}

  virtual void   Exec(Option_t * option) ;
          void   FillOneModule() ;       // Finds range in which RecPoints belonging current PHOS module are

          void   MakeLinks() ;           //Evaluates distances(links) between EMC and PPSD
          void   MakePairs() ;           //Finds pairs(triplets) with smallest link
  virtual void   Print(Option_t * option) const ;
  virtual Bool_t ReadRecPoints() ;
  virtual void   SetMaxEmcPpsdDistance(Float_t r){ fR0 = r ;}
  virtual void   SetRecPointsBranch(const char * title) ; 
  virtual void   SetTrackSegmentsBranch(const char * title) ; 
  virtual void   WriteTrackSegments() ;

  AliPHOSTrackSegmentMakerv1 & operator = (const AliPHOSTrackSegmentMakerv1 & )  {
    // assignement operator requested by coding convention
    // but not needed
    abort() ;
    return *this ; 
  }

private:
  Float_t GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * EmcClu , AliPHOSRecPoint * Ppsd , Bool_t & TooFar ) ; // see R0
  void    Init() ;
  void    PrintTrackSegments(Option_t *option) ;

private:  

  TString fHeaderFileName ;          // name of the file which contains gAlice, Tree headers etc.
  TString fRecPointsBranchTitle ; // name of the file, where RecPoints branchs are stored
  TString fTSBranchTitle ;        // name of the file, where TrackSegment branchs is stored
  AliPHOSClusterizer * fClusterizer ; // !  
  Int_t                  fNTrackSegments ; // number of track segments found 
  AliPHOSGeometry      * fGeom ;           //! pointer to PHOS geometry  
  Int_t          fEvent ;            // ! event being precessed
  TObjArray    * fEmcRecPoints ;     // ! List of EMC Rec Points
  TObjArray    * fCpvRecPoints ;     // ! List of CPV/PPSD recPoints
  TClonesArray * fTrackSegments;     // ! list of final track segments


  Bool_t  fIsInitialized ; //

  Float_t fR0 ;        // Maximum distance between a EMC RecPoint and a PPSD RecPoint   

  TClonesArray * fLinkLowArray ;  //!
  TClonesArray * fLinkUpArray  ;  //!


  Int_t fEmcFirst;     //! Index of first EMC RecPoint belonging to currect PHOS module
  Int_t fEmcLast ;     //!
  Int_t fCpvFirst;     //! Cpv upper layer     
  Int_t fCpvLast;      //! 
  Int_t fPpsdFirst;    //! Cpv low layer     
  Int_t fPpsdLast;     //!
  Int_t fModule ;      //! number of module being processed

  ClassDef( AliPHOSTrackSegmentMakerv1,1)  // Implementation version 1 of algorithm class to make PHOS track segments 

};

#endif // AliPHOSTRACKSEGMENTMAKERV1_H
