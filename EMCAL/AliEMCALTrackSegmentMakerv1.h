#ifndef ALIEMCALTRACKSEGMENTMAKERV1_H
#define ALIEMCALTRACKSEGMENTMAKERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Implementation version 1 of algorithm class to construct EMCAL track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH) & Yves Schutz (SUBATECH) 

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALTrackSegmentMaker.h"

class AliEMCALRecPoint ;

class  AliEMCALTrackSegmentMakerv1 : public AliEMCALTrackSegmentMaker {

public:

  AliEMCALTrackSegmentMakerv1() ;                     
  AliEMCALTrackSegmentMakerv1(const TString alirunFileNameFile, const TString eventFolderName = AliConfig::GetDefaultEventFolderName());                  
  AliEMCALTrackSegmentMakerv1(const AliEMCALTrackSegmentMakerv1 & tsm):AliEMCALTrackSegmentMaker(tsm) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
}
   
  virtual ~ AliEMCALTrackSegmentMakerv1() ; // dtor
  
  virtual Int_t GetTrackSegmentsInRun()const {return fTrackSegmentsInRun ;}  

  virtual void   Exec(Option_t * option) ;
  Float_t HowClose(AliEMCALRecPoint * ec, AliEMCALRecPoint * rp, Bool_t &toofar) const ;
          void   MakeLinks() const;      //Evaluates distances(links) between recpoints
          void   MakePairs() ;           //Finds pairs(triplets) with smallest link
  virtual void   Print(Option_t * option) const ;
  virtual const char * Version() const { return "tsm-v1" ; }  

  AliEMCALTrackSegmentMakerv1 & operator = (const AliEMCALTrackSegmentMakerv1 & )  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }


private:

  const TString BranchName() const ; 
  void    Init() ;
  void    InitParameters() ;
  void    PrintTrackSegments(Option_t *option) ;
  void    Unload() ;
  virtual void   WriteTrackSegments() ;

private:  

  Float_t fClose ;               // Spread within which 2 recpoints are declared to have the same direction 
  Bool_t  fDefaultInit ;         //! Says if the task was created by defaut ctor (only parameters are initialized)
  Int_t fNTrackSegments ;        // number of track segments found 
  Int_t fTrackSegmentsInRun ;    //! Total number of track segments in one run

  ClassDef( AliEMCALTrackSegmentMakerv1,4)  // Implementation version 1 of algorithm class to make EMCAL track segments 

};

#endif // AliEMCALTRACKSEGMENTMAKERV1_H
