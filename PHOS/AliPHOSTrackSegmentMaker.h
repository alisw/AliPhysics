#ifndef ALIPHOSTRACKSEGMENTMAKER_H
#define ALIPHOSTRACKSEGMENTMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Algorithm Base class to construct PHOS track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Kurchatov Institute  & SUBATECH)

// --- ROOT system ---
#include "TTask.h"


// --- Standard library ---

// --- AliRoot header files ---


class AliPHOSClusterizer ;
class AliPHOSGeometry ;

class  AliPHOSTrackSegmentMaker : public TTask {

public:

  AliPHOSTrackSegmentMaker() ;                     
  AliPHOSTrackSegmentMaker(char* headerFile, char* branchTitle = 0) ;                     
  
  virtual ~ AliPHOSTrackSegmentMaker(){
    // dtor 
  } 

  virtual void    Exec(Option_t * option) = 0 ;
  virtual char*   GetRecPointsBranch ()const = 0 ;
  virtual char*   GetTrackSegmentsBranch ()const = 0 ;

  virtual void    Print(Option_t * option)const = 0;
  //  virtual void Set...   // method to choose recPoints: along z only, along x ...???
  //  virtual void SetChoosingAlgirithm() = 0 ;
  //  virtual void SetMaxEmcCpvDistance(Float_t r) = 0 ; 
  virtual Bool_t ReadRecPoints() = 0 ; 
  virtual void SetRecPointsBranch(const char * title) = 0 ;
  virtual void SetTrackSegmentsBranch(const char * title) = 0 ;
  virtual void WriteTrackSegments() = 0 ;
  
 protected:

  ClassDef( AliPHOSTrackSegmentMaker,1)    // Algorithm class to make PHOS track segments (Base Class)

};

#endif // ALIPHOSTRACKSEGMENTMAKER_H
