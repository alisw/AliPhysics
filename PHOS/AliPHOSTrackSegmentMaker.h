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
class TFile ;

// --- Standard library ---

// --- AliRoot header files ---


class AliPHOSClusterizer ;
class AliPHOSGeometry ;

class  AliPHOSTrackSegmentMaker : public TTask {

public:

  AliPHOSTrackSegmentMaker() ;                     
  AliPHOSTrackSegmentMaker(const char* headerFile, const char* name, const Bool_t toSplit) ;                     
  
  virtual ~ AliPHOSTrackSegmentMaker() ;

  virtual void    Exec(Option_t * option){Warning("Exec", "Not Defined" ) ; } 
  //  virtual char*   GetRecPointsBranch ()const{Warning("GetRecPointsBranch", "Not Defined" ) ; return 0 ; } 
  //  virtual char*   GetTrackSegmentsBranch ()const{Warning(" GetTrackSegmentsBranch", "Not Defined" ) ; return 0 ; } 
  virtual const Int_t GetTrackSegmentsInRun()  const {Warning("GetTrackSegmentsInRun", "Not Defined" ) ; return 0 ; } 

  virtual void    Print(Option_t * option)const {Warning("Print", "Not Defined" ) ; }  
  //  virtual void Set...   // method to choose recPoints: along z only, along x ...???
  //  virtual void SetChoosingAlgirithm(){Warning("SetChoosingAlgirithm", "Not Defined" ) ; return 0 ; } 
  //  virtual void SetMaxEmcCpvDistance(Float_t r) {Warning("SetMaxEmcCpvDistance", "Not Defined" ) ; return 0 ; } 
  //  virtual void SetRecPointsBranch(const char * title){Warning("SetRecPointsBranch", "Not Defined" ) ; } 
  //  virtual void SetTrackSegmentsBranch(const char * title){Warning("SetTrackSegmentsBranch", "Not Defined" ) ; } 
  //  virtual void SetSplitFile(const TString splitFileName = "PHOS.RecData.root") const ; 
  virtual const char * Version() const {Warning("Version", "Not Defined" ) ; return 0 ; }   
  virtual void WriteTrackSegments(Int_t event){Warning("WriteTrackSegments", "Not Defined" ) ; } 
  
protected:
  
  TFile * fSplitFile ;             //! file in which TrackSegments will eventually be stored
  Bool_t  fToSplit ;               //! Do we work in the split mode

  ClassDef( AliPHOSTrackSegmentMaker,1)    // Algorithm class to make PHOS track segments (Base Class)

};

#endif // ALIPHOSTRACKSEGMENTMAKER_H
