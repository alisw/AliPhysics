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
#include <iostream.h>

// --- AliRoot header files ---


class AliPHOSClusterizer ;
class AliPHOSGeometry ;

class  AliPHOSTrackSegmentMaker : public TTask {

public:

  AliPHOSTrackSegmentMaker() ;                     
  AliPHOSTrackSegmentMaker(const char* headerFile, const char* name, const Bool_t toSplit) ;                     
  
  virtual ~ AliPHOSTrackSegmentMaker() ;

  virtual void    Exec(Option_t * option){cout << "Not Defined" << endl ; } 
  //  virtual char*   GetRecPointsBranch ()const{cout << "Not Defined" << endl ; return 0 ; } 
  //  virtual char*   GetTrackSegmentsBranch ()const{cout << "Not Defined" << endl ; return 0 ; } 
  virtual const Int_t GetTrackSegmentsInRun()  const {cout << "Not Defined" << endl ; return 0 ; } 

  virtual void    Print(Option_t * option)const {cout << "Not Defined" << endl ; }  
  //  virtual void Set...   // method to choose recPoints: along z only, along x ...???
  //  virtual void SetChoosingAlgirithm(){cout << "Not Defined" << endl ; return 0 ; } 
  //  virtual void SetMaxEmcCpvDistance(Float_t r) {cout << "Not Defined" << endl ; return 0 ; } 
  //  virtual void SetRecPointsBranch(const char * title){cout << "Not Defined" << endl ; } 
  //  virtual void SetTrackSegmentsBranch(const char * title){cout << "Not Defined" << endl ; } 
  //  virtual void SetSplitFile(const TString splitFileName = "PHOS.RecData.root") const ; 
  virtual const char * Version() const {cout << "Not Defined" << endl ; return 0 ; }   
  virtual void WriteTrackSegments(Int_t event){cout << "Not Defined" << endl ; } 
  
protected:
  
  TFile * fSplitFile ;             //! file in which TrackSegments will eventually be stored
  Bool_t  fToSplit ;               //! Do we work in the split mode

  ClassDef( AliPHOSTrackSegmentMaker,1)    // Algorithm class to make PHOS track segments (Base Class)

};

#endif // ALIPHOSTRACKSEGMENTMAKER_H
