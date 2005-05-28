#ifndef ALIPHOSTRACKSEGMENTMAKER_H
#define ALIPHOSTRACKSEGMENTMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

//_________________________________________________________________________
// Algorithm Base class to construct PHOS track segments
// Associates EMC and CPV clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Kurchatov Institute  & SUBATECH)

// --- ROOT system ---
#include "TTask.h"
#include "AliConfig.h"
class TFile ; 

// --- Standard library ---
//#include <iostream>

// --- AliRoot header files ---


class AliPHOSClusterizer ;
class AliPHOSGeometry ;
class AliESD ;

class  AliPHOSTrackSegmentMaker : public TTask {

public:

  AliPHOSTrackSegmentMaker();
  AliPHOSTrackSegmentMaker(const TString alirunFileName, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;                     
  AliPHOSTrackSegmentMaker(const AliPHOSTrackSegmentMaker & tsmaker) : TTask(tsmaker) { ; } 
  virtual ~ AliPHOSTrackSegmentMaker() ;

  virtual Int_t GetTrackSegmentsInRun()  const {Warning("GetTrackSegmentsInRun", "Not Defined" ) ; return 0 ; } 

  virtual void    Print(const Option_t * = "")const {Warning("Print", "Not Defined" ) ; }
  void SetEventRange(Int_t first=0, Int_t last=-1) {fFirstEvent=first; fLastEvent=last; }
  void SetEventFolderName(TString name) { fEventFolderName = name ; }
  void SetESD(AliESD *esd) { fESD = esd; }

  virtual void WriteTrackSegments() = 0;
  
protected:
  TString fEventFolderName ;  // event folder name
  Int_t   fFirstEvent;        // first event to process
  Int_t   fLastEvent;         // last  event to process
  AliESD * fESD;              //! ESD object

  ClassDef( AliPHOSTrackSegmentMaker,4)  // Algorithm class to make PHOS track segments (Base Class)
};

#endif // ALIPHOSTRACKSEGMENTMAKER_H
