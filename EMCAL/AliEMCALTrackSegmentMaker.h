#ifndef ALIEMCALTRACKSEGMENTMAKER_H
#define ALIEMCALTRACKSEGMENTMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Algorithm Base class to construct EMCAL track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Kurchatov Institute  & SUBATECH)
//             Adapted from PHOS by Y. Schutz (SUBATECH)

// --- ROOT system ---
#include "TTask.h"
class TFile ;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h"
class AliEMCALClusterizer ;
class AliEMCALGeometry ;

class  AliEMCALTrackSegmentMaker : public TTask {

public:

  AliEMCALTrackSegmentMaker() ;                     
  AliEMCALTrackSegmentMaker(const TString alirunFileName, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;                       
  AliEMCALTrackSegmentMaker(const AliEMCALTrackSegmentMaker & tsmaker):TTask(tsmaker) { ; } 
  virtual ~ AliEMCALTrackSegmentMaker() ;

  virtual void    Exec(Option_t *) = 0;
  virtual Int_t GetTrackSegmentsInRun()  const {Warning("GetTrackSegmentsInRun", "Not Defined" ) ; return 0 ; } 

  virtual void    Print(Option_t *) = 0;
  virtual const char * Version() const {Warning("Version", "Not Defined" ) ; return 0 ; }   
  virtual void WriteTrackSegments() = 0;
  
protected:
  TString fEventFolderName ;  // event folder name
  
  ClassDef( AliEMCALTrackSegmentMaker,2)    // Algorithm class to make EMCAL track segments (Base Class)

};

#endif // ALIEMCALTRACKSEGMENTMAKER_H
