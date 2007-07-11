#ifndef ALIPHOSTRACKSEGMENTMAKER_H
#define ALIPHOSTRACKSEGMENTMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.40  2006/08/29 11:41:19  kharlov
 * Missing implementation of ctors and = operator are added
 *
 * Revision 1.39  2006/08/25 16:00:53  kharlov
 * Compliance with Effective C++AliPHOSHit.cxx
 *
 * Revision 1.38  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
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
class AliESDEvent ;

class  AliPHOSTrackSegmentMaker : public TTask {

public:

  AliPHOSTrackSegmentMaker();
  AliPHOSTrackSegmentMaker(const TString alirunFileName, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;                     
  AliPHOSTrackSegmentMaker(const AliPHOSTrackSegmentMaker & tsmaker) ;
  virtual ~ AliPHOSTrackSegmentMaker() ;
  AliPHOSTrackSegmentMaker & operator = (const AliPHOSTrackSegmentMaker & /*rvalue*/)  {
    Fatal("operator =", "not implemented") ; return *this ; }

  virtual Int_t GetTrackSegmentsInRun()  const {Warning("GetTrackSegmentsInRun", "Not Defined" ) ; return 0 ; } 

  virtual void    Print(const Option_t * = "")const {Warning("Print", "Not Defined" ) ; }

  void SetEventRange(Int_t first=0, Int_t last=-1) {fFirstEvent=first; fLastEvent=last; }
  void SetEventFolderName(TString name) { fEventFolderName = name ; }
  void SetESD(AliESDEvent *esd) { fESD = esd; }

  TString GetEventFolderName() const {return fEventFolderName;}
  Int_t   GetFirstEvent()      const {return fFirstEvent;     }
  Int_t   GetLastEvent()       const {return fLastEvent;      }
  AliESDEvent *GetESD()             const {return fESD;            }

  virtual void WriteTrackSegments() = 0;
  
protected:
  TString fEventFolderName ;  // event folder name
  Int_t   fFirstEvent;        // first event to process
  Int_t   fLastEvent;         // last  event to process
  AliESDEvent * fESD;              //! ESD object

  ClassDef( AliPHOSTrackSegmentMaker,5)  // Algorithm class to make PHOS track segments (Base Class)
};

#endif // ALIPHOSTRACKSEGMENTMAKER_H
