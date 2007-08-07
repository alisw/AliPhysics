#ifndef AliPHOSTracker_h
#define AliPHOSTracker_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.4  2007/08/03 13:52:16  kharlov
 * Working skeleton of matching the ESD tracks and ESD clusters (Iouri Belikov)
 *
 */

//-------------------------------------------------------------------------
//                          PHOS tracker.
// Matches ESD tracks with the PHOS and makes the PID.  
// 
//-------------------------------------------------------------------------

#include <AliTracker.h>

class AliRunLoader;   // Bad !

class TClonesArray;
class TTree;

class AliCluster;
class AliESDEvent;
class AliPHOSTrackSegmentMaker ; 
class AliPHOSPID ; 

class AliPHOSTracker : public AliTracker
{
public:
  AliPHOSTracker();
  AliPHOSTracker(AliRunLoader *loader);  // Bad !
  virtual ~AliPHOSTracker();
  
  Int_t LoadClusters(TTree *ct);
  Int_t PropagateBack(AliESDEvent *ev);
  AliCluster *GetCluster(Int_t idx) const;
  void UnloadClusters();

  Int_t Clusters2Tracks(AliESDEvent *) {return 0;}
  Int_t RefitInward(AliESDEvent *)     {return 0;}

  static void                SetDebug()   { fgDebug = kTRUE ; }
  static void                ResetDebug() { fgDebug = kFALSE ; }
  static Bool_t              Debug() { return fgDebug ; }

protected:
  AliPHOSTracker(const AliPHOSTracker & rhs): AliTracker(rhs){}

private:
  Int_t PropagateBackOld(AliESDEvent *ev); //Bad function: uses RunLoader ;(

  AliPHOSTracker &operator=(const AliPHOSTracker &rhs);

  AliRunLoader *fRunLoader;  //! Bad !

  static Bool_t fgDebug ;    //! Verbosity controller

  TClonesArray *fModules[5];
  
  AliPHOSTrackSegmentMaker * fTSM ; //! the track segment maker 
  AliPHOSPID * fPID ;               //! the pid maker 
  ClassDef(AliPHOSTracker,1)
};

#endif
