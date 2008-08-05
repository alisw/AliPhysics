#ifndef ALITRDTRACKINGEFFICIENCY_H
#define ALITRDTRACKINGEFFICIENCY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingEfficiency.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"

class TObjArray;
class TList;
class TClonesArray;
class TTreeSRedirector;
class AliTRDtrackingEfficiency : public AliAnalysisTask
{
public:
  AliTRDtrackingEfficiency(const Char_t *name = "TRD Tracking efficiency");
  ~AliTRDtrackingEfficiency(){};
  void  ConnectInputData(Option_t *);
  void  CreateOutputObjects();
  Int_t GetDebugLevel() const {return fDebugLevel;} 
  void  Exec(Option_t *);
  void  SetDebugLevel(Int_t debug){fDebugLevel = debug;}
  void  Terminate(Option_t *);

private:
  AliTRDtrackingEfficiency(const AliTRDtrackingEfficiency&);
  AliTRDtrackingEfficiency& operator=(const AliTRDtrackingEfficiency&);

private:
  TList        *fObjectContainer;       // Container
  TObjArray        *fTracks;            // Array of tracks
  TClonesArray     *fMissed;            // Missed ?

  Int_t            fDebugLevel;         // Debug level
  TTreeSRedirector *fDebugStream;       // Debug stream

  ClassDef(AliTRDtrackingEfficiency, 1) // TRD tracking efficiency
};

#endif

