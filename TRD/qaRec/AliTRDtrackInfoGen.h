#ifndef ALITRDTRACKINFOGEN_H
#define ALITRDTRACKINFOGEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackInfoGen.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"

class AliESDEvent;
class AliMCEvent;
class AliESDfriend;
class AliTRDtrackInfo;
class TObjArray;
class TTreeSRedirector;

class AliTRDtrackInfoGen : public AliAnalysisTask{
public:

  AliTRDtrackInfoGen(const Char_t *name = "TRD Track Info");
  ~AliTRDtrackInfoGen(){};
  
  void  ConnectInputData(Option_t *);
  void  CreateOutputObjects();
  Int_t GetDebugLevel() const {return fDebugLevel;} 
  Bool_t HasMCdata() const { return fHasMCdata; }
  void  Exec(Option_t *);
  void  SetDebugLevel(Int_t level);
  void  SetMCdata(Bool_t mcdata) { fHasMCdata = mcdata; };
  void  Terminate(Option_t *);

private:

  AliTRDtrackInfoGen(const AliTRDtrackInfoGen&);
  AliTRDtrackInfoGen& operator=(const AliTRDtrackInfoGen&);

private:

  AliESDEvent      *fESD;                  // ESD event
  AliMCEvent       *fMC;                   // MC event
  AliESDfriend     *fESDfriend;            // ESD friends
  AliTRDtrackInfo  *fTrackInfo;            // Track info
  TObjArray        *fObjectContainer;      // Object container

	Bool_t           fHasMCdata;             // Contains MC information
  Int_t            fDebugLevel;            // Debug level
  TTreeSRedirector *fDebugStream;          // Debug stream

  ClassDef(AliTRDtrackInfoGen, 1)          // entry to TRD analysis
};
#endif
