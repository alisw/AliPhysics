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

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class AliESDEvent;
class AliMCEvent;
class AliESDfriend;
class AliTRDtrackInfo;
class TObjArray;
class TTreeSRedirector;

class AliTRDtrackInfoGen : public AliTRDrecoTask{
public:

  AliTRDtrackInfoGen();
  virtual ~AliTRDtrackInfoGen();
  
  void  ConnectInputData(Option_t *);
  void  CreateOutputObjects();
  void  Exec(Option_t *);
  void  Terminate(Option_t *);

private:
  AliTRDtrackInfoGen(const AliTRDtrackInfoGen&);
  AliTRDtrackInfoGen& operator=(const AliTRDtrackInfoGen&);

private:

  AliESDEvent      *fESD;                  // ESD event
  AliMCEvent       *fMC;                   // MC event
  AliESDfriend     *fESDfriend;            // ESD friends
  AliTRDtrackInfo  *fTrackInfo;            // Track info

  ClassDef(AliTRDtrackInfoGen, 1)          // entry to TRD analysis
};
#endif
