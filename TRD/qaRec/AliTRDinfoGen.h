#ifndef ALITRDINFOGEN_H
#define ALITRDINFOGEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDinfoGen.h 27496 2008-07-22 08:35:45Z cblume $ */

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
class AliTRDeventInfo;
class TObjArray;
class TTreeSRedirector;
class AliTRDinfoGen : public AliTRDrecoTask{
public:

  AliTRDinfoGen();
  virtual ~AliTRDinfoGen();
  
  void  ConnectInputData(Option_t *);
  void  CreateOutputObjects();
  void  Exec(Option_t *);
  void  Terminate(Option_t *);

  static const Float_t xTPC;
  static const Float_t xTOF;

private:
  AliTRDinfoGen(const AliTRDinfoGen&);
  AliTRDinfoGen& operator=(const AliTRDinfoGen&);

  AliESDEvent      *fESD;                  // ESD event
  AliMCEvent       *fMC;                   // MC event
  AliESDfriend     *fESDfriend;            // ESD friends
  AliTRDtrackInfo  *fTrackInfo;            // Track info
  AliTRDeventInfo  *fEventInfo;		   // Event info

  ClassDef(AliTRDinfoGen, 1)          // entry to TRD analysis
};
#endif
