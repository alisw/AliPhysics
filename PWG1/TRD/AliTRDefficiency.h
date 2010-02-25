#ifndef ALITRDEFFICIENCY_H
#define ALITRDEFFICIENCY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDefficiency.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TList;
class TClonesArray;
class TTreeSRedirector;
class AliTRDefficiency : public AliTRDrecoTask
{
public:
  AliTRDefficiency();
  AliTRDefficiency(char* name);
  virtual ~AliTRDefficiency();
  void    UserCreateOutputObjects();
  void    UserExec(Option_t *);
  Bool_t  GetRefFigure(Int_t ifig);
  Bool_t  PostProcess();

private:
  AliTRDefficiency(const AliTRDefficiency&);
  AliTRDefficiency& operator=(const AliTRDefficiency&);

private:
  TClonesArray     *fMissed;            // Missed ?

  ClassDef(AliTRDefficiency, 1) // TRD tracking efficiency
};

#endif

