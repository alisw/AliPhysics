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

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TList;
class TClonesArray;
class TTreeSRedirector;
class AliTRDtrackingEfficiency : public AliTRDrecoTask
{
public:
  AliTRDtrackingEfficiency();
  virtual ~AliTRDtrackingEfficiency();
  void    CreateOutputObjects();
  void    Exec(Option_t *);
  Bool_t  GetRefFigure(Int_t ifig);
  Bool_t  PostProcess();
  void    Terminate(Option_t *);

private:
  AliTRDtrackingEfficiency(const AliTRDtrackingEfficiency&);
  AliTRDtrackingEfficiency& operator=(const AliTRDtrackingEfficiency&);

private:
  TClonesArray     *fMissed;            // Missed ?

  ClassDef(AliTRDtrackingEfficiency, 1) // TRD tracking efficiency
};

#endif

