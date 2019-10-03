#ifndef ALITRDMULTIPLICITY_H
#define ALITRDMULTIPLICITY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD multiplicity                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class AliTRDmultiplicity : public AliTRDrecoTask
{
public:
  AliTRDmultiplicity();
  AliTRDmultiplicity(char* name);
  virtual ~AliTRDmultiplicity();
  void    UserCreateOutputObjects();
  void    UserExec(Option_t *);

private:
  AliTRDmultiplicity(const AliTRDmultiplicity&);
  AliTRDmultiplicity& operator=(const AliTRDmultiplicity&);

private:
  Int_t fEventCounter;            // Event counter


  ClassDef(AliTRDmultiplicity, 1) // TRD tracking multiplicity
};

#endif

