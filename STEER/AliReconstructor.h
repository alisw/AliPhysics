#ifndef ALIRECONSTRUCTOR_H
#define ALIRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// base class for reconstruction algorithm
// Derived classes should implement a default constructor and
// the virtual methods
//

#include <TObject.h>

class AliRunLoader;
class AliVertexer;
class AliTracker;
class AliESD;


class AliReconstructor: public TObject {
public:
  AliReconstructor(): TObject() {};
  virtual ~AliReconstructor() {};

  virtual void         Reconstruct(AliRunLoader* runLoader) const = 0;
  virtual AliVertexer* CreateVertexer(AliRunLoader* /*runLoader*/) const 
    {return NULL;}
  virtual AliTracker*  CreateTracker(AliRunLoader* /*runLoader*/) const 
    {return NULL;}
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const = 0;

  virtual const char*  GetDetectorName() const;

  ClassDef(AliReconstructor, 0)   // base class for reconstruction algorithms
};

#endif
