#ifndef ALIVMC_H
#define ALIVMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//    Generic interface to MC for AliRoot                                    //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class AliVMC : public TNamed 
{

private:
  static AliVMC* fgVMC;

public:
  AliVMC(const char *name, const char *title);
  AliVMC();
  virtual ~AliVMC() {}
  //Generic access functions
  static inline AliVMC* GetVMC() {return fgVMC;}
  //Generic Alice MonteCarlo Functions
  virtual void FinishGeometry() = 0;
  virtual void BuildPhysics() = 0;
  virtual void ProcessEvent() = 0;

  //
  ClassDef(AliVMC,1) //Generic MonteCarlo Class

};

R__EXTERN AliVMC *gVMC;

#endif

