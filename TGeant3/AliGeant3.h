#ifndef ALIGEANT3_H
#define ALIGEANT3_H
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

#include <AliVMC.h>

class AliGeant3 : public AliVMC
{

private:

public:
  AliGeant3(const char *title);
  AliGeant3();
  virtual ~AliGeant3() {}
  //
  //
  void FinishGeometry();
  void BuildPhysics();
  void ProcessEvent();
  ClassDef(AliGeant3,1) //Generic MonteCarlo Class

};

#endif

