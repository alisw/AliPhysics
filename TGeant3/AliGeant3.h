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

#include <TGeant3.h>

class AliGeant3 : public TGeant3
{

private:

public:
  AliGeant3(const char *title);
  AliGeant3() {}
  virtual ~AliGeant3() {}

  void   SetColors();

  //
  //
  // Control Methods

  virtual void Init();
  virtual void FinishGeometry();
  virtual void ProcessEvent();
  virtual void ProcessRun(Int_t nevent);

  ClassDef(AliGeant3,1) //Generic MonteCarlo Class

};

#endif

