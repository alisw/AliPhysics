#ifndef ALIPIPEVTEMP_H
#define ALIPIPEVTEMP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and class for detector: PIPE  version 0    //
/////////////////////////////////////////////////////////
 
#include "AliPIPE.h"

class AliPIPEvTemp : public AliPIPE {
  
public:
  AliPIPEvTemp();
  AliPIPEvTemp(const char *name, const char *title);
  virtual       ~AliPIPEvTemp() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void  DrawModule();
 private:
  virtual void  Undulation(char *undul, Float_t pitch, Float_t thick, Float_t zundul, Float_t rundul,
                           char (*cone)[5]);
   ClassDef(AliPIPEvTemp,1)  //Class for PIPE temporary version
};
 
#endif
