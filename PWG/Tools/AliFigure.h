/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Compatibility wrapper to provide AliCanvas as AliFigure
//
// Author: Jochen Klein <jochen.klein@cern.ch>

#ifndef ALIFIGURE_H
#define ALIFIGURE_H

#include "AliCanvas.h"

class AliFigure : public AliCanvas
{
 public:
  AliFigure(const char* name = "", const char* title = "", Int_t ww = 800, Int_t wh = 600);

  ClassDef(AliFigure,1);
};

#endif
