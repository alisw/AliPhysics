// $Header$

#ifndef REVE_Pad_H
#define REVE_Pad_H

#include <Reve/Reve.h>

#include <TPad.h>

namespace Reve {

class Pad : public TPad 
{
public:
  Pad();
  Pad(const char* name, const char* title,
      Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
      Color_t color = -1, Short_t bordersize = -1, Short_t bordermode = -2);
  virtual ~Pad() {}

  ClassDef(Pad, 1); // Wrapper for TPad
};

}

#endif
