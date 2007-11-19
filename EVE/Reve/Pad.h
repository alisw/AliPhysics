// $Header$

#ifndef REVE_Pad_H
#define REVE_Pad_H

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

  virtual Bool_t    IsBatch() const { return kTRUE; }

  virtual void      Update() { PaintModified(); }

  virtual TVirtualViewer3D *GetViewer3D(Option_t * /*type*/ = "")
  { return fViewer3D; }

  ClassDef(Pad, 1); // Internal Reve pad (sub-class of TPad).
};

}

#endif
