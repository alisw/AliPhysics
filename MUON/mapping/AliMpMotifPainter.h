// $Id$
// Category: graphics
//
// Class AliMpMotifPainter
// -----------------------
// Class for drawing a motif into canvas
//
// Authors: David Guez, IPN Orsay

#ifndef ALI_MP_MOTIF_PAINTER_H
#define ALI_MP_MOTIF_PAINTER_H

#include "AliMpVPainter.h"

class AliMpMotifPosition;

class AliMpMotifPainter : public AliMpVPainter
{
 public:
  AliMpMotifPainter();
  AliMpMotifPainter(AliMpMotifPosition *motifPos);
  virtual ~AliMpMotifPainter();
  
  virtual void DumpObject(); //-MENU-
  virtual void Paint(Option_t *option);
  virtual TVector2 GetPosition() const;
  virtual TVector2 GetDimensions() const;

 private:
  AliMpMotifPosition *fMotifPos;          // the motif to draw

  ClassDef(AliMpMotifPainter,1) // Motif painter
};
#endif //ALI_MP_MOTIF_PAINTER_H
