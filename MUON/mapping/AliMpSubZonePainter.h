// $Id$
// Category: graphics
//
// Class AliMpSubZonePainter
// -------------------------
// Class for drawing a subzone into canvas
//
// Authors: David Guez, IPN Orsay

#ifndef ALI_MP_SUBZONE_PAINTER_H
#define ALI_MP_SUBZONE_PAINTER_H

#include "AliMpVPainter.h"

class AliMpSubZone;

class AliMpSubZonePainter : public AliMpVPainter
{
 public:
  AliMpSubZonePainter();
  AliMpSubZonePainter(AliMpSubZone *subZone);
  virtual ~AliMpSubZonePainter();
  
  virtual void DumpObject(); //-MENU-
  virtual void Draw(Option_t *option);
  virtual void Paint(Option_t *option);
  // get/set methods
  virtual TVector2 GetPosition() const;
  virtual TVector2 GetDimensions() const;
  virtual Int_t DistancetoPrimitive(Int_t x, Int_t y);

 private: 
  AliMpSubZone *fSubZone;         // the subzone to draw
  ClassDef(AliMpSubZonePainter,1) // SubZone painter
};
#endif //ALI_MP_SUBZONE_PAINTER_H
