// $Id$
// Category: graphics
//
// Class AliMpZonePainter
// ----------------------
// Class for drawing a zone into canvas
//
// Authors: David Guez, IPN Orsay

#ifndef ALI_MP_ZONE_PAINTER_H
#define ALI_MP_ZONE_PAINTER_H

#include "AliMpVPainter.h"

class AliMpZone;

class AliMpZonePainter : public AliMpVPainter
{
 public:
  AliMpZonePainter();
  AliMpZonePainter(AliMpZone *zone);
  virtual ~AliMpZonePainter();
  
  virtual void DumpObject(); //*MENU*
  virtual void Draw(Option_t *option);
  virtual void Paint(Option_t *option);
  // get/set methods

  virtual TVector2 GetPosition() const;
  virtual TVector2 GetDimensions() const;
  virtual Int_t DistancetoPrimitive(Int_t x, Int_t y);
 private: 
  AliMpZone *fZone;            // the zone to draw
  ClassDef(AliMpZonePainter,1) // Zone painter
};
#endif //ALI_MP_ZONE_PAINTER_H
