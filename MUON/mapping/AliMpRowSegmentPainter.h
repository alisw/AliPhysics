// $Id$
// Category: graphics
//
// Class AliMpRowSegmentPainter
// ----------------------------
// Class for drawing a motif into canvas
//
// Authors: David Guez, IPN Orsay

#ifndef M_ROW_SEGMENT_PAINTER_H
#define M_ROW_SEGMENT_PAINTER_H

#include "AliMpVPainter.h"

class AliMpVRowSegment;

class AliMpRowSegmentPainter : public AliMpVPainter
{
 public:
  AliMpRowSegmentPainter();
  AliMpRowSegmentPainter(AliMpVRowSegment *rowSegment);
  virtual void DumpObject(); //*MENU*
  virtual void Draw(Option_t *option);
  virtual void Paint(Option_t *option);
  virtual TVector2 GetPosition() const;
  virtual TVector2 GetDimensions() const;
 private: 
  AliMpVRowSegment *fRowSegment;      // the row segment to draw
  ClassDef(AliMpRowSegmentPainter,1) // Row Segment painter
};
#endif //M_ROW_SEGMENT_PAINTER_H
