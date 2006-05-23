/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpRowSegmentPainter.h,v 1.7 2006/05/23 13:07:35 ivana Exp $

/// \ingroup graphics
/// \class AliMpRowSegmentPainter
/// \brief Class for drawing a motif into canvas
///
/// Authors: David Guez, IPN Orsay

#ifndef ALI_MP_ROW_SEGMENT_PAINTER_H
#define ALI_MP_ROW_SEGMENT_PAINTER_H

#include "AliMpVPainter.h"

class AliMpVRowSegment;

class AliMpRowSegmentPainter : public AliMpVPainter
{
 public:
  AliMpRowSegmentPainter();
  AliMpRowSegmentPainter(AliMpVRowSegment *rowSegment);
  virtual ~AliMpRowSegmentPainter();
  
  virtual void DumpObject(); //-MENU-
  virtual void Draw(Option_t* option);
  virtual void Paint(Option_t* /*option*/);
  virtual TVector2 GetPosition() const;
  virtual TVector2 GetDimensions() const;

 protected:
  AliMpRowSegmentPainter(const AliMpRowSegmentPainter& right);
  AliMpRowSegmentPainter&  operator = (const AliMpRowSegmentPainter& right);

 private: 
  AliMpVRowSegment *fRowSegment; ///< the row segment to draw

  ClassDef(AliMpRowSegmentPainter,1) // Row Segment painter
};
#endif //ALI_MP_ROW_SEGMENT_PAINTER_H
