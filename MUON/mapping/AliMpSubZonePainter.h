/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSubZonePainter.h,v 1.8 2006/05/24 13:58:13 ivana Exp $

/// \ingroup mpgraphics
/// \class AliMpSubZonePainter
/// \brief Class for drawing a subzone into canvas
///
/// \author David Guez, IPN Orsay

#ifndef ALI_MP_SUB_ZONE_PAINTER_H
#define ALI_MP_SUB_ZONE_PAINTER_H

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
  /// Not implemented
  AliMpSubZonePainter(const AliMpSubZonePainter& right);
  /// Not implemented
  AliMpSubZonePainter&  operator = (const AliMpSubZonePainter& right);

  AliMpSubZone *fSubZone; ///< the subzone to draw

  ClassDef(AliMpSubZonePainter,1) // SubZone painter
};
#endif //ALI_MP_SUB_ZONE_PAINTER_H
