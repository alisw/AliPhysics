/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifPainter.h,v 1.8 2006/05/24 13:58:13 ivana Exp $

/// \ingroup graphics
/// \class AliMpMotifPainter
/// \brief Class for drawing a motif into canvas
///
/// \author David Guez, IPN Orsay

#ifndef ALI_MP_MOTIF_PAINTER_H
#define ALI_MP_MOTIF_PAINTER_H

#include "AliMpVPainter.h"

class AliMpMotifPosition;
class AliMpMotifType;

class AliMpMotifPainter : public AliMpVPainter
{
 public:
  AliMpMotifPainter();
  AliMpMotifPainter(AliMpMotifType* motifType);
  AliMpMotifPainter(AliMpMotifPosition *motifPos);
  virtual ~AliMpMotifPainter();
  
  virtual void DumpObject(); //-MENU-
  virtual void Paint(Option_t *option);
  virtual TVector2 GetPosition() const;
  virtual TVector2 GetDimensions() const;

 private:
  AliMpMotifPainter(const AliMpMotifPainter& right);
  AliMpMotifPainter&  operator = (const AliMpMotifPainter& right);

  AliMpMotifPosition *fMotifPos; ///< the motif to draw
  
  ClassDef(AliMpMotifPainter,1) // Motif painter
};
#endif //ALI_MP_MOTIF_PAINTER_H
