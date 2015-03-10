/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSlatPainter.h,v 1.10 2006/05/24 13:58:13 ivana Exp $

/// \ingroup mpgraphics
/// \class AliMpSlatPainter
/// \brief Class for drawing a slat into canvas
///
//  Author: Laurent Aphecetche

#ifndef ALIMPSLATPAINTER_H
#define ALIMPSLATPAINTER_H

#ifndef ALI_MP_V_PAINTER_H
#  include "AliMpVPainter.h"
#endif

class AliMpSlat;

class AliMpSlatPainter : public AliMpVPainter
{
 public:
  AliMpSlatPainter();
  AliMpSlatPainter(const AliMpSlat* slat);
  virtual ~AliMpSlatPainter();

  TVector2 GetDimensions() const;

  TVector2 GetPosition() const;

  void Draw(Option_t* option);

  void Paint(Option_t* option);

 private:
  /// Not implemented
  AliMpSlatPainter(const AliMpSlatPainter& right);
  /// Not implemented
  AliMpSlatPainter&  operator = (const AliMpSlatPainter& right);

  const AliMpSlat* fkSlat; //!<! pointer to the slat to be drawn

  ClassDef(AliMpSlatPainter,1) // A painter for a slat of stations 3,4,5
};

#endif
