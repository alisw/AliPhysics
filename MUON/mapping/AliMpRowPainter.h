/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpRowPainter.h,v 1.7 2006/05/23 13:07:35 ivana Exp $

/// \ingroup graphics
/// \class AliMpRowPainter
/// \brief Class for drawing a row into canvas
///
/// Authors: David Guez, IPN Orsay

#ifndef ALI_MP_ROW_PAINTER_H
#define ALI_MP_ROW_PAINTER_H

#include "AliMpVPainter.h"

class AliMpRow;

class AliMpRowPainter : public AliMpVPainter
{
 public:
  AliMpRowPainter();
  AliMpRowPainter(AliMpRow *row);
  virtual ~AliMpRowPainter();
  
  virtual void DumpObject(); //-MENU-
  virtual void Draw(Option_t *option);
  virtual void Paint(Option_t *option);
  virtual TVector2 GetPosition() const;
  virtual TVector2 GetDimensions() const;

 protected:
  AliMpRowPainter(const AliMpRowPainter& right);
  AliMpRowPainter&  operator = (const AliMpRowPainter& right);

 private: 
  AliMpRow *fRow;             ///< the row to paint
  
  ClassDef(AliMpRowPainter,1) // Row painter
};
#endif //ALI_MP_ROW_PAINTER_H
