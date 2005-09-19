/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPCBPainter.h,v 1.3 2005/08/26 15:43:36 ivana Exp $

/// \ingroup graphics
/// \class AliMpPCBPainter
/// \brief Class for drawing a PCB into canvas
///
/// Authors: Laurent Aphecetche

#ifndef ALIMPPCBPAINTER_H
#define ALIMPPCBPAINTER_H

#include "AliMpVPainter.h"

class AliMpPCB;

class AliMpPCBPainter : public AliMpVPainter
{
public:
  AliMpPCBPainter(AliMpPCB* pcb);
  virtual ~AliMpPCBPainter();

  void Draw(Option_t* option);

  void Paint(Option_t* option);

  TVector2 GetDimensions() const;
  TVector2 GetPosition() const;

 private:
  AliMpPCB* fPCB;

  ClassDef(AliMpPCBPainter,1) // A painter for a PCB of stations 3,4,5
};

#endif
