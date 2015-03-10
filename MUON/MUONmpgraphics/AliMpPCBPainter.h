/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPCBPainter.h,v 1.7 2006/05/24 13:58:13 ivana Exp $

/// \ingroup mpgraphics
/// \class AliMpPCBPainter
/// \brief Class for drawing a PCB into canvas
///
//  Author: Laurent Aphecetche

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
  /// Not implemented
  AliMpPCBPainter(const AliMpPCBPainter& right);
  /// Not implemented
  AliMpPCBPainter&  operator = (const AliMpPCBPainter& right);

  AliMpPCB* fPCB; //!<! PCB to be plotted.

  ClassDef(AliMpPCBPainter,1) // A painter for a PCB of stations 3,4,5
};

#endif
