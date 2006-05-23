/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPCBPainter.h,v 1.6 2006/05/23 13:07:35 ivana Exp $

/// \ingroup graphics
/// \class AliMpPCBPainter
/// \brief Class for drawing a PCB into canvas
///
/// \author Laurent Aphecetche

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

 protected:
  AliMpPCBPainter(const AliMpPCBPainter& right);
  AliMpPCBPainter&  operator = (const AliMpPCBPainter& right);
     
 private:
  AliMpPCB* fPCB; //!< PCB to be plotted.

  ClassDef(AliMpPCBPainter,1) // A painter for a PCB of stations 3,4,5
};

#endif
