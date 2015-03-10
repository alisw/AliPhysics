#ifndef ALIMPITERATORPAINTER_H
#define ALIMPITERATORPAINTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup mpgraphics
/// \class AliMpIteratorPainter
/// \brief Painter for a group of pads defined by an iterator
/// 
//  Author Laurent Aphecetche

#ifndef ALI_MP_V_PAINTER_H
#  include "AliMpVPainter.h"
#endif

#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif

class TObjArray;
class AliMpVPadIterator;

class AliMpIteratorPainter : public AliMpVPainter
{
public:
  AliMpIteratorPainter(AliMpVPadIterator* it);
  virtual ~AliMpIteratorPainter();
  
  void Draw(Option_t* option);
  void Paint(Option_t* option);
  
  TVector2 GetDimensions() const { return fDimensions; }
  TVector2 GetPosition() const { return fPosition; }

private:
  /// Not implemented
  AliMpIteratorPainter();
  /// Not implemented
  AliMpIteratorPainter(const AliMpIteratorPainter&);
  /// Not implemented
  AliMpIteratorPainter& operator=(const AliMpIteratorPainter&);
  
  TObjArray* fPads; //!<! pads of the iterator
  TVector2 fPosition; //!<! position
  TVector2 fDimensions; //!<! dimension
  
  ClassDef(AliMpIteratorPainter,1) // Painter for a group of pads
};

#endif
