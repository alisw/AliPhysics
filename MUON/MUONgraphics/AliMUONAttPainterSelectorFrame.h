#ifndef ALIMUONATTPAINTERSELECTORFRAME_H
#define ALIMUONATTPAINTERSELECTORFRAME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONAttPainterSelectorFrame
/// \brief Widget to select the painter(s) view type
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include "TGFrame.h"
#endif
#ifndef ALIMUONATTPAINTER_H
#  include "AliMUONAttPainter.h"
#endif

class TGButtonGroup;

class AliMUONAttPainterSelectorFrame : public TGHorizontalFrame
{
public:
  AliMUONAttPainterSelectorFrame(TGWindow* p=0x0, UInt_t w=1, UInt_t h=1);
  virtual ~AliMUONAttPainterSelectorFrame();
  
  void Update(const AliMUONAttPainter& att);
  
  void Clicked(const AliMUONAttPainter* newValues); // *SIGNAL*
  
  void CathodeClicked(Int_t buttonId);
  
  void PlaneClicked(Int_t buttonId);
  
  void ViewClicked(Int_t buttonId);
  
private:
  /// Not implemented
  AliMUONAttPainterSelectorFrame(const AliMUONAttPainterSelectorFrame& rhs);
  /// Not implemented
  AliMUONAttPainterSelectorFrame& operator=(const AliMUONAttPainterSelectorFrame& rhs);
  
private:
  
  TGButtonGroup* fCathode; ///< cathode selection buttons
  TGButtonGroup* fPlane;   ///< plane selection buttons
  TGButtonGroup* fViewPoint; ///< viewpoint selection buttons
  
  AliMUONAttPainter fAttributes; ///< attributes
  
  ClassDef(AliMUONAttPainterSelectorFrame,1) // Widget to select painter view type
};

#endif
