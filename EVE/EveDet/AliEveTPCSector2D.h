// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCSector2D_H
#define AliEveTPCSector2D_H

#include "AliEveTPCSectorViz.h"

class AliEveTPCSector2DEditor;
class AliEveTPCSector2DGL;

//------------------------------------------------------------------------------
// AliEveTPCSector2D
//
// Visualization of TPC raw-data in 2D.

class AliEveTPCSector2D : public AliEveTPCSectorViz
{
  friend class AliEveTPCSector2DGL;
  friend class AliEveTPCSector2DEditor;

public:
  AliEveTPCSector2D(const Text_t* n="AliEveTPCSector2D", const Text_t* t=0);
  virtual ~AliEveTPCSector2D() {}

  void SetShowMax(Bool_t sm)  { fShowMax  = sm;  IncRTS(); }
  void SetAverage(Bool_t avg) { fAverage  = avg; IncRTS(); }

  Int_t GetPickMode() const     { return fPickMode; }
  void  SetPickMode(Int_t mode) { fPickMode = mode; }

  void MakeSector3D(); // *MENU*

  virtual void ComputeBBox();

  virtual void PadSelected(Int_t row, Int_t pad);

  virtual void Paint(Option_t* option="");

protected:
  Bool_t      fShowMax;    // Show maximum signal-value in time range.
  Bool_t      fAverage;    // Show average signal value in time range.

  Bool_t      fUseTexture; // Use texture to draw each segment.
  Bool_t      fPickEmpty;  // Pick also empty pads.
  Int_t       fPickMode;   // Pick mode: 0-print, 1-1dhisto of pad, 2-2dhisto of padrow.

  ClassDef(AliEveTPCSector2D, 0); // Visualization of TPC raw-data in 2D.
};

#endif
