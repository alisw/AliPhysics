#ifndef ALIMUONPOINTS_H
#define ALIMUONPOINTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliMUONDigit;
class AliMUONHit;

class TMatrix;

#include <TMarker3DBox.h>
#include "AliPoints.h"

class AliMUONPoints : public AliPoints {
    public:
  AliMUONPoints();
  AliMUONPoints(Int_t npoints);
  AliMUONPoints(const AliMUONPoints& points);  
  virtual ~AliMUONPoints();

  Int_t                 GetHitIndex() {return fHitIndex;}
  Int_t                 GetTrackIndex(); // *MENU*
  Int_t                 GetDigitIndex() {return fDigitIndex;}
  TMarker3DBox         *GetMarker(Int_t i) {return fMarker[i];}
  AliMUONHit           *GetHit() const;
  AliMUONDigit         *GetDigit() const;
  virtual void          InspectHit(); // *MENU*
  virtual void          DumpHit(); // *MENU*
  virtual void          InspectDigit(); // *MENU*
  virtual void          DumpDigit(); // *MENU*
  virtual void          SetHitIndex(Int_t hitindex) {fHitIndex = hitindex;}
  virtual void          SetTrackIndex(Int_t trackindex) {fTrackIndex = trackindex;}
  virtual void          SetDigitIndex(Int_t digitindex) {fDigitIndex = digitindex;}
  virtual void          Set3DMarker(Int_t i,TMarker3DBox *marker) {fMarker[i] = marker;}
  virtual void          SetMatrix(TMatrix *matrix) {fMatrix = matrix;}
  AliMUONPoints& operator = (const AliMUONPoints& rhs);
  
protected:
   Int_t            fHitIndex;         // Link to hit number 
   Int_t            fTrackIndex;       // Link to track number 
   Int_t            fDigitIndex;       // Link to digit 
  TMarker3DBox     *fMarker[3];        // pointer to  associated 3D-marker
  TMatrix          *fMatrix;           // test
  
  ClassDef(AliMUONPoints,1) //Class to draw detector clusters (is PolyMarker3D) for MUON
};
#endif


