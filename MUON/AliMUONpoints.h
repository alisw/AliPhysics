#ifndef AliMUONpoints_H
#define AliMUONpoints_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TPolyMarker3D.h"
#include "TMarker3DBox.h"
#include "TMatrix.h"
#include "AliMUON.h"
#include "AliPoints.h"

class AliMUONpoints : public AliPoints {
protected:
   Int_t            fHitIndex;         // Link to hit number 
   Int_t            fTrackIndex;       // Link to track number 
   Int_t            fDigitIndex;       // Link to digit 
  TMarker3DBox     *fMarker[3];           // pointer to  associated 3D-marker
  TMatrix          *fMatrix;           // test
public:
  AliMUONpoints();
  AliMUONpoints(Int_t npoints);
  virtual ~AliMUONpoints();

  Int_t                 GetHitIndex() {return fHitIndex;}
  Int_t                 GetTrackIndex(); // *MENU*
  Int_t                 GetDigitIndex() {return fDigitIndex;}
  TMarker3DBox         *GetMarker(Int_t i) {return fMarker[i];}
  AliMUONhit           *GetHit() const;
  AliMUONdigit         *GetDigit() const;
  virtual void          InspectHit(); // *MENU*
  virtual void          DumpHit(); // *MENU*
  virtual void          InspectDigit(); // *MENU*
  virtual void          DumpDigit(); // *MENU*
  virtual void          GetCenterOfGravity(); // *MENU*
  virtual void          SetHitIndex(Int_t hitindex) {fHitIndex = hitindex;}
  virtual void          SetTrackIndex(Int_t trackindex) {fTrackIndex = trackindex;}
  virtual void          SetDigitIndex(Int_t digitindex) {fDigitIndex = digitindex;}
  virtual void          Set3DMarker(Int_t i,TMarker3DBox *marker) {fMarker[i] = marker;}
  virtual void          SetMatrix(TMatrix *matrix) {fMatrix = matrix;}
  //  virtual void          ExecuteEvent(Int_t event, Int_t px, Int_t py);
  
  ClassDef(AliMUONpoints,1) //Class to draw detector clusters (is PolyMarker3D)
};
#endif

