#ifndef AliITSpoints_H
#define AliITSpoints_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TPolyMarker3D.h"
#include "TMarker3DBox.h"
#include "TMatrix.h"
#include "AliITS.h"
#include "AliITSdigitNew.h"
#include "AliPoints.h"

class AliITSpoints : public AliPoints {
protected:
   Int_t            fModuleIndex;      // Link to module number 
   Int_t            fRecHitIndex;      // Link to reconstructed hit number 
   Int_t            fHitIndex;         // Link to hit number 
   Int_t            fTrackIndex;       // Link to track number 
   Int_t            fDigitIndex;       // Link to digit 
  TMarker3DBox     *fMarker[3];           // pointer to  associated 3D-marker
  TMatrix          *fMatrix;           // test


  //Bool_t            fConnect;         
public:
  AliITSpoints();
  AliITSpoints(Int_t npoints);
  virtual ~AliITSpoints();

  virtual void          ExecuteEvent(Int_t event, Int_t px, Int_t py);
  virtual void          AnodeProjection(Int_t px, Int_t py);
  virtual void          TimeProjection(Int_t px, Int_t py);
  Int_t                 GetModuleIndex() {return fModuleIndex;} // *MENU* 
  Int_t                 GetHitIndex() {return fHitIndex;}
  Int_t                 GetTrackIndex(); // *MENU*
  Int_t                 GetDigitIndex() {return fDigitIndex;}
  TMarker3DBox         *GetMarker(Int_t i) {return fMarker[i];}
  AliITShit            *GetHit() const;
  AliITSdigitSDD       *GetDigit() const;
  virtual const Text_t *GetName() const;
  virtual Text_t       *GetObjectInfo(Int_t px, Int_t py);
  virtual void          InspectHit(); // *MENU*
  virtual void          DumpHit(); // *MENU*
  virtual void          InspectDigit(); // *MENU*
  virtual void          DumpDigit(); // *MENU*
  virtual void          GetCenterOfGravity(); // *MENU*
  virtual void          FindLocalMaxima(); // *MENU*
  virtual void          DisplayModule(); // *MENU*
  //virtual void          SetConnectOpt(Int_t draw=kFALSE) {fConnect = draw;}// *MENU*
  virtual void          SetModuleIndex(Int_t module) {fModuleIndex = module;}
  virtual void          SetHitIndex(Int_t hitindex) {fHitIndex = hitindex;}
  virtual void          SetTrackIndex(Int_t trackindex) {fTrackIndex = trackindex;}
  virtual void          SetDigitIndex(Int_t digitindex) {fDigitIndex = digitindex;}
  virtual void          Set3DMarker(Int_t i,TMarker3DBox *marker) {fMarker[i] = marker;}
  virtual void          SetMatrix(TMatrix *matrix) {fMatrix = matrix;}
  //  virtual void          ExecuteEvent(Int_t event, Int_t px, Int_t py);
  
  virtual void Neighbours
       (Int_t ix, Int_t iy, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]); 
  
  ClassDef(AliITSpoints,1) //Class to draw detector clusters (is PolyMarker3D)
};
#endif

