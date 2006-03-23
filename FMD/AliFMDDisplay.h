#ifndef AliFMDDISPLAY_H
#define AliFMDDISPLAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//___________________________________________________________________
//
// The classes defined here, are utility classes for reading in data
// for the FMD.  They are  put in a seperate library to not polute the
// normal libraries.  The classes are intended to be used as base
// classes for customized class that do some sort of analysis on the
// various types of data produced by the FMD. 
//
#include "AliFMDInput.h"
class TObjArray;
class TCanvas;
class TPad;
class TButton;

//___________________________________________________________________
class AliFMDDisplay : public AliFMDInput
{
public:
  AliFMDDisplay(const char* gAliceFile="galice.root");
  virtual ~AliFMDDisplay() {}
  static AliFMDDisplay* Instance();

  void  Continue() { fWait = kFALSE; }
  void  Zoom() { fZoomMode = kTRUE; }
  void  Pick() { fZoomMode = kFALSE; }
  void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
  Int_t DistancetoPrimitive(Int_t px, Int_t py);
  void  Paint(Option_t* option="") { (void)option; }

  virtual Bool_t Init();
  virtual Bool_t Begin(Int_t event);
  virtual Bool_t End();
  virtual Bool_t ProcessHit(AliFMDHit* hit, TParticle* p);
  virtual Bool_t ProcessDigit(AliFMDDigit* digit);
  virtual Bool_t ProcessRaw(AliFMDDigit* digit);
  virtual Bool_t ProcessRecPoint(AliFMDRecPoint* recpoint);
  virtual Int_t  LookupColor(Float_t x, Float_t max)  const;
protected:
  static AliFMDDisplay* fgInstance; // Static instance 
  Bool_t                fWait;      // Wait until user presses `Continue'
  TObjArray*            fMarkers;   // Cache of markers
  TObjArray*            fHits;      // Cache of `hits'
  TCanvas*              fCanvas;    // Canvas to draw in 
  TPad*                 fPad;       // View pad. 
  TButton*              fButton;    // Continue button
  TButton*              fZoom;      // Zoom button
  TButton*              fPick;      // Pick button
  Bool_t                fZoomMode;  // Whether we're in Zoom mode
  Float_t               fX0;        // X at lower left corner or range 
  Float_t               fY0;        // Y at lower left corner or range 
  Float_t               fX1;        // X at upper right corner or range 
  Float_t               fY1;        // Y at upper right corner or range 
  Int_t                 fXPixel;    // X pixel of mark
  Int_t                 fYPixel;    // Y pixel of mark
  Int_t                 fOldXPixel; // Old x pixel of mark
  Int_t                 fOldYPixel; // Old y pixel of mark
  Bool_t                fLineDrawn; // Whether we're drawing a box
  ClassDef(AliFMDDisplay,0)  // FMD specialised event display
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
