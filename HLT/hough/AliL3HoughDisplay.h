// @(#) $Id$

#ifndef ALIL3HOUGHDISPLAY_H
#define ALIL3HOUGHDISPLAY_H

#include "AliL3RootTypes.h"

class TGeometry;
class AliL3TrackArray;
class AliL3DigitRowData;
class TPolyMarker3D;
class AliL3Track;

class AliL3HoughDisplay {

 public:
  AliL3HoughDisplay();
  virtual ~AliL3HoughDisplay();
  
  void Init(Char_t *trackfile, Char_t *gfile="$(LEVEL3)/GEO/alice.geom");
  void DisplayEvent();
  void ShowData(AliL3DigitRowData *data,UInt_t size,Int_t slice,Int_t patch);

 private:
  
  TGeometry *fGeom; //!
  AliL3TrackArray *fTracks; //!
  AliL3DigitRowData *fDigitRowData;  //!
  UInt_t fNDigitRowData; //!
  Int_t fShowSlice; //Which slice to show
  Int_t fPatch;//Which patch to show
  
  void GenerateHits(AliL3Track *track,Float_t *x,Float_t *y,Float_t *z,Int_t &n);
  
  TPolyMarker3D *LoadDigits();

  ClassDef(AliL3HoughDisplay,1) 
};

inline void AliL3HoughDisplay::ShowData(AliL3DigitRowData *data,UInt_t size,Int_t slice,Int_t patch)
{
  fShowSlice = slice;
  fPatch = patch;
  fDigitRowData = data;
  fNDigitRowData = size;
}

#endif
