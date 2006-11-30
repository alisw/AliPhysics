// @(#) $Id$

#ifndef ALIL3HOUGHDISPLAY_H
#define ALIL3HOUGHDISPLAY_H

class TGeometry;
class AliHLTTrackArray;
class AliHLTDigitRowData;
class TPolyMarker3D;
class AliHLTTrack;

class AliHLTHoughDisplay {

 public:
  AliHLTHoughDisplay();
  virtual ~AliHLTHoughDisplay();
  
  void Init(Char_t *trackfile, Char_t *gfile="$(LEVEL3)/GEO/alice.geom");
  void DisplayEvent();
  void ShowData(AliHLTDigitRowData *data,UInt_t size,Int_t slice,Int_t patch);

 private:
  
  TGeometry *fGeom; //!
  AliHLTTrackArray *fTracks; //!
  AliHLTDigitRowData *fDigitRowData;  //!
  UInt_t fNDigitRowData; //!
  Int_t fShowSlice; //Which slice to show
  Int_t fPatch;//Which patch to show
  
  void GenerateHits(AliHLTTrack *track,Float_t *x,Float_t *y,Float_t *z,Int_t &n);
  
  TPolyMarker3D *LoadDigits();

  ClassDef(AliHLTHoughDisplay,1) 
};

typedef AliHLTHoughDisplay AliL3HoughDisplay; // for backward comaptibility

inline void AliHLTHoughDisplay::ShowData(AliHLTDigitRowData *data,UInt_t size,Int_t slice,Int_t patch)
{
  fShowSlice = slice;
  fPatch = patch;
  fDigitRowData = data;
  fNDigitRowData = size;
}

#endif
