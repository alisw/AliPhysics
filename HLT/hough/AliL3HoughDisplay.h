#ifndef ALIL3HoughDisplay
#define ALIL3HoughDisplay

#include "AliL3RootTypes.h"

class TGeometry;
class AliL3TrackArray;
class AliL3Transform;
class AliL3DigitRowData;
class TPolyMarker3D;
class AliL3HoughTrack;

class AliL3HoughDisplay {

 private:
  
  TGeometry *fGeom; //!
  AliL3TrackArray *fTracks; //!
  AliL3Transform *fTransform; //!
  AliL3DigitRowData *fDigitRowData;  //!
  UInt_t fNDigitRowData; //!
  Int_t fShowSlice; 
  Int_t fPatch;
  
  void GenerateHits(AliL3HoughTrack *track,Float_t *x,Float_t *y,Float_t *z,Int_t &n);
  void Init();
  TPolyMarker3D *LoadDigits();

 public:
  AliL3HoughDisplay();
  virtual ~AliL3HoughDisplay();

  void DisplayEvent();
  void SetTracks(AliL3TrackArray *tracks) {fTracks = tracks;}
  void ShowData(AliL3DigitRowData *data,UInt_t size,Int_t slice,Int_t patch);

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
