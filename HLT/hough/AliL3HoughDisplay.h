#ifndef ALIL3HoughDisplay
#define ALIL3HoughDisplay

#include "AliL3RootTypes.h"

class TGeometry;
class AliL3TrackArray;
class AliL3Transform;

class AliL3HoughDisplay {

 private:
  
  TGeometry *fGeom; //!
  AliL3TrackArray *fTracks; //!
  AliL3Transform *fTransform; //!

  void GenerateHits(AliL3HoughTrack *track,Float_t *x,Float_t *y,Float_t *z,Int_t &n);
  void Init();
  
 public:
  AliL3HoughDisplay();
  virtual ~AliL3HoughDisplay();

  void DisplayTracks();
  void SetTracks(AliL3TrackArray *tracks) {fTracks = tracks;}

  ClassDef(AliL3HoughDisplay,1) 
};

#endif
