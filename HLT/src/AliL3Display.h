// @(#) $Id$

#ifndef ALIL3DISPLAY_H
#define ALIL3DISPLAY_H

/** \class AliL3Display
<pre>
//_____________________________________________________________
// AliL3Display
//
// Simple display class for the HLT tracks.
</pre>
*/
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group 

#include <TObject.h>
class TGeometry;
class AliL3SpacePointData;
class AliL3TrackArray;

class AliL3Display : public TObject {

 public:
  AliL3Display();
  AliL3Display(Int_t *slice, Char_t *gfile="$(ALIHLT_BASEDIR)/geo/alice.geom");
  virtual ~AliL3Display();

  void Setup(Char_t *trackfile,Char_t *path,Int_t event=-1,Bool_t sp=kFALSE);
  void DisplayTracks(Int_t min_hits=10,Bool_t x3don=kTRUE,Float_t thr=0.);
  void DisplayAll(Int_t min_hits=10,Bool_t x3don=kTRUE);
  void DisplayClusters(Bool_t x3don=kTRUE);

  void DisplayClusterRow(Int_t slice,Int_t padrow,Char_t *digitsFile,Char_t *type="hist");
  void SetTracks(AliL3TrackArray *tracks) {fTracks=tracks;}

 private:
  AliL3Display(const AliL3Display &/*d*/):TObject(){;}
  AliL3Display& operator=(const AliL3Display &/*d*/){return *this;}

  TGeometry *fGeom; //!
  AliL3SpacePointData *fClusters[36][6]; //!
  AliL3TrackArray *fTracks; //!
  UInt_t fNcl[36][6]; //number of cluster
  Int_t fMinSlice; //min slice
  Int_t fMaxSlice; //max slice
  
  ClassDef(AliL3Display,1) //Display class
};

#endif
