// @(#) $Id$

#ifndef ALIL3DISPLAY_H
#define ALIL3DISPLAY_H

/** \class AliHLTDisplay
<pre>
//_____________________________________________________________
// AliHLTDisplay
//
// Simple display class for the HLT tracks.
</pre>
*/
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group 

#include <TObject.h>
class TGeometry;
class AliHLTSpacePointData;
class AliHLTTrackArray;

class AliHLTDisplay : public TObject {

 public:
  AliHLTDisplay();
  AliHLTDisplay(Int_t *slice, Char_t *gfile="$(ALIHLT_BASEDIR)/geo/alice.geom");
  virtual ~AliHLTDisplay();

  void Setup(Char_t *trackfile,Char_t *path,Int_t event=-1,Bool_t sp=kFALSE);
  void DisplayTracks(Int_t min_hits=10,Bool_t x3don=kTRUE,Float_t thr=0.);
  void DisplayAll(Int_t min_hits=10,Bool_t x3don=kTRUE);
  void DisplayClusters(Bool_t x3don=kTRUE);

  void DisplayClusterRow(Int_t slice,Int_t padrow,Char_t *digitsFile,Char_t *type="hist");
  void SetTracks(AliHLTTrackArray *tracks) {fTracks=tracks;}

 private:
  AliHLTDisplay(const AliHLTDisplay &/*d*/):TObject(){;}
  AliHLTDisplay& operator=(const AliHLTDisplay &/*d*/){return *this;}

  TGeometry *fGeom; //!
  AliHLTSpacePointData *fClusters[36][6]; //!
  AliHLTTrackArray *fTracks; //!
  UInt_t fNcl[36][6]; //number of cluster
  Int_t fMinSlice; //min slice
  Int_t fMaxSlice; //max slice
  
  ClassDef(AliHLTDisplay,1) //Display class
};

typedef AliHLTDisplay AliL3Display; // for backward compatibility

#endif
