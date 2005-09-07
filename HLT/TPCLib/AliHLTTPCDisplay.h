// @(#) $Id$

#ifndef ALIHLTTPCDISPLAY_H
#define ALIHLTTPCDISPLAY_H

/** \class AliHLTTPCDisplay
<pre>
//_____________________________________________________________
// AliHLTTPCDisplay
//
// Simple display class for the HLT tracks.
</pre>
*/
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group 

#include <TObject.h>
#include <TCanvas.h>
class TGeometry;
class AliHLTTPCSpacePointData;
class AliHLTTPCTrackArray;

class AliHLTTPCDisplay : public TObject {

 public:
  AliHLTTPCDisplay();
  AliHLTTPCDisplay(Int_t *slice, Char_t *gfile="$(ALIHLT_BASEDIR)/geo/alice.geom");
  virtual ~AliHLTTPCDisplay();

  void Setup(Char_t *trackfile,Char_t *path,Int_t event=-1,Bool_t sp=kFALSE);
  void SetupClusterDataForPatch(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data);
  Bool_t SetSlices(Int_t minslice, Int_t maxslice);
  Bool_t LoadGeometrie(Char_t *gfile);
  void DisplayTracks(Int_t min_hits=10,Bool_t x3don=kTRUE,Float_t thr=0.);
  void DisplayAll(Int_t min_hits=10,Bool_t x3don=kTRUE,Float_t* etaRange=NULL);
  void DisplayClusters(Bool_t x3don=kTRUE,Float_t* etaRange=NULL);

  void DisplayClusterRow(Int_t slice,Int_t padrow,Char_t *digitsFile,Char_t *type="hist");
  void SetTracks(AliHLTTPCTrackArray *tracks) {fTracks=tracks;}

 private:
  AliHLTTPCDisplay(const AliHLTTPCDisplay &/*d*/):TObject(){;}
  AliHLTTPCDisplay& operator=(const AliHLTTPCDisplay &/*d*/){return *this;}

  TGeometry *fGeom; //!
  AliHLTTPCSpacePointData *fClusters[36][6]; //!
  AliHLTTPCTrackArray *fTracks; //!
  UInt_t fNcl[36][6]; //number of cluster
  Int_t fMinSlice; //min slice
  Int_t fMaxSlice; //max slice
 
  TCanvas *fc1;
  ClassDef(AliHLTTPCDisplay,1) //Display class
};

#endif
