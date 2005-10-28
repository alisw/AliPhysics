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
//         Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TObject.h>
#include <TCanvas.h>
#include <TH2.h>

class TGeometry;
class AliHLTTPCSpacePointData;
class AliHLTTPCTrackArray;

class AliHLTTPCDisplay : public TObject {

 public:
  AliHLTTPCDisplay();
  AliHLTTPCDisplay(Int_t *slice, Char_t *gfile="$(ALIHLT_BASEDIR)/geo/alice.geom");
  virtual ~AliHLTTPCDisplay();

  // SETUP
  void Setup(Char_t *trackfile,Char_t *path,Int_t event=-1,Bool_t sp=kFALSE);
  void SetupClusterDataForPatch(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data);
  void SetupPadRow(Int_t histSwitch, Int_t slice, Int_t padrow);

  // SETTER
  void SetSlices(Int_t minslice, Int_t maxslice);
  void SetSlices(Int_t slice);
  void SetSlices();
  void SetSlicesPair(Int_t slice);
  void SetSlicesPair(Int_t minslice, Int_t maxslice);
  void SetInvert(Bool_t invert=kTRUE);
  void SetDrawGeo(Bool_t drawgeo=kTRUE);
  void SetTracks(AliHLTTPCTrackArray *tracks) {fTracks=tracks;}

  // GETTER
  Int_t GetPadrow(){return fPadRow;}
  Int_t GetSlice(){return fSlice;}

  void DrawGeom(Int_t slice);
  Bool_t LoadGeometrie(Char_t *gfile);
  void DisplayGeom(Bool_t x3don=kTRUE);
  void ResetDisplay(){ fc1->Clear(); }
  void DisplayAll(Int_t minhits=10,Bool_t clusterswitch=kTRUE,Bool_t trackswitch=kTRUE,Bool_t x3don=kTRUE, Float_t thr=0., Float_t* etaRange=NULL);


  // PADROW
  void FillPadRow(Int_t patch, ULong_t dataBlock, ULong_t dataLen, UInt_t nofClusters, AliHLTTPCSpacePointData* data);
  void DrawPadRow(Bool_t x3don=kTRUE);

 private:
  AliHLTTPCDisplay(const AliHLTTPCDisplay &/*d*/):TObject(){;}
  AliHLTTPCDisplay& operator=(const AliHLTTPCDisplay &/*d*/){return *this;}

  TCanvas *fc1;
  TGeometry *fGeom; //!

  AliHLTTPCSpacePointData *fClusters[36][6]; //!
  AliHLTTPCTrackArray *fTracks; //!

  UInt_t fNcl[36][6]; //number of cluster
  Int_t fMinSlice; //min slice
  Int_t fMaxSlice; //max slice
  Int_t fBackColor; //Background color
  Int_t fLineColor; //Line color
  Int_t fSlicePair; //draw pair of slices;
  Float_t fTheta;
  Int_t fSlicePairMax;
  Int_t fSlicePairMin;
  Bool_t fDrawGeo;

  // PADROW
  TH1F *fHistrawcl;  // Histogram for cluster in padrow
  TH2F *fHistraw;    // Histogram for signals in padrow
  Int_t fSlice;      // slice
  Int_t fPadRow;     // padrow
  Int_t fNPads;      // number of pads in padrow
  Int_t fNTimes;     // number of timebins
  Int_t fhistSwitch; // switch between histogram and geometrie
  Int_t fcolorbin[20]; // number of entries per colorbin
  Int_t fbinct[20];    // index of colorbin
  Float_t *fpmarr[20]; // contains point data

  ClassDef(AliHLTTPCDisplay,1) //Display class
};

#endif
