// @(#) $Id$

#ifndef ALIL3HOUGHEVAL_H
#define ALIL3HOUGHEVAL_H

#include "AliL3RootTypes.h"

class AliL3TrackArray;
class AliL3HoughBaseTransformer;
class AliL3HoughTrack;
class AliL3DigitRowData;
class AliL3Histogram;
class AliL3Histogram1D;

class AliL3HoughEval {
  
 public:
  AliL3HoughEval(); 
  virtual ~AliL3HoughEval();
  
  void InitTransformer(AliL3HoughBaseTransformer *transformer);
  void GenerateLUT();
  void DisplayEtaSlice(Int_t etaindex,AliL3Histogram *hist);
  Bool_t LookInsideRoad(AliL3HoughTrack *track,Int_t &nrowscrossed,Int_t *rowrange,Bool_t remove=kFALSE);
#ifdef use_root
  void CompareMC(AliL3TrackArray *tracks,Char_t *goodtracks="good_tracks",Int_t treshold=0);
#endif
  void FindEta(AliL3TrackArray *tracks);
  
  //Getters
  AliL3Histogram1D *GetEtaHisto(Int_t i) {if(!fEtaHistos) return 0; if(!fEtaHistos[i]) return 0; return fEtaHistos[i];}

  //Setters:
  void RemoveFoundTracks() {fRemoveFoundTracks = kTRUE;}
  void SetNumOfRowsToMiss(Int_t i) {fNumOfRowsToMiss = i;}
  void SetNumOfPadsToLook(Int_t i) {fNumOfPadsToLook = i;}
  void SetSlice(Int_t i) {fSlice=i;}
  void SetZVertex(Float_t zvertex) {fZVertex=zvertex;}

 private:

  Int_t fSlice;//Index of the slice being processed
  Int_t fPatch;//Index of the patch being processed
  Int_t fNrows;//Number of padrows inside the patch
  Int_t fNEtaSegments;//Number of eta slices
  Double_t fEtaMin;//Minimum allowed eta
  Double_t fEtaMax;//Maximum allowed eta
  Int_t fNumOfPadsToLook;//Padrow search window
  Int_t fNumOfRowsToMiss;//Maximum numbers of padrow which could be missed
  AliL3Histogram1D **fEtaHistos; //!
  Float_t fZVertex;//Z position of the primary vertex

  //Flags
  Bool_t fRemoveFoundTracks;//Remove the found tracks or not?
  
  AliL3HoughBaseTransformer *fHoughTransformer; //!
  AliL3DigitRowData **fRowPointers; //!
  
  ClassDef(AliL3HoughEval,1) //Hough transform verfication class

};

#endif
