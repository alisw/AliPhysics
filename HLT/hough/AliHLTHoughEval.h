// @(#) $Id$

#ifndef ALIL3HOUGHEVAL_H
#define ALIL3HOUGHEVAL_H

#include "AliHLTRootTypes.h"

class AliHLTTrackArray;
class AliHLTHoughBaseTransformer;
class AliHLTHoughTrack;
class AliHLTDigitRowData;
class AliHLTHistogram;
class AliHLTHistogram1D;

class AliHLTHoughEval {
  
 public:
  AliHLTHoughEval(); 
  virtual ~AliHLTHoughEval();
  
  void InitTransformer(AliHLTHoughBaseTransformer *transformer);
  void GenerateLUT();
  void DisplayEtaSlice(Int_t etaindex,AliHLTHistogram *hist);
  Bool_t LookInsideRoad(AliHLTHoughTrack *track,Int_t &nrowscrossed,Int_t *rowrange,Bool_t remove=kFALSE);
#ifdef use_root
  void CompareMC(AliHLTTrackArray *tracks,Char_t *goodtracks="good_tracks",Int_t treshold=0);
#endif
  void FindEta(AliHLTTrackArray *tracks);
  
  //Getters
  AliHLTHistogram1D *GetEtaHisto(Int_t i) {if(!fEtaHistos) return 0; if(!fEtaHistos[i]) return 0; return fEtaHistos[i];}

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
  AliHLTHistogram1D **fEtaHistos; //!
  Float_t fZVertex;//Z position of the primary vertex

  //Flags
  Bool_t fRemoveFoundTracks;//Remove the found tracks or not?
  
  AliHLTHoughBaseTransformer *fHoughTransformer; //!
  AliHLTDigitRowData **fRowPointers; //!
  
  ClassDef(AliHLTHoughEval,1) //Hough transform verfication class

};

typedef AliHLTHoughEval AliL3HoughEval; // for backward comaptibility

#endif
