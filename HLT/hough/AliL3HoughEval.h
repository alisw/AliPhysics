// @(#) $Id$

#ifndef ALIL3_HOUGH_Eval
#define ALIL3_HOUGH_Eval

#include "AliL3RootTypes.h"

class AliL3TrackArray;
class AliL3HoughBaseTransformer;
class AliL3HoughTrack;
class AliL3DigitRowData;
class AliL3Histogram;
class AliL3Histogram1D;

class AliL3HoughEval {
  
 private:

  Int_t fSlice;
  Int_t fPatch;
  Int_t fNrows;
  Int_t fNEtaSegments;
  Double_t fEtaMin;
  Double_t fEtaMax;
  Int_t fNumOfPadsToLook;
  Int_t fNumOfRowsToMiss;
  AliL3Histogram1D **fEtaHistos; //!
  Float_t fZVertex;

  //Flags
  Bool_t fRemoveFoundTracks;
  
  AliL3HoughBaseTransformer *fHoughTransformer; //!
  AliL3DigitRowData **fRowPointers; //!
  
 public:
  AliL3HoughEval(); 
  virtual ~AliL3HoughEval();
  
  void InitTransformer(AliL3HoughBaseTransformer *transformer);
  void GenerateLUT();
  void DisplayEtaSlice(Int_t eta_index,AliL3Histogram *hist);
  Bool_t LookInsideRoad(AliL3HoughTrack *track,Int_t &nrows_crossed,Int_t *rowrange,Bool_t remove=kFALSE);
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

  ClassDef(AliL3HoughEval,1) //Hough transform verfication class

};

#endif
