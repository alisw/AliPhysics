#ifndef ALIL3_HOUGH
#define ALIL3_HOUGH

#include "AliL3RootTypes.h"

class AliL3HoughMaxFinder;
class AliL3HoughTransformer;
class AliL3Histogram;
class AliL3FileHandler;
class AliL3HoughEval;
class AliL3Transform;
class AliL3TrackArray;
class TFile;

class AliL3Hough : public TObject {
  
 private:
  Char_t fPath[256];
  Bool_t fBinary;
  Bool_t fAddHistograms;
  Bool_t fRemoveFoundTracks;
  Bool_t fWriteDigits;
  Int_t fNEtaSegments;
  AliL3FileHandler **fMemHandler; //!
  AliL3HoughTransformer **fHoughTransformer; //!
  AliL3HoughEval **fEval; //!
  AliL3HoughMaxFinder *fPeakFinder; //!
  AliL3TrackArray *fTracks; //!
  TFile *fRootFile; //!
  
  void DeleteEval();
  void DeleteTransformers();
  void DeleteMemory();
  void Init();
  
 public:
  
  AliL3Hough(); 
  AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments=100);
  virtual ~AliL3Hough();
  
  void Process(Int_t minslice,Int_t maxslice);
  void TransformSlice(Int_t slice);
  void FindTrackCandidates();
  AliL3Histogram *AddHistograms(Int_t eta_index);
  void AddAllHistograms();
  void Evaluate(Bool_t remove=kFALSE);
  void WriteDigits(Char_t *outfile="output_digits.root");
  
  //Setters
  void SetNEtaSegments(Int_t i) {fNEtaSegments = i;}
  void SetAddHistograms() {fAddHistograms = kTRUE;}
  void SetRemoveFoundTracks() {fRemoveFoundTracks = kTRUE;}
  void SetWriteDigits() {fWriteDigits = kTRUE;}
  
  //Getters
  AliL3HoughTransformer *GetTransformer(Int_t i) {if(!fHoughTransformer[i]) return 0; return fHoughTransformer[i];}
  AliL3TrackArray *GetTracks() {if(!fTracks) return 0; return fTracks;}

  ClassDef(AliL3Hough,1)

};

#endif
