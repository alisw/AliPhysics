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
class AliL3HoughMerger;
class AliL3HoughIntMerger;

class AliL3Hough {
  
 private:
  Char_t fPath[256];
  Bool_t fBinary;
  Bool_t fAddHistograms;
  Bool_t fDoIterative;
  Bool_t fWriteDigits;
  Int_t fNEtaSegments;
  AliL3FileHandler **fMemHandler; //!
  AliL3HoughTransformer **fHoughTransformer; //!
  AliL3HoughEval **fEval; //!
  AliL3HoughMaxFinder *fPeakFinder; //!
  AliL3TrackArray **fTracks; //!
  AliL3HoughMerger *fMerger; //!
  AliL3HoughIntMerger *fInterMerger; //!
  
  void CleanUp();
  void Init();
  
 public:
  
  AliL3Hough(); 
  AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments=100);
  virtual ~AliL3Hough();
  
  void Process(Int_t minslice,Int_t maxslice);
  void ReadData(Int_t slice);
  void Transform();
  void ProcessSliceIter();
  void ProcessPatchIter(Int_t patch);
  void MergePatches();
  void MergeInternally();

  void FindTrackCandidates();
  AliL3Histogram *AddHistograms(Int_t eta_index);
  void AddAllHistograms();
  void Evaluate(Int_t road_width=1);
  void EvaluateWithEta();
  void WriteDigits(Char_t *outfile="output_digits.root");
  
  //Setters
  void SetNEtaSegments(Int_t i) {fNEtaSegments = i;}
  void SetAddHistograms() {fAddHistograms = kTRUE;}
  void DoIterative() {fDoIterative = kTRUE;}
  void SetWriteDigits() {fWriteDigits = kTRUE;}
  
  //Getters
  AliL3HoughTransformer *GetTransformer(Int_t i) {if(!fHoughTransformer[i]) return 0; return fHoughTransformer[i];}
  AliL3TrackArray *GetTracks(Int_t i) {if(!fTracks[i]) return 0; return fTracks[i];}
  AliL3HoughEval *GetEval(Int_t i) {if(!fEval[i]) return 0; return fEval[i];}
  AliL3HoughMerger *GetMerger() {if(!fMerger) return 0; return fMerger;}
  AliL3HoughIntMerger *GetInterMerger() {if(!fInterMerger) return 0; return fInterMerger;}
  AliL3FileHandler *GetMemHandler(Int_t i) {if(!fMemHandler[i]) return 0; return fMemHandler[i];}

  ClassDef(AliL3Hough,1) //Hough transform base class

};

#endif
