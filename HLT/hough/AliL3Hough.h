#ifndef ALIL3_HOUGH
#define ALIL3_HOUGH

#include "AliL3RootTypes.h"

class AliL3HoughMaxFinder;
class AliL3HoughBaseTransformer;
class AliL3Histogram;
class AliL3MemHandler;
class AliL3FileHandler;
class AliL3HoughEval;
class AliL3TrackArray;
class AliL3HoughMerger;
class AliL3HoughIntMerger;
class AliL3HoughGlobalMerger;

class AliL3Hough {
  
 private:
  Char_t fPath[256];
  Bool_t fBinary;
  Bool_t fAddHistograms;
  Bool_t fDoIterative;
  Bool_t fWriteDigits;
  Bool_t fUse8bits;
  Int_t fNEtaSegments;
  Int_t fNPatches;
  Int_t fversion;
  Int_t fCurrentSlice;
  
  AliL3MemHandler **fMemHandler; //!
  AliL3HoughBaseTransformer **fHoughTransformer; //!
  AliL3HoughEval **fEval; //!
  AliL3HoughMaxFinder *fPeakFinder; //!
  AliL3TrackArray **fTracks; //!
  AliL3HoughMerger *fMerger; //!
  AliL3HoughIntMerger *fInterMerger; //!
  AliL3HoughGlobalMerger *fGlobalMerger; //!

  void CleanUp();
  Double_t GetCpuTime();
  
 public:
  
  AliL3Hough(); 
  AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments=100,Int_t tversion=0);
  virtual ~AliL3Hough();
  
  void Init(Char_t *path,Bool_t binary,Int_t n_eta_segments,Bool_t bit8=kFALSE,Int_t tv=0);
  void Process(Int_t minslice,Int_t maxslice);
  void ReadData(Int_t slice,Int_t eventnr=0);
  void Transform(Int_t row_range = -1);
  void ProcessSliceIter();
  void ProcessPatchIter(Int_t patch);
  void MergePatches();
  void MergeInternally();
  
  void FindTrackCandidates();
  void AddAllHistograms();
  Int_t Evaluate(Int_t road_width=1,Int_t nrowstomiss=1);
  void EvaluateWithEta();
  void WriteTracks(Int_t slice,Char_t *path="./");
  void WriteDigits(Char_t *outfile="output_digits.root");
  void InitEvaluate();
  
  //Setters
  void SetNEtaSegments(Int_t i) {fNEtaSegments = i;}
  void SetAddHistograms() {fAddHistograms = kTRUE;}
  void DoIterative() {fDoIterative = kTRUE;}
  void SetWriteDigits() {fWriteDigits = kTRUE;}
  
  //Getters
  AliL3HoughBaseTransformer *GetTransformer(Int_t i) {if(!fHoughTransformer[i]) return 0; return fHoughTransformer[i];}
  AliL3TrackArray *GetTracks(Int_t i) {if(!fTracks[i]) return 0; return fTracks[i];}
  AliL3HoughEval *GetEval(Int_t i) {if(!fEval[i]) return 0; return fEval[i];}
  AliL3HoughMerger *GetMerger() {if(!fMerger) return 0; return fMerger;}
  AliL3HoughIntMerger *GetInterMerger() {if(!fInterMerger) return 0; return fInterMerger;}
  AliL3MemHandler *GetMemHandler(Int_t i) {if(!fMemHandler[i]) return 0; return fMemHandler[i];}
  AliL3HoughMaxFinder *GetMaxFinder() {return fPeakFinder;}

  ClassDef(AliL3Hough,1) //Hough transform base class
};

#endif
