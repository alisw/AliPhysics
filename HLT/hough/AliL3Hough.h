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
  Int_t fNEtaSegments;
  Int_t fNPatches;
  Int_t fPeakThreshold;
#ifdef use_aliroot
  AliL3FileHandler **fMemHandler; //!
#else
  AliL3MemHandler **fMemHandler; //!
#endif
  AliL3HoughBaseTransformer **fHoughTransformer; //!
  AliL3HoughEval **fEval; //!
  AliL3HoughMaxFinder *fPeakFinder; //!
  AliL3TrackArray **fTracks; //!
  AliL3HoughMerger *fMerger; //!
  AliL3HoughIntMerger *fInterMerger; //!
  AliL3HoughGlobalMerger *fGlobalMerger; //!

  void CleanUp();
  void Init();
  
 public:
  
  AliL3Hough(); 
  AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments=100);
  virtual ~AliL3Hough();
  
  void Process(Int_t minslice,Int_t maxslice);
  void ReadData(Int_t slice);
  void Transform(Int_t row_range = -1);
  void ProcessSliceIter();
  void ProcessPatchIter(Int_t patch);
  void MergePatches();
  void MergeInternally();

  void FindTrackCandidates();
  void AddAllHistograms();
  void Evaluate(Int_t road_width=1);
  void EvaluateWithEta();
  void WriteTracks(Char_t *path="./");
#ifdef use_aliroot
  void WriteDigits(Char_t *outfile="output_digits.root");
#endif

  //Setters
  void SetNEtaSegments(Int_t i) {fNEtaSegments = i;}
  void SetAddHistograms() {fAddHistograms = kTRUE;}
  void DoIterative() {fDoIterative = kTRUE;}
  void SetWriteDigits() {fWriteDigits = kTRUE;}
  void SetPeakThreshold(Int_t i) {fPeakThreshold = i;}
  
  //Getters
  AliL3HoughBaseTransformer *GetTransformer(Int_t i) {if(!fHoughTransformer[i]) return 0; return fHoughTransformer[i];}
  AliL3TrackArray *GetTracks(Int_t i) {if(!fTracks[i]) return 0; return fTracks[i];}
  AliL3HoughEval *GetEval(Int_t i) {if(!fEval[i]) return 0; return fEval[i];}
  AliL3HoughMerger *GetMerger() {if(!fMerger) return 0; return fMerger;}
  AliL3HoughIntMerger *GetInterMerger() {if(!fInterMerger) return 0; return fInterMerger;}
#ifdef use_aliroot
  AliL3FileHandler *GetMemHandler(Int_t i) {if(!fMemHandler[i]) return 0; return fMemHandler[i];}
#else
  AliL3MemHandler *GetMemHandler(Int_t i) {if(!fMemHandler[i]) return 0; return fMemHandler[i];}
#endif

  ClassDef(AliL3Hough,1) //Hough transform base class

};

#endif
