// @(#) $Id$

#ifndef ALIL3HOUGH_H
#define ALIL3HOUGH_H

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
class AliL3Benchmark;

#include "TThread.h"
#ifdef use_newio
#include <AliRunLoader.h>
#include <../RAW/AliRawEvent.h>
#endif
#ifdef use_aliroot
#include <AliESD.h>
#include <AliESDHLTtrack.h>
#endif

class AliL3Hough {
 public:
  
  AliL3Hough(); 
  AliL3Hough(Char_t *path,Bool_t binary,Int_t netasegments=100,Bool_t bit8=kFALSE,Int_t tv=0,Char_t *infile=0,Char_t *ptr=0);
  virtual ~AliL3Hough();

#ifdef use_newio  
  void SetRunLoader(AliRunLoader *runloader) {fRunLoader = runloader;}
#endif

  void Init(Int_t netasegments,Int_t tv,AliRawEvent *rawevent,Float_t zvertex=0.0);
  void Init(Char_t *path,Bool_t binary,Int_t netasegments=100,Bool_t bit8=kFALSE,Int_t tv=0,Char_t *infile=0,Char_t *ptr=0,Float_t zvertex=0.0);
  void Init(Bool_t doit=kFALSE, Bool_t addhists=kFALSE);

  void Process(Int_t minslice,Int_t maxslice);
  void ReadData(Int_t slice,Int_t eventnr=0);
  void Transform(Int_t *rowrange = 0);
  void ProcessSliceIter();
  void ProcessPatchIter(Int_t patch);
  void MergePatches();
  void MergeInternally();
  void MergeEtaSlices();

  void FindTrackCandidates();
  void FindTrackCandidatesRow();
  void AddAllHistograms();
  void AddAllHistogramsRows();
  void PrepareForNextPatch(Int_t nextpatch);
  Int_t Evaluate(Int_t roadwidth=1,Int_t nrowstomiss=1);
  void EvaluatePatch(Int_t i,Int_t roadwidth,Int_t nrowstomiss);
  void WriteTracks(Int_t slice,Char_t *path="./");
  void WriteTracks(Char_t *path);
#ifdef use_aliroot
  Int_t FillESD(AliESD *esd);
#endif
  void WriteDigits(Char_t *outfile="output_digits.root");
  void InitEvaluate();
  void DoBench(Char_t *filename);
  void AddTracks();
  
  //Setters
  void SetNEtaSegments(Int_t i) {fNEtaSegments = i;}
  void SetAddHistograms() {fAddHistograms = kTRUE;}
  void DoIterative() {fDoIterative = kTRUE;}
  void SetWriteDigits() {fWriteDigits = kTRUE;}
  void SetTransformerParams(Float_t ptres=0,Float_t ptmin=0,Float_t ptmax=0,Int_t ny=0,Int_t patch=-1);
  //{fPtRes=ptres;fNBinY=ny;fLowPt=ptmin;fUpperPt=ptmax;fPhi=psi;}
  void SetTransformerParams(Int_t nx,Int_t ny,Float_t lpt,Int_t patch);
  void SetTransformerParams(Int_t nx,Int_t ny,Float_t lpt,Float_t phi);
  //{fNBinX=nx;fNBinY=ny;fLowPt=lpt;fPhi=phi;}
  void SetThreshold(Int_t t=3,Int_t patch=-1);
  void SetNSaveIterations(Int_t t=10) {fNSaveIterations=t;}
  void SetPeakThreshold(Int_t threshold=0,Int_t patch=-1);
  
  void SetPeakParameters(Int_t kspread,Float_t pratio) {fKappaSpread=kspread; fPeakRatio=pratio;}
  
  //Getters
  AliL3HoughBaseTransformer *GetTransformer(Int_t i) {if(!fHoughTransformer[i]) return 0; return fHoughTransformer[i];}
  AliL3TrackArray *GetTracks(Int_t i) {if(!fTracks[i]) return 0; return fTracks[i];}
  AliL3HoughEval *GetEval(Int_t i) {if(!fEval[i]) return 0; return fEval[i];}
  AliL3HoughMerger *GetMerger() {if(!fMerger) return 0; return fMerger;}
  AliL3HoughIntMerger *GetInterMerger() {if(!fInterMerger) return 0; return fInterMerger;}
  AliL3MemHandler *GetMemHandler(Int_t i) {if(!fMemHandler[i]) return 0; return fMemHandler[i];}
  AliL3HoughMaxFinder *GetMaxFinder() {return fPeakFinder;}

  //Special methods for executing Hough Transform as a thread
  static void *ProcessInThread(void *args);
  void StartProcessInThread(Int_t minslice,Int_t maxslice);
  Int_t WaitForThreadFinish();
  void SetMinMaxSlices(Int_t minslice,Int_t maxslice) {fMinSlice = minslice; fMaxSlice = maxslice;} 
  Int_t GetMinSlice() {return fMinSlice;}
  Int_t GetMaxSlice() {return fMaxSlice;}
  
 private:
  Char_t *fInputFile;//!
  Char_t *fInputPtr;//!
  AliRawEvent *fRawEvent;//!
  Char_t fPath[1024]; // Path to the files
  Bool_t fBinary; // Is input binary
  Bool_t fAddHistograms; // Add all patch histograms at the end or not
  Bool_t fDoIterative; // Iterative or not
  Bool_t fWriteDigits; // Write Digits or not
  Bool_t fUse8bits; // Use 8 bits or not
  Int_t fNEtaSegments; // Number of eta slices
  Int_t fNPatches; // Number of patches
  Int_t fLastPatch; //The index of the last processed patch
  Int_t fVersion; //which HoughTransformer to use
  Int_t fCurrentSlice; // Current TPC slice (sector)
  Int_t fEvent; // Current event number

  Int_t fPeakThreshold[6]; // Threshold for the peak finder
  Float_t fLowPt[6]; // Lower limit on Pt
  Float_t fUpperPt[6]; // Upper limit on Pt
  Float_t fPtRes[6]; // Desired Pt resolution
  Float_t fPhi[6]; // Limit on the emission angle
  Int_t fNBinX[6]; // Number of bins in the Hough space
  Int_t fNBinY[6]; // Number of bins in the Hough space
  Int_t fThreshold[6]; // Threshold for digits
  Int_t fNSaveIterations; //for HoughtransformerVhdl
  
  //parameters for the peak finder:
  Int_t fKappaSpread; // Kappa spread
  Float_t fPeakRatio; // Peak ratio

  Float_t fZVertex; // Z position of the primary vertex

  Int_t fMinSlice;
  Int_t fMaxSlice;

  AliL3MemHandler **fMemHandler; //!
  AliL3HoughBaseTransformer **fHoughTransformer; //!
  AliL3HoughEval **fEval; //!
  AliL3HoughMaxFinder *fPeakFinder; //!
  AliL3TrackArray **fTracks; //!
  AliL3TrackArray *fGlobalTracks; //!
  AliL3HoughMerger *fMerger; //!
  AliL3HoughIntMerger *fInterMerger; //!
  AliL3HoughGlobalMerger *fGlobalMerger; //!
  AliL3Benchmark *fBenchmark; //!

#ifdef use_newio
  AliRunLoader *fRunLoader; // Run Loader
#endif

  void CleanUp();
  Double_t GetCpuTime();

  TThread *fThread;

  ClassDef(AliL3Hough,1) //Hough transform base class
};

#endif







