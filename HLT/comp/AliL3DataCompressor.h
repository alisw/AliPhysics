// @(#) $Id$

#ifndef AliL3_DataCompressor
#define AliL3_DataCompressor

#include "AliL3RootTypes.h"

class AliL3SpacePointData;
class AliL3Benchmark;
class AliL3TrackArray;
class AliL3Track;

#ifdef use_root
class TH2F;
class TFile;
#endif

struct TempCluster {
  Float_t pad;
  Float_t time;
  Float_t sigmaY2;
  Float_t sigmaZ2;
  Int_t charge;
  Int_t padrow;
};

class AliL3DataCompressor {
  
 private:
  AliL3Benchmark *fBenchmark;    //!
  AliL3TrackArray *fInputTracks; //!
  AliL3SpacePointData *fClusters[36][6]; //!
  ofstream *fCompRatioFile;      //!
#ifdef use_root
  TFile *fOutputFile;            //!
#else
  FILE *fOutputFile;
#endif
  UInt_t fNcl[36][6];
   
  void SelectRemainingClusters();
  void ExpandTrackData(AliL3TrackArray *tracks);
  void ReadUncompressedData(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints);
  void ReadRemaining(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints);
  void QSort(TempCluster **a, Int_t first, Int_t last);
  Int_t Compare(TempCluster *a,TempCluster *b);
  void OpenOutputFile();
  void CloseOutputFile();
  
 protected:
  Char_t fPath[1024];   //!
  Int_t fEvent;
  Int_t fNusedClusters;
  Int_t fNunusedClusters;
  
  Bool_t fWriteClusterShape;
  Bool_t fKeepRemaining;
  Bool_t fSinglePatch;
  Bool_t fWriteIdsToFile;
  Bool_t fNoCompression; //Just process the data through the chain, but do not compress. (input=output). Mostly for debugging...

 public:
  AliL3DataCompressor();
  AliL3DataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape);
  virtual ~AliL3DataCompressor();
  
  virtual void LoadData(Int_t event,Bool_t sp=kTRUE);
  virtual void FillData(Int_t minhits,Bool_t expand);
  virtual void WriteRemaining(Bool_t select);
  void DetermineMinBits();
  void CompressAndExpand(Bool_t arithmetic_coding=kTRUE);
  void RestoreData(Bool_t remaining_only=kFALSE);
  void DoBench(Char_t *fname="benchmark");
  void DoNotCompress() {fNoCompression=kTRUE;}

  Int_t GetNusedClusters() {return fNusedClusters;}
  Int_t GetNunusedClusters() {return fNunusedClusters;}

  ClassDef(AliL3DataCompressor,1) 

};

#endif
