// @(#) $Id$

#ifndef AliHLT_DataCompressor
#define AliHLT_DataCompressor

#include "AliHLTRootTypes.h"

class AliHLTSpacePointData;
class AliHLTBenchmark;
class AliHLTTrackArray;
class AliHLTTrack;

#ifdef use_root
class TH2F;
class TFile;
#endif

struct TempCluster {
  Float_t fPad; // Pad
  Float_t fTime; // Time
  Float_t fSigmaY2; // SigmaY2
  Float_t fSigmaZ2; // SigmaZ2
  Int_t fCharge; // Charge
  Int_t fPadrow; // Pad row
};

class AliHLTDataCompressor {
  
 public:
  AliHLTDataCompressor();
  AliHLTDataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape);
  virtual ~AliHLTDataCompressor();
  
  virtual void LoadData(Int_t event,Bool_t sp=kTRUE);
  virtual void FillData(Int_t minhits,Bool_t expand);
  virtual void WriteRemaining(Bool_t select);
  void DetermineMinBits();
  void CompressAndExpand(Bool_t arithmeticCoding=kTRUE);
  void RestoreData(Bool_t remainingOnly=kFALSE);
  void DoBench(Char_t *fname="benchmark");
  void DoNotCompress() {fNoCompression=kTRUE;}

  Int_t GetNusedClusters() const {return fNusedClusters;}
  Int_t GetNunusedClusters() const {return fNunusedClusters;}

 protected:
  Char_t fPath[1024];   //! Path to the files
  Int_t fEvent; // Current event
  Int_t fNusedClusters; // Number of used clusters
  Int_t fNunusedClusters; // Number of unused clusters
  
  Bool_t fWriteClusterShape; // Flag to write the cluster's shape
  Bool_t fKeepRemaining; // Flag to keep the remaining clusters
  Bool_t fSinglePatch; // Flag to run over single patch (?)
  Bool_t fWriteIdsToFile; // Flag (not used?)
  Bool_t fNoCompression; //Just process the data through the chain, but do not compress. (input=output). Mostly for debugging...

 private:
  AliHLTBenchmark *fBenchmark;    //! Benchmark
  AliHLTTrackArray *fInputTracks; //! Array of input tracks
  AliHLTSpacePointData *fClusters[36][6]; //! Array of pointers to clusters
  ofstream *fCompRatioFile;      //! Stream to write the ration between use and unused clusters
#ifdef use_root
  TFile *fOutputFile;            //! Output file
#else
  FILE *fOutputFile; // Output file
#endif
  UInt_t fNcl[36][6]; // Array containing numbers of clusters
   
  void SelectRemainingClusters();
  void ExpandTrackData(AliHLTTrackArray *tracks);
  void ReadUncompressedData(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints);
  void ReadRemaining(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints);
  void QSort(TempCluster **a, Int_t first, Int_t last);
  Int_t Compare(TempCluster *a,TempCluster *b);
  void OpenOutputFile();
  void CloseOutputFile();
  
  ClassDef(AliHLTDataCompressor,1) 

};

typedef AliHLTDataCompressor AliL3DataCompressor; // for backward compatibility

#endif
