// @(#) $Id$

#ifndef AliL3_DataCompressor
#define AliL3_DataCompressor

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
  Char_t fPath[1024];            //!
  
  UInt_t fNcl[36][6];
  Bool_t fKeepRemaining;
  Bool_t fWriteClusterShape;
  Int_t fEvent;
  Bool_t fSinglePatch;
  
  static Int_t fNumPadBits;
  static Int_t fNumTimeBits;
  static Int_t fNumChargeBits;
  static Int_t fNumShapeBits;
  
  static Float_t fXYResidualStep;
  static Float_t fZResidualStep;
  static Float_t fXYWidthStep;
  static Float_t fZWidthStep;
  static Int_t fClusterCharge;
  
  void SelectRemainingClusters();
  void ExpandTrackData(AliL3TrackArray *tracks);
  void ReadUncompressedData(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints);
  void ReadRemaining(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints);
  void QSort(TempCluster **a, Int_t first, Int_t last);
  Int_t Compare(TempCluster *a,TempCluster *b);
  void OpenOutputFile();
  void CloseOutputFile();

 public:
  AliL3DataCompressor();
  AliL3DataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape);
  virtual ~AliL3DataCompressor();
  
  void LoadData(Int_t event,Bool_t sp=kTRUE);
  void FillData(Int_t minhits,Bool_t expand);
  void LoadOfflineData(Int_t event);
  void CompressAndExpand();
  void WriteRemaining(Bool_t select);
  void RestoreData();
  void DoBench(Char_t *fname="benchmark");
  
  void SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape);
  void SetResolutions(Float_t xyresidual,Float_t zresidual,Int_t clustercharge,Float_t xywidth=0.005,Float_t zwidth=0.005);
  
  static const Int_t GetNPadBits() {return fNumPadBits;}
  static const Int_t GetNTimeBits() {return fNumTimeBits;}
  static const Int_t GetNChargeBits() {return fNumChargeBits;}
  static const Int_t GetNShapeBits() {return fNumShapeBits;}
  static const Float_t GetXYResidualStep() {return fXYResidualStep;}
  static const Float_t GetZResidualStep() {return fZResidualStep;}
  static const Float_t GetXYWidthStep() {return fXYWidthStep;}
  static const Float_t GetZWidthStep() {return fZWidthStep;}
  static const Int_t GetClusterCharge() {return fClusterCharge;}

  ClassDef(AliL3DataCompressor,1) 

};

#endif
