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
  
  UInt_t fNcl[36][6];
  
  static Int_t fNumPadBits;
  static Int_t fNumTimeBits;
  static Int_t fNumChargeBits;
  static Int_t fNumPadShapeBits;
  static Int_t fNumTimeShapeBits;
  
  static Float_t fPadResidualStep1;
  static Float_t fPadResidualStep2;
  static Float_t fPadResidualStep3;
  static Float_t fTimeResidualStep1;
  static Float_t fTimeResidualStep2;
  static Float_t fTimeResidualStep3;
  static Float_t fPadSigma2Step1;
  static Float_t fPadSigma2Step2;
  static Float_t fTimeSigma2Step;
  static Int_t fClusterCharge;
  
  void SelectRemainingClusters();
  void ExpandTrackData(AliL3TrackArray *tracks);
  void ReadUncompressedData(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints);
  void ReadRemaining(TempCluster **clusters,Int_t *ncl,const Int_t maxpoints);
  void QSort(TempCluster **a, Int_t first, Int_t last);
  Int_t Compare(TempCluster *a,TempCluster *b);
  void OpenOutputFile();
  void CloseOutputFile();
  
 protected:
  Char_t fPath[1024];            //!
  Int_t fEvent;
  Int_t fNusedClusters;
  Int_t fNunusedClusters;
  
  Bool_t fWriteClusterShape;
  Bool_t fKeepRemaining;
  Bool_t fSinglePatch;
  Bool_t fWriteIdsToFile;
  
 public:
  AliL3DataCompressor();
  AliL3DataCompressor(Char_t *path,Bool_t keep,Bool_t writeshape);
  virtual ~AliL3DataCompressor();
  
  virtual void LoadData(Int_t event,Bool_t sp=kTRUE);
  virtual void FillData(Int_t minhits,Bool_t expand);
  virtual void WriteRemaining(Bool_t select);
  void CompressAndExpand();
  void RestoreData(Bool_t remaining_only=kFALSE);
  void DoBench(Char_t *fname="benchmark");
  
  void SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shapepad,Int_t shapetime);
  void SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width=0.005);
  void SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width=0.005);
  
  Int_t GetNusedClusters() {return fNusedClusters;}
  Int_t GetNunusedClusters() {return fNunusedClusters;}

  static const Int_t GetNPadBits() {return fNumPadBits;}
  static const Int_t GetNTimeBits() {return fNumTimeBits;}
  static const Int_t GetNChargeBits() {return fNumChargeBits;}
  static const Int_t GetNPadShapeBits() {return fNumPadShapeBits;}
  static const Int_t GetNTimeShapeBits() {return fNumTimeShapeBits;}
  static const Float_t GetPadSigma2Step(Int_t patch) {return patch < 2 ? fPadSigma2Step1 : fPadSigma2Step2;}
  static const Float_t GetTimeSigma2Step() {return fTimeSigma2Step;}
  static const Int_t GetClusterCharge() {return fClusterCharge;}
  static const Float_t GetPadResidualStep(Int_t row);
  static const Float_t GetTimeResidualStep(Int_t row);


  ClassDef(AliL3DataCompressor,1) 

};

#endif
