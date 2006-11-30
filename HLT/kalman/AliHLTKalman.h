// @(#) $Id$

#ifndef ALIL3_KALMAN
#define ALIL3_KALMAN


class AliHLTSpacePointData;
class AliHLTTrackArray;
class AliHLTBenchmark;
class AliHLTTrackSegmentData;
class AliHLTKalmanTrack;
class AliHLTTrack;

class AliHLTKalman {

 private:

  Int_t fMinSlice;
  Int_t fMaxSlice;
  AliHLTSpacePointData *fClusters[36][6];
  Char_t fPath[1024];
  UInt_t fNcl[36][6];
  AliHLTTrackArray *fTracks;
  AliHLTTrackArray *fKalmanTracks;
  AliHLTTrackArray *fSeeds;

  AliHLTBenchmark *fBenchmark;
  Int_t fMinPointsOnTrack;
  Int_t fRow[6][2];
  Char_t fWriteOutPath[256];
  Bool_t fWriteOut;
  Int_t fEvent;

 public:

  AliHLTKalman(Char_t *datapath, Int_t *slice, Int_t min_clusters);
  virtual ~AliHLTKalman();
  void Init();
  void LoadTracks(Int_t event, Bool_t sp);
  void ProcessTracks();
  Int_t MakeKalmanSeed(AliHLTKalmanTrack *kalmantrack, AliHLTTrack *track);
  Int_t InitKalmanTrack(AliHLTKalmanTrack *kalmantrack, AliHLTTrack *track);
  Int_t Propagate(AliHLTKalmanTrack *kalmantrack, AliHLTTrack *track);
  Int_t Update(AliHLTSpacePointData *points, UInt_t pos, AliHLTKalmanTrack *kalmantrack);
  void WriteFiles(Char_t *path="data"){fWriteOut = kTRUE; sprintf(fWriteOutPath,"%s",path);}
  Double_t GetCpuTime();
  AliHLTTrackArray *GetTracks() {return fKalmanTracks;}
};

typedef AliHLTKalman AliL3Kalman; // for backward compatibility

#endif
