// @(#) $Id$

#ifndef ALIL3_KALMAN
#define ALIL3_KALMAN


class AliL3SpacePointData;
class AliL3TrackArray;
class AliL3Benchmark;
class AliL3TrackSegmentData;
class AliL3KalmanTrack;
class AliL3Track;

class AliL3Kalman {

 private:

  Int_t fMinSlice;
  Int_t fMaxSlice;
  AliL3SpacePointData *fClusters[36][6];
  Char_t fPath[1024];
  UInt_t fNcl[36][6];
  AliL3TrackArray *fTracks;
  AliL3TrackArray *fKalmanTracks;
  AliL3TrackArray *fSeeds;

  AliL3Benchmark *fBenchmark;
  Int_t fMinPointsOnTrack;
  Int_t fRow[6][2];
  Char_t fWriteOutPath[256];
  Bool_t fWriteOut;
  Int_t fEvent;

 public:

  AliL3Kalman(Char_t *datapath, Int_t *slice, Int_t min_clusters);
  virtual ~AliL3Kalman();
  void Init();
  void LoadTracks(Int_t event, Bool_t sp);
  void ProcessTracks();
  Int_t MakeKalmanSeed(AliL3KalmanTrack *kalmantrack, AliL3Track *track);
  Int_t InitKalmanTrack(AliL3KalmanTrack *kalmantrack, AliL3Track *track);
  Int_t Propagate(AliL3KalmanTrack *kalmantrack, AliL3Track *track);
  Int_t Update(AliL3SpacePointData *points, UInt_t pos, AliL3KalmanTrack *kalmantrack);
  void WriteFiles(Char_t *path="data"){fWriteOut = kTRUE; sprintf(fWriteOutPath,"%s",path);}
  Double_t GetCpuTime();
  AliL3TrackArray *GetTracks() {return fKalmanTracks;}
};

#endif
