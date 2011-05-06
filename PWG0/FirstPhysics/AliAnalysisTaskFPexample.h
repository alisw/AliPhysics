
#ifndef ALIANALYSISTASKFPEXAMPLE_H
#define ALIANALYSISTASKFPEXAMPLE_H

class AliAnalysisHistosVertex;
#include "AliAnalysisTaskFirstPhysics.h"

class AliAnalysisTaskFPexample : public AliAnalysisTaskFirstPhysics {
 public:
  AliAnalysisTaskFPexample(const char *name = "You should have given a name to this analysis");
  virtual ~AliAnalysisTaskFPexample();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

 private:

  TH1D *fhTrackPt; // pT spectrum of tracks
  TH2D *fh2TrackPhiEta; // ESD tracks
  TH1D *fhMulITSTPC; // multiplicity distribution from TPC tracks and ITS tracks
  TH1D *fhMulITSSA; // multiplicity distribution from ITS standalone
  TH1D *fhMulSPD; // multiplicity distribution from SPD tracklets
  TH2D *fh2TrackletsPhiEta; // tracklet distribution
  TH2D *fh2TracksPhiTPCchi2; // track chi2 from TPC vs. azimuthal angle

  AliAnalysisTaskFPexample(const AliAnalysisTaskFPexample&); // not implemented
  AliAnalysisTaskFPexample& operator=(const AliAnalysisTaskFPexample&); // not implemented

  ClassDef(AliAnalysisTaskFPexample, 1);
};

#endif

