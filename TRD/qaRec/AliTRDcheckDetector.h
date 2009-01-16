#ifndef ALITRDCHECKDETECTOR_H
#define ALITRDCHECKDETECTOR_H

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TH1;
class TMap;
class AliESDHeader;
class AliTRDcluster;
class AliTRDseedV1;
class AliTRDgeometry;
class AliTRDReconstructor;
class AliTRDrecoParam;
class AliTRDeventInfo;

class AliTRDcheckDetector : public AliTRDrecoTask{
public:
  // The Histogram number
  enum  HistType_t {
    kNclustersTrack     = 0,
    kNclustersTracklet   = 1,
    kNtrackletsTrack    = 2,
    kNtrackletsCross    = 3,
    kNtrackletsFindable = 4,
    kNtracksEvent       = 5,
    kNtracksSector      = 6,
    kPH                 = 7,
    kChi2               = 8,
    kChargeCluster      = 9,
    kChargeTracklet     = 10,
    kNeventsTrigger     = 11,
    kNeventsTriggerTracks=12,
    kTriggerPurity      = 13
  };

  AliTRDcheckDetector();
  virtual ~AliTRDcheckDetector();

  virtual void ConnectInputData(const Option_t *);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *);
  virtual void Terminate(Option_t *);

  virtual TObjArray *Histos();

  // Plotting Functions:
  TH1 *PlotNClustersTracklet(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotNClustersTrack(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotNTrackletsTrack(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotNTrackletsRowCross(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotFindableTracklets(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNTracksSector(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotPHt(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotPHx(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotChi2(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotChi2Norm(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotChargeCluster(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotChargeTracklet(const AliTRDtrackV1 *t = 0x0);

  virtual Bool_t PostProcess();
  virtual Bool_t GetRefFigure(Int_t ifig);
  
  void SetRecoParam(AliTRDrecoParam *r);

private:
  AliTRDcheckDetector(const AliTRDcheckDetector &);
  AliTRDcheckDetector& operator=(const AliTRDcheckDetector &);
  void GetDistanceToTracklet(Double_t *dist, AliTRDseedV1 *tracklet, AliTRDcluster *c);
  AliTRDeventInfo *fEventInfo;         //! ESD Header
  TMap *fTriggerNames;                 //! Containing trigger class names
  AliTRDReconstructor *fReconstructor; // TRD Reconstructor
  AliTRDgeometry *fGeo;                // TRD Geometry object
    
  ClassDef(AliTRDcheckDetector, 1)
};
#endif

