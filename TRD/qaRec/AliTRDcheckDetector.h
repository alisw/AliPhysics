#ifndef __ALITRDCHECKDETECTOR_H__
#define __ALITRDCHECKDETECTOR_H__

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TH1;
class TMap;
class AliESDHeader;
class AliTRDcheckDetector : public AliTRDrecoTask{
// common constants
enum{
  kNTimeBins = 30
};
// The Histogram number
enum{
  kNTracksEventHist=0,
  kNEventsTriggerTracks=1,
  kNclustersHist=2,
  kNtrackletsHist=3,
  kNclusterTrackletHist=4,
  kChi2=5, 
  kChi2Normalized=6,
  kNTracksSectorHist=7,
  kPulseHeight=8,
  kClusterCharge=9,
  kChargeDeposit=10,
  kNEventsTrigger=11,
  kPurity = 12
};
public:
  AliTRDcheckDetector();
  virtual ~AliTRDcheckDetector();
  
  virtual void ConnectInputData(const Option_t *);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *);
  virtual void Terminate(Option_t *);
  
  virtual TObjArray *Histos();
  
  // Plotting Functions:
  TH1 *PlotMeanNClusters(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotNClusters(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotNTracklets(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotTracksSector(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotPulseHeight(const AliTRDtrackV1 *t = 0x0);
	TH1 *PlotChi2(const AliTRDtrackV1 *t = 0x0);
	TH1 *PlotNormalizedChi2(const AliTRDtrackV1 *t = 0x0);
	TH1 *PlotClusterCharge(const AliTRDtrackV1 *t = 0x0);
	TH1 *PlotChargeDeposit(const AliTRDtrackV1 *t = 0x0);

  virtual Bool_t PostProcess();
  virtual void  GetRefFigure(Int_t ifig);
  
private:
  AliTRDcheckDetector(const AliTRDcheckDetector &);
  AliTRDcheckDetector& operator=(const AliTRDcheckDetector &);
  AliTRDeventInfo *fEventInfo;						//! ESD Header
  TMap *fTriggerNames;										//! Containing trigger class names
  ClassDef(AliTRDcheckDetector, 1)
};
#endif

