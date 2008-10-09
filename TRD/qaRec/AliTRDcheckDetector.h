#ifndef __ALITRDCHECKDETECTOR_H__
#define __ALITRDCHECKDETECTOR_H__

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TMap;
class AliESDHeader;
class AliTRDcheckDetector : public AliTRDrecoTask{
// common constants
enum{
  kNDetectors = 540,
  kNDetectorsSector = 30,
  kNSectors = 18,
  kNLayers = 6,
  kNTimebins = 30
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
  
  virtual Bool_t PostProcess();
  virtual void  GetRefFigure(Int_t ifig, Int_t &first, Int_t &last, Option_t *opt = "lp");
  
private:
  AliTRDcheckDetector(const AliTRDcheckDetector &);
  AliTRDcheckDetector& operator=(const AliTRDcheckDetector &);
  AliTRDeventInfo *fEventInfo;						//! ESD Header
  TMap *fTriggerNames;										//! Containing trigger class names
  ClassDef(AliTRDcheckDetector, 1)
};
#endif

