#ifndef __ALITRDCHECKDETECTOR_H__
#define __ALITRDCHECKDETECTOR_H__

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
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
  kNclustersHist=1,
  kNtrackletsHist=2,
  kNclusterTrackletHist = 3,
  kChi2 = 4, 
  kChi2Normalized = 5,
  kNTracksSectorHist = 6,
  kPulseHeight = 7,
  kClusterCharge = 8,
  kChargeDeposit = 9
};
public:
  AliTRDcheckDetector();
  virtual ~AliTRDcheckDetector();
  
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *);
  virtual void Terminate(Option_t *);
  
  virtual Bool_t PostProcess();
  virtual void  GetRefFigure(Int_t ifig, Int_t &first, Int_t &last);
  
private:
  AliTRDcheckDetector(const AliTRDcheckDetector &);
  AliTRDcheckDetector& operator=(const AliTRDcheckDetector &);
  TObjArray *fPHSdetector;						//! PHS container for single Detectors
  TObjArray *fPHSsector;							//! PHS container for whole sector
  TObjArray *fQCLdetector;						//!	Cluster Charge Container for single detector
  TObjArray *fQCLsector;							//! Cluster Charge Container for whole sector
  TObjArray *fQTdetector;							//! Total charge Deposit for single detector 
  TObjArray *fQTsector;								//! Total charge Deposit for whole sector
  
  ClassDef(AliTRDcheckDetector, 1)
};
#endif

