#ifndef ALITRDCHECKDETECTOR_H
#define ALITRDCHECKDETECTOR_H

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TH1;
class TMap;
class AliESDHeader;
class AliTRDgeometry;
class AliTRDReconstructor;
class AliTRDrecoParam;
class AliTRDeventInfo;
class AliTRDcheckDetector : public AliTRDrecoTask{
  // The Histogram number
  typedef enum{
    kNTracksEventHist=0,
    kNEventsTriggerTracks=1,
    kNclustersHist=2,
    kNtrackletsHist=3,
    kNTrackletsVsFindable = 4,
    kNclusterTrackletHist=5,
    kChi2=6, 
    kChi2Normalized=7,
    kNTracksSectorHist=8,
    kPulseHeight=9,
    kClusterCharge=10,
    kChargeDeposit=11,
    kNEventsTrigger=12,
    kPurity = 12
  }HistType_t;
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
    TH1 *PlotTrackletsVsFindable(const AliTRDtrackV1 *track = 0x0);
    TH1 *PlotTracksSector(const AliTRDtrackV1 *t = 0x0);
    TH1 *PlotPulseHeight(const AliTRDtrackV1 *t = 0x0);
	  TH1 *PlotChi2(const AliTRDtrackV1 *t = 0x0);
	  TH1 *PlotNormalizedChi2(const AliTRDtrackV1 *t = 0x0);
	  TH1 *PlotClusterCharge(const AliTRDtrackV1 *t = 0x0);
	  TH1 *PlotChargeDeposit(const AliTRDtrackV1 *t = 0x0);
  
    virtual Bool_t PostProcess();
    virtual void  GetRefFigure(Int_t ifig);
    
    void SetRecoParam(AliTRDrecoParam *r);
  
  private:
    AliTRDcheckDetector(const AliTRDcheckDetector &);
    AliTRDcheckDetector& operator=(const AliTRDcheckDetector &);
    Int_t GetNTracklets(const AliTRDtrackV1 *track);
    AliTRDeventInfo *fEventInfo;						//! ESD Header
    TMap *fTriggerNames;										//! Containing trigger class names
    AliTRDReconstructor *fReconstructor;    // TRD Reconstructor
    AliTRDgeometry *fGeo;                   // TRD Geometry object
    
  ClassDef(AliTRDcheckDetector, 1)
};
#endif

