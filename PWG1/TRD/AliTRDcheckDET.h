#ifndef ALITRDCHECKDET_H
#define ALITRDCHECKDET_H

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif
////////////////////////////////////////////////////////////////////////////
//  Basic checks for tracking and detector performance                    //
//                                                                        //
//  Authors:                                                              //
//    Anton Andronic <A.Andronic@gsi.de>                                  //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
////////////////////////////////////////////////////////////////////////////


class TObjArray;
class TH1;
class TMap;
class AliESDHeader;
class AliExternalTrackParam;
class AliTRDcluster;
class AliTRDseedV1;
class AliTRDgeometry;
class AliTRDcheckDET : public AliTRDrecoTask{
public:
  // The Histogram number
  enum  HistType_t {
    kNclustersTrack     = 0,
    kNclustersTracklet  = 1,
    kNtrackletsTrack    = 2,
    kNtrackletsSTA      = 3,
    kNtrackletsBAR      = 4,
    kNtrackletsCross    = 5,
    kNtrackletsFindable = 6,
    kNtracksEvent       = 7,
    kNtracksSector      = 8,
    kPH                 = 9,
    kChi2               = 10,
    kChargeCluster      = 11,
    kChargeTracklet     = 12,
    kNeventsTrigger     = 13,
    kNeventsTriggerTracks=14,
    kTriggerPurity      = 15,
    kTrackStatus        = 16,
    kTrackletStatus     = 17,
    kNTrackletsP        = 18,
    kNclustersLayer     = 19
  };
  enum FigureType_t{
    kFigNclustersTrack,
    kFigNclustersTracklet,
    kFigNtrackletsTrack,
    kFigNTrackletsP,
    kFigNtrackletsCross,
    kFigNtrackletsFindable,
    kFigNtracksEvent,
    kFigNtracksSector,
    kFigTrackStatus,
    kFigTrackletStatus,
    kFigChi2,
    kFigPH,
    kFigChargeCluster,
    kFigChargeTracklet,
    kFigNeventsTrigger,
    kFigNeventsTriggerTracks,
    kFigTriggerPurity
  };
 
  AliTRDcheckDET();
  AliTRDcheckDET(char* name);
  virtual ~AliTRDcheckDET();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  virtual TObjArray *Histos();

  // Plotting Functions:
  TH1 *PlotTrackStatus(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotTrackletStatus(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNClustersTracklet(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotNClustersTrack(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotNTrackletsTrack(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotNTrackletsRowCross(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotFindableTracklets(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotNTracksSector(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotPHt(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotPHx(const AliTRDtrackV1 *track = 0x0);
  TH1 *PlotChi2(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotChargeCluster(const AliTRDtrackV1 *t = 0x0);
  TH1 *PlotChargeTracklet(const AliTRDtrackV1 *t = 0x0);

  virtual Bool_t PostProcess();
  virtual Bool_t GetRefFigure(Int_t ifig);
  virtual void MakeSummary();
  
  Bool_t IsUsingClustersOutsideChamber() const {return TESTBIT(fFlags, kUseClustersOutsideChamber);}
  void UseClustersOutsideChamber(Bool_t b = kTRUE) {if(b) SETBIT(fFlags, kUseClustersOutsideChamber); else CLRBIT(fFlags, kUseClustersOutsideChamber);}

private:
  enum{
    kUseClustersOutsideChamber = 0 
  };
  AliTRDcheckDET(const AliTRDcheckDET &);
  AliTRDcheckDET& operator=(const AliTRDcheckDET &);
  void GetDistanceToTracklet(Double_t *dist, AliTRDseedV1 * const tracklet, AliTRDcluster * const c);
  //----------------------------------------------------
  // Functions creating the reference figures
  TH1* MakePlotChi2() const;
  TH1* MakePlotNTracklets();
  Bool_t MakePlotPulseHeight();
  void MakePlotnTrackletsVsP();
  void MakePlotMeanClustersLayer();
  void MakePlotNclustersTrack();
  void MakePlotNclustersTracklet();
  void MakePlotTrackletCharge();
  Bool_t MakeBarPlot(TH1 *histo, Int_t Color);
  //----------------------------------------------------
  void GetEtaPhiAt(const AliExternalTrackParam *track, Double_t x, Double_t &eta, Double_t &phi);
  TH1 *ProjectCentrality(TH2 *h2d, Int_t centralityBin = -1);

  Int_t fCentralityClass;              // Centrality Class
  TMap *fTriggerNames;                 //! Containing trigger class names
  UChar_t fFlags;                      // Flags for setting

  static const Color_t fkColorsCentrality[AliTRDeventInfo::kCentralityClasses];		// Colors for the different centrality classes in the Ref Figures
  ClassDef(AliTRDcheckDET, 2)
};
#endif
