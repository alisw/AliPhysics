#ifndef ALITRDRESOLUTION_H
#define ALITRDRESOLUTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDresolution.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD Resolution performance                                            //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TH1;
class TH2;
class TF1;
class TGraphErrors;
class TObjArray;
class TDatabasePDG;
class AliTRDrecoParam;
class AliTRDseedV1;
class AliTRDtrackInfo;
class AliTrackPoint;
class AliTRDresolution : public AliTRDrecoTask
{
public:
  enum ETRDresolutionPlot {
     kCharge     =  0 // charge resolution
    ,kCluster    =  1 // cluster - track
    ,kTrack      =  2 // tracklet - track residuals/pulls
    ,kTrackIn    =  3 // tracklet - track residuals/pulls at lower TRD entrance 
    ,kTrackOut   =  4 // tracklet - track residuals/pulls at lower TRD entrance during refit
    ,kMCcluster  =  5 // cluster-mc resolution/pulls
    ,kMCtracklet =  6 // tracklet-mc resolution/pulls
    ,kMCtrackIn  =  7 // TPC track monitor
    ,kMCtrackOut =  8 // TOF/HMPID track monitor
    ,kMCtrack    =  9 // TRD track monitor
    ,kNviews     = 10 // total number of resolution views
    ,kNprojs     = 70 // total number of projections for all views
  };
  enum ETRDresolutionSteer {
     kVerbose    = BIT(18) // switch verbosity
    ,kVisual     = BIT(19) // show partial results during processing
    ,kTrackRefit = BIT(20) // steer track refit
    ,kTrackSelect= BIT(21) // steer track selection
    ,kLoadCorr   = BIT(23) // steer load corrections
  };
  enum ETRDresolutionOutSlots {
     kClToTrk    = 2
    ,kClToMC     = 3
    ,kTrkltToTrk = 4
    ,kTrkltToMC  = 5
    ,kNOutSlots  = 4
  };
  enum ETRDresolutionSegmentation {
     kSector    = 0
    ,kStack
    ,kDetector
  };

  AliTRDresolution();
  AliTRDresolution(char* name);
  virtual ~AliTRDresolution();
  
  static Bool_t  FitTrack(const Int_t np, AliTrackPoint *points, Float_t params[10]);
  static Bool_t  FitTracklet(const Int_t ly, const Int_t np, const AliTrackPoint *points, const Float_t trackPars[10], Float_t trackletPars[3]);
  void    UserCreateOutputObjects();
  Float_t GetCorrectionX(Int_t det, Int_t tb) const {return fXcorr[det][tb];}
  Float_t GetDyRange() const {return fDyRange;}
  Float_t GetPtThreshold() const {return fPtThreshold;}
  Float_t GetSegmentationLevel() const {return fSegmentLevel;}
  Bool_t  GetRefFigure(Int_t ifig);
  TObjArray*  Histos(); 
  Bool_t  Load(const Char_t *file = "AnalysisResults.root", const Char_t *dir="TRD_Performance");
  Bool_t  LoadCorrection(const Char_t *file=NULL);
  void    MakeSummary();

  TObjArray*  Results(Int_t i=0) const {return i ? fGraphS : fGraphM;} 
  void    UserExec(Option_t * opt);
  void    InitExchangeContainers();
  Bool_t  HasLoadCorrection() const {return TestBit(kLoadCorr);}
  Bool_t  HasTrackRefit() const {return TestBit(kTrackRefit);}
  Bool_t  HasTrackSelection() const {return TestBit(kTrackSelect);}
  Bool_t  IsVerbose() const {return TestBit(kVerbose);}
  Bool_t  IsVisual() const {return TestBit(kVisual);}
  Bool_t  PostProcess();

  TH1*    PlotCharge(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotCluster(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotTracklet(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotTrackIn(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotTrackOut(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotMC(const AliTRDtrackV1 *t=NULL);

  static Bool_t  Process(TH2* const h2, TGraphErrors **g, Int_t stat=100);
  void    SetDyRange(Float_t dy) {fDyRange = dy;}
  void    SetSegmentationLevel(Int_t l=0);
  void    SetPtThreshold(Float_t pt) {fPtThreshold = pt;}
  void    SetLoadCorrection(Bool_t v = kTRUE) {SetBit(kLoadCorr, v);}
  void    SetVerbose(Bool_t v = kTRUE) {SetBit(kVerbose, v);}
  void    SetVisual(Bool_t v = kTRUE) {SetBit(kVisual, v);}
  void    SetTrackRefit(Bool_t v = kTRUE) {SetBit(kTrackRefit, v);}
  void    SetTrackSelection(Bool_t v = kTRUE) {SetBit(kTrackSelect, v);}

  void    Terminate(Option_t * opt);
  Bool_t  GetGraph(Float_t *bb, ETRDresolutionPlot ip, Int_t idx=-1, Bool_t kLEG=kTRUE, const Char_t *explain=NULL);
  Bool_t  GetGraphArray(Float_t *bb, ETRDresolutionPlot ip, Int_t idx, Bool_t kLEG=kTRUE, Int_t n=0, Int_t *sel=NULL, const Char_t *explain=NULL);
  static Bool_t  UseTrack(const Int_t np, const AliTrackPoint *points, Float_t params[10]);
private:
  AliTRDresolution(const AliTRDresolution&);
  AliTRDresolution& operator=(const AliTRDresolution&);

  void    AdjustF1(TH1 *h, TF1 *f);
  TObjArray*  BuildMonitorContainerCluster(const char* name, Bool_t expand=kFALSE, Float_t range=-1.);
  TObjArray*  BuildMonitorContainerTracklet(const char* name, Bool_t expand=kFALSE);
  TObjArray*  BuildMonitorContainerTrack(const char* name);
  void    GetLandauMpvFwhm(TF1 * const f, Float_t &mpv, Float_t &xm, Float_t &xM);
  void    GetRange(TH2 *h2, Char_t mod, Float_t *range);
  void    MakeSummaryPlot(TObjArray *a, TH2 *h2);
  Bool_t  Process(TH2* const h2, TF1 *f, Float_t k, TGraphErrors **g);
  Bool_t  Process2D(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=NULL,  Float_t scale=1., Int_t gidx=-1);
  Bool_t  Process2Darray(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=NULL,  Float_t scale=1.);
  Bool_t  Process3D(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=NULL,  Float_t scale=1.);
  Bool_t  Process3Dlinked(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=NULL,  Float_t scale=1.);
  Bool_t  Process3DL(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=NULL,  Float_t scale=1.);
  Bool_t  Process3Darray(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=NULL,  Float_t scale=1.);
  Bool_t  Process3DlinkedArray(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=NULL,  Float_t scale=1.);
  Bool_t  Pulls(Double_t dyz[2], Double_t cc[3], Double_t tilt) const;

  UChar_t             fSegmentLevel;    // segmentation level [sector/stack/chamber]
  UShort_t            fIdxPlot;         // plot counter (internal)
  UShort_t            fIdxFrame;        // frame counter (internal)
  UShort_t            fNcomp[kNprojs];  // number of projections per task
  Char_t              *fAxTitle[kNprojs][4]; //! Title for all ref histos
  Float_t             fPtThreshold;     // pt threshold for some performance plots
  Float_t             fDyRange;         // min/max dy
  Float_t             fXcorr[540][30];  // min/max dy
  static Char_t const *fgPerformanceName[kNviews]; //! name of performance plot
  static Char_t const *fgParticle[11];    //! latex name of particles/sign 
  static UChar_t const fgNproj[kNviews]; //! number of projections per task
  static Int_t const  fgkNresYsegm[3];   //! number of segments for saving y resolution
  static Char_t const *fgkResYsegmName[3];//! name of segment for saving y resolution
  TDatabasePDG        *fDBPDG;          // PDG database
  TObjArray           *fGraphS;         //! result holder - sigma values
  TObjArray           *fGraphM;         //! result holder - mean values

  // calibration containers
  TObjArray           *fCl;     //! cluster2track calib
  TObjArray           *fMCcl;   //! cluster2mc calib
/*  TObjArray           *fTrklt;  //! tracklet2track calib
  TObjArray           *fMCtrklt;//! tracklet2mc calib*/
  
  ClassDef(AliTRDresolution, 9) // TRD tracking resolution task
};
#endif
