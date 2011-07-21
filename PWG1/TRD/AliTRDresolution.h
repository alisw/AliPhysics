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
  enum ETRDresolutionSteer {
     kVerbose    = BIT(18) // switch verbosity
    ,kVisual     = BIT(19) // show partial results during processing
    ,kTrackRefit = BIT(20) // steer track refit
    ,kTrackSelect= BIT(21) // steer track selection
//    ,kLoadCorr   = BIT(23) // steer load corrections
  };
  enum ETRDresolutionSlots {
     kClToTrk    = 2
    ,kClToMC
    ,kTrkltToTrk
    ,kTrkltToMC
    ,kNSlots     = 4
  };
  enum ETRDresolutionClass {
    kCluster=0        // cluster - track
    ,kTracklet        // tracklet - track residuals/pulls
    ,kTrackIn         // tracklet - track residuals/pulls at lower TRD entrance
    ,kTrackOut        // tracklet - track residuals/pulls at lower TRD entrance during refit
    ,kMCcluster       // cluster-mc resolution/pulls
    ,kMCtracklet      // tracklet-mc resolution/pulls
    ,kMCtrackIn       // TPC track monitor
    ,kMCtrackOut      // TOF/HMPID track monitor
    ,kMCtrack         // TRD track monitor
    ,kNclasses        // total number of resolution classes
  };
  enum ETRDresolutionProjs {
    kBC    = 0 // bunch cross
    ,kPhi
    ,kEta
    ,kSpeciesChgRC
    ,kPt
    ,kYrez
    ,kZrez
    ,kPrez
    ,kNdim  // no of dimensions in the THnSparse
  };
  enum ETRDresolutionSegmentation {
     kSector    = 0
    ,kStack
    ,kDetector
    ,kNbunchCross = 3  // no of classes for bunch crossing
    ,kNpt         = 24 // no of log bins in pt spectrum
    ,kNcharge     = 2
  };

  AliTRDresolution();
  AliTRDresolution(char* name);
  virtual ~AliTRDresolution();
  
  static Bool_t  FitTrack(const Int_t np, AliTrackPoint *points, Float_t params[10]);
  static Bool_t  FitTracklet(const Int_t ly, const Int_t np, const AliTrackPoint *points, const Float_t trackPars[10], Float_t trackletPars[3]);
  void    UserCreateOutputObjects();
//  Float_t GetCorrectionX(Int_t det, Int_t tb) const {return fXcorr[det][tb];}
  Float_t GetDyRange() const {return fDyRange;}
  Float_t GetMeanWithBoundary(TH1 *h, Float_t zm, Float_t zM, Float_t dz) const;
  Float_t GetPtThreshold() const {return fPtThreshold;}
  static Int_t GetPtBin(Float_t pt);
  Bool_t  GetRefFigure(Int_t ifig);

  TObjArray*  Histos(); 
//  Bool_t  Load(const Char_t *file = "AnalysisResults.root", const Char_t *dir="TRD_Performance");
//  Bool_t  LoadCorrection(const Char_t *file=NULL);
  void    MakeSummary();
  static void MakePtSegmentation(Float_t pt0=0.2, Float_t dpt=0.002);

  TObjArray*  Results(ETRDresolutionClass c) const {if(!fProj) return NULL; return (TObjArray*)fProj->At(c);}
  void    UserExec(Option_t * opt);
  void    InitExchangeContainers();
//  Bool_t  HasLoadCorrection() const {return TestBit(kLoadCorr);}
  Bool_t  HasTrackRefit() const {return TestBit(kTrackRefit);}
  Bool_t  HasTrackSelection() const {return TestBit(kTrackSelect);}
  Bool_t  IsVerbose() const {return TestBit(kVerbose);}
  Bool_t  IsVisual() const {return TestBit(kVisual);}
  Bool_t  PostProcess();

  TH1*    PlotCluster(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotTracklet(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotTrackIn(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotTrackOut(const AliTRDtrackV1 *t=NULL);
  TH1*    PlotMC(const AliTRDtrackV1 *t=NULL);

  static Bool_t  Process(TH2* const h2, TGraphErrors **g, Int_t stat=100);
  void    SetDyRange(Float_t dy) {fDyRange = dy;}
//  void    SetSegmentationLevel(Int_t l=0);
  void    SetPtThreshold(Float_t pt) {fPtThreshold = pt;}
//  void    SetLoadCorrection(Bool_t v = kTRUE) {SetBit(kLoadCorr, v);}
  void    SetVerbose(Bool_t v = kTRUE) {SetBit(kVerbose, v);}
  void    SetVisual(Bool_t v = kTRUE) {SetBit(kVisual, v);}
  void    SetTrackRefit(Bool_t v = kTRUE) {SetBit(kTrackRefit, v);}
  void    SetTrackSelection(Bool_t v = kTRUE) {SetBit(kTrackSelect, v);}

  void    Terminate(Option_t * opt);
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
  Bool_t  MakeProjectionCluster();
  Bool_t  MakeProjectionTracklet();
  Bool_t  MakeProjectionTrackIn();
  Bool_t  Process(TH2* const h2, TF1 *f, Float_t k, TGraphErrors **g);
  Bool_t  Pulls(Double_t dyz[2], Double_t cc[3], Double_t tilt) const;

  UShort_t              fIdxPlot;         // plot counter (internal)
  UShort_t              fIdxFrame;        // frame counter (internal)
  Float_t               fPtThreshold;     // pt threshold for some performance plots
  Float_t               fDyRange;         // min/max dy
  static Char_t const  *fgPerformanceName[kNclasses]; //! name of performance plot
  static UChar_t const  fgNproj[kNclasses];//! number of projections per task
  static Int_t const    fgkNbins[kNdim];  //! no of bins/projection
  static Double_t const fgkMin[kNdim];    //! low limits for projections
  static Double_t const fgkMax[kNdim];    //! high limits for projections
  static Char_t const  *fgkTitle[kNdim];  //! title of projection 
  static Float_t        fgPtBin[kNpt+1];  //! pt segmentation
  TObjArray            *fProj;            //! result holder - sigma values
  TDatabasePDG         *fDBPDG;           //! PDG database

  // calibration containers
  TObjArray           *fCl;               //! cluster2track calib
  TObjArray           *fMCcl;             //! cluster2mc calib
  
  ClassDef(AliTRDresolution, 10) // TRD tracking resolution task
};
#endif
