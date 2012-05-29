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

#ifndef Root_TNamed
#include "TNamed.h"
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
    ,kXchange    = BIT(22) // use exchange containers
  };
  enum ETRDresolutionSlots {
     kClToTrk    = 2
    ,kClToMC
    ,kTrkltToTrk
    ,kTrkltToMC
    ,kNSlots     = 4
  };
  enum ETRDresolutionClass {
     kDetector=0      // cluster - detector
    ,kCluster         // cluster - track
    ,kTracklet        // tracklet - track residuals/pulls
    ,kTrackIn         // tracklet - track residuals/pulls at lower TRD entrance
    ,kMCcluster       // cluster-mc resolution/pulls
    ,kMCtracklet      // tracklet-mc resolution/pulls
    ,kMCtrackIn       // TPC track monitor
    ,kMCtrack         // TRD track monitor
/*    ,kTrackOut        // tracklet - track residuals/pulls at lower TRD entrance during refit
    ,kMCtrackOut      // TOF/HMPID track monitor*/
    ,kNclasses        // total number of resolution classes
  };
  enum ETRDresolutionClassProjs {
    kDetNproj=460      // detector projections
    ,kClNproj=2800     // cluster projections
    ,kTrkltNproj=8000  // tracklet projections
    ,kTrkInNproj=300   // trackIn projections
    ,kTrkNproj=2500    // track projections
    ,kMCTrkInNproj=260 // trackIn projections
  };
  enum ETRDresolutionProjs {
    kBC    = 0 // bunch cross
    ,kPhi
    ,kEta
    ,kYrez
    ,kPrez
    ,kZrez
    ,kSpeciesChgRC
    ,kPt
    ,kNdim  // no of dimensions in the THnSparse
    ,kNdimDet     = 4
    ,kNdimCl      = 4
    ,kNdimTrklt   = 4
    ,kNdimTrkIn   = 7
    ,kNbunchCross = 3  // no of classes for bunch crossing
    ,kNpt         = 3  // no of log bins in pt spectrum
    ,kNspc        = 3  // no of species e, mu+pi, K+p
    ,kNcharge     = 2  // no of charges
    ,kNpads       = 4  // no of charges
  };

  AliTRDresolution();
  AliTRDresolution(char* name, Bool_t xchange=kTRUE);
  AliTRDresolution(const AliTRDresolution&);
  AliTRDresolution& operator=(const AliTRDresolution&);
  virtual ~AliTRDresolution();
  
  static Bool_t   FitTrack(const Int_t np, AliTrackPoint *points, Float_t params[10]);
  static Bool_t   FitTracklet(const Int_t ly, const Int_t np, const AliTrackPoint *points, const Float_t trackPars[10], Float_t trackletPars[3]);
  void            UserCreateOutputObjects();
//  Float_t GetCorrectionX(Int_t det, Int_t tb) const {return fXcorr[det][tb];}
  static void     GetRangeZ(TH2 *h2, Float_t &m, Float_t &M);
  Float_t         GetDyRange() const {return fDyRange;}
  Float_t         GetPtThreshold() const {return fPtThreshold;}
  static Int_t    GetPtBin(Float_t pt);
  Bool_t          GetRefFigure(Int_t ifig);

  virtual TObjArray*  Histos(); 
//  Bool_t  Load(const Char_t *file = "AnalysisResults.root", const Char_t *dir="TRD_Performance");
//  Bool_t  LoadCorrection(const Char_t *file=NULL);
  void            MakeSummary();
  static void     MakePtSegmentation(Float_t pt0=0.5, Float_t dpt=0.002);

  TObjArray*      Results(ETRDresolutionClass c) const  { if(!fProj) return NULL; return (TObjArray*)fProj->At(c);}
  void            UserExec(Option_t * opt);
  void            InitExchangeContainers();
  Bool_t          HasTrackRefit() const                 { return TestBit(kTrackRefit);}
  Bool_t          HasTrackSelection() const             { return TestBit(kTrackSelect);}
  Bool_t          IsVerbose() const                     { return TestBit(kVerbose);}
  Bool_t          IsVisual() const                      { return TestBit(kVisual);}
  Bool_t          UseBCselectTOF() const                { return fBCbinTOF>0;}
  Bool_t          UseBCselectFill() const               { return fBCbinFill>0;}
  Bool_t          UseExchangeContainers() const         { return TestBit(kXchange);}
  Bool_t          PostProcess();

  TH1*            DetCluster(const TObjArray *cl=NULL);
  TH1*            PlotCluster(const AliTRDtrackV1 *t=NULL);
  TH1*            PlotTracklet(const AliTRDtrackV1 *t=NULL);
  TH1*            PlotTrackIn(const AliTRDtrackV1 *t=NULL);
//  TH1*            PlotTrackOut(const AliTRDtrackV1 *t=NULL);
  TH1*            PlotMC(const AliTRDtrackV1 *t=NULL);

  static Bool_t   Process(TH2* const /*h2*/, TGraphErrors **/*g*/, Int_t stat=100){ return Bool_t(stat);}
  void            SetDyRange(Float_t dy) {fDyRange = dy;}
  void            SetPtThreshold(Float_t pt)            { fPtThreshold = pt;}
  void            SetBCselectTOF(Int_t b=0)             { fBCbinTOF = b==0?2:(b<0?1:3);}
  void            SetBCselectFill(Int_t b=0)            { fBCbinFill = b<0||b>3499?1:b+1;}
  void            SetBsign(Int_t b=0)                   { fBsign = Bool_t(b);}
  void            SetProcesses(Bool_t det, Bool_t cl, Bool_t trklt, Bool_t trkin);
  void            SetDump3D(Bool_t det, Bool_t cl, Bool_t trklt, Bool_t trkin);
  void            SetVerbose(Bool_t v = kTRUE)          { SetBit(kVerbose, v);}
  void            SetVisual(Bool_t v = kTRUE)           { SetBit(kVisual, v);}
  void            SetTrackRefit(Bool_t v = kTRUE)       { SetBit(kTrackRefit, v);}
  void            SetTrackSelection(Bool_t v = kTRUE)   { SetBit(kTrackSelect, v);}
  void            SetUseExchangeContainers(Bool_t v = kTRUE) { SetBit(kXchange, v);}

  void            Terminate(Option_t * opt);
  static Bool_t   UseTrack(const Int_t np, const AliTrackPoint *points, Float_t params[10]);

  void        AdjustF1(TH1 *h, TF1 *f);
  TObjArray*  BuildMonitorContainerCluster(const char* name, Bool_t expand=kFALSE, Float_t range=-1.);
  TObjArray*  BuildMonitorContainerTracklet(const char* name, Bool_t expand=kFALSE);
  TObjArray*  BuildMonitorContainerTrack(const char* name);
  void        DrawSigma(TH2 *h2, const Char_t *t, Float_t m=0., Float_t M=-1., Float_t scale=1);
  void        GetLandauMpvFwhm(TF1 * const f, Float_t &mpv, Float_t &xm, Float_t &xM);
  void        GetRange(TH2 *h2, Char_t mod, Float_t *range);

protected:
  Bool_t      HasDump3DFor(ETRDresolutionClass cls) const { return TESTBIT(fSteer, 4+cls);}
  Bool_t      HasProcess(ETRDresolutionClass cls) const   { return TESTBIT(fSteer, cls);}
  Bool_t      MakeProjectionDetector();
  Bool_t      MakeProjectionCluster(Bool_t mc=kFALSE);
  Bool_t      MakeProjectionTracklet(Bool_t mc=kFALSE);
  Bool_t      MakeProjectionTrackIn(Bool_t mc=kFALSE);
  Bool_t      MakeProjectionTrack();
  Bool_t      Process(TH2* const /*h2*/, TF1 */*f*/, Float_t /*k*/, TGraphErrors **/*g*/) { return kTRUE;}
  Bool_t      Pulls(Double_t dyz[2], Double_t cc[3], Double_t tilt) const;

  UShort_t              fSteer;           // bit map to steer internal behaviour of class
                                          // MakeProjection [kTrackIn kTracklet kCluster kDetector]
                                          // Dump3D [4+kTrackIn 4+kTracklet 4+kCluster 4+kDetector]
  UShort_t              fIdxFrame;        // frame counter (internal)
  Float_t               fPtThreshold;     // pt threshold for some performance plots
  Float_t               fDyRange;         // min/max dy
  Int_t                 fBCbinTOF;        // set/select by TOF BC index
  Int_t                 fBCbinFill;       // set/select by Bunch Fill index
  Bool_t                fBsign;           // sign of magnetic field (kFALSE[-] kTRUE[+])
  static Char_t const  *fgPerformanceName[kNclasses]; //! name of performance plot
  static Int_t const    fgkNbins[kNdim];  //! no of bins/projection
  static Double_t const fgkMin[kNdim];    //! low limits for projections
  static Double_t const fgkMax[kNdim];    //! high limits for projections
  static Char_t const  *fgkTitle[kNdim];  //! title of projection
  static Float_t        fgPtBin[25];      //! pt segmentation
  TObjArray            *fProj;            //! result holder - sigma values
  TDatabasePDG         *fDBPDG;           //! PDG database

  // calibration containers
  TObjArray            *fCl;              //! cluster2track calib
  TObjArray            *fMCcl;            //! cluster2mc calib
  
  ClassDef(AliTRDresolution, 10) // TRD tracking resolution task
};
#endif
