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

class TAxis;
class TH1;
class TH2;
class TH3;
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
    kCluster=0        // cluster - track
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
    kClNproj=60        // cluster projections
    ,kTrkltNproj=2580  // tracklet projections
    ,kTrkInNproj=142   // trackIn projections
    ,kTrkNproj=2500    // track projections
    ,kMCTrkInNproj=162 // trackIn projections
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
  };
  enum ETRDresolutionSize {
     kNdimCl      = 4
    ,kNdimTrklt   = 4
    ,kNdimTrkIn   = 7
    ,kNbunchCross = 3  // no of classes for bunch crossing
    ,kNpt         = 3  // no of log bins in pt spectrum
    ,kNcharge     = 2
  };

  AliTRDresolution();
  AliTRDresolution(char* name, Bool_t xchange=kTRUE);
  virtual ~AliTRDresolution();
  
  static Bool_t   FitTrack(const Int_t np, AliTrackPoint *points, Float_t params[10]);
  static Bool_t   FitTracklet(const Int_t ly, const Int_t np, const AliTrackPoint *points, const Float_t trackPars[10], Float_t trackletPars[3]);
  void            UserCreateOutputObjects();
//  Float_t GetCorrectionX(Int_t det, Int_t tb) const {return fXcorr[det][tb];}
  Float_t         GetDyRange() const {return fDyRange;}
  static Float_t  GetMeanStat(TH1 *h, Float_t cut=0., Option_t *opt="");
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

  TH1*            PlotCluster(const AliTRDtrackV1 *t=NULL);
  TH1*            PlotTracklet(const AliTRDtrackV1 *t=NULL);
  TH1*            PlotTrackIn(const AliTRDtrackV1 *t=NULL);
//  TH1*            PlotTrackOut(const AliTRDtrackV1 *t=NULL);
  TH1*            PlotMC(const AliTRDtrackV1 *t=NULL);

  static Bool_t   Process(TH2* const h2, TGraphErrors **g, Int_t stat=100);
  void            SetDyRange(Float_t dy) {fDyRange = dy;}
  void            SetPtThreshold(Float_t pt)            { fPtThreshold = pt;}
  void            SetBCselectTOF(Int_t b=0)             { fBCbinTOF = b==0?2:(b<0?1:3);}
  void            SetBCselectFill(Int_t b=0)            { fBCbinFill = b<0||b>3499?1:b+1;}
  void            SetVerbose(Bool_t v = kTRUE)          { SetBit(kVerbose, v);}
  void            SetVisual(Bool_t v = kTRUE)           { SetBit(kVisual, v);}
  void            SetTrackRefit(Bool_t v = kTRUE)       { SetBit(kTrackRefit, v);}
  void            SetTrackSelection(Bool_t v = kTRUE)   { SetBit(kTrackSelect, v);}
  void            SetUseExchangeContainers(Bool_t v = kTRUE) { SetBit(kXchange, v);}

  void            Terminate(Option_t * opt);
  static Bool_t   UseTrack(const Int_t np, const AliTrackPoint *points, Float_t params[10]);

  class AliTRDresolutionProjection : public TNamed
  {
    friend class AliTRDresolution;  // Friend class
  public:
    AliTRDresolutionProjection();
    virtual ~AliTRDresolutionProjection();
    AliTRDresolutionProjection& operator+=(const AliTRDresolutionProjection& other);
    AliTRDresolutionProjection& operator=(const AliTRDresolutionProjection& other);
    void  Build(const Char_t *n, const Char_t *t, Int_t ix, Int_t iy, Int_t iz, TAxis *aa[]);
    void  Increment(Int_t bin[], Double_t v);
    TH2*  Projection2D(const Int_t nstat, const Int_t ncol, const Int_t mid=0);
    void  SetRebinStrategy(Int_t n, Int_t rebx[], Int_t reby[]);
    void  SetShowRange(Float_t zm, Float_t zM, Float_t em=0., Float_t eM=0.) {fRange[0] = zm; fRange[1] = zM; fRange[2] = em; fRange[3] = eM;}
  private:
    AliTRDresolutionProjection(const AliTRDresolutionProjection&);
  protected:
    TH3  *fH;          // data container
    Int_t fAx[3];      // projection axes
    Int_t fNrebin;     // no. of rebinning steps
    Int_t *fRebinX;    //[fNrebin] rebinning of the X axis
    Int_t *fRebinY;    //[fNrebin] rebinning of the Y axis
    Float_t fRange[4]; //show range of the z processed

    ClassDef(AliTRDresolutionProjection, 2)  // wrapper for a projection container THnSparse -> TH3
  };

  AliTRDresolution(const AliTRDresolution&);
  AliTRDresolution& operator=(const AliTRDresolution&);

  void        AdjustF1(TH1 *h, TF1 *f);
  TObjArray*  BuildMonitorContainerCluster(const char* name, Bool_t expand=kFALSE, Float_t range=-1.);
  TObjArray*  BuildMonitorContainerTracklet(const char* name, Bool_t expand=kFALSE);
  TObjArray*  BuildMonitorContainerTrack(const char* name);
  void        DrawSigma(TH2 *h2, const Char_t *t, Float_t m=0., Float_t M=-1., Float_t scale=1);
  void        GetLandauMpvFwhm(TF1 * const f, Float_t &mpv, Float_t &xm, Float_t &xM);
  void        GetRange(TH2 *h2, Char_t mod, Float_t *range);

protected:
  Bool_t      MakeProjectionCluster(Bool_t mc=kFALSE);
  Bool_t      MakeProjectionTracklet(Bool_t mc=kFALSE);
  Bool_t      MakeProjectionTrackIn(Bool_t mc=kFALSE);
  Bool_t      MakeProjectionTrack();
  Bool_t      Process(TH2* const h2, TF1 *f, Float_t k, TGraphErrors **g);
  Bool_t      Pulls(Double_t dyz[2], Double_t cc[3], Double_t tilt) const;

  UShort_t              fIdxPlot;         // plot counter (internal)
  UShort_t              fIdxFrame;        // frame counter (internal)
  Float_t               fPtThreshold;     // pt threshold for some performance plots
  Float_t               fDyRange;         // min/max dy
  Int_t                 fBCbinTOF;        // set/select by TOF BC index
  Int_t                 fBCbinFill;       // set/select by Bunch Fill index
  static Char_t const  *fgPerformanceName[kNclasses]; //! name of performance plot
  static Int_t const    fgkNbins[kNdim];  //! no of bins/projection
  static Double_t const fgkMin[kNdim];    //! low limits for projections
  static Double_t const fgkMax[kNdim];    //! high limits for projections
  static Char_t const  *fgkTitle[kNdim];  //! title of projection 
  static Float_t        fgPtBin[kNpt+1];  //! pt segmentation
  TObjArray            *fProj;            //! result holder - sigma values
  TDatabasePDG         *fDBPDG;           //! PDG database

  // calibration containers
  TObjArray            *fCl;              //! cluster2track calib
  TObjArray            *fMCcl;            //! cluster2mc calib
  
  ClassDef(AliTRDresolution, 10) // TRD tracking resolution task
};
#endif
