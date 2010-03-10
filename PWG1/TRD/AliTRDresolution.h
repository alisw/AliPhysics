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
class AliTRDReconstructor;
class AliTRDgeometry;
class AliTRDrecoParam;
class AliTRDseedV1;
class AliTRDtrackInfo;
class AliTRDresolution : public AliTRDrecoTask
{
public:
  enum ETRDresolutionPlot {
     kCharge     =  0 // charge resolution
    ,kCluster    =  1 // cluster - track
    ,kTrackTRD   =  2 // tracklet - track residuals/pulls
    ,kTrackTPC   =  3 // tracklet - track residuals/pulls at lower TRD entrance 
    ,kMCcluster  =  4 // cluster-mc resolution/pulls
    ,kMCtracklet =  5 // tracklet-mc resolution/pulls
    ,kMCtrackTPC =  6 // TPC track monitor
    ,kMCtrackTOF =  7 // TOF/HMPID track monitor
    ,kMCtrackTRD =  8 // TRD track monitor
    ,kNviews     =  9 // total number of resolution views
    ,kNprojs     = 53 // total number of projections for all views
  };
  enum ETRDresolutionSteer {
    kVerbose  = 0
    ,kVisual  = 1
  };
  enum ETRDresolutionOutSlots {
     kClToTrk    = 2
    ,kTrkltToTrk = 3
    ,kClToMC     = 4
    ,kTrkltToMC  = 5
    ,kNOutSlots  = 4
  };

  AliTRDresolution();
  AliTRDresolution(char* name);
  virtual ~AliTRDresolution();
  
  void    UserCreateOutputObjects();
  Bool_t  GetRefFigure(Int_t ifig);
  TObjArray*  Histos(); 
  TObjArray*  Results(Int_t i=0) const {return i ? fGraphS : fGraphM;} 
  void    UserExec(Option_t * opt);
  Bool_t  IsVerbose() const {return TESTBIT(fStatus, kVerbose);}
  Bool_t  IsVisual() const {return TESTBIT(fStatus, kVisual);}
  Bool_t  PostProcess();

  TH1*    PlotCharge(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotCluster(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotTracklet(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotTrackTPC(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotMC(const AliTRDtrackV1 *t=0x0);

  void    SetRecoParam(AliTRDrecoParam *r);
  void    SetVerbose(Bool_t v = kTRUE) {v ? SETBIT(fStatus ,kVerbose): CLRBIT(fStatus ,kVerbose);}
  void    SetVisual(Bool_t v = kTRUE) {v ? SETBIT(fStatus, kVisual) : CLRBIT(fStatus, kVisual);}

  void    Terminate(Option_t * opt);
  Bool_t  GetGraphPlot(Float_t *bb, ETRDresolutionPlot ip, Int_t idx=-1);
  Bool_t  GetGraphTrack(Float_t *bb, Int_t idx, Int_t ist, Int_t n=123456789, Bool_t kLEG=kFALSE);
  Bool_t  GetGraphTrackTPC(Float_t *bb, Int_t idx, Int_t ist=0, Int_t n=123456789, Bool_t kLEG=kFALSE);
  
private:
  AliTRDresolution(const AliTRDresolution&);
  AliTRDresolution& operator=(const AliTRDresolution&);

  void    AdjustF1(TH1 *h, TF1 *f);
  TObjArray*  BuildMonitorContainerCluster(const char* name);
  TObjArray*  BuildMonitorContainerTracklet(const char* name);
  TObjArray*  BuildMonitorContainerTrack(const char* name);
  void    GetLandauMpvFwhm(TF1 * const f, Float_t &mpv, Float_t &xm, Float_t &xM);
  Bool_t  Process(TH2* const h2, TF1 *f, Float_t k, TGraphErrors **g);
  Bool_t  Process2D(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=0x0,  Float_t scale=1., Int_t gidx=-1);
  Bool_t  Process3D(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=0x0,  Float_t scale=1.);
  Bool_t  Process3Drange(ETRDresolutionPlot ip, Int_t idx=-1, Int_t gidx=-1, TF1 *f=0x0,  Float_t scale=1., Int_t zbin0=0, Int_t zbin1=0);
  Bool_t  Process3DL(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=0x0,  Float_t scale=1.);
  Bool_t  Process4D(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=0x0,  Float_t scale=1., Int_t n=-1);

  UChar_t             fStatus;          // steer parameter of the task
  UShort_t            fIdxPlot;         //! plot counter (internal)
  UShort_t            fIdxFrame;        //! frame counter (internal)
  static Char_t const *fgPerformanceName[kNviews]; // name of performance plot
  static UChar_t const fgNhistos[kNviews]; // number of histos per task
  static UChar_t const fgNproj[kNviews]; // number of projections per task
  static UChar_t const fgNcomp[kNprojs]; // number of projections per task
  static Char_t const *fgAxTitle[kNprojs][4]; // Title for all ref histos
  AliTRDReconstructor *fReconstructor;  //! local reconstructor
  AliTRDgeometry      *fGeo;            //! TRD geometry
  TDatabasePDG        *fDBPDG;          //! PDG database
  TObjArray           *fGraphS;         //! result holder - sigma values
  TObjArray           *fGraphM;         //! result holder - mean values

  // calibration containers
  TObjArray           *fCl;     //! cluster2track calib
  TObjArray           *fTrklt;  //! tracklet2track calib
  TObjArray           *fMCcl;   //! cluster2mc calib
  TObjArray           *fMCtrklt;//! tracklet2mc calib
  
  ClassDef(AliTRDresolution, 4) // TRD tracking resolution task
};
#endif
