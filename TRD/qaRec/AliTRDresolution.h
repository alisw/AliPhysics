#ifndef ALITRDRESOLUTION_H
#define ALITRDRESOLUTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDresolution.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Resolution QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TH1;
class TF1;
class TObjArray;
class AliTRDReconstructor;
class AliTRDgeometry;
class AliTRDrecoParam;
class AliTRDseedV1;
class AliTRDtrackInfo;
class AliTRDresolution : public AliTRDrecoTask
{
public:
  enum ETRDresolutionPlot {
    kCluster          =  0 // cluster - track
    ,kTracklet        =  1 // tracklet - track y pulls
    ,kTrackletPhi     =  2 // tracklet - track angular pulls residuals
    ,kMCcluster       =  3 // cluster-mc resolution/pulls
    ,kMCtracklet      =  4 // tracklet-mc resolution/pulls
    ,kMCtrackTPC      =  5 // TPC track monitor
    ,kMCtrackHMPID    =  6 // TOF/HMPID track monitor
    ,kMCtrack         =  7 // TRD track monitor
    ,kNhistos         =  8
  };
  enum ETRDresolutionSteer {
    kVerbose  = 0
    ,kVisual  = 1
  };

  AliTRDresolution();
  virtual ~AliTRDresolution();
  
  void    CreateOutputObjects();
  Bool_t  GetRefFigure(Int_t ifig);
  TObjArray*  Histos(); 
  void    Exec(Option_t *);
  Bool_t  IsVerbose() const {return TESTBIT(fStatus, kVerbose);}
  Bool_t  IsVisual() const {return TESTBIT(fStatus, kVisual);}
  Bool_t  PostProcess();

  TH1*    PlotCluster(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotTracklet(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotTrackletPhi(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotTrackIn(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotMC(const AliTRDtrackV1 *t=0x0);
  void    DumpAxTitle(Int_t plot, Int_t fig=-1);

  void    SetRecoParam(AliTRDrecoParam *r);
  void    SetVerbose(Bool_t v = kTRUE) {v ? SETBIT(fStatus ,kVerbose): CLRBIT(fStatus ,kVerbose);}
  void    SetVisual(Bool_t v = kTRUE) {v ? SETBIT(fStatus, kVisual) : CLRBIT(fStatus, kVisual);}

  void    Terminate(Option_t *);
  Bool_t  GetGraphPlot(Float_t *bb, ETRDresolutionPlot ip, Int_t idx=-1);
  
private:
  AliTRDresolution(const AliTRDresolution&);
  AliTRDresolution& operator=(const AliTRDresolution&);
  void    AdjustF1(TH1 *h, TF1 *f);
  Bool_t  Process(ETRDresolutionPlot ip, Int_t idx=-1, TF1 *f=0x0,  Float_t scale=1.);
  Bool_t  Process3D(ETRDresolutionPlot ip, TF1 *f=0x0,  Float_t scale=1.);

  UChar_t             fStatus;          // steer parameter of the task
  static UChar_t      fNElements[kNhistos]; // number of componets per task
  static Char_t       *fAxTitle[32][3]; // axis title for all ref histos
  AliTRDReconstructor *fReconstructor;  //! local reconstructor
  AliTRDgeometry      *fGeo;            //! TRD geometry
  TObjArray           *fGraphS;         //! result holder - sigma values
  TObjArray           *fGraphM;         //! result holder - mean values

  // calibration containers
  TObjArray           *fCl;     //! cluster2track calib
  TObjArray           *fTrklt;  //! tracklet2track calib
  TObjArray           *fMCcl;   //! cluster2mc calib
  TObjArray           *fMCtrklt;//! tracklet2mc calib
  
  ClassDef(AliTRDresolution, 2) // TRD tracking resolution task
};
#endif
