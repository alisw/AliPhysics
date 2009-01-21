#ifndef ALITRDTRACKINGRESOLUTION_H
#define ALITRDTRACKINGRESOLUTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingResolution.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
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
class AliTRDtrackingResolution : public AliTRDrecoTask
{
public:
  enum{
    kCluster        = 0 // cluster - track
    ,kTrackletY     = 1 // tracklet - track y pools
    ,kTrackletPhi   = 2 // tracklet - track angular pools residuals
    ,kMCcluster     = 3 // cluster - mc residuals/systematics
    ,kMCtrackletY   = 4 // tracklet - mc y resolution/systematics
    ,kMCtrackletZ   = 5 // tracklet - mc z resolution/systematics (pad row cross)
    ,kMCtrackletPhi = 6 // tracklet - mc phi resolution/systematics
    ,kMCtrackY      = 7 // Kalman Y resolution
    ,kMCtrackZ      = 8 // Kalman Z resolution
    ,kMCtrackPt     = 9 // Kalman Pt resolution
  };
  enum{
    kVerbose  = 0
    ,kVisual  = 1
  };

  AliTRDtrackingResolution();
  virtual ~AliTRDtrackingResolution();
  
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
  TH1*    PlotResolution(const AliTRDtrackV1 *t=0x0);

  void    SetRecoParam(AliTRDrecoParam *r);
  void    SetVerbose(Bool_t v = kTRUE) {v ? SETBIT(fStatus ,kVerbose): CLRBIT(fStatus ,kVerbose);}
  void    SetVisual(Bool_t v = kTRUE) {v ? SETBIT(fStatus, kVisual) : CLRBIT(fStatus, kVisual);}

  void    Terminate(Option_t *);
  
private:
  AliTRDtrackingResolution(const AliTRDtrackingResolution&);
  AliTRDtrackingResolution& operator=(const AliTRDtrackingResolution&);
  void        AdjustF1(TH1 *h, TF1 *f);

private:
  UChar_t             fStatus;          // steer parameter of the task
  AliTRDReconstructor *fReconstructor;  //! local reconstructor
  AliTRDgeometry      *fGeo;            //! TRD geometry
  TObjArray           *fGraphS;         //! result holder - sigma values
  TObjArray           *fGraphM;         //! result holder - mean values

  // calibration containers
  TObjArray           *fCl;     //! cluster2track calib
  TObjArray           *fTrklt;  //! tracklet2track calib
  TObjArray           *fMCcl;   //! cluster2mc calib
  TObjArray           *fMCtrklt;//! tracklet2mc calib
  
  ClassDef(AliTRDtrackingResolution, 1) // TRD tracking resolution task
};
#endif
