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
    kCluster        = 0
    ,kTrackletY     = 1 // tracklet - track y pools
    ,kTrackletPhi   = 2 // tracklet - track angular pools residuals
    ,kMCcluster     = 3/*5*/
    ,kMCtrackletY   = 4/*6*/
    ,kMCtrackletZ   = 5/*6*/
    ,kMCtrackletPhi = 6/*7*/
    ,kMCtrackY      = 7 // Kalman Y resolution
    ,kMCtrackZ      = 8 // Kalman Z resolution
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
  TObjArray           *fClResiduals;    //!
  TObjArray           *fTrkltResiduals;    //!
  TObjArray           *fTrkltPhiResiduals;    //!
  TObjArray           *fClResolution;   //!
  TObjArray           *fTrkltResolution;//!
  
  ClassDef(AliTRDtrackingResolution, 1) // tracking resolution task
};
#endif
