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
    kClusterResidual         = 0
    ,kTrackletYResidual      = 1 // tracklet - track residuals
    ,kTrackletPhiResidual    = 2 // tracklet - track angular residuals
//     ,kTrackletKalmanYResidual = 3 // Kalman track model
//     ,kTrackletKalmanAngleResidual = 4
    ,kClusterResolution      = 3/*5*/
    ,kTrackletYResolution    = 4/*6*/
    ,kTrackletZResolution    = 5/*6*/
    ,kTrackletAngleResolution= 6/*7*/
//     ,kTrackRYResolution       = 8 // Riemann track model
//     ,kTrackRZResolution       = 9
//     ,kTrackRAngleResolution   = 10
//     ,kTrackKYResolution       = 11 // Kalman track model
//     ,kTrackKZResolution       = 12
//     ,kTrackKAngleResolution   = 13
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

  TH1*    PlotClusterResiduals(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotTrackletResiduals(const AliTRDtrackV1 *t=0x0);
  TH1*    PlotTrackletPhiResiduals(const AliTRDtrackV1 *t=0x0);
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
