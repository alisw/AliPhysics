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

class AliTRDReconstructor;
class AliTRDrecoParam;
class AliTRDseedV1;
class AliTRDtrackInfo;
class AliTRDtrackingResolution : public AliTRDrecoTask
{
public:
  enum{
    kClusterYResidual         = 0
    ,kTrackletRiemanYResidual = 1 // Riemann track model
    ,kTrackletRiemanAngleResidual = 2
    ,kTrackletKalmanYResidual = 3 // Kalman track model
    ,kTrackletKalmanAngleResidual = 4
    ,kTrackletYResolution     = 5
    ,kTrackletAngleResolution = 6
    ,kTrackRYResolution       = 7 // Riemann track model
    ,kTrackRZResolution       = 8
    ,kTrackRAngleResolution   = 9
    ,kTrackKYResolution       = 10 // Kalman track model
    ,kTrackKZResolution       = 11
    ,kTrackKAngleResolution   = 12
    ,kGraphStart              = 13 // First graph
  };

  AliTRDtrackingResolution();
  virtual ~AliTRDtrackingResolution();
  
  void    CreateOutputObjects();
  void    Exec(Option_t *);
  void    GetRefFigure(Int_t ifig, Int_t &first, Int_t &last);  
  void    SetRecoParam(AliTRDrecoParam *r);
  Bool_t  PostProcess();
  void    Terminate(Option_t *);
  
private:
  AliTRDtrackingResolution(const AliTRDtrackingResolution&);
  AliTRDtrackingResolution& operator=(const AliTRDtrackingResolution&);
  TObjArray*  Histos(); 
  Bool_t      Resolution(AliTRDseedV1 *tracklet, AliTRDtrackInfo *info, Double_t &p, Double_t &y, Double_t &z, Double_t &phi, Double_t &theta);

private:
  enum{
    kNLayers = 6
  };
  
  AliTRDReconstructor   *fReconstructor;  //! local reconstructor
  
  ClassDef(AliTRDtrackingResolution, 1) // tracking resolution task
};
#endif
