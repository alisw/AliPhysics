#ifndef ALITRDCHECKTRK_H
#define ALITRDCHECKTRK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD tracker systematic                                                //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

template <typename Value> class TVectorT;
typedef struct TVectorT<Double_t> TVectorD;
class TH1;
class TObjArray;
class AliExternalTrackParam;
class AliTRDtrackV1;
class AliTRDcheckTRK : public AliTRDrecoTask
{
public:
  enum ETRDcheckTRKconst {
     kNbunchCross = 3  // no of classes for bunch crossing
    ,kNpt         = 24 // no of log bins in pt spectrum
    ,kNcharge     = 2
  };
  enum ETRDcheckTRKclasses {
     kEntry = 0
    ,kPropagation
    ,kNclasses
  };
  enum ETRDcheckTRKprojs {
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
  AliTRDcheckTRK();
  AliTRDcheckTRK(char* name);
  virtual ~AliTRDcheckTRK();
  static Float_t  GetKalmanStep()                      { return fgKalmanStep;}
  Int_t           GetSpeciesByMass(Float_t m);
  Int_t           GetPtBin(Float_t pt);
  Bool_t          GetRefFigure(Int_t ifig);
  static Bool_t   HasClRecalibrate()                   { return fgClRecalibrate;}
  static Bool_t   HasKalmanUpdate()                    { return fgKalmanUpdate;}
  TObjArray*      Histos();
  TH1*            PlotEntry(const AliTRDtrackV1 *t=NULL);
  TH1*            PlotPropagation(const AliTRDtrackV1 *t=NULL);
  static Bool_t   PropagateKalman(const AliTRDtrackV1 *t, AliExternalTrackParam *ref,
                    TVectorD *dx, TVectorD *dy, TVectorD *dz, TVectorD *dphi,
                    TVectorD *pt, TVectorD *phi, TVectorD *eta,
                    TVectorD *budget=NULL, TVectorD *c=NULL, Option_t *opt="");
  static void     SetKalmanStep(Float_t step)          { fgKalmanStep=step;}
  static void     SetClRecalibrate(Bool_t set=kTRUE)   { fgClRecalibrate=set;}
  static void     SetKalmanUpdate(Bool_t set=kTRUE)    { fgKalmanUpdate=set;}

private:
  AliTRDcheckTRK(const AliTRDcheckTRK&);
  AliTRDcheckTRK& operator=(const AliTRDcheckTRK&);
  Bool_t          MakeProjectionEtaPhi();

  // kalman related settings
  static Bool_t  fgKalmanUpdate;  // update Kalman with TRD point
  static Bool_t  fgClRecalibrate; // recalibrate clusters and recalculate tracklet fit
  static Float_t fgKalmanStep;    // Kalman stepping

  Double_t       fPtBins[kNpt+1]; // discretization of pt range
  TH1            *fProj[10];      //! array of histo projections

  ClassDef(AliTRDcheckTRK, 1) // TRD tracker systematic
};

#endif