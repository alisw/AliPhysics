#ifndef ALITRDTRACKLETOFLHELPER_H
#define ALITRDTRACKLETOFLHELPER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackletOflHelper.h 45666 2010-11-24 17:00:44Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// \class AliTRDtrackletOflHelper                                         //
// \brief The TRD offline tracklet helper class                           //
// \author Alexandru Bercuci                                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ALITRDTRACKERV1_H
#include "AliTRDtrackerV1.h"
#endif

class TVirtualPad;
class TObjArray;
class AliTRDpadPlane;
class AliTRDtrackletOflHelper:public TObject
{
public:
  enum ETRDtrackletOflHelper{
    kNcls     = 200
   ,kNormal   = 0
   ,kDeltaRay = 1
   ,kSecond   = 2
   ,kElephant = 3
   ,kUnknown  = 4
  };
  AliTRDtrackletOflHelper();
  AliTRDtrackletOflHelper(const AliTRDtrackletOflHelper &ref);
  AliTRDtrackletOflHelper& operator=(const AliTRDtrackletOflHelper &rhs);
  virtual ~AliTRDtrackletOflHelper();
  
  Int_t         ClassifyTopology();
  Int_t         Expand(TObjArray *cls, Int_t *mark=NULL, Int_t groupId=-1);
  static void   FindSolidCls(Bool_t *mark, Int_t *q);
  static Bool_t Fit(Int_t n, Double_t *dx, Double_t *dy, Double_t *s2, Double_t *par, Double_t sCut=1.5, Double_t *cov=NULL);
  Bool_t        Fit(Double_t *par, Double_t sCut=1.5) const;
  Bool_t        FitPSR(Double_t dy[200], Bool_t useSolid=kFALSE);
  TObjArray*    GetClusters() const {return fClusters;};
  void          GetColSignals(Int_t col, Int_t adc[32], Bool_t mainRow=kTRUE) const;
  Int_t         GetColSpread() const {return fCol[1]-fCol[0]+1;}
  Int_t         GetColStart() const {return fCol[0];}
  Int_t         GetColStop() const {return fCol[1];}
  Int_t         GetRow() const {return fRow;}
  Int_t         GetRMS(Double_t &r, Double_t &m, Double_t &s, Double_t &xm) const;
  Double_t      GetSyMean() const;
  Int_t         GetTbStart() const {return fTBrange[0];}
  Int_t         GetTbStop() const {return fTBrange[1];}
  Int_t         Init(AliTRDpadPlane *p, TObjArray *cls, Int_t *mark=NULL, Int_t groupId=-1);
  static Int_t  Segmentation(Int_t n, Double_t *x, Double_t *y, Int_t *Index);
  void          SetTbRange(Float_t t0, Float_t vd);
  void          View(TVirtualPad *pad=NULL);

protected:
private:
  static AliTRDtrackerV1::AliTRDLeastSquare& Fitter();

  Int_t                              fCol[2];     //! pad col start/stop
  Int_t                              fRow;        //! main pad row
  Int_t                              fTBrange[2]; //! start/stop time bin (t0, vd)
  TObjArray                          *fClusters;  // cluster array
  AliTRDpadPlane                     *fPadPlane;  //! pad plane for detector

  ClassDef(AliTRDtrackletOflHelper, 1)            // Offline tracklet helper
};

#endif
