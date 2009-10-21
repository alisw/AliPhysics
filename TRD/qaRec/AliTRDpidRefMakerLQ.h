#ifndef ALITRDPIDREFMAKERLQ_H
#define ALITRDPIDREFMAKERLQ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDpidRefMakerLQ.h 34125 2009-08-06 09:35:40Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for building reference data for PID                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ALIPID_H
#include "AliPID.h"
#endif

#ifndef ALITRDPIDREFMAKER_H
#include "AliTRDpidRefMaker.h"
#endif

#ifndef ALITRDPIDUTIL_H
#include "AliTRDpidUtil.h"
#endif

class TObjArray;
class AliTRDpidRefMakerLQ : public AliTRDpidRefMaker {
public:
  enum ETRDpidRefMakerLQsteer{
    kMinStat = 50    // minimum statistics/bucket 14%
   ,kMinBuckets = 450 // minimum number of buckets [lambda(6)*alpha(1.5)*regions(50)]
  };
  AliTRDpidRefMakerLQ();
  ~AliTRDpidRefMakerLQ();
 
  void      CreateOutputObjects();
  TObject*  GetOCDBEntry(Option_t *opt);
  Bool_t    GetRefFigure(Int_t ifig);
  Bool_t    PostProcess();

protected:
  Float_t*  CookdEdx(AliTRDseedV1 *trklt);
  Int_t     GetNslices() const { return 2;}
  void      Fill();

private:
  AliTRDpidRefMakerLQ(const AliTRDpidRefMakerLQ &ref);
  AliTRDpidRefMakerLQ& operator=(const AliTRDpidRefMakerLQ &ref);
 
private:
  UChar_t     fPIDbin;       //! momentum&species bin

  ClassDef(AliTRDpidRefMakerLQ, 5)  // Reference builder for Multidim-LQ TRD-PID

};

#endif

