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

class TKDNodeInfo;
class TKDInterpolator;
class TObjArray;
class AliTRDpidRefMakerLQ : public AliTRDpidRefMaker {
public:
  enum ETRDpidRefMakerLQsteer{
    kMaxStat    = 40000 // maximum statistics/PID bin
   ,kMinStat    = 50     // minimum statistics/bucket 14%
   ,kMinBuckets = 100    // minimum number of buckets [lambda(6)*alpha(1.5)*regions(50)]
  };
  AliTRDpidRefMakerLQ();
  AliTRDpidRefMakerLQ(const char *n);
  ~AliTRDpidRefMakerLQ();

  TObject*    GetOCDBEntry(Option_t *opt);
  Bool_t      GetRefFigure(Int_t ifig);
  Bool_t      HasOnlineMonitor() const {return kTRUE;}
  TObjArray*  Histos();
  Bool_t      Load(const Char_t *file = "AnalysisResults.root", const Char_t *dir = "TRD.CalibPIDrefMaker");
  Bool_t      PostProcess();
  void        UserCreateOutputObjects();
  void        UserExec(Option_t *opt);

private:
  AliTRDpidRefMakerLQ(const AliTRDpidRefMakerLQ &ref);
  AliTRDpidRefMakerLQ& operator=(const AliTRDpidRefMakerLQ &ref);
  void        SetZeroes(TKDInterpolator *in, TKDNodeInfo *node, Int_t n0, Int_t& idx, Float_t x, Float_t dx, Float_t y, Float_t dy, const Char_t opt='x');

  TObjArray   *fPDF;          // list of PDF estimations

  ClassDef(AliTRDpidRefMakerLQ, 6)  // Reference builder for Multidim-LQ TRD-PID

};

#endif

