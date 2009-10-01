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

class TH2;
class AliTRDpidRefMakerLQ : public AliTRDpidRefMaker {
public:
  enum ETRDpidRefMakerLQsteer{
    kMaxStat = 100000 // maximum statistics/species/momentum 
  };
  AliTRDpidRefMakerLQ();
  ~AliTRDpidRefMakerLQ();
 
  void      CreateOutputObjects();
  Bool_t    PostProcess();

protected:
  Float_t*  GetdEdx(AliTRDseedV1*);
  Int_t     GetNslices() { return 2;}
  void      Fill();

private:
  AliTRDpidRefMakerLQ(const AliTRDpidRefMakerLQ &ref);
  AliTRDpidRefMakerLQ& operator=(const AliTRDpidRefMakerLQ &ref);
  void   Reset();
  void   SaveReferences(const Int_t mom, const char *fn);
 
private:
  UChar_t   fPbin;        //! momentum bin
  UChar_t   fSbin;        //! species bin
  TH2       *fH2dEdx[5];  //! dE/dx data holders

  ClassDef(AliTRDpidRefMakerLQ, 4)  // Reference builder for Multidim-LQ TRD-PID

};

#endif

