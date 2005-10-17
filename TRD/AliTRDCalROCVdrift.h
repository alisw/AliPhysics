#ifndef ALITRDCALROCVDRIFT_H
#define ALITRDCALROCVDRIFT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalROCVdrift.h,v */

///////////////////////////////////////////////////
//                                               //
//  TRD calibration class for Vdrift in one ROC  //
//                                               //
///////////////////////////////////////////////////

#include "AliTRDCalROC.h"

//_____________________________________________________________________________
class AliTRDCalROCVdrift : public AliTRDCalROC {

 public:

  AliTRDCalROCVdrift();
  AliTRDCalROCVdrift(Int_t p, Int_t c);
  AliTRDCalROCVdrift(const AliTRDCalROCVdrift &c);
  virtual             ~AliTRDCalROCVdrift();
  AliTRDCalROCVdrift  &operator=(const AliTRDCalROCVdrift &c);
  virtual void         Copy(TObject &c) const;

  Int_t        GetChannel(Int_t c, Int_t r)     { return r+c*fNrows; };
  Int_t        GetNchannels()       const       { return fNchannels;   };
  Float_t      GetVdrift(Int_t ich) const       { return fVdrift[ich]; };
  Float_t      GetVdrift(Int_t col, Int_t row)  { return fVdrift[GetChannel(col,row)]; };

  void         SetVdrift(Int_t ich, Float_t vd) { fVdrift[ich] = vd;   };
  void         SetVdrift(Int_t col, Int_t row, Float_t vd) 
                                                { fVdrift[GetChannel(col,row)] = vd; };

 protected:

  Int_t     fNchannels;             //  Number of channels
  Float_t  *fVdrift;                //[fNchannels] Drift velocities

  ClassDef(AliTRDCalROCVdrift,1)    //  TRD ROC calibration class for Vdrift

};

#endif
