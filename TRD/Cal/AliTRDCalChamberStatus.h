#ifndef ALITRDCalChamberStatus_H
#define ALITRDCalChamberStatus_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for the status of a readout chamber                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>
#include "TNamed.h"

class TH2D;
class AliTRDCalChamberStatus : public TNamed {

 public:

  enum { kNdet = 540, kNstacks = 90, kNcham = 5, kNsect = 18 };
  enum { kGood = 0, kNoData = 1, kNoDataHalfChamberSideA = 2, kNoDataHalfChamberSideB = 3, kBadCalibrated = 4};
  
  AliTRDCalChamberStatus();
  AliTRDCalChamberStatus(const Text_t* name, const Text_t* title);

  Char_t GetStatus(Int_t det) const          { return fStatus[det];   };
  void   SetStatus(Int_t det, Char_t status);
  void   UnsetStatusBit(Int_t det, Char_t status);

  Bool_t IsGood(Int_t det) const             { return TESTBIT(fStatus[det], kGood);                   }
  Bool_t IsNoData(Int_t det) const           { return TESTBIT(fStatus[det], kNoData);                 }
  Bool_t IsNoDataSideA(Int_t det) const      { return TESTBIT(fStatus[det], kNoDataHalfChamberSideA); }
  Bool_t IsNoDataSideB(Int_t det) const      { return TESTBIT(fStatus[det], kNoDataHalfChamberSideB); }
  Bool_t IsBadCalibrated(Int_t det) const    { return TESTBIT(fStatus[det], kBadCalibrated);          }

 TH2D *Plot(Int_t sm, Int_t rphi);           // Plot fStatus for sm and halfchamberside  
 TH2D *PlotNoData(Int_t sm, Int_t rphi);     // Plot data status for sm and halfchamberside  
 TH2D *PlotBadCalibrated(Int_t sm, Int_t rphi); // Plot calibration status for sm and halfchamberside  
 TH2D *Plot(Int_t sm);                       // Plot fStatus for sm 


 protected:

  Char_t fStatus[kNdet];                    //  Status byte

  ClassDef(AliTRDCalChamberStatus,1)        //  Defines the status of a single readout chamber

};

#endif
