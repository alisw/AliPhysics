#ifndef ALITRDCALPADSTATUS_H
#define ALITRDCALPADSTATUS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for the single pad status                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalSingleChamberStatus;

class AliTRDCalPadStatus : public TNamed {

 public:
 
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };
  enum { kMasked = 2, kPadBridgedLeft = 4, kPadBridgedRight = 8, kReadSecond = 16 };

  AliTRDCalPadStatus();
  AliTRDCalPadStatus(const Text_t* name, const Text_t* title);
  AliTRDCalPadStatus(const AliTRDCalPadStatus &c);   
  virtual            ~AliTRDCalPadStatus();
  AliTRDCalPadStatus &operator=(const AliTRDCalPadStatus &c);

  virtual void        Copy(TObject &c) const;

          Bool_t      IsMasked(Int_t d, Int_t col, Int_t row) const 
                                               { return CheckStatus(d, col, row, kMasked);          };
          Bool_t      IsBridgedLeft(Int_t d, Int_t col, Int_t row) const 
                                               { return CheckStatus(d, col, row, kPadBridgedLeft);  };
          Bool_t      IsBridgedRight(Int_t d, Int_t col, Int_t row) const 
                                               { return CheckStatus(d, col, row, kPadBridgedRight); };
	  Bool_t      IsReadSecond(Int_t d, Int_t col, Int_t row) const 
	                                       { return CheckStatus(d, col, row, kReadSecond); };
          Bool_t      CheckStatus(Int_t d, Int_t col, Int_t row, Int_t bitMask) const;

  AliTRDCalSingleChamberStatus *GetCalROC(Int_t d) const { return fROC[d]; };
  AliTRDCalSingleChamberStatus *GetCalROC(Int_t p, Int_t c, Int_t s) const;

  // Plot functions
  TH1F    *MakeHisto1D();
  TH2F    *MakeHisto2DSmPl(Int_t sm, Int_t pl);
  void     PlotHistos2DSm(Int_t sm, const Char_t *name);

 protected:

  AliTRDCalSingleChamberStatus *fROC[kNdet];          //  Array of ROC objects which contain the values per pad

  ClassDef(AliTRDCalPadStatus,1)                      //  TRD calibration class for the single pad status

};
          
#endif
