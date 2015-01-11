#ifndef ALITRDCALSINGLECHAMBERSTATUS_H
#define ALITRDCALSINGLECHAMBERSTATUS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalSingleChamberStatus.h,v */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD calibration base class containing status values for one ROC       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

//_____________________________________________________________________________
class AliTRDCalSingleChamberStatus : public TObject {

 public:

  enum { kMasked          = 2
       , kPadBridgedLeft  = 4
       , kPadBridgedRight = 8
       , kReadSecond      = 16 
       , kNotConnected    = 32};

  AliTRDCalSingleChamberStatus();
  AliTRDCalSingleChamberStatus(Int_t p, Int_t c, Int_t cols);
  AliTRDCalSingleChamberStatus(const AliTRDCalSingleChamberStatus &c);
  virtual                      ~AliTRDCalSingleChamberStatus();
  AliTRDCalSingleChamberStatus &operator=(const AliTRDCalSingleChamberStatus &c);

  virtual void    Copy(TObject &c) const;

          Bool_t  IsMasked(Int_t col, Int_t row) const       { return ((GetStatus(col,row) & kMasked) 
                                                                       ? kTRUE 
                                                                       : kFALSE);                 };
	  Bool_t  IsBridgedLeft(Int_t col, Int_t row) const  { return ((GetStatus(col, row) & kPadBridgedLeft)                                            ? kTRUE                                                                         : kFALSE);                 };
	  Bool_t  IsBridgedRight(Int_t col, Int_t row) const { return ((GetStatus(col, row) & kPadBridgedRight)                                           ? kTRUE                                                                         : kFALSE);                 };
	  Bool_t  IsNotConnected(Int_t col, Int_t row) const { return ((GetStatus(col, row) & kNotConnected)                                           ? kTRUE                                                                         : kFALSE);                 };
          Int_t   GetNrows() const                           { return fNrows;                     };
          Int_t   GetNcols() const                           { return fNcols;                     };

          Int_t   GetChannel(Int_t col, Int_t row) const     { return row+col*fNrows;             };
          Int_t   GetNchannels() const                       { return fNchannels;                 };
          Char_t  GetStatus(Int_t ich) const                 { return fData[ich];                 };
          Char_t  GetStatus(Int_t col, Int_t row) const      { return fData[GetChannel(col,row)]; };

          void    SetStatus(Int_t ich, Char_t vd)            { fData[ich] = vd;                   };
          void    SetStatus(Int_t col, Int_t row, Char_t vd) { fData[GetChannel(col,row)] = vd;   };

 protected:

          Int_t   fPla;                    //  Plane number
          Int_t   fCha;                    //  Chamber number

          Int_t   fNrows;                  //  Number of rows
          Int_t   fNcols;                  //  Number of columns

          Int_t   fNchannels;              //  Number of channels
          Char_t *fData;                   //[fNchannels] Data

  ClassDef(AliTRDCalSingleChamberStatus,1) //  TRD ROC calibration class

};

#endif
