#ifndef ALITRDRAWDATA_H
#define ALITRDRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Converts TRD digits into a raw data stream                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TTree;

class AliRawReader;

class AliTRDdigitsManager;
class AliTRDgeometry;
class AliTRDfeeParam;
class AliTRDdataArrayS;

class AliTRDrawData : public TObject {

 public:

  AliTRDrawData();
  AliTRDrawData(const AliTRDrawData &r);
  virtual ~AliTRDrawData();

  AliTRDrawData &operator=(const AliTRDrawData &/*r*/) { return *this; }

  virtual Bool_t       Digits2Raw(TTree *digits, TTree *tracks = NULL);

  virtual AliTRDdigitsManager *Raw2Digits(AliRawReader *rawReader);
  virtual AliTRDdigitsManager *Raw2DigitsOLD(AliRawReader *rawReader);

 protected:

  virtual Bool_t       Digits2Raw(AliTRDdigitsManager* digitsManager); // for fRawVersion > 0
  virtual Int_t        ProduceHcDataV1andV2(AliTRDdataArrayS *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize);
  virtual Int_t        ProduceHcDataV3(AliTRDdataArrayS *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize);
  
  AliTRDgeometry      *fGeo;            //! Geometry
  AliTRDfeeParam      *fFee;            //! Fee Parameters
  Int_t                fNumberOfDDLs;   //  Number of DDLs

  ClassDef(AliTRDrawData,5)             //  TRD raw data class

};
#endif
