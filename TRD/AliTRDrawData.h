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
class AliTRDCommonParam;
class AliTRDcalibDB;
class AliTRDgeometry;
class AliTRDdataArrayI;

//
// The following are defined in AliTRDRawStream.h:
//const Int_t end_of_tracklet_marker = 0xAAAAAAAA; /*This marks the end of tracklet data words*/
//const Int_t end_of_event_marker    = 0x00000000; /*This marks the end of half-chamber-data*/
//

class AliTRDrawData : public TObject {

 public:

  AliTRDrawData();
  AliTRDrawData(const AliTRDrawData &r);
  virtual ~AliTRDrawData();

  AliTRDrawData &operator=(const AliTRDrawData &/*r*/) { return *this; }

  virtual Bool_t       Digits2Raw(TTree *digits, TTree *tracks = NULL);
  virtual Bool_t       SetRawVersion(Int_t v);

  virtual AliTRDdigitsManager* Raw2Digits(AliRawReader *rawReader);

 protected:

  virtual Bool_t       Digits2RawV0(AliTRDdigitsManager* digitsManager); // for fRawVersion == 0
  virtual Bool_t       Digits2RawVx(AliTRDdigitsManager* digitsManager); // for fRawVersion > 0
  virtual Int_t        ProduceHcDataV1andV2(AliTRDdataArrayI *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize);
  //virtual Int_t        ProduceHcDataV3(AliTRDdataArrayI *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize);
  //virtual Int_t        ProduceHcDataV4(AliTRDdataArrayI *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize);

  virtual AliTRDdigitsManager* Raw2DigitsV0(AliRawReader* rawReader);
  virtual AliTRDdigitsManager* Raw2DigitsVx(AliRawReader* rawReader);

  Int_t                fRawVersion;     //  Which version of raw simulator is used
  AliTRDCommonParam   *fCommonParam;    //! Common parameters
  AliTRDcalibDB       *fCalibration;    //! Offline database interface
  AliTRDgeometry      *fGeo;            //! Geometry
  Int_t                fNumberOfDDLs;   //  Number of DDLs

  ClassDef(AliTRDrawData,3)             //  TRD raw data class

};
#endif
