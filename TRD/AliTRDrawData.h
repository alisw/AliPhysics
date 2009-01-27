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
class AliTRDarrayADC;

class AliTRDrawData : public TObject {

 public:

  AliTRDrawData();
  AliTRDrawData(const AliTRDrawData &r);
  virtual ~AliTRDrawData();

  AliTRDrawData &operator=(const AliTRDrawData &/*r*/) { return *this; }

  virtual Bool_t       Digits2Raw(TTree *digits, TTree *tracks = NULL);

  virtual AliTRDdigitsManager *Raw2Digits(AliRawReader *rawReader);
  virtual AliTRDdigitsManager *Raw2DigitsOLD(AliRawReader *rawReader);
  static void SetRawFormatVersion(Int_t iver){ fgRawFormatVersion=iver; };
  static void SetSuppressionLevel(Int_t ilevel){ fgDataSuppressionLevel=ilevel; };

  enum FORMATTYPE
    {
      kRawOldFormat  =  0,
      kRawNewFormat  =  1
    };

 protected:

  virtual Bool_t       Digits2Raw(AliTRDdigitsManager* digitsManager); // for fRawVersion > 0
  virtual Int_t        ProduceHcData(AliTRDarrayADC *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize, Bool_t newEvent, Bool_t newSM);
  virtual Int_t        ProduceHcDataV1andV2(AliTRDarrayADC *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize);
  virtual Int_t        ProduceHcDataV3(AliTRDarrayADC *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize, Bool_t newEvent);
  //virtual Int_t      ProduceHcDataV3(AliTRDarrayADC *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize);
  		  void 	       ProduceSMIndexData(UInt_t *buf, Int_t& nw);				// SM index words and header - real data format
          void         WriteIntermediateWords(UInt_t *buf, Int_t& nw, Int_t& of, const Int_t& maxSize, const Int_t& det, const Int_t& side); // writes tracklet-endmarker and additional words between tracklet and raw-data
          void   	   WriteIntermediateWordsV2(UInt_t *buf, Int_t& nw, Int_t& of, const Int_t& maxSize, const Int_t& det, const Int_t& side); // real data format
		  void         AssignStackMask(UInt_t *buf, Int_t nStack);  // re-assignment of stack mask in the SM index word
          void         AssignLinkMask(UInt_t *buf, Int_t nLayer);   // re-assignment of link mask in the stack index word
          Int_t        AddStackIndexWords(UInt_t *buf, Int_t nStack, Int_t nMax);   // add stack index words and stack header when there is no data for the stack 
          Bool_t       ShiftWords(UInt_t *buf, Int_t nStart, Int_t nWords, Int_t nMax); // shifts n words
  
  AliTRDgeometry      *fGeo;            //! Geometry
  AliTRDfeeParam      *fFee;            //! Fee Parameters
  Int_t                fNumberOfDDLs;   //  Number of DDLs

  ClassDef(AliTRDrawData,5)             //  TRD raw data class

 private:

	static       Int_t  fgRawFormatVersion;           	  // simulation raw data version - 0:old , 1:new(real data format)
	static       Int_t  fgDataSuppressionLevel;                 // Data suppression level - 0:no su, 1: su, 2: deep suppression 
 	static const UInt_t fgkEndOfTrackletMarker  = 0x10001000; // This marks the end of tracklet data words

	Int_t   fSMindexPos;                // Position of SM index word
    Int_t   fStackindexPos;             // Position of SM index word
    UInt_t  fEventCounter;              // Event counter(starting from 1)




};
#endif






