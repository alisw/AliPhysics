//-*- Mode: C++ -*-
// $Id$

#ifndef ALIRAWREADERHLT_H
#define ALIRAWREADERHLT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliRawReaderHLT.h
    @author Matthias Richter
    @date   
    @brief  AliRawReader implementation which replaces original input of
            detectors with the appropriate HLT output.                    */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliRawReader.h"
#include "TString.h"

/**
 * @class AliRawReaderHLT
 * Handler of HLTOUT data for AliRawReader input.
 */
class AliRawReaderHLT : public AliRawReader {
 public:
  /** constructor */
  AliRawReaderHLT(AliRawReader* pParentReader, const char* options=NULL);
  /** destructor */
  virtual ~AliRawReaderHLT();

  // interface methods of AliRawReader
  void     Select(Int_t detectorID, 
		  Int_t minDDLID = -1, Int_t maxDDLID = -1);
//   void     Select(const char *detectorName, 
// 		  Int_t minDDLID = -1, Int_t maxDDLID = -1);
  void     SelectEquipment(Int_t equipmentType, 
			   Int_t minEquipmentId = -1, 
			   Int_t maxEquipmentId = -1);
  void     SkipInvalid(Bool_t skip = kTRUE);
  void     SelectEvents(Int_t type);

  UInt_t   GetType() const;
  UInt_t   GetRunNumber() const;
  const UInt_t* GetEventId() const;
  const UInt_t* GetTriggerPattern() const;
  const UInt_t* GetDetectorPattern() const;
  const UInt_t* GetAttributes() const;
  const UInt_t* GetSubEventAttributes() const;
  UInt_t   GetLDCId() const;
  UInt_t   GetGDCId() const;
  UInt_t   GetTimestamp() const;

  const UInt_t* GetEquipmentAttributes() const;
  Int_t    GetEquipmentElementSize() const;
  Int_t    GetEquipmentHeaderSize() const;

  Int_t    GetEquipmentSize() const;
  Int_t    GetEquipmentType() const;
  Int_t    GetEquipmentId() const;
  Bool_t   ReadHeader();
  Bool_t   ReadNextData(UChar_t*& data);
  Bool_t   ReadNextInt(UInt_t& data);
  Bool_t   ReadNextShort(UShort_t& data);
  Bool_t   ReadNextChar(UChar_t& data);
  Bool_t   ReadNext(UChar_t* data, Int_t size);

  Bool_t   Reset();

  Bool_t   NextEvent();
  Bool_t   RewindEvents();

 protected:

 private:
  /** standard constructor prohibited */
  AliRawReaderHLT();
  /** copy constructor prohibited */
  AliRawReaderHLT(const AliRawReaderHLT&);
  /** assignment operator prohibited */
  AliRawReaderHLT& operator=(const AliRawReaderHLT&);

  /** the rawreader */
  AliRawReader* fpParentReader; //!transient

  /** options */
  TString fOptions; //!transient
  
  ClassDef(AliRawReaderHLT, 0)
};

#define ALIHLTREC_LIBRARY                   "libHLTrec.so"
#define ALIHLTREC_LIBRARY_VERSION           0
#define ALIRAWREADERHLT_CREATE_INSTANCE     "AliRawReaderHLTCreateInstance"

#ifdef __cplusplus
extern "C" {
#endif
  typedef AliRawReader* (*AliRawReaderHLTCreateInstance_t)(AliRawReader* pParentReader, const char* options);

  /**
   * Create an instance of the AliRawReader class
   */
  AliRawReader* AliRawReaderHLTCreateInstance(AliRawReader* pParentReader, const char* options);
#ifdef __cplusplus
}
#endif
#endif
