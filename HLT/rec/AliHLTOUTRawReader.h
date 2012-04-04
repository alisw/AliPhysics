//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTRAWREADER_H
#define ALIHLTOUTRAWREADER_H
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTOUTRawReader.h
/// @author Matthias Richter
/// @date   
/// @brief  HLTOUT data wrapper for AliRawReader.
///

#include "AliHLTOUTHomerCollection.h"

class AliRawReader;
class AliHLTHOMERReader;

/**
 * @class AliHLTOUTRawReader
 * Handler of HLTOUT data for AliRawReader input.
 */
class AliHLTOUTRawReader : public AliHLTOUTHomerCollection {
 public:
  /** standard constructor */
  AliHLTOUTRawReader();
  /** constructor */
  AliHLTOUTRawReader(AliRawReader* pRawReader, int event=-1, AliHLTEsdManager* pEsdManager=NULL);
  /** destructor */
  virtual ~AliHLTOUTRawReader();

 protected:
  // interface functions of AliHLTOUTHomerCollection
  Bool_t ReadNextData(UChar_t*& data);
  int Reset();
  int GetDataSize();
  const AliRawDataHeader* GetDataHeader();
  void SelectEquipment(int equipmentType, int minEquipmentId = -1, int maxEquipmentId = -1);
  int GetEquipmentId();

 private:
  /** copy constructor prohibited */
  AliHLTOUTRawReader(const AliHLTOUTRawReader&);
  /** assignment operator prohibited */
  AliHLTOUTRawReader& operator=(const AliHLTOUTRawReader&);

  /**
   * Set the RawReader as parameter.
   * The function is for internal use only in conjunction with the
   * AliHLTOUT::New() functions.
   */
  void SetParam(AliRawReader* pRawReader) {fpRawreader=pRawReader;}

  /** the rawreader */
  AliRawReader* fpRawreader; //!transient

  ClassDef(AliHLTOUTRawReader, 1)
};
#endif
