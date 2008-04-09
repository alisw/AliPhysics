//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTOUTRAWREADER_H
#define ALIHLTOUTRAWREADER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTOUTRawReader.h
    @author Matthias Richter
    @date   
    @brief  HLTOUT data wrapper for AliRawReader.
*/

#include "AliHLTOUTHomerCollection.h"

class AliRawReader;
class AliHLTHOMERReader;

/**
 * @class AliHLTOUTRawReader
 * Handler of HLTOUT data for AliRawReader input.
 */
class AliHLTOUTRawReader : public AliHLTOUTHomerCollection {
 public:
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
  /** standard constructor prohibited */
  AliHLTOUTRawReader();
  /** copy constructor prohibited */
  AliHLTOUTRawReader(const AliHLTOUTRawReader&);
  /** assignment operator prohibited */
  AliHLTOUTRawReader& operator=(const AliHLTOUTRawReader&);

  /** the rawreader */
  AliRawReader* fpRawreader; //!transient

  ClassDef(AliHLTOUTRawReader, 1)
};
#endif
