//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTOUTRAWREADER_H
#define ALIHLTOUTRAWREADER_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOUTRawReader.h
    @author Matthias Richter
    @date   
    @brief  HLTOUT data wrapper for AliRawReader.

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
                                                                          */
#include "AliHLTOUTHomerBuffer.h"

class AliRawReader;
class AliHLTHOMERReader;

/**
 * @class AliHLTOUTRawReader
 * Handler of HLTOUT data for AliRawReader input.
 */
class AliHLTOUTRawReader : public AliHLTOUTHomerBuffer {
 public:
  /** constructor */
  AliHLTOUTRawReader(AliRawReader* pRawReader);
  /** destructor */
  virtual ~AliHLTOUTRawReader();

 protected:

 private:
  /** standard constructor prohibited */
  AliHLTOUTRawReader();
  /** copy constructor prohibited */
  AliHLTOUTRawReader(const AliHLTOUTRawReader&);
  /** assignment operator prohibited */
  AliHLTOUTRawReader& operator=(const AliHLTOUTRawReader&);

  /**
   * Generate the index of the HLTOUT data from the data buffer.
   */
  int GenerateIndex();

  /**
   * Get the data buffer
   * @param index   [in]  index of the block
   * @param pBuffer [out] buffer of the selected data block
   * @param size    [out] size of the selected data block
   */
  int GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
		    AliHLTUInt32_t& size);

  /** the rawreader */
  AliRawReader* fpRawreader; //!transient

  /** current instance of the HOMER reader */
  AliHLTHOMERReader* fpCurrent;  //!transient

  /** DDL id offset shift for index
   *  bit 16-31: DDL id, bit 0-15 block no
   */
  static const int fgkIdShift; //!transient
  
  ClassDef(AliHLTOUTRawReader, 0)
};
#endif
