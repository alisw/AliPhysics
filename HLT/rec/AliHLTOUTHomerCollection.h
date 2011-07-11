//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTHOMERCOLLECTION_H
#define ALIHLTOUTHOMERCOLLECTION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTOUTHomerCollection.h
    @author Matthias Richter
    @date   
    @brief  General collection for HLTOUT data in DDL format.
*/
#include "AliHLTOUTHomerBuffer.h"

class AliHLTHOMERReader;
class AliRawDataHeader;
class AliHLTEsdManager;

/**
 * @class AliHLTOUTHomerCollection
 * Handler of HLTOUT data in DDL format, base class for specific
 * handlers for RawReader or Digit data.
 *
 * The class expects the data to be in the DDL data format.
 * In contrast to the AliHLTOUTHomerBuffer, it also takes the
 * additional CDH and HLT headers and optional HLT decision into
 * account when opening the HOMER reader (see OpenReader()).
 *
 * The data access must be provided by the child class, the
 * interface is pretty much like the AliRawReader interface.
 */
class AliHLTOUTHomerCollection : public AliHLTOUTHomerBuffer {
 public:
  /** constructor */
  AliHLTOUTHomerCollection(int event=-1, AliHLTEsdManager* pEsdManager=NULL);
  /** destructor */
  virtual ~AliHLTOUTHomerCollection();

 protected:
  /**
   * Read next data form the data source.
   */
  virtual Bool_t ReadNextData(UChar_t*& data)=0;

  /**
   * Reset data stream position.
   * Reeading of data starts with the first data block.
   */
  virtual int Reset()=0;

  /**
   * Get size of the current data block.
   */
  virtual int GetDataSize()=0;

  /**
   * Get the header of the current data block
   */
  virtual const AliRawDataHeader* GetDataHeader()=0;

  /**
   * Select equipment for data readout.
   */
  virtual void SelectEquipment(int equipmentType, 
			       int minEquipmentId = -1, 
			       int maxEquipmentId = -1)=0;

  /**
   * Get equipment id of the current data block.
   */
  virtual int GetEquipmentId()=0;

  /**
   * Get the current event no.
   * The event no is set during creation of the HLTOUT object.
   * For the reconstruction it is taken from the ESD provided by the
   * AliReconstruction framework AliESDEvent::GetEventNumberInFile
   * @return event no or -1 if not set
   */
  int GetCurrentEventNo() const {return fEvent;}

  // interface function of AliHLTOUT
  int WriteESD(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size, AliHLTComponentDataType dt, AliESDEvent* tgtesd=NULL) const;

 private:
  /** copy constructor prohibited */
  AliHLTOUTHomerCollection(const AliHLTOUTHomerCollection&);
  /** assignment operator prohibited */
  AliHLTOUTHomerCollection& operator=(const AliHLTOUTHomerCollection&);

  /**
   * Generate the index of the HLTOUT data from the data buffer.
   */
  int GenerateIndex();

  /**
   * Get the data buffer
   * @param [in]  index   index of the block
   * @param [out] pBuffer buffer of the selected data block
   * @param [out] size    size of the selected data block
   */
  int GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
		    AliHLTUInt32_t& size);

  /**
   * Open HOMER reader for the data buffer.
   * The function expects the data buffer including all headers (CDH
   * and HLTOUT header). The offset for the HLT payload is determined from
   * the headers and the optional HLT decision data.
   * @param pSrc    data buffer
   * @param size    size of the buffer in byte
   * @return instance of HOMER reader
   */
  AliHLTHOMERReader* OpenReader(UChar_t* pSrc, unsigned int size);

  /** current instance of the HOMER reader */
  AliHLTHOMERReader* fpCurrent;  //!transient

protected:
  /** current event no */
  int fEvent; //!transient

private:
  /** instance of the ESD manager for writing and merging */
  AliHLTEsdManager* fpEsdManager; //!transient

  /** DDL id offset shift for index
   *  bit 16-31: DDL id, bit 0-15 block no
   */
  static const int fgkIdShift; //!transient
  
  ClassDef(AliHLTOUTHomerCollection, 0)
};
#endif
