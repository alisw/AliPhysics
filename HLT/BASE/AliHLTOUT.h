//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTOUT_H
#define ALIHLTOUT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOUT.h
    @author Matthias Richter
    @date   
    @brief  The control class for HLTOUT data.

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

                                                                          */
#include <vector>
#include "AliHLTLogging.h"

/**
 * @class AliHLTOUT
 * The control class for HLTOUT data.
 * The output of the HLT, either from the HLTOUT nodes or simulated output,
 * is transferred and stored in the HOMER format. The AliHLTOUT class 
 * implements scanning of the HOMER data for all HLTOUT DDL links and
 * abstracts access to the complete HLTOUT data.
 * 
 */
class AliHLTOUT : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTOUT();
  /** standard destructor */
  virtual ~AliHLTOUT();

  /**
   * Get number of data blocks in the HLTOUT data
   */
  int GetNofDataBlocks();

  /**
   * Select the first data block of a certain data type and specification.
   * The selection criteria can be of kAliHLTAnyDataType and/or
   * kAliHLTVoidDataSpec in order to be ignored and just match any data block.
   * @param dt    [in]  data type to match                                <br>
   * @param spec  [in]  data specification to match                       <br>
   * @return identifier >=0 if success, neg. error code if failed         <br>
   *                        -ENOENT if no block found                     <br>
   *                        -EPERM if access denied (object locked)
   */
  int SelectFirstDataBlock(AliHLTComponentDataType dt, AliHLTUInt32_t spec);

  /**
   * Select the next data block of data type and specification of the previous
   * call to @ref SelectFirstDataBlock.
   * @return identifier >=0 if success, neg. error code if failed         <br>
   *                        -ENOENT if no block found                     <br>
   *                        -EPERM if access denied (object locked)
   */
  int SelectNextDataBlock();

  /**
   * Get properties of the selected data block.
   * @param dt    [out] data type of the selected block
   * @param spec  [out] data specification of the selected block
   */
  int GetDataBlockDescription(AliHLTComponentDataType& dt, AliHLTUInt32_t& spec);

  /**
   * Get buffer of the selected data block.
   * @param pBuffer [out] buffer of the selected data block
   * @param size    [out] size of the selected data block
   */
  int GetDataBuffer(const AliHLTUInt8_t* &pBuffer, AliHLTUInt32_t& size);

  /**
   * Release buffer after use.
   * @param pBuffer [in]  buffer of the selected data block
   */
  int ReleaseDataBuffer(const AliHLTUInt8_t* pBuffer);

  /**
   * Locking guard for the AliHLTOUT object.
   * If the object is locked, the selection of data blocks can not be changed.
   */
  class AliHLTOUTLockGuard {
  public:
    /** constructor */
    AliHLTOUTLockGuard(AliHLTOUT* pInstance) : fpInstance(pInstance)
    {if (fpInstance) fpInstance->fbLocked=1;}
    /** destructor */
    ~AliHLTOUTLockGuard()
    {if (fpInstance) fpInstance->fbLocked=0;}

  private:
    /** standard constructor prohibited */
    AliHLTOUTLockGuard();
    /** copy constructor prohibited */
    AliHLTOUTLockGuard(const AliHLTOUTLockGuard&);
    /** assignment operator prohibited */
    AliHLTOUTLockGuard& operator=(const AliHLTOUTLockGuard&);

    /** the AliHLTOUT instance the guard is locking */
    AliHLTOUT* fpInstance; //!transient
  };

  /**
   * Block descriptor.
   */
  class AliHLTOUTBlockDescriptor {
  public:
    AliHLTOUTBlockDescriptor(AliHLTComponentDataType dt, AliHLTUInt32_t spec, AliHLTUInt32_t index)
      : fDataType(dt), fSpecification(spec), fIndex(index) {};
    ~AliHLTOUTBlockDescriptor() {}

    operator AliHLTComponentDataType() const {return fDataType;}
    operator AliHLTUInt32_t() const {return fSpecification;}
    int operator==(AliHLTComponentDataType dt) const {return dt==fDataType;}
    int operator==(AliHLTUInt32_t spec) const {return spec==fSpecification;}

    AliHLTUInt32_t GetIndex() const {return fIndex;}
  private:
    /** data type of the block */
    AliHLTComponentDataType fDataType; //!transient
    /** data specification of the block */
    AliHLTUInt32_t          fSpecification; //!transient
    /** index in the data stream */
    AliHLTUInt32_t          fIndex; //!transient
  };

 protected:
  /**
   * Add a block descriptor.
   * This is done by the child classes generating the index. The AliHLTOUT
   * object must be locked for index generation.
   * @param desc    the block descriptor
   * @return 0 if success, -EPERM if access denied
   */
  int AddBlockDescriptor(const AliHLTOUTBlockDescriptor desc);

 private:
  /** copy constructor prohibited */
  AliHLTOUT(const AliHLTOUT&);
  /** assignment operator prohibited */
  AliHLTOUT& operator=(const AliHLTOUT&);

  /**
   * Generate the index of the HLTOUT data.
   * Must be implemented by the child classes.
   */
  virtual int GenerateIndex()=0;

  /**
   * Get the data buffer
   * @param index   [in]  index of the block
   * @param pBuffer [out] buffer of the selected data block
   * @param size    [out] size of the selected data block
   */
  virtual int GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
			    AliHLTUInt32_t& size)=0;

  /** data type for the current block search, set from @ref SelectFirstDataBlock */
  AliHLTComponentDataType fSearchDataType; //!transient

  /** data specification for the current block search */
  AliHLTUInt32_t fSearchSpecification; //!transient

  /** instance locked or not */
  int fbLocked; //!transient

  /** list of block descriptors */
  vector<AliHLTOUTBlockDescriptor> fBlockDescList; //!transient

  /** current position in the list */
  vector<AliHLTOUTBlockDescriptor>::iterator fCurrent; //!transient

  /** data buffer under processing */
  const AliHLTUInt8_t* fpBuffer; //!transient

  ClassDef(AliHLTOUT, 0)
};
#endif
