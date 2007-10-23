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
    {if (fpInstance) fpInstance->SetStatusFlag(kLocked);}
    /** destructor */
    ~AliHLTOUTLockGuard()
    {if (fpInstance) fpInstance->ClearStatusFlag(kLocked);}

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

  enum AliHLTOUTByteOrder_t {
    /** no data block selected */
    kInvalidByteOrder=-1,
    kUnknownByteOrder=0,
    kLittleEndian,
    kBigEndian
  };

  enum AliHLTOUTDataType_t {
    kUint64 = 0,
    kUint32 = 1,
    kUint16 = 2,
    kUint8  = 3,
    kDouble = 4,
    kFloat  = 5
  };

  /**
   * Check byte order of selected block
   */
  AliHLTOUTByteOrder_t CheckByteOrder();

  /**
   * Check alignment of selected block
   */
  int CheckAlignment(AliHLTOUT::AliHLTOUTDataType_t type);

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
   * Internal status flags
   */
  enum {
    /** the HLTOUT object is locked with the current data block */
    kLocked = 0x1,
    /** childs can add block descriptors */
    kCollecting = 0x2,
    /** user of the data block has checked the byte order */
    kByteOrderChecked = 0x4,
    /** warning on byte order missmatch has been printed */
    kByteOrderWarning = 0x8,
    /** user of the data block has checked the alignment */
    kAlignmentChecked = 0x10,
    /** warning on alignment missmatch has been printed */
    kAlignmentWarning = 0x20
  };

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

  /**
   * Check byte order of data block
   */
  virtual AliHLTOUTByteOrder_t CheckBlockByteOrder(AliHLTUInt32_t index)=0;

  /**
   * Check alignment of data block
   */
  virtual int CheckBlockAlignment(AliHLTUInt32_t index, AliHLTOUT::AliHLTOUTDataType_t type)=0;

  /**
   * Select the data block of data type and specification of the previous
   * call to @ref SelectFirstDataBlock. Core function of @ref SelectFirstDataBlock
   * and @ref SelectNextDataBlock, starts to find a block at the current list
   * position. 
   * @return identifier >=0 if success, neg. error code if failed         <br>
   *                        -ENOENT if no block found                     <br>
   *                        -EPERM if access denied (object locked)
   */
  int FindAndSelectDataBlock();

  /**
   * Set status flag.
   * @param flag     flag to set
   * @return current status flags
   */
  unsigned int SetStatusFlag(unsigned int flag) {return fFlags|=flag;}

  /**
   * Clear status flag.
   * @param flag     flag to clear
   * @return current status flags
   */
  unsigned int ClearStatusFlag(unsigned int flag) {return fFlags&=~flag;}

  /**
   * Check status flag.
   * @param flag     flag to check
   * @return 1 if flag is set
   */
  int CheckStatusFlag(unsigned int flag) const {return (fFlags&flag)==flag;}

  /** data type for the current block search, set from @ref SelectFirstDataBlock */
  AliHLTComponentDataType fSearchDataType; //!transient

  /** data specification for the current block search */
  AliHLTUInt32_t fSearchSpecification; //!transient

  /** instance flags: locked, collecting, ... */
  unsigned int fFlags; //!transient

  /** list of block descriptors */
  vector<AliHLTOUTBlockDescriptor> fBlockDescList; //!transient

  /** current position in the list */
  vector<AliHLTOUTBlockDescriptor>::iterator fCurrent; //!transient

  /** data buffer under processing */
  const AliHLTUInt8_t* fpBuffer; //!transient

  ClassDef(AliHLTOUT, 0)
};
#endif
