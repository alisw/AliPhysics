//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUT_H
#define ALIHLTOUT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTOUT.h
/// @author Matthias Richter
/// @date   
/// @brief  The control class for HLTOUT data.
///

#include <vector>
#include "AliHLTLogging.h"
#include "AliHLTModuleAgent.h"

class AliHLTOUTHandler;
class AliHLTOUTHandlerDesc; // AliHLTModuleAgent.h
class AliESDEvent;
class AliHLTReconstructor;
class AliRawReader;
class TTree;

#define AliHLTOUTInvalidIndex (~(AliHLTUInt32_t)0)

typedef vector<AliHLTUInt32_t> AliHLTOUTIndexList;

/**
 * @class AliHLTOUT
 * The control class for HLTOUT data.
 * The output of the HLT, either from the HLTOUT nodes or simulated output,
 * is transferred and stored in the HOMER format. The AliHLTOUT class 
 * implements scanning of the HOMER data for all HLTOUT DDL links and
 * abstracts access to the complete HLTOUT data.
 * 
 */
class AliHLTOUT {
 public:
  /** standard constructor */
  AliHLTOUT();
  /** standard destructor */
  virtual ~AliHLTOUT();

  /**
   * Create an AliHLTOUTRawReader instance.
   * Helper function to transparently access classes from the
   * libHLTrec library.
   */
  static AliHLTOUT* New(AliRawReader* pRawReader);

  /**
   * Create an AliHLTOUTDigitReader instance
   * Helper function to transparently access classes from the
   * libHLTrec library.
   */
  static AliHLTOUT* New(TTree* pDigitTree, int event=-1);

  /**
   * Create an AliHLTOUTDigitReader instance
   * Helper function to transparently access classes from the
   * libHLTrec library.
   */
  static AliHLTOUT* New(const char* filename, int event=-1);

  /**
   * Get the global instance.
   * The global instance is set temporarily by the AliHLTOUTGlobalInstanceGuard
   * mainly for the sake of data input to an analysis chain. The
   * AliHLTOUTPublisherComponent objects can access the global instance.
   */
  static AliHLTOUT* GetGlobalInstance() {return fgGlobalInstance;}

  /**
   * Delete an instance of the HLTOUT.
   * Helper function to transparently access classes from the
   * libHLTrec library. Before deleting, the availability of the
   * library is checked.
   */
  static void Delete(AliHLTOUT* pInstance);

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
   * Guard for handling of global AliHLTOUT instance.
   * The global HLTOUT instance can be set for certain steps of the
   * processing. The initial objective is to support input to the
   * AliHLTOUTPublisherComponent running in a kChain handler.
   *
   * The Guard restores the original instance when closed.
   */
  class AliHLTOUTGlobalInstanceGuard {
  public:
    /** constructor */
    AliHLTOUTGlobalInstanceGuard(AliHLTOUT* pInstance) : fpLastInstance(AliHLTOUT::fgGlobalInstance)
    {AliHLTOUT::fgGlobalInstance=pInstance;}
    /** destructor */
    ~AliHLTOUTGlobalInstanceGuard()
    {AliHLTOUT::fgGlobalInstance=fpLastInstance;}

  private:
    /** standard constructor prohibited */
    AliHLTOUTGlobalInstanceGuard();
    /** copy constructor prohibited */
    AliHLTOUTGlobalInstanceGuard(const AliHLTOUTGlobalInstanceGuard&);
    /** assignment operator prohibited */
    AliHLTOUTGlobalInstanceGuard& operator=(const AliHLTOUTGlobalInstanceGuard&);

    /** the AliHLTOUT instance the guard is locking */
    AliHLTOUT* fpLastInstance; //!transient
  };

  /**
   * The HLT OUT Event Header.
   * Defined between HLT and DAQ.
   */
  struct AliHLTOUTEventHeader {
    /**Total length of the data in bytes, including HLT event header, excluding CDH. */
    AliHLTUInt32_t fLength; //! see above
    /** version of the header */
    AliHLTUInt32_t fVersion; //! see above
    /** High 32 bit word of event id */
    AliHLTUInt32_t fEventIDHigh; //! see above
    /** Low 32 bit word of event id */
    AliHLTUInt32_t fEventIDLow; //! see above
  };

  enum {
    /// versions 1 of the HLT header
    kVersion1 = 1,
    /// versions 2 of the HLT header
    kVersion2 = 2
  };

  enum {
    /// size of HLT decision in data format version 1: 29x4
    kSizeDecisionVersion1 = 116,
    /// size of HLT decision in data format version 2: 30x4
    kSizeDecisionVersion2 = 120
  };

  // definitions from ALICE internal notes ALICE-INT-2002-010 and
  // ALICE-INT-2006-XXX
  enum {
    /** the 32bit word in the CDH containing the status flags */
    kCDHStatusWord=4, //! see above
    /** start of the flags in word fgkCDHStatusWord */
    kCDHStatusFlagsOffset=12, //! see above
    /** bit indicating HLT decision in the HLTOUT*/
    kCDHFlagsHLTDecision=6, //! see above
    /** bit indicating HLT payload in the HLTOUT*/
    kCDHFlagsHLTPayload=7 //! see above
  };

  class AliHLTOUTHandlerListEntry;
  /**
   * Block descriptor.
   */
  class AliHLTOUTBlockDescriptor {
  public:
    AliHLTOUTBlockDescriptor(AliHLTComponentDataType dt, AliHLTUInt32_t spec, AliHLTUInt32_t index, AliHLTOUT* pCollection)
      : fDataType(dt), fSpecification(spec), fIndex(index), fSelected(false), fProcessed(false), fpCollection(pCollection) {};
    AliHLTOUTBlockDescriptor(const AliHLTOUTBlockDescriptor& src)
      : fDataType(src.fDataType), fSpecification(src.fSpecification), fIndex(src.fIndex), fSelected(false), fProcessed(false), fpCollection(src.fpCollection) {}
    AliHLTOUTBlockDescriptor& operator=(const AliHLTOUTBlockDescriptor& src)
    { fDataType=src.fDataType; fSpecification=src.fSpecification; fIndex=src.fIndex; fSelected=false; fProcessed=false; fpCollection=src.fpCollection; return *this; }
    ~AliHLTOUTBlockDescriptor() {}

    operator AliHLTComponentDataType() const {return fDataType;}
    operator AliHLTUInt32_t() const {return fSpecification;}
    int operator==(AliHLTComponentDataType dt) const {return dt==fDataType;}
    int operator==(AliHLTUInt32_t spec) const {return spec==fSpecification;}
    int operator==(AliHLTOUT* collection) const {return collection==fpCollection;}

    AliHLTUInt32_t GetIndex() const {return fIndex;}
    AliHLTOUT* GetCollection() const {return fpCollection;}

    bool IsSelected() const {return fSelected;}
    void Select(bool selected=true) {fSelected=selected;}
    bool IsProcessed() const {return fProcessed;}
    void MarkProcessed() {fProcessed=true;}

    /**
     * Get the data buffer
     * @param [out] pBuffer buffer of the selected data block
     * @param [out] size    size of the selected data block
     */
    int GetDataBuffer(const AliHLTUInt8_t* &pBuffer, AliHLTUInt32_t& size) {
      if (fpCollection) return fpCollection->GetDataBuffer(GetIndex(), pBuffer, size);
      return -ENODEV;
    }

    /**
     * Get the handler descriptor
     */
    const AliHLTOUTHandlerListEntry& GetHandlerDesc();

    /// print info
    void Print(const char* option="") const;

  private:      
    /** data type of the block */
    AliHLTComponentDataType fDataType; //!transient
    /** data specification of the block */
    AliHLTUInt32_t          fSpecification; //!transient
    /** index in the data stream */
    AliHLTUInt32_t          fIndex; //!transient
    /** selection flag */
    bool                    fSelected; //!transient
    /** processed flag */
    bool                    fProcessed; //!transient
    /** the collection */
    AliHLTOUT*              fpCollection; //!transient
  };

  class AliHLTOUTHandlerListEntry {
  public:
    AliHLTOUTHandlerListEntry(AliHLTOUTHandler* pHandler, 
			      AliHLTModuleAgent::AliHLTOUTHandlerDesc& handlerDesc,
			      AliHLTModuleAgent* pAgent,
			      AliHLTUInt32_t index);

    /** copy constructor for vector handling */
    AliHLTOUTHandlerListEntry(const AliHLTOUTHandlerListEntry& src);

    /** assignment operator for vector handling */
    AliHLTOUTHandlerListEntry& operator=(const AliHLTOUTHandlerListEntry& src);

    ~AliHLTOUTHandlerListEntry();

    static const AliHLTOUTHandlerListEntry fgkVoidHandlerListEntry; //! initializer

    operator AliHLTOUTHandler*() const {return fpHandler;}

    // please note that fpHandlerDesc is really a pointer and is created
    // in the constructor. Thats why it is dereferenced here. The pointer
    // type is on purpose, even though it is a bit confusing with the 
    // argument by reference in the AliHLTOUTHandlerListEntry constructor.
    operator const AliHLTModuleAgent::AliHLTOUTHandlerDesc&() const 
    {return fpHandlerDesc?*fpHandlerDesc:AliHLTModuleAgent::fgkVoidHandlerDesc;}
    operator AliHLTModuleAgent*() const {return fpAgent;}

    /**
     * Two list entries are considered to be equal if the handlers
     * are equal.
     */
    bool operator==(const AliHLTOUTHandlerListEntry& entry) const;

    bool operator==(const AliHLTModuleAgent::AliHLTOUTHandlerType handlerType) const;

    /**
     * Compare the handler descriptor of this list entry with another
     * descriptor.
     * @return true if both handler and agent match
     */
    bool operator==(const AliHLTModuleAgent::AliHLTOUTHandlerDesc desc) const; 

    AliHLTUInt32_t operator[](int i) const;

    /**
     * Add a block index to this descriptor.
     * One descriptor can serve multiple blocks if the agent returns the
     * same handler for all of the blocks. Instead of creating a new entry
     * the block index ist just added
     */
    void AddIndex(AliHLTUInt32_t index);

    /**
     * Add all indices of the descriptor.
     */
    void AddIndex(AliHLTOUTHandlerListEntry &desc);

    /**
     * Check if an index is served by this descriptor.
     * @return true if the index is in the table
     */
    bool HasIndex(AliHLTUInt32_t index) const;

    /**
     * Invalidate all block indices
     */
    void InvalidateBlocks() {fBlocks.clear();}

    /**
     * Check whether the entry has valid blocks.
     */
    bool IsEmpty() {return fBlocks.size()==0;}

    /// print info
    void Print(const char* option="") const;

  private:
    /** standard constructor prohibited */
    AliHLTOUTHandlerListEntry();

    /** pointer to the handler */
    AliHLTOUTHandler* fpHandler; //! transient

    /** pointer to handler description */
    AliHLTModuleAgent::AliHLTOUTHandlerDesc* fpHandlerDesc; //! transient

    /** pointer to module agent */
    AliHLTModuleAgent* fpAgent; //! transient

    /** list of block indices */
    AliHLTOUTIndexList fBlocks; //!transient
  };

  typedef vector<AliHLTOUTHandlerListEntry> AliHLTOUTHandlerListEntryVector;
  typedef vector<AliHLTOUTBlockDescriptor>  AliHLTOUTBlockDescriptorVector;
  typedef vector<AliHLTOUT*>                AliHLTOUTPVector;

  /**
   * Selection guard for the AliHLTOUT object.
   * If the object is locked, the selection of data blocks can not be changed.
   */
  class AliHLTOUTSelectionGuard {
  public:
    /** constructor */
    AliHLTOUTSelectionGuard(AliHLTOUT* pInstance) : fpInstance(pInstance)
    {if (fpInstance) fpInstance->SelectDataBlock();}
    /** constructor */
    AliHLTOUTSelectionGuard(AliHLTOUT* pInstance, const AliHLTOUTHandlerListEntry* pHandlerDesc) : fpInstance(pInstance)
    {if (fpInstance) fpInstance->SelectDataBlocks(pHandlerDesc);}
    /** destructor */
    ~AliHLTOUTSelectionGuard()
    {if (fpInstance) fpInstance->DisableBlockSelection();}

  private:
    /** standard constructor prohibited */
    AliHLTOUTSelectionGuard();
    /** copy constructor prohibited */
    AliHLTOUTSelectionGuard(const AliHLTOUTSelectionGuard&);
    /** assignment operator prohibited */
    AliHLTOUTSelectionGuard& operator=(const AliHLTOUTSelectionGuard&);

    /** the AliHLTOUT instance the guard is locking */
    AliHLTOUT* fpInstance; //!transient
  };

  /**
   * Init for processing.
   * The HLTOUT is scanned for all available data blocks and the
   * AliHLTOUTHandler objects for the data blocks are created according
   * to the module agents (see AliHLTModuleAgent).
   */
  int Init();

  /**
   * Reset and clear all data block descriptors.
   * @note Since sub-collections are only referred via the block
   * descriptors, all information on sub-collections is also lost.
   * Child classes should implement ResetInput in order to cleanup
   * the input devices.
   */
  int Reset();

  /**
   * Get the current event id
   */
  AliHLTUInt64_t EventId() {return fCurrentEventId;}

  /**
   * Get number of data blocks in the HLTOUT data
   */
  int GetNofDataBlocks();

  /**
   * Select the first data block of a certain data type and specification.
   * The selection criteria can be of kAliHLTAnyDataType and/or
   * kAliHLTVoidDataSpec in order to be ignored and just match any data block.
   *
   * The search criteria can be combined with a handler type (e.g. kRawReader)
   * @param [in] dt     data type to match                                <br>
   * @param [in] spec   data specification to match                       <br>
   * @param [in] handlerType  type of the handler
   * @param [in] skipProcessed skip all block marked processed
   * @return block index: >= 0 if success, -ENOENT if no block found      <br>
   *         neg. error code if failed                                    <br>
   *                        -EPERM if access denied (object locked)
   */
  int SelectFirstDataBlock(AliHLTComponentDataType dt=kAliHLTAnyDataType,
			   AliHLTUInt32_t spec=kAliHLTVoidDataSpec,
			   AliHLTModuleAgent::AliHLTOUTHandlerType handlerType=AliHLTModuleAgent::kUnknownOutput,
			   bool skipProcessed=true);

  /**
   * Select the next data block of data type and specification of the previous
   * call to @ref SelectFirstDataBlock.
   * @return block index: >= 0 if success, -ENOENT if no block found      <br>
   *         neg. error code if failed                                    <br>
   *                        -EPERM if access denied (object locked)
   */
  int SelectNextDataBlock();

  /**
   * Get properties of the selected data block.
   * @param [out] dt    data type of the selected block
   * @param [out] spec  data specification of the selected block
   */
  int GetDataBlockDescription(AliHLTComponentDataType& dt, AliHLTUInt32_t& spec);

  /**
   * Get handler description of the current data block.
   */
  const AliHLTOUTHandlerListEntry& GetDataBlockHandlerDesc();

  /**
   * Get handler type of the selected data block.
   * @return handler type for the selected data block
   */
  AliHLTModuleAgent::AliHLTOUTHandlerType GetDataBlockHandlerType();

  /**
   * Get the index of the current data block.
   * @return index, AliHLTOUTInvalidIndex if no block selected
   */
  AliHLTUInt32_t GetDataBlockIndex();

  /**
   * Get buffer of the selected data block.
   * @param [out] pBuffer buffer of the selected data block
   * @param [out] size    size of the selected data block
   */
  int GetDataBuffer(const AliHLTUInt8_t* &pBuffer, AliHLTUInt32_t& size);

  /**
   * Release buffer after use.
   * @param [in] pBuffer  buffer of the selected data block
   */
  int ReleaseDataBuffer(const AliHLTUInt8_t* pBuffer);

  /**
   * Get a TObject from the data buffer
   * @return TObject pointer if data block is a streamed object
   */
  TObject* GetDataObject();

  /**
   * Release data object
   */
  int ReleaseDataObject(TObject* pObject);

  /**
   * Add the current data block to the selection.
   * Note: enables also the block selection
   */
  int SelectDataBlock();

  /**
   * Add the all data blocks of a certain handler to the selection.
   * Note: enables also the block selection
   */
  int SelectDataBlocks(const AliHLTOUTHandlerListEntry* pHandlerDesc);

  /**
   * Enable the selection of data blocks.
   */
  int EnableBlockSelection();

  /**
   * Disable the selection of data blocks.
   */
  int DisableBlockSelection();

  /**
   * Reset the data block selection.
   * Resets the selection list, none of the blocks is selected.
   */
  int ResetBlockSelection();

  /**
   * Mark the current block as processed.
   */
  int MarkDataBlockProcessed();

  /**
   * Mark all data blocks of a certain handler processed.
   */
  int MarkDataBlocksProcessed(const AliHLTOUTHandlerListEntry* pHandlerDesc);

  /**
   * Add a sub collection to the HLTOUT.
   */
  int AddSubCollection(AliHLTOUT* pCollection);

  /**
   * Release a previously added sub collection.
   */
  int ReleaseSubCollection(AliHLTOUT* pCollection);

  /**
   * Get module agent for the selected data block.
   */
  AliHLTModuleAgent* GetAgent();

  /**
   * Get handler for the selected data block.
   */
  AliHLTOUTHandler* GetHandler();

  /**
   * Convert data buffer to ESD.
   * The buffer is supposed to describe a streamed AliESDEvent object.
   * If no target object is specified, the ESD is written to a file AliHLTdetESDs.root,
   * where 'det' is derived from the data type origin. Each time the function is invoked
   * a new event is created. Dummy events are added if the previous events did not contain
   *
   * The function needs AliRoot and might not be implemented by all AliHLTOUT
   * implementations.
   * a data block of this specification.
   * @param [in]  pBuffer  the data buffer
   * @param [in]  size     data buffer size
   * @param [in]  dt       data type of the block
   * @param [out] tgtesd   optional target
   */
  virtual int WriteESD(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size, AliHLTComponentDataType dt, AliESDEvent* tgtesd=NULL) const;

  enum AliHLTOUTByteOrder {
    /** no data block selected */
    kInvalidByteOrder=-1,
    kUnknownByteOrder=0,
    kLittleEndian,
    kBigEndian
  };

  enum AliHLTOUTDataType {
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
  AliHLTOUTByteOrder CheckByteOrder();

  /**
   * Check alignment of selected block
   */
  int CheckAlignment(AliHLTOUT::AliHLTOUTDataType type);

  /**
   * Helper function to byte swap a 64 bit value.
   */
  static AliHLTUInt64_t ByteSwap64(AliHLTUInt64_t src);

  /**
   * Helper function to byte swap a 32 bit value.
   */
  static AliHLTUInt32_t ByteSwap32(AliHLTUInt32_t src);

  /**
   * Insert a handler item.
   * The current list entries are checked if the handler is already in
   * the list. It is added if not in the list, otherwise the block index
   * is added to the existing entry.
   * @param list     the handler list
   * @param entry    handler list entry
   * @return 0 if added, EEXIST (non negative!) if merged with existing entry <br>
   *         neg. error code if failed
   */
  static int InsertHandler(AliHLTOUTHandlerListEntryVector& list, const AliHLTOUTHandlerListEntry &entry);
  
  /**
   * Insert all handlers of the specified type to the list.
   */
  int FillHandlerList(AliHLTOUTHandlerListEntryVector& list, AliHLTModuleAgent::AliHLTOUTHandlerType handlerType);

  /**
   * Remove empty items which have a duplicate in the list.
   */
  static int RemoveEmptyDuplicateHandlers(AliHLTOUTHandlerListEntryVector& list);

  /**
   * Find an entry of a certain description in the list.
   * @return index of the entry if found, -ENOENT if not found
   */
  static int FindHandler(AliHLTOUTHandlerListEntryVector& list, const AliHLTModuleAgent::AliHLTOUTHandlerDesc desc);

  /**
   * Invalidate all blocks of the handlers in the list
   */
  static int InvalidateBlocks(AliHLTOUTHandlerListEntryVector& list);

 protected:
  /**
   * Add a block descriptor.
   * This is done by the child classes generating the index. The AliHLTOUT
   * object must be locked for index generation.
   * @param desc    the block descriptor
   * @return 0 if success, -EPERM if access denied
   */
  int AddBlockDescriptor(const AliHLTOUTBlockDescriptor desc);

  /**
   * Set the event id, only for child classes
   */
  void SetEventId(AliHLTUInt64_t id);

  /**
   * Print output or suppress.
   */
  bool BeVerbose() {return fbVerbose;}

  /**
   * Switch output.
   */
  void SwitchVerbosity(bool verbose) {fbVerbose=verbose;}

  /// print info
  void Print(const char* option="") const;

 private:
  /** copy constructor prohibited */
  AliHLTOUT(const AliHLTOUT&);
  /** assignment operator prohibited */
  AliHLTOUT& operator=(const AliHLTOUT&);

  /**
   * Internal status flags
   */
  enum {
    /** the HLTOUT object is locked with the current data block selection */
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
    kAlignmentWarning = 0x20,
    /** enable block selection list */
    kBlockSelection = 0x40,
    /** skip processed data blocks */
    kSkipProcessed = 0x80,
    /** marked as sub collection */
    kIsSubCollection = 0x100
  };

  /**
   * Generate the index of the HLTOUT data.
   * Must be implemented by the child classes.
   */
  virtual int GenerateIndex()=0;

  /**
   * Find AliHLTOUTHandler objects for the data blocks.
   * The available AliHLTModuleAgents are probed whether they provide
   * handlers for data processing.
   */
  int InitHandlers();

  /**
   * Cleanup and reset the data input.
   */
  virtual int ResetInput();

  /**
   * Get the data buffer
   * @param [in]  index   index of the block
   * @param [out] pBuffer buffer of the selected data block
   * @param [out] size    size of the selected data block
   */
  virtual int GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
			    AliHLTUInt32_t& size)=0;

  /**
   * Check byte order of data block
   */
  virtual AliHLTOUTByteOrder CheckBlockByteOrder(AliHLTUInt32_t index)=0;

  /**
   * Check alignment of data block
   */
  virtual int CheckBlockAlignment(AliHLTUInt32_t index, AliHLTOUT::AliHLTOUTDataType type)=0;

  /**
   * Select the data block of data type and specification of the previous
   * call to @ref SelectFirstDataBlock. Core function of @ref SelectFirstDataBlock
   * and @ref SelectNextDataBlock, starts to find a block at the current list
   * position.
   * 
   * The data block is searched from the conditions of fSearchDataType,
   * fSearchSpecification, fSearchHandlerType and the selection list.
   *
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

  /**
   * Find handler description for a certain block index.
   */
  const AliHLTOUTHandlerListEntry& FindHandlerDesc(AliHLTUInt32_t blockIndex);

  /**
   * Set the RawReader as parameter.
   * The function is for internal use only in conjunction with the
   * New() functions.
   */
  virtual void SetParam(AliRawReader* pRawReader);

  /**
   * Set the RunLoader as parameter
   * The function is for internal use only in conjunction with the
   * New() functions.
   */
  virtual void SetParam(TTree* pDigitTree, int event=-1);

  /**
   * Set name of the digit file as parameter
   * The function is for internal use only in conjunction with the
   * New() functions.
   */
  virtual void SetParam(const char* filename, int event=-1);

  /** data type for the current block search, set from @ref SelectFirstDataBlock */
  AliHLTComponentDataType fSearchDataType; //!transient

  /** data specification for the current block search */
  AliHLTUInt32_t fSearchSpecification; //!transient

  /** handler type for the current block search */
  AliHLTModuleAgent::AliHLTOUTHandlerType fSearchHandlerType; // !transient

  /** instance flags: locked, collecting, ... */
  unsigned int fFlags; //!transient

  /** list of block descriptors */
  AliHLTOUTBlockDescriptorVector fBlockDescList; //!transient

  /** current position in the list */
  unsigned int fCurrent; //!transient

  /** data buffer under processing */
  const AliHLTUInt8_t* fpBuffer; //!transient

  /** list of AliHLTOUTHandlers */
  AliHLTOUTHandlerListEntryVector fDataHandlers; // !transient

  /** verbose or silent output */
  bool fbVerbose; //!transient

  /** gobal instance set for certain steps of the analysis */
  static AliHLTOUT* fgGlobalInstance; //! transient

  /** logging methods */
  AliHLTLogging fLog; //! transient

  /** current buffer converted to a TObject */
  TObject* fpDataObject; //!
  const AliHLTUInt8_t* fpObjectBuffer; //!
  AliHLTUInt32_t fObjectBufferSize; //!

  /** current event id */
  AliHLTUInt64_t fCurrentEventId; //!

  ClassDef(AliHLTOUT, 6)
};
#endif
