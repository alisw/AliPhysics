//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTDATABUFFER_H
#define ALIHLTDATABUFFER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTDataBuffer.h
//  @author Matthias Richter
//  @date   
//  @brief  Handling of Data Buffers for HLT components.
//  @note   The class is used in Offline (AliRoot) context

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <vector>
#include "TObject.h"
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"

class AliHLTConsumerDescriptor;
class AliHLTTask;

/** list of AliHLTConsumerDescriptor pointers */
typedef vector<AliHLTConsumerDescriptor*> AliHLTConsumerDescriptorPList;

typedef AliHLTUInt8_t* AliHLTUInt8Pointer_t;

/**
 * @class AliHLTDataBuffer
 * @brief  Handling of data buffers for the HLT.
 * 
 * The class provides handling of data buffers for HLT tasks. Each task gets
 * its own Data Buffer instance. The buffer is grouped into different data
 * segments according to the output of the component.<br>
 * The Data Buffer keeps control over the data requests of the 'child'
 * components. Each component can subscribe to a certain segment of the data
 * buffer. It's state is then changed from 'reserved' to 'active'. After the
 * data processing, the component has to release the segment and it's state is
 * set to 'processed'. If all components have requested and released their data,
 * the Raw Buffer is released and pushed back in the list of available buffers.
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTDataBuffer : public TObject, public AliHLTLogging 
{
 public:
  //////////////////////////////////////////////////////////////////////////////
  // constructors and destructors

  /* standard constructor
   */
  AliHLTDataBuffer();
  /** destructor */
  virtual ~AliHLTDataBuffer();

  //////////////////////////////////////////////////////////////////////////////
  // initialization

  /**
   * Add component to the list of consumers
   * @param pConsumer - a consumer of type AliHLTComponent
   */
  int SetConsumer(AliHLTComponent* pConsumer);

  //////////////////////////////////////////////////////////////////////////////
  // component to component communication

  /**
   * Determine the number of matching data blocks for the component and a
   * consumer component. <br>
   * The first approach will support only one output data type for processing
   * components.
   * @param pConsumer       the component which subscribes to the buffer
   * @param tgtList         (optional) the list to receive the data types
   * @return: number of data blocks which match the input data types 
   *          of the consumer, neg. error code if failed <br>
   *          -EINVAL       invalid parameter <br>
   */
  int FindMatchingDataBlocks(const AliHLTComponent* pConsumer,
			     AliHLTComponentDataTypeList* tgtList=NULL);

  /**
   * Subscribe to a segment of the data buffer.
   * The function prepares the block descriptor for subsequent use with the
   * AliHLTComponent::ProcessEvent method, the method can prepare several block
   * descriptors up to the array size specified by iArraySize. The return value
   * is independent from the array size the number of block descriptors which
   * would have been prepared if there was enough space in the array<br>
   * The method is used by the consumer component.
   * @param pConsumer       the component which subscribes to the buffer
   * @param blockDescList   block descriptor vector to be filled
   * @return: number of matching data blocks, neg. error code if failed<br>
   *          -EACCESS      the consumer state can't be changed (activated)
   *          -EBADF        unresolved data segments <br>
   *          -ENOENT       consumer component not found <br>
   *          -ENODATA      data buffer does not have raw data <br>
   *          -EINVAL       invalid parameter <br>
   */
  int Subscribe(const AliHLTComponent* pConsumer,
		AliHLTComponentBlockDataList& blockDescList);

  /**
   * Release an instance of the data buffer.
   * Resets the variables of the block descriptor.
   * If all buffer segments are released, the Data Buffer is reseted
   * and the Raw Buffer released.<br>
   * The method is used by the consumer component.
   * @param pBlockDesc      descriptor of the data segment
   * @param pConsumer       the component which subscribes to the buffer
   * @param pOwnerTask      task owning this buffer
   * @return: >0 if success, negative error code if failed <br>
   *          -EACCESS      the consumer state can not be changed (de-activated)
   *          -ENOENT       consumer has not subscribed to the buffer <br>
   *          -EINVAL       invalid parameter <br>
   */
  int Release(AliHLTComponentBlockData* pBlockDesc, const AliHLTComponent* pConsumer,
	      const AliHLTTask* pOwnerTask);

  /**
   * Release a forwarded data block.
   */
  int ReleaseForwardedBlock(AliHLTComponentBlockData* pBlockDesc,
			     const AliHLTTask* pOwnerTask);

  /**
   * Register an input data block for forwarding.
   * Consumer of this data buffer subscribe to forwarded data blocks in te same way.
   * Forwarded data blocks are released when the last consumer has released the
   * blocks.
   * @param pSrcTask        original source task of the data block
   * @param pBlockDesc      descriptor of the data segment
   */
  int Forward(AliHLTTask* pSrcTask, AliHLTComponentBlockData* pBlockDesc);

  /**
   * Get a target buffer of minimum size iMinSize.
   * The method is used by the component which owns the Data Buffer to 
   * allocate a buffer for the data it is going to produce.
   * @param iMinSize        minumum size of the requested buffer
   * @return: pointer to target buffer if 
   */
  AliHLTUInt8_t* GetTargetBuffer(int iMinSize);

  /**
   * Set the segments for the data buffer.
   * This is usually done after the component has written the data to the buffer, 
   * which was requested by the @ref GetTargetBuffer method. The component might
   * produce different types of data, for each type a segment has to be defined
   * which describes the data inside the buffer.<br>
   * The @ref AliHLTComponentBlockData segment descriptor comes directly from
   * the @ref AliHLTComponent::ProcessEvent method.
   * @param pTgt            the target buffer which the segments refer to
   * @param arraySegments   the output block descriptors of the component
   * @param iSize           size of the array
   */
  int SetSegments(AliHLTUInt8_t* pTgt, AliHLTComponentBlockData* arraySegments, int iSize);

  /**
   * Check if the data buffer is empty.
   * @return 1 if empty, 0 if not
   */
  int IsEmpty();

  /**
   * Get the total and maximum size of the buffer.
   * Lets see if this is needed later
   */
  //int GetTotalSize();

  /**
   * Get the number of segments including the forwarded data blocks.
   * @return number of segments
   */
  int GetNofSegments() const;

  /**
   * Get the total number of consumers.
   * This gives the number of consumers regardless of their state.
   * @return number of consumers
   */
  int GetNofConsumers() const;

  /**
   * Get the number of consumers which still need to be processed during
   * the current event.
   * @return number of consumers
   */
  int GetNofPendingConsumers() const;

  /**
   * Get the number of consumers currently under processing.
   * @return number of active consumers
   */
  int GetNofActiveConsumers() const;

  /**
   * Check if a consumer is already in the list
   * @param pConsumer   pointer to consumer component
   * @param bAllLists   search in all lists if 1
   *                    search only in fConsumer list if 0
   * @return 1 if found, 0 if not
   */
  int FindConsumer(const AliHLTComponent* pConsumer, int bAllLists=1);

  /**
   * Public method to reset the buffer.
   * Eventually with some additional checks. In normal operation,
   * an external reset should not be necessary.
   */
  int Reset();

  /**
   * Print info about the buffer
   */
  virtual void Print(const char* option) const;

  /**
   * Set local logging level
   * logging filter for individual object
   */
  void SetLocalLoggingLevel(AliHLTComponentLogSeverity level)
  {fgLogging.SetLocalLoggingLevel(level); AliHLTLogging::SetLocalLoggingLevel(level);}

  /**
   * Print summary of the global buffer management.
   */
  static int PrintStatistics();

  /**
   * Set the global event count.
   * The event count is deployed to find buffers which have not been used
   * for a while. In such a case to policy to find an appropriate buffer is
   * adjusted.
   */
  static int SetGlobalEventCount(AliHLTUInt32_t eventCount) {fgEventCount=eventCount; return 0;}

  /**
   * @class AliHLTDataSegment
   * @brief  Descriptor of a data segment within the buffer.
   */
  class AliHLTDataSegment {
    friend class AliHLTDataBuffer; // TODO: implement some getters/setters
  public:
    AliHLTDataSegment()
      :
      fDataType(kAliHLTVoidDataType),
      fPtr(NULL),
      fSegmentOffset(0),
      fSegmentSize(0),
      fSpecification(0)
    {
    }

    AliHLTDataSegment(AliHLTUInt8_t* ptr, AliHLTUInt32_t offset, AliHLTUInt32_t size) 
      :
      fDataType(kAliHLTVoidDataType),
      fPtr(ptr),
      fSegmentOffset(offset),
      fSegmentSize(size),
      fSpecification(0)
    {
    }

    AliHLTDataSegment(void* ptr, AliHLTUInt32_t offset, AliHLTUInt32_t size) 
      :
      fDataType(kAliHLTVoidDataType),
      fPtr(reinterpret_cast<AliHLTUInt8_t*>(ptr)),
      fSegmentOffset(offset),
      fSegmentSize(size),
      fSpecification(0)
    {
    }

    AliHLTDataSegment(void* ptr, AliHLTUInt32_t offset, AliHLTUInt32_t size, AliHLTComponentDataType dt, AliHLTUInt32_t spec)
      :
      fDataType(dt),
      fPtr(reinterpret_cast<AliHLTUInt8_t*>(ptr)),
      fSegmentOffset(offset),
      fSegmentSize(size),
      fSpecification(spec)
    {
    }

    AliHLTDataSegment(const AliHLTDataSegment& src)
      :
      fDataType(src.fDataType),
      fPtr(src.fPtr),
      fSegmentOffset(src.fSegmentOffset),
      fSegmentSize(src.fSegmentSize),
      fSpecification(src.fSpecification)
    {
      // AliHLTDataSegment just stores external pointers and properties
    }

    AliHLTDataSegment& operator=(const AliHLTDataSegment& src)
    {
      // AliHLTDataSegment just stores external pointers and properties
      if (this==&src) return *this;
      fDataType=src.fDataType;
      fPtr=src.fPtr;
      fSegmentOffset=src.fSegmentOffset;
      fSegmentSize=src.fSegmentSize;
      fSpecification=src.fSpecification;
      return *this;
    }

    virtual ~AliHLTDataSegment() {}

    AliHLTUInt8_t* GetPtr() const {return (AliHLTUInt8_t*)*this;}

    AliHLTUInt32_t GetSize() const {return fSegmentSize;}
    
    int operator==(const AliHLTDataSegment& seg) const
    {
      return (fPtr+fSegmentOffset==seg.fPtr+seg.fSegmentOffset) && (fSegmentSize==seg.fSegmentSize);
    }
    operator AliHLTUInt8_t*() const {return fPtr+fSegmentOffset;}

    virtual void Print(const char* option) const;

  private:
    /** the data type of this segment */
    AliHLTComponentDataType fDataType;                             // see above
    /** pointer to the buffer */
    AliHLTUInt8Pointer_t fPtr;                                     //!transient
    /** offset in byte within the data buffer */
    AliHLTUInt32_t fSegmentOffset;                                 // see above
    /** size of the actual content */
    AliHLTUInt32_t fSegmentSize;                                   // see above
    /** data specification */
    AliHLTUInt32_t fSpecification;                                 // see above

  };

  /**
   * @class AliHLTForwardedDataSegment
   * @brief  Descriptor of a forwarded data segment.
   * Contains in addition information about the parent of this forwarded
   * block and the original data type and specification
   */
  class AliHLTForwardedDataSegment : public AliHLTDataSegment {
    friend class AliHLTDataBuffer; // TODO: implement some getters/setters
  public:
    AliHLTForwardedDataSegment()
      : AliHLTDataSegment()
      , fParentSegment()
      , fParentTask(NULL)
    {
    }

    AliHLTForwardedDataSegment(AliHLTDataSegment& mySegment, AliHLTDataSegment& parentSegment, AliHLTTask* parentTask)
      : AliHLTDataSegment(mySegment)
      , fParentSegment(parentSegment)
      , fParentTask(parentTask)
    {
    }

    AliHLTForwardedDataSegment(const AliHLTForwardedDataSegment& src)
      : AliHLTDataSegment(src),
      fParentSegment(src.fParentSegment),
      fParentTask(src.fParentTask)
    {
      // AliHLTForwardedDataSegment just stores external pointers and properties
    }

    AliHLTForwardedDataSegment& operator=(const AliHLTForwardedDataSegment& src)
    {
      // AliHLTForwardedDataSegment just stores external pointers and properties
      AliHLTDataSegment::operator=(src);
      fParentSegment=src.fParentSegment;
      fParentTask=src.fParentTask;
      return *this;
    }

    virtual ~AliHLTForwardedDataSegment() {}

    virtual void Print(const char* option) const;

  private:
    /// description of the original segment
    AliHLTDataSegment fParentSegment;                              // see above
    /// the parent task
    AliHLTTask* fParentTask;                                       //!transient
  };

  class AliHLTRawBuffer;
  typedef vector<AliHLTRawBuffer*>  AliHLTRawBufferPList;

  /**
   * @class AliHLTRawPage
   * Memory allocation is organized in pages of a fixed size. Within a
   * page, AliHLTRawBuffer chunks are created.
   */
  class AliHLTRawPage : public AliHLTLogging {
  public:
    /** standard constructor */
  AliHLTRawPage() : fSize(0), fPtr(NULL), fFreeBuffers(), fUsedBuffers() {}
    /** constructor */
    AliHLTRawPage(AliHLTUInt32_t pagesize);
    /** destructor */
    virtual ~AliHLTRawPage();

    /** alloc a buffer of specified size from the global pages*/
    static AliHLTRawBuffer* GlobalAlloc(AliHLTUInt32_t size, int verbosity=0);
    /** find buffer in the global pages */
    static AliHLTRawPage* FindPage(AliHLTRawBuffer* buffer);
    /** cleanup the global pages */
    static int GlobalClean();
    /** adjust global page size */
    static void SetGlobalPageSize(AliHLTUInt32_t size) {fgGlobalPageSize=size;}
    /** find next page after prev, or first page */
    static AliHLTRawPage* NextPage(const AliHLTRawPage* prev=NULL);

    /** alloc a buffer of specified size */
    AliHLTRawBuffer* Alloc(AliHLTUInt32_t size);
    /** free a buffer and merge consecutive free buffers */
    int Free(AliHLTRawBuffer* pBuffer);
    /** set the size of a raw buffer and release the remaining part */
    int SetSize(const AliHLTRawBuffer* pBuffer, AliHLTUInt32_t size);
    /// check if the buffer is in this page
    bool HasBuffer(const AliHLTRawBuffer* pBuffer);

    AliHLTUInt32_t Size() const {return fSize;}
    AliHLTUInt32_t Capacity() const;
    bool IsUsed() const {return fUsedBuffers.size()>0;}
    bool IsFragmented() const {return (fFreeBuffers.size()+fUsedBuffers.size())>1;}

    /**
     * Print page information
     */
    virtual void Print(const char* option);

  private:
    /** copy constructor prohibited */
    AliHLTRawPage(const AliHLTRawPage&);
    /** assignment operator prohibited */
    AliHLTRawPage& operator=(const AliHLTRawPage&);

    /// list of global pages
    static vector<AliHLTDataBuffer::AliHLTRawPage*> fgGlobalPages; //! transient
    /// pages size of global pages
    static AliHLTUInt32_t fgGlobalPageSize;                        //! transient

    /** page size */
    AliHLTUInt32_t fSize;                                          // see above
    /** the memory segment */
    AliHLTUInt8_t* fPtr;                                           //! transient

    /** list of free buffers */
    AliHLTRawBufferPList fFreeBuffers;                             //! transient
    /** list of used buffers */
    AliHLTRawBufferPList fUsedBuffers;                             //! transient
  };

  /**
   * @class AliHLTRawBuffer
   * @brief  Descriptor of the raw data buffer which can host several segments.
   */
  class AliHLTRawBuffer {
  public:
    /** standard constructor */
  AliHLTRawBuffer() : fSize(0), fTotalSize(0), fExternalPtr(NULL), fPtr(NULL), fLastEventCount(0) {}
    /** constructor */
    AliHLTRawBuffer(AliHLTUInt32_t size);
    /** constructor */
    AliHLTRawBuffer(AliHLTUInt32_t size, AliHLTUInt8_t* buffer);
    /** destructor */
    virtual ~AliHLTRawBuffer();

    /**
     * Use a fraction of the buffer.
     * @param size    size in bytes to be used
     * @return pointer to buffer
     */
    AliHLTUInt8_t* UseBuffer(AliHLTUInt32_t size);

    /**
     * split a buffer at specified size
     * only possible for buffers with external memory
     */
    AliHLTRawBuffer* Split(AliHLTUInt32_t size);

    /**
     * Check whether buffer fits for a request.
     * A buffer fits if it is at least of the requested size and at most
     * the requested size plus a margin. The margin increases with the
     * number of events the buffer has not been used.
     * @param size    size of the request in bytes
     * @return 1 if buffer is big enough, 0 if not
     */
    int CheckSize(AliHLTUInt32_t size) const;

    /**
     * Get used size of the buffer
     */
    AliHLTUInt32_t GetUsedSize() const {return fSize;}

    /**
     * Get total size of the buffer
     */
    AliHLTUInt32_t GetTotalSize() const {return fTotalSize;}

    /**
     * Get pointer of data buffer
     */
    AliHLTUInt8_t* GetPointer() const {return fPtr;}

    /**
     * Write check pattern
     */
    int WritePattern(const char* pattern, int size);

    /**
     * Check pattern
     */
    int CheckPattern(const char* pattern, int size) const;

    /**
     * Reset buffer.
     * Data buffer remains allocated, used size set to 0
     */
    int Reset();

    /*
     * Merge buffer with succeeding buffer.
     * Only possible if the buffers are consecutive with out any gap.
     */
    int Merge(const AliHLTRawBuffer& succ);

    /**
     * Print buffer information
     */
    virtual void Print(const char* option) const;

    int operator==(void* ptr) const;
    int operator==(AliHLTUInt8_t* ptr) const {return fPtr==ptr;}
    int operator<(void* ptr) const;
    int operator<=(void* ptr) const;
    int operator>(void* ptr) const;
    int operator-(void* ptr) const;
    int operator<(const AliHLTRawBuffer& op) const;
    int operator<=(const AliHLTRawBuffer& op) const;
    int operator>(const AliHLTRawBuffer& op) const;

    operator void*() const {return fPtr;}
    operator AliHLTUInt8_t*() const {return fPtr;}

  private:
    /** copy constructor prohibited */
    AliHLTRawBuffer(const AliHLTRawBuffer&);
    /** assignment operator prohibited */
    AliHLTRawBuffer& operator=(const AliHLTRawBuffer&);

    /** size of the currently occupied partition of the buffer */
    AliHLTUInt32_t fSize;                                          // see above
    /** total size of the buffer, including safety margin */
    AliHLTUInt32_t fTotalSize;                                     // see above
    /** optional external buffer */
    AliHLTUInt8_t* fExternalPtr;                                   //! transient
    /** the buffer, external or allocated */
    AliHLTUInt8_t* fPtr;                                           //! transient
    /** last event count where the buffer has been used */
    AliHLTUInt32_t fLastEventCount;                                //! transient
  };

 private:
  /** copy constructor prohibited */
  AliHLTDataBuffer(const AliHLTDataBuffer&);
  /** assignment operator prohibited */
  AliHLTDataBuffer& operator=(const AliHLTDataBuffer&);

  /* lets see if this is needed
     AliHLTDataSegment* FindDataSegment(AliHLTComponentDataType datatype);
  */

  /**
   * Find those data segments which match the input types of a component.
   * @param pConsumer       the component which subscribes to the buffer
   * @param tgtList         the list to receive the data segment descriptors
   * @return: number of data blocks which match the input data types 
   *          of the consumer, neg. error code if failed <br>
   *          -EINVAL       invalid parameter <br>
   */
  int FindMatchingDataSegments(const AliHLTComponent* pConsumer, 
			       vector<AliHLTDataBuffer::AliHLTDataSegment>& tgtList);

 protected:
  // 2010-02-01 make function protected in order to be used from unit test
  /**
   * Reset the data buffer.
   * Removes all consumers back to the @ref fConsumers list, deletes
   * segments and releases the Raw Buffer.
   */
  int ResetDataBuffer();
 private:

  //////////////////////////////////////////////////////////////////////////////

  // the data description

  // the data segments within this buffer
  vector<AliHLTDataSegment> fSegments;                             // see above

  // the list of all consumers which are going to subscribe to the buffer
  AliHLTConsumerDescriptorPList fConsumers;                         // see above
  // the list of all consumers which are currently subscribed to the buffer
  AliHLTConsumerDescriptorPList fActiveConsumers;                   // see above
  // the list of all consumers which are already released for the current event
  AliHLTConsumerDescriptorPList fReleasedConsumers;                 // see above

  // the buffer instance
  AliHLTRawBuffer* fpBuffer;                                       //! transient

  // flags indicating the state of the buffer
  AliHLTUInt32_t fFlags;                                           // see above

  /** list of tasks with forwarded data blocks */
  vector<AliHLTTask*> fForwardedSegmentSources;                    //! transient

  /** list of forwarded block descriptors */
  vector<AliHLTDataSegment> fForwardedSegments;                    //! transient

  //////////////////////////////////////////////////////////////////////////////
  // global buffer handling, internal use only

  /**
   * Create a raw buffer of a certain size.
   * The function tries to find a buffer of the given size (or a bit bigger by a 
   * certain margin @ref fgMargin) from the list of free buffers.
   * If no buffer is available, a new one is created and added to the buffer handling.
   * @param size            min. size of the requested buffer
   * @return pointer to raw buffer
   */
  static AliHLTRawBuffer* CreateRawBuffer(AliHLTUInt32_t size);

  /**
   * Set the data size of a raw buffer after it has been filled by
   * the component.
   */
  int SetRawBufferDataSize(AliHLTRawBuffer* pBuffer, AliHLTUInt32_t size) const;

  /**
   * Mark a buffer as free.
   * After the Data Buffer has finnished using the raw buffer, it is released
   * and added to the list of available buffers.
   * @param pBuffer         the raw buffer to release
   * @return >=0 if succeeded, neg. error code if failed
   */
  static int ReleaseRawBuffer(AliHLTRawBuffer* pBuffer);

  /**
   * Deletes all the raw buffers.
   * When the last Data Buffer object is destructed, all raw data buffers are
   * relesed.
   */
  static int DeleteRawBuffers();

  /**
   * Number of instances of AliHLTDataBuffer.
   * The statice variable is incremented and decremented in the constructor/
   * destructor. All internal data structures are cleaned up when the last
   * instance is exiting.
   */
  static int fgNofInstances;                                       // see above
  /** global list of free raw buffers */
  static vector<AliHLTRawBuffer*> fgFreeBuffers;                   // see above
  /** global list of currently active raw buffers */
  static vector<AliHLTRawBuffer*> fgActiveBuffers;                 // see above
  /** determines the raw buffer size margin at buffer requests */
  static AliHLTUInt32_t fgMargin;                                  // see above

  /** global instance to HLT logging class for static methods */
  static AliHLTLogging fgLogging;                                  // see above

  /** size of the safety pattern */
  static const Int_t fgkSafetyPatternSize;                         // see above

  /** the safety pattern */
  static const char fgkSafetyPattern[];                            //!transient

  static AliHLTUInt32_t fgEventCount;                              //!transient

  //////////////////////////////////////////////////////////////////////////////
  // internal helper functions

  /**
   * Find the consumer descriptor for a certain component and data type in 
   * a list of consumers.<br>
   * <b>Note:</b> There are three lists which contain the consumers in the
   * different states.
   * @param pConsumer       pointer to consumer component
   * @param list            list where to search for the consumer
   */
  AliHLTConsumerDescriptor* FindConsumer(const AliHLTComponent* pConsumer,
					 AliHLTConsumerDescriptorPList &list) const;

  /**
   * Change the state of a consumer.
   * The state of a consumer is determined by the list it is strored in, the
   * method moves a consumer from the source to the target list.
   * @param pDesc           pointer to consumer descriptor
   * @param srcList         list where the consumer is currently to be found
   * @param tgtList         list where to move the consumer
   */
  int ChangeConsumerState(AliHLTConsumerDescriptor* pDesc,
			  AliHLTConsumerDescriptorPList &srcList,
			  AliHLTConsumerDescriptorPList &tgtList);

  /**
   * Cleanup a consumer list.
   * Release all allocated data structures. <b>Note:</b> Not the component itself!
   */
  int CleanupConsumerList();

  ClassDef(AliHLTDataBuffer, 1)
};

#endif // ALIHLTDATABUFFER_H
