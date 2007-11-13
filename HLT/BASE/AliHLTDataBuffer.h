// @(#) $Id$

#ifndef ALIHLTDATABUFFER_H
#define ALIHLTDATABUFFER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTDataBuffer.h
    @author Matthias Richter
    @date   
    @brief  Handling of Data Buffers for HLT components.
    @note   The class is used in Offline (AliRoot) context
*/

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
   * @param arrayBlockDesc  pointer to block descriptor to be filled
   * @param iArraySize      size of the block descriptor array
   * @return: number of matching data blocks, neg. error code if failed<br>
   *          -EACCESS      the consumer state can't be changed (activated)
   *          -EBADF        unresolved data segments <br>
   *          -ENOENT       consumer component not found <br>
   *          -ENODATA      data buffer does not have raw data <br>
   *          -EINVAL       invalid parameter <br>
   */
  int Subscribe(const AliHLTComponent* pConsumer,
		AliHLTComponentBlockData* arrayBlockDesc,
		int iArraySize);

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
   * Get the number of segments
   * @return number of segments
   */
  int GetNofSegments();

  /**
   * Get the total number of consumers.
   * This gives the number of consumers regardless of their state.
   * @return number of consumers
   */
  int GetNofConsumers();

  /**
   * Get the number of consumers which still need to be processed during
   * the current event.
   * @return number of consumers
   */
  int GetNofPendingConsumers();

  /**
   * Get the number of consumers currently under processing.
   * @return number of active consumers
   */
  int GetNofActiveConsumers();

  /**
   * Check if a consumer is already in the list
   * @param pConsumer   pointer to consumer component
   * @param bAllLists   search in all lists if 1
   *                    search only in fConsumer list if 0
   * @return 1 if found, 0 if not
   */
  int FindConsumer(AliHLTComponent* pConsumer, int bAllLists=1);

  /**
   * Public method to reset the buffer.
   * Eventually with some additional checks. In normal operation,
   * an external reset should not be necessary.
   */
  int Reset();

  /**
   * Set local logging level
   * logging filter for individual object
   */
  void SetLocalLoggingLevel(AliHLTComponentLogSeverity level)
  {fgLogging.SetLocalLoggingLevel(level); AliHLTLogging::SetLocalLoggingLevel(level);}

  /**
   * @class AliHLTDataSegment
   * @brief  Descriptor of a data segment within the buffer.
   */
  class AliHLTDataSegment {
    friend class AliHLTDataBuffer;
    friend class AliHLTConsumerDescriptor;
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
      fPtr((AliHLTUInt8_t*)ptr),
      fSegmentOffset(offset),
      fSegmentSize(size),
      fSpecification(0)
    {
    }

    AliHLTDataSegment(void* ptr, AliHLTUInt32_t offset, AliHLTUInt32_t size, AliHLTComponentDataType dt, AliHLTUInt32_t spec)
      :
      fDataType(dt),
      fPtr((AliHLTUInt8_t*)ptr),
      fSegmentOffset(offset),
      fSegmentSize(size),
      fSpecification(spec)
    {
    }

    AliHLTUInt8_t* GetPtr() const {return (AliHLTUInt8_t*)*this;}

    AliHLTUInt32_t GetSize() const {return fSegmentSize;}
    
    int operator==(const AliHLTDataSegment& seg) const
    {
      return (fPtr+fSegmentOffset==seg.fPtr+seg.fSegmentOffset) && (fSegmentSize==seg.fSegmentSize);
    }
    operator AliHLTUInt8_t*() const {return fPtr+fSegmentOffset;}

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
   * @class AliHLTRawBuffer
   * @brief  Descriptor of the raw data buffer which can host several segments.
   */
  class AliHLTRawBuffer {
  public:
    /** standard constructor */
    AliHLTRawBuffer() : fSize(0), fTotalSize(0), fPtr(NULL) {}
    /** constructor */
    AliHLTRawBuffer(AliHLTUInt32_t size);
    /** destructor */
    virtual ~AliHLTRawBuffer();

    /**
     * Use a fraction of the buffer.
     * @param size    size in bytes to be used
     * @return pointer to buffer
     */
    AliHLTUInt8_t* UseBuffer(AliHLTUInt32_t size);

    /**
     * Check whether buffer fits for a request.
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

    int operator==(void*) const;
    int operator==(AliHLTUInt8_t* ptr) const {return fPtr==ptr;}
    int operator<=(void*) const;
    int operator>(void*) const;
    int operator-(void*) const;

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
    /** the buffer */
    AliHLTUInt8_t* fPtr;                                           //! transient
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

  /**
   * Reset the data buffer.
   * Removes all consumers back to the @ref fConsumers list, deletes
   * segments and releases the Raw Buffer.
   */
  int ResetDataBuffer();

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

  ClassDef(AliHLTDataBuffer, 0)
};

#endif // ALIHLTDATABUFFER_H
