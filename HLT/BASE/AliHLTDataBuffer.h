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

#include <cerrno>
#include <vector>
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTDefinitions.h"
#include "TObject.h"
//#include "TList.h"

class AliHLTComponent;
/* @name internal data structures
 */

/**
 * @struct AliHLTDataSegment
 * @brief  Descriptor of a data segment within the buffer.
 * @ingroup alihlt_system
 */
struct AliHLTDataSegment {
  AliHLTDataSegment()
    :
    fDataType(),
    fSegmentOffset(0),
    fSegmentSize(0),
    fSpecification(0)
  {
    memset(&fDataType, 0, sizeof(AliHLTComponentDataType));
  }
  AliHLTDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size) 
    :
    fDataType(),
    fSegmentOffset(offset),
    fSegmentSize(size),
    fSpecification(0)
  {
    memset(&fDataType, 0, sizeof(AliHLTComponentDataType));
  }
  /** the data type of this segment */
  AliHLTComponentDataType fDataType;                               // see above
  /** offset in byte within the data buffer */
  AliHLTUInt32_t fSegmentOffset;                                   // see above
  /** size of the actual content */
  AliHLTUInt32_t fSegmentSize;                                     // see above
  /** data specification */
  AliHLTUInt32_t fSpecification;                                   // see above
};

/**
 * @struct AliHLTRawBuffer
 * @brief  Descriptor of the raw data buffer which can host several segments.
 * @ingroup alihlt_system
 */
struct AliHLTRawBuffer {
  /** size of the currently occupied partition of the buffer */
  AliHLTUInt32_t fSize;                                            // see above
  /** total size of the buffer, including safety margin */
  AliHLTUInt32_t fTotalSize;                                       // see above
  /** the buffer */
  void* fPtr;                                                      //! transient
};

/**
 * @class AliHLTConsumerDescriptor
 * @brief Helper class to describe a consumer component.
 *
 * There is unfortunately no unique determination of the data type from the
 * component itself possible, thats why both component and data type have to
 * be initialized and are stored in a compound. The class is intended to make
 * bookkeeping easier.
 *
 * @note This class is only used for the @ref alihlt_system.
 *
 * @ingroup alihlt_system
 */
class AliHLTConsumerDescriptor : public TObject, public AliHLTLogging {
 public:
  /** standard constructur */
  AliHLTConsumerDescriptor();
  /** constructur 
   * @param pConsumer pointer to the consumer component
   */
  AliHLTConsumerDescriptor(AliHLTComponent* pConsumer);
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTConsumerDescriptor(const AliHLTConsumerDescriptor&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTConsumerDescriptor& operator=(const AliHLTConsumerDescriptor&);
  /** destructor */
  ~AliHLTConsumerDescriptor();

  /**
   * Get the component of this descriptor.
   * @return pointer to the component
   */
  AliHLTComponent* GetComponent() {return fpConsumer;}

  /**
   * Set an active data segment.
   * the pointer will be handled in a container, no allocation, copy or
   * cleanup.
   * @param offset  offset of the segment in the buffer
   * @param size    size of the segment in the buffer
   * @return >=0 if succeeded
   */
  int SetActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size);

  /**
   * Check whether there is an active data segment of certain size with
   * certain offset.
   * @param offset  offset of the data segment in the data buffer
   * @param size    size of the data segment in the data buffer
   * @return > if existend, 0 if not
   */
  int CheckActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size);

  /** find an active data segment of certain size with certain offset
   * will see if this is necessary
   * @param offset  offset of the data segment in the data buffer
   * @param size    size of the data segment in the data buffer
   * @return offset of the data segment
   */
  //AliHLTUInt32_t FindActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size);

  /** get the number of active segments for this consumer
   * @return number of active segments
   */
  int GetNofActiveSegments() {return fSegments.size();};

  /**
   */
  int ReleaseActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size);

 private:
  /** consumer object */
  AliHLTComponent* fpConsumer;                                     //! transient

  /** list of data segments */
  vector<AliHLTDataSegment> fSegments;                             // see above

  //ClassDef(AliHLTConsumerDescriptor, 0)
};

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
class AliHLTDataBuffer : public TObject, public AliHLTLogging {
 public:
  //////////////////////////////////////////////////////////////////////////////
  // condtructors and destructors

  /* standard constructor
   */
  AliHLTDataBuffer();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTDataBuffer(const AliHLTDataBuffer&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTDataBuffer& operator=(const AliHLTDataBuffer&);
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
			     vector<AliHLTComponentDataType>* tgtList=NULL);

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
   * @return: >0 if success, negative error code if failed <br>
   *          -EACCESS      the consumer state can not be changed (de-activated)
   *          -ENOENT       consumer has not subscribed to the buffer <br>
   *          -EINVAL       invalid parameter <br>
   */
  int Release(AliHLTComponentBlockData* pBlockDesc, const AliHLTComponent* pConsumer);

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
   * Get the number of consumers
   * @return number of consumers
   */
  int GetNofConsumers();

  /**
   * Get the number of active consumers
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

 private:
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
  int FindMatchingDataSegments(const AliHLTComponent* pConsumer, vector<AliHLTDataSegment>& tgtList);

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
  vector<AliHLTConsumerDescriptor*> fConsumers;                    // see above
  // the list of all consumers which are currently subscribed to the buffer
  vector<AliHLTConsumerDescriptor*> fActiveConsumers;              // see above
  // the list of all consumers which are already released for the current event
  vector<AliHLTConsumerDescriptor*> fReleasedConsumers;            // see above

  // the buffer instance
  AliHLTRawBuffer* fpBuffer;                                       //! transient

  // flags indicating the state of the buffer
  AliHLTUInt32_t fFlags;                                           // see above

  //////////////////////////////////////////////////////////////////////////////
  // global buffer handling, internal use only

  /**
   * Create a raw buffer of a certain size.
   * The function tries to find a buffer of the given size (or a bit bigger by a 
   * certain margin @ref fMargin) from the list of free buffers.
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
					 vector<AliHLTConsumerDescriptor*> &list) const;

  /**
   * Change the state of a consumer.
   * The state of a consumer is determined by the list it is strored in, the
   * method moves a consumer from the source to the target list.
   * @param pDesc           pointer to consumer descriptor
   * @param srcList         list where the consumer is currently to be found
   * @param tgtList         list where to move the consumer
   */
  int ChangeConsumerState(AliHLTConsumerDescriptor* pDesc,
			  vector<AliHLTConsumerDescriptor*> &srcList,
			  vector<AliHLTConsumerDescriptor*> &tgtList);

  /**
   * Cleanup a consumer list.
   * Release all allocated data structures. <b>Note:</b> Not the component itself!
   */
  int CleanupConsumerList();

  ClassDef(AliHLTDataBuffer, 0)
};
#endif // ALIHLTDATABUFFER_H
