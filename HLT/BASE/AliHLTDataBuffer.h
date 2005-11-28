// @(#) $Id$

#ifndef ALIHLTDATABUFFER_H
#define ALIHLTDATABUFFER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTDataBuffer
   handling of data buffers for the HLT
 */

#include <cerrno>
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTDefinitions.h"
#include "AliHLTComponent.h"
#include "TObject.h"
#include "TList.h"

/* internal data structure
 */
struct AliHLTDataSegment {
  AliHLTComponent_DataType fDataType; // the data type of this buffer
  Int_t fSegmentOffset;               // offset in byte within the data buffer
  Int_t fSegmentSize;                 // size of the actual content
  AliHLTUInt32_t fSpecification;      // data specification
};

/* internal data structure
 */
struct AliHLTRawBuffer {
  AliHLTUInt32_t fSize;                        // size of the buffer
  AliHLTUInt32_t fTotalSize;                   // total size of the buffer, including safety margin
  void* fPtr;                         // the buffer
};

/* internal data structure
 * there is unfortunately no unique determination of the data type from the component
 * itself possible, thats way both component and data type have to be initialized
 * and are stored in a compound
 */
class AliHLTConsumerDescriptor : public TObject, public AliHLTLogging {
 private:
  AliHLTComponent* fpConsumer;
  AliHLTComponent_DataType fDataType;
  AliHLTDataSegment* fpSegment;

 public:
  AliHLTConsumerDescriptor();
  AliHLTConsumerDescriptor(AliHLTComponent* pConsumer, AliHLTComponent_DataType datatype);
  ~AliHLTConsumerDescriptor();

  AliHLTComponent* GetComponent() {return fpConsumer;}
  AliHLTComponent_DataType GetDataType() {return fDataType;}

  int SetActiveDataSegment(AliHLTDataSegment* pSegment);
  int CheckActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size);
  int ReleaseActiveDataSegment();
};

class AliHLTDataBuffer : public AliHLTLogging, public TObject {
 public:
  AliHLTDataBuffer();
  virtual ~AliHLTDataBuffer();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // initialization

  /* add component to the list of consumers
   * parameter:
   *   pConsumer - a consumer of type AliHLTComponent
   *   datatype - data type of the segement, the consumer is registered for
   */
  int SetConsumer(AliHLTComponent* pConsumer, AliHLTComponent_DataType datatype);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // component to component communication

  /* subscribe to a segment of the data buffer
   * the function prepares the block descriptor for subsequent use with the AliHLTComponent::ProcessEvent
   * method
   * parameter:
   *   datatype - type of the data segment
   *   pConsumer - the component which subscribes to the buffer
   *   pBlockDesc - pointer to receive the prepared block descriptor
   * return: >0 if success, negative error code if failed
   */
  int Subscribe(AliHLTComponent_DataType datatype, const AliHLTComponent* pConsumer, AliHLTComponent_BlockData* pBlockDesc);

  /* release an instance of the data buffer
   * resets the variables of the block descriptor
   * parameter:
   *   pBlockDesc - descriptor of the data segment
   *   pConsumer - the component which subscribes to the buffer
   * return: >0 if success, negative error code if failed
   */
  int Release(AliHLTComponent_BlockData* pBlockDesc, const AliHLTComponent* pConsumer);

  /* get a target buffer if minimum size iMinSize
   */
  AliHLTUInt8_t* GetTargetBuffer(int iMinSize);

  /* set the segments for the data buffer
   * this is usually done after the component has written the data to the buffer
   * parameter:
   *   pTgt - the target buffer the segments refer to
   *   arraySegments - the output block descriptors of the component
   *   iSize - size of the array
   */
  int SetSegments(AliHLTUInt8_t* pTgt, AliHLTComponent_BlockData* arraySegments, int iSize);

  /* check if the data buffer is empty
   */
  int IsEmpty();

  /* get the total and maximum size of the buffer
   * lets see if this is needed later
   */
  //int GetTotalSize();

  /* get the number of segments
   */
  int GetNofSegments();

  /* get the number of consumers
   */
  int GetNofConsumers();

  /* get the number of consumers
   */
  int GetNofActiveConsumers();

 private:
  AliHLTDataSegment* FindDataSegment(AliHLTComponent_DataType datatype);

  /* reset the data buffer
   * removes all condumers back to the fConsumers list
   */
  int ResetDataBuffer();


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // the data description
  vector<AliHLTDataSegment> fSegments;// the data segments within this buffer

  vector<AliHLTConsumerDescriptor*> fConsumers;                   // the list of all consumers which are going to subscribe to the buffer
  vector<AliHLTConsumerDescriptor*> fActiveConsumers;             // the list of all consumers which are currently subscribed to the buffer
  vector<AliHLTConsumerDescriptor*> fReleasedConsumers;           // the list of all consumers which are already released for the current event

  AliHLTRawBuffer* fpBuffer;           // the buffer instance

  AliHLTUInt32_t fFlags;                // flags indicating the state of the buffer

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // global buffer handling

  /* create a raw buffer of a certain size
   * the function tries to find a buffer of the given size (or a little bit bigger) from the list of free buffers
   * if no buffer is available, a new one is created
   */
  static AliHLTRawBuffer* CreateRawBuffer(AliHLTUInt32_t size);

  /* mark a buffer as free
   */
  static int ReleaseRawBuffer(AliHLTRawBuffer* pBuffer);

  /* deletes all the raw buffers
   */
  static int DeleteRawBuffers();

  static int fNofInstances;
  static vector<AliHLTRawBuffer*> fFreeBuffers;
  static vector<AliHLTRawBuffer*> fActiveBuffers;
  static AliHLTUInt32_t fMargin;

  /*
   */
  AliHLTConsumerDescriptor* FindConsumer(const AliHLTComponent* pConsumer, AliHLTComponent_DataType datatype, vector<AliHLTConsumerDescriptor*> &pList);

  int ChangeConsumerState(AliHLTConsumerDescriptor* pDesc, vector<AliHLTConsumerDescriptor*> &srcList, vector<AliHLTConsumerDescriptor*> &tgtList);

  int CleanupConsumerList();

  ClassDef(AliHLTDataBuffer, 0)
};
#endif // ALIHLTDATABUFFER_H
