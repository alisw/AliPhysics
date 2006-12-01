// @(#) $Id$

#ifndef ALIHLTCOMPONENT_H
#define ALIHLTCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTComponent.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Base class declaration for HLT components. 
    @note   The class is both used in Online (PubSub) and Offline (AliRoot)
            context
                                                                          */
/**
 * @defgroup alihlt_component Component handling of the HLT module
 * This section describes the the component handling for the HLT module.
 */
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTDefinitions.h"
#include "TObject.h"

class AliHLTComponentHandler;

/**
 * @class AliHLTComponent
 * Base class of HLT data processing components.
 * The class provides a common interface for HLT data processing components.
 * The interface can be accessed from the online HLT framework or the AliRoot
 * offline analysis framework.
 * Components can be of type 
 * - @ref kSource:    components which only produce data 
 * - @ref kProcessor: components which consume and produce data
 * - @ref kSink:      components which only consume data
 *
 * where data production and consumption refer to the analysis data stream.<br>
 *
 * In order to adapt to different environments (on-line/off-line), the component
 * gets an environment structure with function pointers. The base class provides
 * member functions for those environment dependend functions. The member 
 * functions are used by the component implementation and are re-mapped to the
 * corresponding functions.
 * @ingroup alihlt_component
 */
class AliHLTComponent : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTComponent();
  /** standard destructor */
  virtual ~AliHLTComponent();

  /** component type definitions */
  enum TComponentType { kUnknown=0, kSource=1, kProcessor=2, kSink=3 };

  /**
   * Init function to prepare data processing.
   * Initialization of common data structures for a sequence of events.
   * The call is redirected to the internal method @ref DoInit which can be
   * overridden by the child class.<br>
   * During Init also the environment structure is passed to the component.
   * @param environ        environment pointer with environment dependend function
   *                       calls
   * @param environ_param  additionel parameter for function calls, the pointer
   *                       is passed as it is
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   */
  virtual int Init( AliHLTComponentEnvironment* environ, void* environ_param, int argc, const char** argv );

  /**
   * Clean-up function to terminate data processing.
   * Clean-up of common data structures after data processing.
   * The call is redirected to the internal method @ref DoDeinit which can be
   * overridden by the child class.
   */
  virtual int Deinit();

  /**
   * Processing of one event.
   * The method is pure virtual and implemented by the child classes 
   * - @ref AliHLTProcessor
   * - @ref AliHLTDataSource
   * - @ref AliHLTDataSink
   *
   * @param evtData
   * @param blocks
   * @param trigData
   * @param outputPtr
   * @param size
   * @param outputBlockCnt
   * @param outputBlocks
   * @param edd
   * @return neg. error code if failed
   */
  virtual int ProcessEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
			    AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
			    AliHLTUInt32_t& size, AliHLTUInt32_t& outputBlockCnt, 
			    AliHLTComponent_BlockData*& outputBlocks,
			    AliHLTComponent_EventDoneData*& edd ) = 0;

  // Information member functions for registration.

  /**
   * Get the type of the component.
   * The function is pure virtual and must be implemented by the child class.
   * @return component type id
   */
  virtual TComponentType GetComponentType() = 0; // Source, sink, or processor

  /**
   * Get the id of the component.
   * Each component is identified by a unique id.
   * The function is pure virtual and must be implemented by the child class.
   * @return component id (string)
   */
  virtual const char* GetComponentID() = 0;

  /**
   * Get the input data types of the component.
   * The function is pure virtual and must be implemented by the child class.
   * @return list of data types in the vector reference
   */
  virtual void GetInputDataTypes( vector<AliHLTComponent_DataType>& ) = 0;

  /**
   * Get the output data type of the component.
   * The function is pure virtual and must be implemented by the child class.
   * @return output data type
   */
  virtual AliHLTComponent_DataType GetOutputDataType() = 0;

  /**
   * Get a ratio by how much the data volume is shrinked or enhanced.
   * The function is pure virtual and must be implemented by the child class.
   * @param constBase        <i>return</i>: additive part, independent of the
   *                                   input data volume  
   * @param inputMultiplier  <i>return</i>: multiplication ratio
   * @return values in the reference variables
   */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) = 0;

  /**
   * Spawn function.
   * Each component must implement a spawn function to create a new instance of 
   * the class. Basically the function must return <i>new <b>my_class_name</b></i>.
   * @return new class instance
   */
  virtual AliHLTComponent* Spawn() = 0;

  /**
   * Find matching data types between this component and a consumer component.
   * Currently, a component can produce only one type of data. This restriction is most
   * likely to be abolished in the future.
   * @param pConsumer a component and consumer of the data produced by this component
   * @param tgtList   reference to a vector list to receive the matching data types.
   * @return >= 0 success, neg. error code if failed
   */ 
  int FindMatchingDataTypes(AliHLTComponent* pConsumer, vector<AliHLTComponent_DataType>* tgtList);
 
  /**
   * Set the global component handler.
   * The static method is needed for the automatic registration of components. 
   */
  static int SetGlobalComponentHandler(AliHLTComponentHandler* pCH, int bOverwrite=0);

  /**
   * Clear the global component handler.
   * The static method is needed for the automatic registration of components. 
   */
  static int UnsetGlobalComponentHandler();

 protected:

  /**
   * Fill AliHLTComponent_BlockData structure with default values.
   * @param blockData   reference to data structure
   */
  void FillBlockData( AliHLTComponent_BlockData& blockData ) {
    blockData.fStructSize = sizeof(blockData);
    FillShmData( blockData.fShmKey );
    blockData.fOffset = ~(AliHLTUInt32_t)0;
    blockData.fPtr = NULL;
    blockData.fSize = 0;
    FillDataType( blockData.fDataType );
    blockData.fSpecification = ~(AliHLTUInt32_t)0;
  }

  /**
   * Fill AliHLTComponent_ShmData structure with default values.
   * @param shmData   reference to data structure
   */
  void FillShmData( AliHLTComponent_ShmData& shmData ) {
    shmData.fStructSize = sizeof(shmData);
    shmData.fShmType = gkAliHLTComponent_InvalidShmType;
    shmData.fShmID = gkAliHLTComponent_InvalidShmID;
  }

  /**
   * Fill AliHLTComponent_DataType structure with default values.
   * @param dataType   reference to data structure
   */
  void FillDataType( AliHLTComponent_DataType& dataType ) {
    dataType.fStructSize = sizeof(dataType);
    memset( dataType.fID, '*', 8 );
    memset( dataType.fOrigin, '*', 4 );
  }
  
  /**
   * Default method for the internal initialization.
   * The method is called by @ref Init
   */
  virtual int DoInit( int argc, const char** argv );

  /**
   * Default method for the internal clean-up.
   * The method is called by @ref Deinit
   */
  virtual int DoDeinit();

  /**
   * General memory allocation method.
   * All memory which is going to be used 'outside' of the interface must
   * be provided by the framework (online or offline).
   * The method is redirected to a function provided by the current
   * framework. Function pointers are transferred via the @ref
   * AliHLTComponentEnvironment structure.
   */
  void* AllocMemory( unsigned long size );

  /**
   * Helper function to create a monolithic BlockData description block out
   * of a list BlockData descriptors.
   * For convenience, inside the interface vector lists are used, to make the
   * interface pure C style, monilithic blocks must be exchanged. 
   * The method is redirected to a function provided by the current
   * framework. Function pointers are transferred via the @ref
   * AliHLTComponentEnvironment structure.
   */
  int MakeOutputDataBlockList( const vector<AliHLTComponent_BlockData>& blocks, AliHLTUInt32_t* blockCount,
			       AliHLTComponent_BlockData** outputBlocks );

  /**
   * Fill the EventDoneData structure.
   * The method is redirected to a function provided by the current
   * framework. Function pointers are transferred via the @ref
   * AliHLTComponentEnvironment structure.
   */
  int GetEventDoneData( unsigned long size, AliHLTComponent_EventDoneData** edd );

  /**
   * Helper function to convert the data type to a string.
   */
  void DataType2Text( const AliHLTComponent_DataType& type, char output[14] );

 private:
  /** The global component handler instance */
  static AliHLTComponentHandler* fpComponentHandler;
  /** The environment where the component is running in */
  AliHLTComponentEnvironment fEnvironment;

  /** 
   * Set by ProcessEvent before the processing starts (e.g. before calling 
   * @ref AliHLTProcessor::DoEvent)
   */
  AliHLTEventID_t fCurrentEvent;

  ClassDef(AliHLTComponent, 0)
};
#endif
