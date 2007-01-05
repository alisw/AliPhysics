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

#include <vector>
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTDefinitions.h"
#include "TObject.h"

/* Matthias Dec 2006
 * The names have been changed for Aliroot's coding conventions sake
 * The old names are defined for backward compatibility with the 
 * stand alone SampleLib package
 */
typedef AliHLTComponentLogSeverity AliHLTComponent_LogSeverity;
typedef AliHLTComponentEventData AliHLTComponent_EventData;
typedef AliHLTComponentShmData AliHLTComponent_ShmData;
typedef AliHLTComponentDataType AliHLTComponent_DataType;
typedef AliHLTComponentBlockData AliHLTComponent_BlockData;
typedef AliHLTComponentTriggerData AliHLTComponent_TriggerData;
typedef AliHLTComponentEventDoneData AliHLTComponent_EventDoneData;

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
   * @param outputBlockCnt  out: size of the output block array, set by the component
   * @param outputBlocks    out: the output block array is allocated internally
   * @param edd
   * @return neg. error code if failed
   */
  virtual int ProcessEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
			    AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
			    AliHLTUInt32_t& size, AliHLTUInt32_t& outputBlockCnt, 
			    AliHLTComponentBlockData*& outputBlocks,
			    AliHLTComponentEventDoneData*& edd ) = 0;

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
  virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& ) = 0;

  /**
   * Get the output data type of the component.
   * The function is pure virtual and must be implemented by the child class.
   * @return output data type
   */
  virtual AliHLTComponentDataType GetOutputDataType() = 0;

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
  int FindMatchingDataTypes(AliHLTComponent* pConsumer, vector<AliHLTComponentDataType>* tgtList);
 
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
   * Fill AliHLTComponentBlockData structure with default values.
   * @param blockData   reference to data structure
   */
  void FillBlockData( AliHLTComponentBlockData& blockData );

  /**
   * Fill AliHLTComponentShmData structure with default values.
   * @param shmData   reference to data structure
   */
  void FillShmData( AliHLTComponentShmData& shmData );

  /**
   * Fill AliHLTComponentDataType structure with default values.
   * @param dataType   reference to data structure
   */
  void FillDataType( AliHLTComponentDataType& dataType );
  
  /**
   * Copy data type structure
   * Copies the value an AliHLTComponentDataType structure to another one
   * @param[out] tgtdt   target structure
   * @param[in] srcdt   source structure
   */
  void CopyDataType(AliHLTComponentDataType& tgtdt, const AliHLTComponentDataType& srcdt);

  /**
   * Set the ID and Origin of an AliHLTComponentDataType structure.
   * The function sets the fStructureSize member and copies the strings
   * to the ID and Origin. Only characters from the valid part of the string
   * are copied, the rest is fille with 0's.
   * Please note that the fID and fOrigin members are not strings, just arrays of
   * chars of size @ref kAliHLTComponentDataTypefIDsize and
   * @ref kAliHLTComponentDataTypefOriginSize respectively and not necessarily with
   * a terminating zero.
   * @param id      ID string
   * @param origin  Origin string
   */
  void SetDataType(AliHLTComponentDataType& tgtdt, const char* id, const char* origin);

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
  int MakeOutputDataBlockList( const vector<AliHLTComponentBlockData>& blocks, AliHLTUInt32_t* blockCount,
			       AliHLTComponentBlockData** outputBlocks );

  /**
   * Fill the EventDoneData structure.
   * The method is redirected to a function provided by the current
   * framework. Function pointers are transferred via the @ref
   * AliHLTComponentEnvironment structure.
   */
  int GetEventDoneData( unsigned long size, AliHLTComponentEventDoneData** edd );

  /**
   * Helper function to convert the data type to a string.
   */
  void DataType2Text( const AliHLTComponentDataType& type, char output[14] );

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
