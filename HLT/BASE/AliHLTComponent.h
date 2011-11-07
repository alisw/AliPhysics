//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTCOMPONENT_H
#define ALIHLTCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTComponent.h
//  @author Matthias Richter, Timm Steinbeck
//  @date   
//  @brief  Base class declaration for HLT components. 
//  @note   The class is both used in Online (PubSub) and Offline (AliRoot)
//          context


/**
 * @defgroup alihlt_component Component handling of the HLT module
 * This section describes the the component base classes and handling for
 * the HLT module.
 *
 * @section alihlt_component_intro General remarks
 * HLT analysis is organized in so called components. Each component can
 * subscribe to the data produced by other components and can from the
 * analysis publish new data for the subsequent components. Only the
 * input data blocks and entries from CDB are available for the analysis. 
 *
 * @section alihlt_component_implementation Component implementation
 * AliHLTComponent provides the interface for all components, see there
 * for details. Three types are provided:
 * - AliHLTProcessor
 * - AliHLTDataSource
 * - AliHLTDataSink
 *
 * The two last represent data sinks and sources for the HLT integration
 * into AliRoot. When running only, only the processors are relevant,
 * sources and sinks are provided by the HLT PubSub framework. Please check
 * AliHLTComponent for detailed description.
 *
 * @section alihlt_component_registration Component registration
 * Components need to be registered with the AliHLTComponentHandler in
 * order to be used with the system. Registration is purely done from the
 * module library. Two methods are possible:
 * - the module library implements an AliHLTModuleAgent and overloads the
 *   AliHLTModuleAgent::RegisterComponents() function
 * - in the implementation file, one object is defined. The global object is
 *   automatically instantiated when the library is loaded for the first
 *   time and the object is used for registration.
 *
 * In both cases, the library must be loaded via the method
 * <pre>
 *  AliHLTComponentHandler::LoadComponentLibraries()
 * </pre>
 * For the global object approach it is important that the library is
 * not loaded elsewhere before (e.g. a gSystem->Load operation in your
 * rootlogon.C).
 *
 *
 */

#include <vector>
#include <string>
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTCommonCDBEntries.h"

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
class TObjArray;
class TMap;
class TStopwatch;
class TUUID;
class AliRawDataHeader;
class AliHLTComponent;
class AliHLTMemoryFile;
class AliHLTCTPData;
class AliHLTReadoutList;

/** list of component data type structures */
typedef vector<AliHLTComponentDataType>   AliHLTComponentDataTypeList;
/** list of component block data structures */
typedef vector<AliHLTComponentBlockData>  AliHLTComponentBlockDataList;
/** list of component statistics struct */
typedef vector<AliHLTComponentStatistics> AliHLTComponentStatisticsList;
/** list of component pointers */
typedef vector<AliHLTComponent*>          AliHLTComponentPList;
/** list of memory file pointers */
typedef vector<AliHLTMemoryFile*>         AliHLTMemoryFilePList;

/**
 * @class AliHLTComponent
 * Base class of HLT data processing components.
 * The class provides a common interface for HLT data processing components.
 * The interface can be accessed from the online HLT framework or the AliRoot
 * offline analysis framework.
 * @section alihltcomponent-properties Component identification and properties
 * Each component must provide a unique ID, input and output data type indications,
 * and a spawn function.
 * @subsection alihltcomponent-req-methods Required property methods
 * - @ref GetComponentID
 * - @ref GetInputDataTypes (see @ref alihltcomponent-type for default
 *   implementations.)
 * - @ref GetOutputDataType (see @ref alihltcomponent-type for default
 *   implementations.)
 * - @ref GetOutputDataSize (see @ref alihltcomponent-type for default
 *   implementations.)
 * - @ref Spawn
 *
 * @subsection alihltcomponent-opt-mehods Optional handlers
 * - @ref DoInit
 * - @ref DoDeinit
 * - @ref GetOutputDataTypes
 *   If the component has multiple output data types @ref GetOutputDataType
 *   should return @ref kAliHLTMultipleDataType. The framework will invoke
 *   @ref GetOutputDataTypes, a list can be filled.
 * - @ref Reconfigure
 *   This function is invoked by the framework on a special event which
 *   triggers the reconfiguration of the component.
 *
 * @subsection alihltcomponent-processing-mehods Data processing
 * 
 * 
 * @subsection alihltcomponent-type Component type
 * Components can be of type
 * - @ref kSource     components which only produce data 
 * - @ref kProcessor  components which consume and produce data
 * - @ref kSink       components which only consume data
 *
 * where data production and consumption refer to the analysis data stream. In
 * order to indicate the type, a child component can overload the
 * @ref GetComponentType function.
 * @subsubsection alihltcomponent-type-std Standard implementations
 * Components in general do not need to implement this function, standard
 * implementations of the 3 types are available:
 * - AliHLTDataSource for components of type @ref kSource <br>
 *   All types of data sources can inherit from AliHLTDataSource and must
 *   implement the @ref AliHLTDataSource::GetEvent method. The class
 *   also implements a standard method for @ref GetInputDataTypes.
 *   
 * - AliHLTProcessor for components of type @ref kProcessor <br>
 *   All types of data processors can inherit from AliHLTProcessor and must
 *   implement the @ref AliHLTProcessor::DoEvent method.
 *
 * - AliHLTDataSink for components of type @ref kSink <br>
 *   All types of data processors can inherit from AliHLTDataSink and must
 *   implement the @ref AliHLTDataSink::DumpEvent method. The class
 *   also implements a standard method for @ref GetOutputDataType and @ref
 *   GetOutputDataSize.
 *
 * @subsection alihltcomponent-environment Running environment
 *
 * In order to adapt to different environments (on-line/off-line), the component
 * gets an environment structure with function pointers. The base class provides
 * member functions for those environment dependend functions. The member 
 * functions are used by the component implementation and are re-mapped to the
 * corresponding functions.
 *
 * @section alihltcomponent-interfaces Component interfaces
 * Each of the 3 standard component base classes AliHLTProcessor, AliHLTDataSource
 * and AliHLTDataSink provides it's own processing method (see
 * @ref alihltcomponent-type-std), which splits into a high and a low-level
 * method. For the @ref alihltcomponent-low-level-interface, all parameters are
 * shipped as function arguments, the component is supposed to write data to the
 * output buffer and handle all block descriptors. 
 * The @ref alihltcomponent-high-level-interface is the standard processing
 * method and will be used whenever the low-level method is not overloaded.
 *
 * In both cases it is necessary to calculate/estimate the size of the output
 * buffer before the processing. Output buffers can never be allocated inside
 * the component because of the push-architecture of the HLT.
 * For that reason the @ref GetOutputDataSize function should return a rough
 * estimatian of the data to be produced by the component. The component is
 * responsible for checking the memory size and must return -ENOSPC if the
 * available buffer is too small, and update the estimator respectively. The
 * framework will allocate a buffer of appropriate size and call the processing
 * again.
 *
 * @subsection alihltcomponent-error-codes Return values/Error codes
 * For return codes, the following scheme applies:
 * - The data processing methods have to indicate error conditions by a negative
 * error/return code. Preferably the system error codes are used like
 * e.g. -EINVAL. This requires to include the header
 * <pre>
 * \#include \<cerrno\>
 * </pre>
 * This schema aplies to all interface functions of the component base class.
 * For data processing it is as follows:
 * - If no suitable input block could be found (e.g. no clusters for the TPC cluster
 * finder) set size to 0, block list is empty, return 0
 * - If no ususable or significant signal could be found in the input blocks
 * return an empty output block, set size accordingly, and return 0. An empty output
 * block here could be either a real empty one of size 0 (in which case size also
 * would have to be set to zero) or a block filled with just the minimum necessary
 * accounting/meta-structures. E.g. in the TPC
 *
 * @subsection alihltcomponent-high-level-interface High-level interface
 * The high-level component interface provides functionality to exchange ROOT
 * structures between components. In contrast to the 
 * @ref alihltcomponent-low-level-interface, a couple of functions can be used
 * to access data blocks of the input stream
 * and send data blocks or ROOT TObject's to the output stream. The functionality
 * is hidden from the user and is implemented by using ROOT's TMessage class.
 *
 * @subsubsection alihltcomponent-high-level-int-methods Interface methods
 * The interface provides a couple of methods in order to get objects from the
 * input, data blocks (non TObject) from the input, and to push back objects and
 * data blocks to the output. For convenience there are several functions of 
 * identical name (and similar behavior) with different parameters defined.
 * Please refer to the function documentation.
 * - @ref GetNumberOfInputBlocks <br>
 *        return the number of data blocks in the input stream
 * - @ref GetFirstInputObject <br>
 *        get the first object of a specific data type
 * - @ref GetNextInputObject <br>
 *        get the next object of same data type as last GetFirstInputObject/Block call
 * - @ref GetFirstInputBlock <br>
 *        get the first block of a specific data type
 * - @ref GetNextInputBlock <br>
 *        get the next block of same data type as last GetFirstInputBlock/Block call
 * - @ref PushBack <br>
 *        insert an object or data buffer into the output
 * - @ref CreateEventDoneData <br>
 *        add event information to the output
 * 
 * In addition, the processing methods are simplified a bit by cutting out most of
 * the parameters.
 * @see 
 * - @ref AliHLTProcessor::DoEvent
 * - @ref AliHLTDataSource::GetEvent
 * - @ref AliHLTDataSink::DumpEvent
 *
 * \em IMPORTANT: objects and block descriptors provided by the high-level interface
 *  <b>MUST NOT BE DELETED</b> by the caller.
 *
 * @subsubsection alihltcomponent-high-level-int-guidelines High-level interface guidelines
 * - Structures must inherit from the ROOT object base class TObject in order be 
 * transported by the transportation framework.
 * - all pointer members must be transient (marked <tt>//!</tt> behind the member
 * definition), i.e. will not be stored/transported, or properly marked
 * (<tt>//-></tt>) in order to call the streamer of the object the member is pointing
 * to. The latter is not recomended. Structures to be transported between components
 * should be streamlined.
 * - no use of stl vectors/strings, use appropriate ROOT classes instead 
 * 
 * @subsection alihltcomponent-low-level-interface Low-level interface
 * The low-level component interface consists of the specific data processing
 * methods for @ref AliHLTProcessor, @ref AliHLTDataSource, and @ref AliHLTDataSink.
 * - @ref AliHLTProcessor::DoEvent
 * - @ref AliHLTDataSource::GetEvent
 * - @ref AliHLTDataSink::DumpEvent
 * 
 * The base class passes all relevant parameters for data access directly on to the
 * component. Input blocks can be accessed by means of the array <tt> blocks </tt>.
 * Output data are written directly to shared memory provided by the pointer
 * <tt> outputPtr </tt> and output block descriptors are inserted directly to the
 * list <tt> outputBlocks </tt>.
 *
 * \b NOTE: The high-level input data access methods can be used also from the low
 * level interface. Also the PushBack functions can be used BUT ONLY if no data is
 * written to the output buffer and no data block descriptors are inserted into the
 * output block list.
 *
 * @section alihltcomponent-initialization Component initialization and configuration
 * The component interface provides two optional methods for component initialization
 * and configuration. The @ref DoInit function is called once before the processing.
 * During the event processing, a special event can trigger a reconfiguration and the
 * @ref Reconfigure method is called. There are three possible options of initialization
 * and configuration:
 * - default values: set directly in the source code
 * - OCDB objects: all necessary information must be loaded from OCDB objects. The
 *   Offline Conditions Data Base stores objects specifically valid for individual runs
 *   or run ranges.
 * - Component arguments: can be specified for every component in the chain
 *   configuration. The arguments can be used to override specific parameters of the
 *   component.
 *
 * As a general rule, the three options should be processed in that sequence, i.e
 * default parameters might be overridden by OCDB configuration, and the latter one
 * by component arguments.
 *
 * @subsection alihltcomponent-initialization-arguments Component arguments
 * In normal operation, components are supposed to run without any additional argument,
 * however such arguments can be useful for testing and debugging. The idea follows
 * the format of command line arguments. A keyword is indicated by a dash and an
 * optional argument might follow, e.g.:
 * <pre>
 * -argument1 0.5 -argument2
 * </pre>
 * In this case argument1 requires an additional parameter whereas argument2 does not.
 * The arguments will be provided as an array of separated arguments.
 *
 * Component arguments can be classified into initialization arguments and configuration
 * arguments. The latter are applicable for both the @ref DoInit and @ref Reconfigure
 * method whereas initialization arguments are not applicable after DoInit.
 *
 * @subsection alihltcomponent-initialization-ocdb OCDB objects
 * OCDB objects are ROOT <tt>TObjects</tt> and can be of any type. This is in particular
 * useful for complex parameter sets. However in most cases, a simple approach of human
 * readable command line arguments is appropriate. Such a string can be simply stored
 * in a TObjString (take note that the TString does not derive from TObject). The
 * same arguments as for the command line can be used. Take note that in the TObjString
 * all arguments are separated by blanks, instead of being in an array of separate
 * strings.
 *
 * The base class provides two functions regarding OCDB objects: 
 * - LoadAndExtractOCDBObject() loads the OCDB entry for the specified path and extracts
 *                              the TObject from it. An optional key allows to access
 *                              a TObject within a TMap
 * - ConfigureFromCDBTObjString() can load a number of OCDB objects and calls the
 *                              argument parsing ConfigureFromArgumentString
 *
 *
 * @subsection alihltcomponent-initialization-sequence Initialization sequence
 * Using the approach of <tt>TObjString</tt>-type configuration objects allows to treat
 * configuration from both @ref DoInit and @ref Reconfigure in the same way.
 *
 * The base class provides the function ConfigureFromArgumentString() which loops over
 * all arguments and calls the child's method ScanConfigurationArgument(). Here the
 * actual treatment of the argument and its parameters needs to be implemented.
 * ConfigureFromArgumentString() can treat both arrays of arguments and arguments in
 * one single string separated by blanks. The two options can be mixed.
 *
 * A second base class function ConfigureFromCDBTObjString() allows to configure
 * directly from a number of OCDB objects. This requires the entries to be of
 * type TObjString and the child implementation of ScanConfigurationArgument().
 * The object can also be of type TMap with TObjStrings as key-value pairs. The
 * key identifier can be chosen by the component implementation. Normally it will
 * be the run type ("p","A-A", "p-A", ...) or e.g. the trigger code secified by
 * ECS.
 *
 * @section alihltcomponent-handling Component handling 
 * The handling of HLT analysis components is carried out by the AliHLTComponentHandler.
 * Component are registered automatically at load-time of the component shared library
 * under the following suppositions:
 * - the component library has to be loaded from the AliHLTComponentHandler using the
 *   @ref AliHLTComponentHandler::LoadLibrary method.
 * - the library defines an AliHLTModuleAgent which registers all components.
 *   See AliHLTModuleAgent::RegisterComponents                               <br>
 *     or                                                                    <br>
 * - the component implementation defines one global object (which is generated
 *   when the library is loaded)                                             <br>
 *
 * @subsection alihltcomponent-design-rules General design considerations
 * The analysis code should be implemented in one or more destict class(es). A 
 * \em component should be implemented which interface the destict analysis code to the
 * component interface. This component generates the analysis object dynamically. <br>
 *
 * Assume you have an implemetation <tt> AliHLTDetMyAnalysis </tt>, another class <tt>
 * AliHLTDetMyAnalysisComponent </tt> contains:
 * <pre>
 * private:
 *   AliHLTDetMyAnalysis* fMyAnalysis;  //!
 * </pre>
 * The object should then be instantiated in the DoInit handler of 
 * <tt>AliHLTDetMyAnalysisComponent </tt>, and cleaned in the DoDeinit handler.
 *
 * Further rules:
 * - avoid big static arrays in the component, allocate the memory at runtime
 * - allocate all kind of complex data members (like classes, ROOT TObjects of
 *   any kind) dynamically in DoInit and clean up in DoDeinit
 *
 * @section alihlt_component_arguments Default arguments
 * The component base class provides some default arguments:
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -loglevel=level     <br>
 * \li -object-compression=level     <br>
 *      compression level for ROOT objects, default is defined by
 *      @ref ALIHLTCOMPONENT_DEFAULT_OBJECT_COMPRESSION
 * \li -pushback-period=period     <br>
 *      scale down for PushBack of objects, shipped only for one event
 *      every <i>period</i> seconds
 *
 * @ingroup alihlt_component
 * @section alihltcomponent-members Class members
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
   * The call is redirected to the internal method DoInit which can be
   * overridden by the child class.
   * During Init also the environment structure is passed to the component.
   * @param comenv         environment pointer with environment dependent function
   *                       calls
   * @param environParam   additional parameter for function calls, the pointer
   *                       is passed as it is
   * @param argc           size of the argument array
   * @param argv           augment array for component initialization
   */
  virtual int Init( const AliHLTAnalysisEnvironment* comenv, void* environParam, int argc, const char** argv );

  /**
   * Clean-up function to terminate data processing.
   * Clean-up of common data structures after data processing.
   * The call is redirected to the internal method @ref DoDeinit which can be
   * overridden by the child class.
   */
  virtual int Deinit();

  /**
   * Processing of one event.
   * The method is the entrance of the event processing. The parameters are
   * cached for uses with the high-level interface and the DoProcessing
   * implementation is called.
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
  int ProcessEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
			    AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
			    AliHLTUInt32_t& size, AliHLTUInt32_t& outputBlockCnt, 
			    AliHLTComponentBlockData*& outputBlocks,
			    AliHLTComponentEventDoneData*& edd );

  /**
   * Internal processing of one event.
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
   * @param outputBlocks    out: the output block array is allocated internally
   * @param edd
   * @return neg. error code if failed
   */
  virtual int DoProcessing( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
			    AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
			    AliHLTUInt32_t& size,
			    AliHLTComponentBlockDataList& outputBlocks,
			    AliHLTComponentEventDoneData*& edd ) = 0;

  /**
   * Init the CDB.
   * The function must not be called when running in AliRoot unless it it
   * really wanted. The CDB path will be set to the specified path, which might
   * override the path initialized at the beginning of the AliRoot reconstruction.
   *
   * The method is used from the external interface in order to set the correct
   * path when running on-line. The function also initializes the function
   * callback for setting the run no during operation.
   *
   * A separation of library and component handling is maybe appropriate in the
   * future. Using the global component handler here is maybe not the cleanest
   * solution.
   * @param cdbPath      path of the CDB
   * @param pHandler     the component handler used for llibrary handling.
   */
  int InitCDB(const char* cdbPath, AliHLTComponentHandler* pHandler);

  /**
   * Set the run no for the CDB.
   * The function must not be called when running in AliRoot unless it it
   * really wanted. The CDB path will be set to the specified path, which might
   * override the run no initialized at the beginning of the AliRoot reconstruction.
   * InitCDB() has to be called before in order to really change the CDB settings.
   *
   * The method is used from the external interface in order to set the correct
   * path when running on-line.
   */
  int SetCDBRunNo(int runNo);

  /**
   * Set the run description.
   * The run description is set before the call of Init() -> DoInit().
   * @note: This functionality has been added in Juli 2008. The transmission of
   * run properties by a special SOR (SOD event in DAQ terminalogy but this was
   * changed after the HLT interface was designed) event is not sufficient because
   * the data might be needed already in the DoInit handler of the component.
   * @param desc    run descriptor, currently only the run no member is used
   * @param runType originally, run type was supposed to be a number and part
   *                of the run descriptor. But it was defined as string later
   */
  int SetRunDescription(const AliHLTRunDesc* desc, const char* runType);

  /**
   * Set the component description.
   * The description string can contain tokens separated by blanks, a token
   * consists of a key and an optional value separated by '='.
   * Possible keys:
   * \li -chainid=id        string id within the chain of the instance
   *
   * @param desc    component description
   */
  int SetComponentDescription(const char* desc);

  /**
   * Set the running environment for the component.
   * Originally, the environment was set in the Init function. However, the setup of
   * the CDB is required before. In order to have proper logging functionality, the
   * environment is required.
   * @param comenv         environment pointer with environment dependent function
   *                       calls
   * @param environParam   additional parameter for function calls, the pointer
   *                       is passed as it is
   */
  int SetComponentEnvironment(const AliHLTAnalysisEnvironment* comenv, void* environParam);

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
  virtual void GetInputDataTypes( AliHLTComponentDataTypeList& ) = 0;

  /**
   * Get the output data type of the component.
   * The function is pure virtual and must be implemented by the child class.
   * @return output data type
   */
  virtual AliHLTComponentDataType GetOutputDataType() = 0;

  /**
   * Get the output data types of the component.
   * The function can be implemented to indicate multiple output data types
   * in the target array.
   * @ref GetOutputDataType must return @ref kAliHLTMultipleDataType in order
   * to invoke this method.
   * @param tgtList          list to receive the data types
   * @return no of output data types, data types in the target list
   */
  virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /**
   * Get a ratio by how much the data volume is shrunken or enhanced.
   * The function is pure virtual and must be implemented by the child class.
   * @param constBase        <i>return</i>: additive part, independent of the
   *                                   input data volume  
   * @param inputMultiplier  <i>return</i>: multiplication ratio
   * @return values in the reference variables
   */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) = 0;

  /**
   * Get a list of OCDB object description.
   * The list of objects is provided in a TMap
   * - key: complete OCDB path, e.g. GRP/GRP/Data
   * - value: short description why the object is needed
   * Key and value objects created inside this class go into ownership of
   * target TMap.
   * @param targetMap   TMap instance receiving the list
   * @return void
   */
  virtual void GetOCDBObjectDescription( TMap* const targetArray);

  /**
   * Spawn function.
   * Each component must implement a spawn function to create a new instance of 
   * the class. Basically the function must return <i>new <b>my_class_name</b></i>.
   * @return new class instance
   */
  virtual AliHLTComponent* Spawn() = 0;

  /**
   * check the availability of the OCDB entry descriptions in the TMap
   *  key : complete OCDB path of the entry
   *  value : auxiliary object - short description
   * if the external map was not provided the function invokes
   * interface function GetOCDBObjectDescription() to retrieve the list.
   * @param externList  map of entries to be tested
   * @result 0 if all found, -ENOENT if objects not found
   */
  int CheckOCDBEntries(const TMap* const externList=NULL);

  /**
   * Find matching data types between this component and a consumer component.
   * Currently, a component can produce only one type of data. This restriction is most
   * likely to be abolished in the future.
   * @param pConsumer a component and consumer of the data produced by this component
   * @param tgtList   reference to a vector list to receive the matching data types.
   * @return >= 0 success, neg. error code if failed
   */ 
  int FindMatchingDataTypes(AliHLTComponent* pConsumer, AliHLTComponentDataTypeList* tgtList);
 
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

  /**
   * Helper function to convert the data type to a string.
   * @param type        data type structure
   * @param mode        0 print string origin:type          <br>
   *                    1 print chars                       <br>
   *                    2 print numbers                     <br>
   *                    3 print 'type' 'origin' 
   */
  static string DataType2Text( const AliHLTComponentDataType& type, int mode=0);

  /**
   * Calculate a CRC checksum of a data buffer.
   * Polynomial for the calculation is 0xD8.
   */
  static AliHLTUInt32_t CalculateChecksum(const AliHLTUInt8_t* buffer, int size);

  /**
   * Helper function to print content of data type.
   */
  static void PrintDataTypeContent(AliHLTComponentDataType& dt, const char* format=NULL);

  /**
   * helper function to initialize AliHLTComponentEventData structure
   */
  static void FillEventData(AliHLTComponentEventData& evtData);

  /**
   * Print info on an AliHLTComponentDataType structure
   * This is just a helper function to examine an @ref AliHLTComponentDataType
   * structur.
   */
  static void PrintComponentDataTypeInfo(const AliHLTComponentDataType& dt);

  /**
   * Fill AliHLTComponentBlockData structure with default values.
   * @param blockData   reference to data structure
   */
  static void FillBlockData( AliHLTComponentBlockData& blockData );

  /**
   * Fill AliHLTComponentShmData structure with default values.
   * @param shmData   reference to data structure
   */
  static void FillShmData( AliHLTComponentShmData& shmData );

  /**
   * Fill AliHLTComponentDataType structure with default values.
   * @param dataType   reference to data structure
   */
  static void FillDataType( AliHLTComponentDataType& dataType );
  
  /**
   * Copy data type structure
   * Copies the value an AliHLTComponentDataType structure to another one
   * @param [out] tgtdt   target structure
   * @param [in] srcdt   source structure
   */
  static void CopyDataType(AliHLTComponentDataType& tgtdt, const AliHLTComponentDataType& srcdt);

  /**
   * Set the ID and Origin of an AliHLTComponentDataType structure.
   * The function sets the fStructureSize member and copies the strings
   * to the ID and Origin. Only characters from the valid part of the string
   * are copied, the rest is filled with 0's. <br>
   * Please note that the fID and fOrigin members are not strings, just arrays of
   * chars of size @ref kAliHLTComponentDataTypefIDsize and
   * @ref kAliHLTComponentDataTypefOriginSize respectively and not necessarily with
   * a terminating zero. <br>
   * It is possible to pass NULL pointers as id or origin argument, in that case they
   * are just ignored.
   * @param tgtdt   target data type structure
   * @param id      ID string
   * @param origin  Origin string
   */
  static void SetDataType(AliHLTComponentDataType& tgtdt, const char* id, const char* origin);

  /**
   * Set the ID and Origin of an AliHLTComponentDataType structure.
   * Given the fact that the data type ID is 64bit wide and origin 32, this helper
   * function sets the data type from those two parameters.
   * @param dt      target data type structure
   * @param id      64bit id
   * @param orig    32bit origin
   */
  static void SetDataType(AliHLTComponentDataType& dt, AliHLTUInt64_t id, AliHLTUInt32_t orig); 

  /**
   * Extract a component table entry from the payload buffer.
   * The entry consists of the AliHLTComponentTableEntry structure, the array of
   * parents and a description string of the format 'chain-id{component-id:component-args}'.
   * The function fills all the variables after a consistency check.
   */
  static int ExtractComponentTableEntry(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size,
					string& chainId, string& compId, string& compParam,
					vector<AliHLTUInt32_t>& parents) {
    int dummy=0;
    return ExtractComponentTableEntry(pBuffer, size, chainId, compId, compParam, parents, dummy);
  }

  static int ExtractComponentTableEntry(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size,
					string& chainId, string& compId, string& compParam,
					vector<AliHLTUInt32_t>& parents, int& level);

  /**
   * Extracts the different data parts from the trigger data structure.
   * @param [in] trigData  The trigger data as passed to the DoProcessing method.
   * @param [out] attributes  The data block attributes given by the HLT framework.
   * @param [out] status  The HLT status bits given by the HLT framework.
   * @param [out] cdh  The common data header received from DDL links.
   * @param [out] readoutlist  The readout list to fill with readout list bits
   *                           passed on by the HLT framework.
   * @param [in] printErrors  If true then error messages are generated as necessary
   *                          and suppressed otherwise.
   * @note If any of the output parameters are set to NULL then the field is not set.
   *   For example, the following line will only fill the CDH pointer.
   *   \code
   *     AliRawDataHeader* cdh;
   *     ExtractTriggerData(trigData, NULL, NULL, &cdh, NULL);
   *   \endcode
   * @return zero on success or one of the following error codes on failure.
   *   if a non-zero error code is returned then none of the output parameters are
   *   modified.
   *    \li -ENOENT  The <i>trigData</i> structure size is wrong.
   *    \li -EBADF   The <i>trigData</i> data size is wrong.
   *    \li -EBADMSG The common data header (CDH) in the trigger data has the wrong
   *                 number of words indicated.
   *    \li -EPROTO  The readout list structure in the trigger data has the wrong
   *                 number of words indicated.
   */
  static int ExtractTriggerData(
      const AliHLTComponentTriggerData& trigData,
      const AliHLTUInt8_t (**attributes)[gkAliHLTBlockDAttributeCount],
      AliHLTUInt64_t* status,
      const AliRawDataHeader** cdh,
      AliHLTReadoutList* readoutlist,
      bool printErrors = false
    );

  /**
   * Extracts the readout list from a trigger data structure.
   * @param [in] trigData  The trigger data as passed to the DoProcessing method.
   * @param [out] list  The output readout list to fill.
   * @param [in] printErrors  If true then error messages are generated as necessary
   *                          and suppressed otherwise.
   * @return zero on success or one of the error codes returned by ExtractTriggerData.
   */
  static int GetReadoutList(
      const AliHLTComponentTriggerData& trigData, AliHLTReadoutList& list,
      bool printErrors = false
    )
  {
    return ExtractTriggerData(trigData, NULL, NULL, NULL, &list, printErrors);
  }

  /**
   * Extracts the event type from the given Common Data Header.
   * @param [in] cdh  The Common Data Header to extract the event type from.
   * @return the event type code from the CDH.
   */
  static AliHLTUInt32_t ExtractEventTypeFromCDH(const AliRawDataHeader* cdh);
  
  /**
   * Stopwatch type for benchmarking.
   */
  enum AliHLTStopwatchType {
    /** total time for event processing */
    kSWBase,
    /** detector algorithm w/o interface callbacks */
    kSWDA,
    /** data sources */
    kSWInput,
    /** data sinks */
    kSWOutput,
    /** number of types */
    kSWTypeCount
  };

  /**
   * Helper class for starting and stopping a stopwatch.
   * The guard can be used by instantiating an object in a function. The
   * specified stopwatch is started and the previous stopwatch put on
   * hold. When the function is terminated, the object is deleted automatically
   * deleted, stopping the stopwatch and starting the one on hold.<br>
   * \em IMPORTANT: never create dynamic objects from this guard as this violates
   * the idea of a guard.
   */
  class AliHLTStopwatchGuard {
  public:
    /** standard constructor (not for use) */
    AliHLTStopwatchGuard();
    /** constructor */
    AliHLTStopwatchGuard(TStopwatch* pStart);
    /** copy constructor (not for use) */
    AliHLTStopwatchGuard(const AliHLTStopwatchGuard&);
    /** assignment operator (not for use) */
    AliHLTStopwatchGuard& operator=(const AliHLTStopwatchGuard&);
    /** destructor */
    ~AliHLTStopwatchGuard();

  private:
    /**
     * Hold the previous guard for the existence of this guard.
     * Checks whether this guard controls a new stopwatch. In that case, the
     * previous guard and its stopwatch are put on hold.
     * @param pSucc        instance of the stopwatch of the new guard
     * @return    1 if pSucc is a different stopwatch which should
     *            be started<br>
     *            0 if it controls the same stopwatch
     */
    int Hold(const TStopwatch* pSucc);

    /**
     * Resume the previous guard.
     * Checks whether the peceeding guard controls a different stopwatch. In that
     * case, the its stopwatch is resumed.
     * @param pSucc        instance of the stopwatch of the new guard
     * @return    1 if pSucc is a different stopwatch which should
     *            be stopped<br>
     *            0 if it controls the same stopwatch
     */
    int Resume(const TStopwatch* pSucc);

    /** the stopwatch controlled by this guard */
    TStopwatch* fpStopwatch;                                                //!transient

    /** previous stopwatch guard, put on hold during existence of the guard */
    AliHLTStopwatchGuard* fpPrec;                                           //!transient

    /** active stopwatch guard */
    static AliHLTStopwatchGuard* fgpCurrent;                                //!transient
  };

  /**
   * Set a stopwatch for a given purpose.
   * @param pSW         stopwatch object
   * @param type        type of the stopwatch
   */
  int SetStopwatch(TObject* pSW, AliHLTStopwatchType type=kSWBase);

  /**
   * Init a set of stopwatches.
   * @param pStopwatches object array of stopwatches
   */
  int SetStopwatches(TObjArray* pStopwatches);

  /**
   * Customized logging function.
   * The chain id, component id and pointer is added at the beginning of each message.
   */
  int LoggingVarargs(AliHLTComponentLogSeverity severity, 
		     const char* originClass, const char* originFunc,
		     const char* file, int line, ... ) const;

  /**
   * Get size of last serialized object.
   * During PushBack, TObjects are serialized in a separate buffer. The
   * size of the last object can be retrieved by this function.
   *
   * This might be especially useful for PushBack failures caused by too
   * small output buffer.
   */
  int GetLastObjectSize() const {return fLastObjectSize;}

  /**
   * This method generates a V4 Globally Unique Identifier (GUID) using the
   * ROOT TRandom3 pseudo-random number generator with the process' UID, GID
   * PID and host address as seeds. For good measure MD5 sum hashing is also
   * applied.
   * @return the newly generated GUID structure.
   */
  static TUUID GenerateGUID();

  /// get the compression level for TObjects
  int GetCompressionLevel() const {return fCompressionLevel;}

 protected:

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
   * Reconfigure the component.
   * The method is called when an event of type @ref kAliHLTDataTypeComConf
   * {COM_CONF:PRIV} is received by the component. If the event is sent as
   * part of a normal event, the component configuration is called first.
   *
   * The CDB path parameter specifies the path in the CDB, i.e. without
   * leading absolute path of the CDB location. The framework might also
   * provide the id of the component in the analysis chain.
   *
   * The actual sequence of configuration depends on the component. As a
   * general rule, the component should load the specific OCDB object if
   * provided as parameter, and load the default objects if the parameter
   * is NULL. However, other schemes are possible. See @ref 
   *
   * \b Note: The CDB will be initialized by the framework, either already set
   * from AliRoot or from the wrapper interface during initialization.
   *
   * @param cdbEntry     path of the cdbEntry
   * @param chainId      the id/name of the component in the current analysis
   *                     chain. This is not necessarily the same as what is
   *                     returned by the GetComponentID() method.
   * @note both parameters can be NULL, check before usage
   */
  virtual int Reconfigure(const char* cdbEntry, const char* chainId);

  /**
   * Read the Preprocessor values.
   * The function is invoked when the component is notified about available/
   * updated data points from the detector Preprocessors. The 'modules'
   * argument contains all detectors for which the Preprocessors have
   * updated data points. The component has to implement the CDB access to
   * get the desired data points.
   * @param modules     detectors for which the Preprocessors have updated
   *                    data points: TPC, TRD, ITS, PHOS, MUON, or ALL if
   *                    no argument was received.
   * @return neg. error code if failed
   */
  virtual int ReadPreprocessorValues(const char* modules);

  /**
   * Child implementation to scan a number of configuration arguments.
   * The method is invoked by the framework in conjunction with the
   * common framework functions ConfigureFromArgumentString and
   * ConfigureFromCDBTObjString.
   * Function needs to scan the argument and optional additional
   * parameters and returns the number of elements in the array which
   * have been treated.
   * @param argc
   * @param argv
   * @return number of arguments which have been scanned or neg error
   *         code if failed                                              <br>
   *         \li -EINVAL      unknown argument
   *         \li -EPROTO      protocol error, e.g. missing parameter
   */
  virtual int ScanConfigurationArgument(int argc, const char** argv);

  /**
   * Custom handler for the SOR event.
   * Is invoked from the base class if an SOR event is in the block list.
   * The handler is called before the processing function. The processing
   * function is skipped if there are no other data blocks available.
   *
   * The SOR event is generated by the PubSub framework in response to
   * the DAQ start of data (SOD - has been renamed after HLT interface
   * was designed). The SOD event consists of 3 blocks:
   * - ::kAliHLTDataTypeEvent block: spec ::gkAliEventTypeStartOfRun
   * - SOD block of type ::kAliHLTDataTypeSOR, payload: AliHLTRunDesc struct
   * - run type block ::kAliHLTDataTypeRunType, payload: run type string 
   *
   * Run properties can be retrieved by getters like GetRunNo().
   * @return neg. error code if failed
   */
  virtual int StartOfRun();

  /**
   * Custom handler for the EOR event.
   * Is invoked from the base class if an EOR event is in the block list.
   * The handler is called before the processing function. The processing
   * function is skipped if there are no other data blocks available.
   *
   * See StartOfRun() for more comments of the sequence of steering events.
   *
   * @return neg. error code if failed
   */
  virtual int EndOfRun();

  /**
   * Check whether a component requires all steering blocks.
   * Childs can overload in order to indicate that they want to
   * receive also the steering data blocks. There is also the
   * possibility to add the required data types to the input
   * data type list in GetInputDataTypes().
   */
  virtual bool RequireSteeringBlocks() const {return false;}

  /**
   * General memory allocation method.
   * All memory which is going to be used 'outside' of the interface must
   * be provided by the framework (online or offline).
   * The method is redirected to a function provided by the current
   * framework. Function pointers are transferred via the @ref
   * AliHLTAnalysisEnvironment structure.
   */
  void* AllocMemory( unsigned long size );

  /**
   * Helper function to create a monolithic BlockData description block out
   * of a list BlockData descriptors.
   * For convenience, inside the interface vector lists are used, to make the
   * interface pure C style, monilithic blocks must be exchanged. 
   * The method is redirected to a function provided by the current
   * framework. Function pointers are transferred via the @ref
   * AliHLTAnalysisEnvironment structure.
   */
  int MakeOutputDataBlockList( const AliHLTComponentBlockDataList& blocks, AliHLTUInt32_t* blockCount,
			       AliHLTComponentBlockData** outputBlocks );

  /**
   * Fill the EventDoneData structure.
   * The method is redirected to a function provided by the current
   * framework. Function pointers are transferred via the @ref
   * AliHLTAnalysisEnvironment structure.
   */
  int GetEventDoneData( unsigned long size, AliHLTComponentEventDoneData** edd ) const;

  /**
   * Allocate an EventDoneData structure for the current event .
   * The method allocates the memory internally and does not interact with the current Framework.
   * The allocated data structure is empty initially and can be filled by calls to the 
   * @ref PushEventDoneData method. The memory will be automatically released after the event has been processed.
   * 
   */
  int ReserveEventDoneData( unsigned long size );

  /**
   * Push a 32 bit word of data into event done data for the current event which
   * has previously been allocated by the @ref ReserveEventDoneData method.
   */
  int PushEventDoneData( AliHLTUInt32_t eddDataWord );

  /**
   * Release event done data previously reserved by @ref ReserveEventDoneData
   */
   void ReleaseEventDoneData();

  /**
   * Get the pointer to the event done data available/built so far for the current event via
   * @ref ReserveEventDoneData and @ref PushEventDoneData
   */
  AliHLTComponentEventDoneData* GetCurrentEventDoneData() const
    {
    return fEventDoneData;
    }

  /**
   * Helper function to convert the data type to a string.
   */
  void DataType2Text(const AliHLTComponentDataType& type, char output[kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2]) const;

  /**
   * Loop through a list of component arguments.
   * The list can be either an array of separated strings or one single
   * string containing blank separated arguments, or both mixed.
   * ScanConfigurationArgument() is called to allow the component to treat
   * the individual arguments.
   * @return neg. error code if failed
   */
  int ConfigureFromArgumentString(int argc, const char** argv);

  /**
   * Read configuration objects from OCDB and configure from
   * the content of TObjString entries.
   * @param entries   blank separated list of OCDB paths
   * @param key       if the entry is a TMap, search for the corresponding object
   * @return neg. error code if failed
   */
  int ConfigureFromCDBTObjString(const char* entries, const char* key=NULL);

  /**
   * Load specified entry from the OCDB and extract the object.
   * The entry is explicitely unloaded from the cache before it is loaded.
   * If parameter key is specified the OCDB object is treated as TMap
   * and the TObject associated with 'key' is loaded.
   * @param path      path of the entry under to root of the OCDB
   * @param version   version of the entry
   * @param subVersion  subversion of the entry
   * @param key       key of the object within TMap
   */
  TObject* LoadAndExtractOCDBObject(const char* path, const char* key=NULL) const;

  /**
   * Get event number.
   * @return value of the internal event counter
   */
  int GetEventCount() const;

  /**
   * Get the number of input blocks.
   * @return number of input blocks
   */
  int GetNumberOfInputBlocks() const;

  /**
   * Get id of the current event
   * @return event id
   */
  AliHLTEventID_t GetEventId() const;

  /**
   * Get the first object of a specific data type from the input data.
   * The High-level methods provide functionality to transfer ROOT data
   * structures which inherit from TObject.
   *
   * The method looks for the first ROOT object of type dt in the input stream.
   * If also the class name is provided, the object is checked for the right
   * class type. The input data block needs a certain structure, namely the 
   * buffer size as first word. If the cross check fails, the retrieval is
   * silently abandoned, unless the \em bForce parameter is set.<br>
   * \b Note: THE OBJECT MUST NOT BE DELETED by the caller.
   *
   * If called without parameters, the function tries to create objects from
   * all available input blocks, also the ones of data type kAliHLTVoidDataType
   * which are not matched by kAliHLTAnyDataType.
   *
   * @param dt          data type of the object
   * @param classname   class name of the object
   * @param bForce      force the retrieval of an object, error messages
   *                    are suppressed if \em bForce is not set
   * @return pointer to @ref TObject, NULL if no objects of specified type
   *         available
   */
  const TObject* GetFirstInputObject(const AliHLTComponentDataType& dt=kAliHLTAllDataTypes,
				     const char* classname=NULL,
				     int bForce=0);

  /**
   * Get the first object of a specific data type from the input data.
   * The High-level methods provide functionality to transfer ROOT data
   * structures which inherit from TObject.
   * The method looks for the first ROOT object of type specified by the ID and 
   * Origin strings in the input stream.
   * If also the class name is provided, the object is checked for the right
   * class type. The input data block needs a certain structure, namely the 
   * buffer size as first word. If the cross check fails, the retrieval is
   * silently abandoned, unless the \em bForce parameter is set.<br>
   * \em Note: THE OBJECT MUST NOT BE DELETED by the caller.
   * @param dtID        data type ID of the object
   * @param dtOrigin    data type origin of the object
   * @param classname   class name of the object
   * @param bForce      force the retrieval of an object, error messages
   *                    are suppressed if \em bForce is not set
   * @return pointer to @ref TObject, NULL if no objects of specified type
   *         available
   */
  const TObject* GetFirstInputObject(const char* dtID, 
				     const char* dtOrigin,
				     const char* classname=NULL,
				     int bForce=0);

  /**
   * Get the next object of a specific data type from the input data.
   * The High-level methods provide functionality to transfer ROOT data
   * structures which inherit from TObject.
   * The method looks for the next ROOT object of type and class specified
   * to the previous @ref GetFirstInputObject call.<br>
   * \em Note: THE OBJECT MUST NOT BE DELETED by the caller.
   * @param bForce      force the retrieval of an object, error messages
   *                    are suppressed if \em bForce is not set
   * @return pointer to @ref TObject, NULL if no more objects available
   */
  const TObject* GetNextInputObject(int bForce=0);

  /**
   * Get data type of an input block.
   * Get data type of the object previously fetched via
   * GetFirstInputObject/NextInputObject or the last one if no object
   * specified.
   * @param pObject     pointer to TObject
   * @return data specification, kAliHLTVoidDataSpec if failed
   */
  AliHLTComponentDataType GetDataType(const TObject* pObject=NULL);

  /**
   * Get data specification of an input block.
   * Get data specification of the object previously fetched via
   * GetFirstInputObject/NextInputObject or the last one if no object
   * specified.
   * @param pObject     pointer to TObject
   * @return data specification, kAliHLTVoidDataSpec if failed
   */
  AliHLTUInt32_t GetSpecification(const TObject* pObject=NULL);

  /**
   * Get the first block of a specific data type from the input data.
   * The method looks for the first block of type dt in the input stream.
   * It is intended to be used within the high-level interface.<br>
   * \em Note: THE BLOCK DESCRIPTOR MUST NOT BE DELETED by the caller.
   *
   * If called without parameters, the function works on all input blocks,
   * also the ones of data type kAliHLTVoidDataType which are not matched by
   * kAliHLTAnyDataType.
   *
   * @param dt          data type of the block
   * @return pointer to @ref AliHLTComponentBlockData
   */
  const AliHLTComponentBlockData* GetFirstInputBlock(const AliHLTComponentDataType& dt=kAliHLTAllDataTypes);

  /**
   * Get the first block of a specific data type from the input data.
   * The method looks for the first block of type specified by the ID and 
   * Origin strings in the input stream.  It is intended
   * to be used within the high-level interface.<br>
   * \em Note: THE BLOCK DESCRIPTOR MUST NOT BE DELETED by the caller.
   * @param dtID        data type ID of the block
   * @param dtOrigin    data type origin of the block
   * @return pointer to @ref AliHLTComponentBlockData
   */
  const AliHLTComponentBlockData* GetFirstInputBlock(const char* dtID, 
						      const char* dtOrigin);

  /**
   * Get input block by index.<br>
   * \em Note: THE BLOCK DESCRIPTOR MUST NOT BE DELETED by the caller.
   * @return pointer to AliHLTComponentBlockData, NULL if index out of range
   */
  const AliHLTComponentBlockData* GetInputBlock(int index) const;

  /**
   * Get the next block of a specific data type from the input data.
   * The method looks for the next block  of type and class specified
   * to the previous @ref GetFirstInputBlock call.
   * To be used within the high-level interface.<br>
   * \em Note: THE BLOCK DESCRIPTOR MUST NOT BE DELETED by the caller.
   */
  const AliHLTComponentBlockData* GetNextInputBlock();

  /**
   * Get data specification of an input block.
   * Get data specification of the input block previously fetched via
   * GetFirstInputObject/NextInputObject or the last one if no block
   * specified.
   * @param pBlock     pointer to input block
   * @return data specification, kAliHLTVoidDataSpec if failed
   */
  AliHLTUInt32_t GetSpecification(const AliHLTComponentBlockData* pBlock);

  /**
   * Forward an input object to the output.
   * Forward the input block of an object previously fetched via
   * GetFirstInputObject/NextInputObject or the last one if no object
   * specified.
   * The block descriptor of the input block is forwarded to the
   * output block list.
   * @param pObject     pointer to TObject
   * @return neg. error code if failed
   */
  int Forward(const TObject* pObject);

  /**
   * Forward an input block to the output.
   * Forward the input block fetched via GetFirstInputObject/
   * NextInputBlock or the last one if no block specified.
   * The block descriptor of the input block is forwarded to the
   * output block list.
   * @param pBlock     pointer to input block
   * @return neg. error code if failed
   */
  int Forward(const AliHLTComponentBlockData* pBlock=NULL);

  /**
   * Insert an object into the output.
   * If header is specified, it will be inserted before the root object,
   * default is no header.
   * The publishing can be downscaled by means of the -pushback-period
   * parameter. This is especially useful for histograms which do not
   * need to be sent for every event.
   * @param pObject     pointer to root object
   * @param dt          data type of the object
   * @param spec        data specification
   * @param pHeader     pointer to header
   * @param headerSize  size of Header
   * @return neg. error code if failed 
   */
  int PushBack(const TObject* pObject, const AliHLTComponentDataType& dt, 
	       AliHLTUInt32_t spec=kAliHLTVoidDataSpec, 
	       void* pHeader=NULL, int headerSize=0);

  /**
   * Insert an object into the output.
   * If header is specified, it will be inserted before the root object,
   * default is no header.
   * The publishing can be downscaled by means of the -pushback-period
   * parameter. This is especially useful for histograms which do not
   * need to be sent for every event.
   * @param pObject     pointer to root object
   * @param dtID        data type ID of the object
   * @param dtOrigin    data type origin of the object
   * @param spec        data specification
   * @param pHeader     pointer to header
   * @param headerSize  size of Header
   * @return neg. error code if failed 
   */
  int PushBack(const TObject* pObject, const char* dtID, const char* dtOrigin,
	       AliHLTUInt32_t spec=kAliHLTVoidDataSpec,
	       void* pHeader=NULL, int headerSize=0);
 
  /**
   * Insert an object into the output.
   * @param pBuffer     pointer to buffer
   * @param iSize       size of the buffer
   * @param dt          data type of the object
   * @param spec        data specification
   * @param pHeader     pointer to header
   * @param headerSize size of Header
   * @return neg. error code if failed 
   */
  int PushBack(const void* pBuffer, int iSize, const AliHLTComponentDataType& dt,
	       AliHLTUInt32_t spec=kAliHLTVoidDataSpec,
	       const void* pHeader=NULL, int headerSize=0);

  /**
   * Insert an object into the output.
   * @param pBuffer     pointer to buffer
   * @param iSize       size of the buffer
   * @param dtID        data type ID of the object
   * @param dtOrigin    data type origin of the object
   * @param spec        data specification
   * @param pHeader     pointer to header
   * @param headerSize size of Header
   * @return neg. error code if failed 
   */
  int PushBack(const void* pBuffer, int iSize, const char* dtID, const char* dtOrigin,
	       AliHLTUInt32_t spec=kAliHLTVoidDataSpec,
	       const void* pHeader=NULL, int headerSize=0);

  /**
   * Estimate size of a TObject
   * @param pObject
   * @return buffer size in byte
   */
  int EstimateObjectSize(const TObject* pObject) const;

  /**
   * Create a memory file in the output stream.
   * This method creates a TFile object which stores all data in
   * memory instead of disk. The TFile object is published as binary data.
   * The instance can be used like a normal TFile object. The TFile::Close
   * or @ref CloseMemoryFile method has to be called in order to flush the
   * output stream.
   *
   * \b Note: The returned object is deleted by the framework.
   * @param capacity    total size reserved for the memory file
   * @param dtID        data type ID of the file
   * @param dtOrigin    data type origin of the file
   * @param spec        data specification
   * @return file handle, NULL if failed 
   */
  AliHLTMemoryFile* CreateMemoryFile(int capacity, const char* dtID, const char* dtOrigin,
				     AliHLTUInt32_t spec=kAliHLTVoidDataSpec);

  /**
   * Create a memory file in the output stream.
   * This method creates a TFile object which stores all data in
   * memory instead of disk. The TFile object is published as binary data.
   * The instance can be used like a normal TFile object. The TFile::Close
   * or @ref CloseMemoryFile method has to be called in order to flush the
   * output stream.
   *
   * \b Note: The returned object is deleted by the framework.
   * @param capacity    total size reserved for the memory file
   * @param dt          data type of the file
   * @param spec        data specification
   * @return file handle, NULL if failed 
   */
  AliHLTMemoryFile* CreateMemoryFile(int capacity, 
				     const AliHLTComponentDataType& dt=kAliHLTAnyDataType,
				     AliHLTUInt32_t spec=kAliHLTVoidDataSpec);

  /**
   * Create a memory file in the output stream.
   * This method creates a TFile object which stores all data in
   * memory instead of disk. The TFile object is published as binary data.
   * The instance can be used like a normal TFile object. The TFile::Close
   * or @ref CloseMemoryFile method has to be called in order to flush the
   * output stream.
   *
   * \b Note: The returned object is deleted by the framework.
   * @param dtID        data type ID of the file
   * @param dtOrigin    data type origin of the file
   * @param spec        data specification
   * @param capacity    fraction of the available output buffer size
   * @return file handle, NULL if failed 
   */
  AliHLTMemoryFile* CreateMemoryFile(const char* dtID, const char* dtOrigin,
				     AliHLTUInt32_t spec=kAliHLTVoidDataSpec,
				     float capacity=1.0);

  /**
   * Create a memory file in the output stream.
   * This method creates a TFile object which stores all data in
   * memory instead of disk. The TFile object is published as binary data.
   * The instance can be used like a normal TFile object. The TFile::Close
   * or @ref CloseMemoryFile method has to be called in order to flush the
   * output stream.
   *
   * \b Note: The returned object is deleted by the framework.
   * @param dt          data type of the file
   * @param spec        data specification
   * @param capacity    fraction of the available output buffer size
   * @return file handle, NULL if failed 
   */
  AliHLTMemoryFile* CreateMemoryFile(const AliHLTComponentDataType& dt=kAliHLTAnyDataType,
				     AliHLTUInt32_t spec=kAliHLTVoidDataSpec,
				     float capacity=1.0);

  /**
   * Write an object to memory file in the output stream.
   * @param pFile       file handle
   * @param pObject     pointer to root object
   * @param key         key in ROOT file
   * @param option      options, see TObject::Write
   * @return neg. error code if failed
   *         - -ENOSPC    no space left
   */
  int Write(AliHLTMemoryFile* pFile, const TObject* pObject, const char* key=NULL, int option=TObject::kOverwrite);

  /**
   * Close object memory file.
   * @param pFile       file handle
   * @return neg. error code if failed
   *         - -ENOSPC    buffer size too small
   */
  int CloseMemoryFile(AliHLTMemoryFile* pFile);

  /**
   * Insert event-done data information into the output.
   * @param edd          event-done data information
   */
  int CreateEventDoneData(AliHLTComponentEventDoneData edd);

  /**
   * Get current run number
   */
  AliHLTUInt32_t GetRunNo() const;

  /**
   * Get the current run type.
   */
  AliHLTUInt32_t GetRunType() const;

  /**
   * Get the chain id of the component.
   */
  const char* GetChainId() const {return fChainId.c_str();}

  /**
   * Get a timestamp of the current event
   * Exact format needs to be documented.
   */
  AliHLTUInt32_t    GetTimeStamp() const;

  /**
   * Get the period number.
   * Upper 28 bits (36 to 63) of the 64-bit event id 
   */
  AliHLTUInt32_t    GetPeriodNumber() const;

  /**
   * Get the period number.
   * 24 bits, 12 to 35 of the 64-bit event id 
   */
  AliHLTUInt32_t    GetOrbitNumber() const;

  /**
   * Get the bunch crossing number.
   * 12 bits, 0 to 12 of the 64-bit event id 
   */
  AliHLTUInt16_t    GetBunchCrossNumber() const;

  /**
   * Setup the CTP accounting functionality of the base class.
   * The method can be invoked from DoInit() for componenets which want to
   * use the CTP functionality of the base class.
   *
   * The AliHLTCTPData is initialized with the trigger classes from the ECS
   * parameters. The base class automatically increments the counters according
   * to the trigger pattern in the CDH before the event processing. 
   */
  int SetupCTPData();

  /**
   * Get the instance of the CTP data.
   */
  const AliHLTCTPData* CTPData() const {return fpCTPData;}

  /**
   * Check whether a combination of trigger classes is fired.
   * The expression can contain trigger class ids and logic operators
   * like &&, ||, !, and ^, and may be grouped by parentheses.
   * @note the function requires the setup of the CTP handling for the component by
   * invoking SetupCTPData() from DoInit()
   * @param expression     a logic expression of trigger class ids
   * @param trigData       the trigger data data
   */
  bool EvaluateCTPTriggerClass(const char* expression, AliHLTComponentTriggerData& trigData) const;

  /**
   * Check state of a trigger class.
   * If the class name is not part of the current trigger setup (i.e. ECS parameter
   * does not contain a trigger definition for this class name) the function
   * returns -1
   * @note the function requires the setup of the CTP handling for the component by
   * invoking SetupCTPData() from DoInit()
   * @return -1 class name not initialized, 
   *          0 trigger not active
   *          1 trigger active
   */
  int CheckCTPTrigger(const char* name) const;

  /**
   * Get the overall solenoid field.
   */
  Double_t GetBz();
  /**
   * Get the solenoid field at point r.
   */
  Double_t GetBz(const Double_t *r);
  /**
   * Get the solenoid field components at point r.
   */
  void GetBxByBz(const Double_t r[3], Double_t b[3]);

  /**
   * Check whether the current event is a valid data event.
   * @param pTgt    optional pointer to get the event type
   * @return true if the current event is a real data event
   */
  bool IsDataEvent(AliHLTUInt32_t* pTgt=NULL) const;

  /**
   * Copy a struct from block data.
   * The function checks for block size and struct size. The least common
   * size will be copied to the target struct, remaining fields are initialized
   * to zero.<br>
   * The target struct must have a 32bit struct size indicator as first member.
   * @param pStruct     target struct
   * @param iStructSize size of the struct
   * @param iBlockNo    index of input block
   * @param structname  name of the struct (log messages)
   * @param eventname   name of the event (log messages)
   * @return size copied, neg. error if failed
   */
  int CopyStruct(void* pStruct, unsigned int iStructSize, unsigned int iBlockNo,
		 const char* structname="", const char* eventname="");

 private:
  /** copy constructor prohibited */
  AliHLTComponent(const AliHLTComponent&);
  /** assignment operator prohibited */
  AliHLTComponent& operator=(const AliHLTComponent&);

  /**
   * Increment the internal event counter.
   * To be used by the friend classes AliHLTProcessor, AliHLTDataSource
   * and AliHLTDataSink.
   * @return new value of the internal event counter
   * @internal
   */
  int IncrementEventCounter();

  /**
   * Find the first input block of specified data type beginning at index.
   * Input blocks containing a TObject have the size of the object as an
   * unsigned 32 bit number in the first 4 bytes. This has to match the block
   * size minus 4.
   *
   * kAliHLTAllDataTypes is a special data type which includes both 
   * kAliHLTVoidDataType and kAliHLTAnyDataType.
   *
   * @param dt          data type
   * @param startIdx    index to start the search
   * @param bObject     check if this is an object
   * @return index of the block, -ENOENT if no block found
   *
   * @internal
   */
  int FindInputBlock(const AliHLTComponentDataType& dt, int startIdx=-1, int bObject=0) const;

  /**
   * Get index in the array of input bocks.
   * Calculate index and check integrety of a block data structure pointer.
   * @param pBlock      pointer to block data
   * @return index of the block, -ENOENT if no block found
   *
   * @internal
   */
  int FindInputBlock(const AliHLTComponentBlockData* pBlock) const;

  /**
   * Create an object from a specified input block.
   * @param idx         index of the input block
   * @param bForce      force the retrieval of an object, error messages
   *                    are suppressed if \em bForce is not set
   * @return pointer to TObject, caller must delete the object after use
   *
   * @internal
   */
  TObject* CreateInputObject(int idx, int bForce=0);

  /**
   * Get input object
   * Get object from the input block list. The methods first checks whether the
   * object was already created. If not, it is created by @ref CreateInputObject
   * and inserted into the list of objects.
   * @param idx         index in the input block list
   * @param classname   name of the class, object is checked for correct class
   *                    name if set
   * @param bForce      force the retrieval of an object, error messages
   *                    are suppressed if \em bForce is not set
   * @return pointer to TObject
   *
   * @internal
   */
  TObject* GetInputObject(int idx, const char* classname=NULL, int bForce=0);

  /**
   * Clean the list of input objects.
   * Cleanup is done at the end of each event processing.
   */
  int CleanupInputObjects();

  /**
   * Insert a buffer into the output block stream.
   * This is the only method to insert blocks into the output stream, called
   * from all types of the Pushback method. The actual data might have been
   * written to the output buffer already. In that case NULL can be provided
   * as buffer, only the block descriptor will be build. If a header is specified, 
   * it will be inserted before the buffer, default is no header.
   * @param pBuffer     pointer to buffer
   * @param iBufferSize size of the buffer in byte
   * @param dt          data type
   * @param spec        data specification
   * @param pHeader     pointer to header
   * @param iHeaderSize size of Header
   * @return neg. error code if failed
   */
  int InsertOutputBlock(const void* pBuffer, int iBufferSize,
			const AliHLTComponentDataType& dt,
			AliHLTUInt32_t spec,
			const void* pHeader=NULL, int iHeaderSize=0);

  /**
   * Add a component statistics block to the output.
   * @return size of the added data
   */
  int AddComponentStatistics(AliHLTComponentBlockDataList& blocks, 
			     AliHLTUInt8_t* buffer,
			     AliHLTUInt32_t bufferSize,
			     AliHLTUInt32_t offset,
			     AliHLTComponentStatisticsList& stats) const;

  /**
   * Add a component table entry (descriptor) to the output
   * This is done at SOR/EOR. The component table is a list of chain ids
   * and 32bit ids calculated by a crc algorithm from the chian id. This
   * allows to tag data blocks with the id number rather than the string.
   *
   * The kAliHLTDataTypeComponentTable data block currently has the string
   * as payload and the crc id as specification.
   * @return size of the added data
   */
  int  AddComponentTableEntry(AliHLTComponentBlockDataList& blocks, 
			      AliHLTUInt8_t* buffer,
			      AliHLTUInt32_t bufferSize,
			      AliHLTUInt32_t offset,
			      const vector<AliHLTUInt32_t>& parents,
			      int processingLevel) const;

  /**
   * Scan the ECS parameter string.
   * The framework provides both the parameters of CONFIGURE and ENGAGE
   * in one string in a special data block kAliHLTDataTypeECSParam
   * {ECSPARAM:PRIV}. The general format is
   * <command>;<parameterkey>=<parametervalue>;<parameterkey>=<parametervalue>;...
   */
  int ScanECSParam(const char* ecsParam);

  /**
   * The trigger classes are determined from the trigger and propagated by
   * ECS as part of the ENGAGE command parameter which is sent through the
   * framework during the SOR event. This function treats the value of the
   * parameter key CTP_TRIGGER_CLASS.
   */
  int InitCTPTriggerClasses(const char* ctpString);

  enum {
    kRequireSteeringBlocks = 0x1,
    kDisableComponentStat = 0x2
  };

  /** The global component handler instance */
  static AliHLTComponentHandler* fgpComponentHandler;              //! transient

  /** The environment where the component is running in */
  AliHLTAnalysisEnvironment fEnvironment;                         // see above

  /** Set by ProcessEvent before the processing starts */
  AliHLTEventID_t fCurrentEvent;                                   // see above

  /** internal event no */
  int fEventCount;                                                 // see above

  /** the number of failed events */
  int fFailedEvents;                                               // see above

  /** event data struct of the current event under processing */
  AliHLTComponentEventData fCurrentEventData;                      // see above

  /** array of input data blocks of the current event */
  const AliHLTComponentBlockData* fpInputBlocks;                   //! transient

  /** index of the current input block */
  int fCurrentInputBlock;                                          // see above

  /** data type of the last block search */
  AliHLTComponentDataType fSearchDataType;                         // see above

  /** name of the class for the object to search for */
  string fClassName;                                               // see above

  /** array of generated input objects */
  TObjArray* fpInputObjects;                                       //! transient
 
  /** the output buffer */
  AliHLTUInt8_t* fpOutputBuffer;                                   //! transient

  /** size of the output buffer */
  AliHLTUInt32_t fOutputBufferSize;                                // see above

  /** size of data written to output buffer */
  AliHLTUInt32_t fOutputBufferFilled;                              // see above

  /** list of ouput block data descriptors */
  AliHLTComponentBlockDataList fOutputBlocks;                      // see above

  /** stopwatch array */
  TObjArray* fpStopwatches;                                        //! transient

  /** array of memory files AliHLTMemoryFile */
  AliHLTMemoryFilePList fMemFiles;                                 //! transient

  /** descriptor of the current run */
  AliHLTRunDesc* fpRunDesc;                                        //! transient

  /** external fct to set CDB run no, indicates external CDB initialization */
  void (*fCDBSetRunNoFunc)();                                      //! transient

  /** id of the component in the analysis chain */
  string fChainId;                                                 //! transient

  /** crc value of the chainid, used as a 32bit id */
  AliHLTUInt32_t fChainIdCrc;                                      //! transient

  /** optional benchmarking for the component statistics */
  TStopwatch* fpBenchmark;                                         //! transient

  /** component flags, cleared in Deinit */
  AliHLTUInt32_t fFlags;                                           //! transient

  /** current event type */
  AliHLTUInt32_t fEventType;                                       //! transient

  /** component arguments */
  string fComponentArgs;                                           //! transient


  /** event done data */
  AliHLTComponentEventDoneData* fEventDoneData;                    //! transient

  /** Reserved size of the memory stored at fEventDoneData */
  unsigned long fEventDoneDataSize;                                //! transient

  /** Comression level for ROOT objects */
  int fCompressionLevel;                                           //! transient

  /** size of last PushBack-serialized object */
  int fLastObjectSize;                                             //! transient

 /**  array of trigger class descriptors */
  AliHLTCTPData* fpCTPData;                                        //! transient

  /// update period for PushBack calls
  int fPushbackPeriod;                                             //! transient
  /// time of last executed PushBack
  int fLastPushBackTime;                                           //! transient

  ClassDef(AliHLTComponent, 0)
};
#endif
