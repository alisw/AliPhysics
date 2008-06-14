// -*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTROOTFILESTREAMERCOMPONENT_H
#define ALIHLTROOTFILESTREAMERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTRootFileStreamerComponent.h
    @author Matthias Richter
    @date   
    @brief  Save objects in a ROOT memory file
*/

#include "AliHLTProcessor.h"

/**
 * @class AliHLTRootFileStreamerComponent
 * The RootFileStreamer provides a stand alone component to write incoming
 * TObject like structures into a ROOT memory file. A ROOT memory file is
 * a ROOT file stored in memory instead on disk (AliHLTMemoryFile) The file
 * is published via the output stream. On the receiver side the file can
 * be directly written to disk and appears like a normal root file.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ROOTFileStreamer                                    <br>
 * Library: \b libAliHLTUtil.so                                         <br>
 * Input Data Types: ::kAliHLTAnyDataType                               <br>
 * Output Data Types: according to component arguments,
 *                    ::kAliHLTVoidDataType by default                  <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -datatype     <i> datatype   dataorigin </i> <br>
 *      data type ID and origin, e.g. <tt>-datatype CLUSTERS TPC </tt>
 * \li -dataspec     <i> specification </i> <br>
 *      data specification treated as decimal number or hex number if
 *      prepended by '0x'
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * The component does not process any event data.
 *
 * <h2>Memory consumption:</h2>
 * The component does not process any event data.
 *
 * <h2>Output size:</h2>
 * No data published (AliHLTDataSink).
 *
 * @ingroup alihlt_util_components
 */
class AliHLTRootFileStreamerComponent : public AliHLTProcessor
{
 public:
  /** standard constructor */
  AliHLTRootFileStreamerComponent();
  /** destructor */
  virtual ~AliHLTRootFileStreamerComponent();

  /**
   * The id of the component.
   * @return component id (string)
   */
  const char* GetComponentID() {return "ROOTFileStreamer";};

  /**
   * Spawn function.
   * @return new class instance
   */
  AliHLTComponent* Spawn() {return new AliHLTRootFileStreamerComponent;}

  /**
   * Get the input data types of the component.
   * The function is pure virtual and must be implemented by the child class.
   * @return list of data types in the vector reference
   */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& );

  /**
   * Get the output data type of the component.
   * The function is pure virtual and must be implemented by the child class.
   * @return output data type
   */
  AliHLTComponentDataType GetOutputDataType();

  /**
   * Get a ratio by how much the data volume is shrinked or enhanced.
   * The function is pure virtual and must be implemented by the child class.
   * @param constBase        <i>return</i>: additive part, independent of the
   *                                   input data volume  
   * @param inputMultiplier  <i>return</i>: multiplication ratio
   * @return values in the reference variables
   */
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

 protected:
  /**
   * Internal initialization.
   * @see @ref AliHLTComponent::DoInit for description and parameters
   */
  virtual int DoInit( int argc, const char** argv );

  /**
   * The high-level data processing method.
   * All incoming objects are saved into a ROOT file in memory.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @return neg. error code if failed
   */
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  
  using AliHLTProcessor::DoEvent;

 private:
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTRootFileStreamerComponent(const AliHLTRootFileStreamerComponent&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTRootFileStreamerComponent& operator=(const AliHLTRootFileStreamerComponent&);

  /** data type */
  AliHLTComponentDataType fDataType;                               // see above
  /** data specification */
  AliHLTUInt32_t          fSpecification;                          // see above

  ClassDef(AliHLTRootFileStreamerComponent, 0)
};
#endif
