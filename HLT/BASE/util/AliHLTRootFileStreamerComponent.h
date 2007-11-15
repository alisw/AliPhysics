// -*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTROOTFILESTREAMERCOMPONENT_H
#define ALIHLTROOTFILESTREAMERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTRootFileStreamerComponent.h
    @author Matthias Richter
    @date   
    @brief  Save objects in a ROOT memory file

                                                                          */
#include "AliHLTProcessor.h"

/**
 * @class AliHLTRootFileStreamerComponent
 * The RootFileStreamer provides a stand alone component to write incoming
 * TObject like structures into a ROOT memory file. The memory file is
 * published via the output stream.
 *
 * Component ID: \b ROOTFileStreamer <br>
 * Library: \b libAliHLTUtil.
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 *
 * Optional arguments:<br>
 * \li -datatype     <i> datatype   dataorigin </i> <br>
 *      data type ID and origin, e.g. <tt>-datatype CLUSTERS TPC </tt>
 * \li -dataspec     <i> specification </i> <br>
 *      data specification treated as decimal number or hex number if
 *      prepended by '0x'
 *
 * @ingroup alihlt_component
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
