//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTOUTPUBLISHERCOMPONENT_H
#define ALIHLTOUTPUBLISHERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTOUTPublisherComponent.h
    @author Matthias Richter
    @date   2008-06-11
    @brief  A data publisher for data block out of the HLTOUT data
*/

#include "AliHLTOfflineDataSource.h"

/**
 * @class AliHLTOUTPublisherComponent
 * A data publisher component for data blocks out of the HLTOUT data.
 * All data blocks forwarded to the HLTOUT (either real or simulated),
 * are encoded in HOMER format and stored in the HLT data links, or
 * eventually the HLT digit file in case of simulation.
 *
 * This component publishes data blocks out of the HLTOUT data. To the
 * subscribing components the data blocks seem to come directly from
 * the producing component in the HLT analysis chain. The HLTOUT is just
 * a transparent transport layer. Filter rules by data type and
 * specification can be applied.
 * 
 * Component ID: \b AliHLTOUTPublisher <br>
 * Library: \b libAliHLTUtil.
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *      
 * Optional arguments:<br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -datatype     <i> id origin      </i>
 *      e.g. <tt> -datatype 'ESD_TREE' 'TPC ' </tt> <br>
 *      \b Note: due to the 4-character data origin it might be necessary to
 *      append a blank to the detectorname, e.g. TPC -> 'TPC '
 * \li -origin  <i> origin      </i>
 *      e.g. -origin 'TPC ', \b Note:  the filter rule has type id 'ANY'
 * \li -typeid  <i> id      </i>
 *      e.g. -typeid ESD_TREE, \b Note: the filter rule has origin 'ANY'
 * \li -dataspec     <i> specification </i> <br>
 *      data specification treated as decimal number or hex number if
 *      prepended by '0x'
 * \li -verbose<br>
 *      print out some more info messages, mainly for the sake of tutorials
 *
 * By default, all blocks will be published. By means of the -datatype,
 * -origin, and -typeid arguments, the blocks can be selected. A list of filter
 * rules can be built up by multiple usage of the arguments. Each time a new
 * filter rule is added.
 *
 * No filtering by the data specification is applied unless the -specification
 * argument is used. The specification applies to to the current filter rule,
 * regardless of the sequence of -datatype/-specification arguments.
 *
 * @ingroup alihlt_system
 */
class AliHLTOUTPublisherComponent : public AliHLTOfflineDataSource {
 public:
  /** standard constructor */
  AliHLTOUTPublisherComponent();
  /** destructor */
  virtual ~AliHLTOUTPublisherComponent();

  const char* GetComponentID();
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  virtual AliHLTComponent* Spawn();

 protected:
  int DoInit( int argc, const char** argv );
  int DoDeinit();

  /**
   * Data source method.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @param outputPtr	  pointer to target buffer
   * @param size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param outputBlocks  list to receive output block descriptors
   * @return neg. error code if failed
   */
  int GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks );

  using AliHLTOfflineDataSource::GetEvent;

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTOUTPublisherComponent(const AliHLTOUTPublisherComponent&);
  /** assignment operator prohibited */
  AliHLTOUTPublisherComponent& operator=(const AliHLTOUTPublisherComponent&);

  /** filtering rules, only the data type and specification members are use */
  AliHLTComponentBlockDataList fFilterRules;                       //! transient

  /** maximum output size */
  int fMaxSize; //!

  ClassDef(AliHLTOUTPublisherComponent, 0);
};

#endif
