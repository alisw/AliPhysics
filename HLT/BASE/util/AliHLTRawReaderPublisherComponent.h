//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTRAWREADERPUBLISHERCOMPONENT_H
#define ALIHLTRAWREADERPUBLISHERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTRawReaderPublisherComponent.h
/// @author Matthias Richter
/// @date   
/// @brief  A general data publisher component for the AliRawReader.
///

#include "AliHLTOfflineDataSource.h"

/**
 * @class AliHLTRawReaderPublisherComponent
 * A general data publisher component for the AliRawReader.
 * The component publishs the data of a given detector and equipment ID.
 * 
 * If no data specification is given, the equipment id is used as default.
 * A child class can implement @ref GetSpecificationFromEquipmentId to
 * provide a different rule.
 *
 * The component publishes one data block for each equipment id in the
 * give range. If the RawReader does not provide any data, an empty data
 * block consisting of the Common Data Header is produced. 
 * 
 * <h2>General properties:</h2>
 *
 * Component ID: \b AliRawReaderPublisher                               <br>
 * Library: \b libAliHLTUtil.so						<br>
 * Input Data Types: none						<br>
 * Output Data Types: according to parameters and available RAW data	<br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * The equipment id(s) must be defined by argument(s) out of:
 * \li -detector     <i> detector name      </i>
 *      e.g. <tt> -detector TPC </tt>
 * \li -equipmentid  <i> id      </i>
 *      the equipmentid within the detector, e.g. TPC 0 is 768
 * \li -minid  <i> id      </i>
 *      the minimum equipmentid including detector offset, e.g. 768 is TPC 0<br>
 *      if the -detector option is used, the id is without detector offset
 * \li -maxid  <i> id      </i>
 *      the maximum equipmentid including detector offset (default = minid)<br>
 *      if the -detector option is used, the id is without detector offset
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -verbose<br>
 *      print out some more info messages, mainly for the sake of tutorials,
 *      repetitive arguments increase the level each
 * \li -silent<br>
 *      suppress all info messages
 * \li -skipempty
 *      skip all empty ddls in the specified range; by default, the component
 *      generates and inserts empty data blocks
 * \li -datatype     <i> datatype   dataorigin </i> <br>
 *      data type ID and origin, e.g. <tt>-datatype DIGITS TPC </tt>
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
 * According to the available data. The component is an AliHLTDataSource
 * and inteded to be used in the AliHLTSystem framework only. The component
 * implements the standard AliHLTSystem adaptive buffer allocation. 
 *
 * @ingroup alihlt_util_components
 */
class AliHLTRawReaderPublisherComponent : public AliHLTOfflineDataSource {
 public:
  /** standard constructor */
  AliHLTRawReaderPublisherComponent();
  /** destructor */
  virtual ~AliHLTRawReaderPublisherComponent();

  /**
   * Get the id of the component.
   * Each component is identified by a unique id.
   * The function is pure virtual and must be implemented by the child class.
   * @return component id (string)
   */
  const char* GetComponentID();

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

  /**
   * Spawn function.
   * Each component must implement a spawn function to create a new instance of 
   * the class. Basically the function must return <i>new <b>my_class_name</b></i>.
   * @return new class instance
   */
  virtual AliHLTComponent* Spawn();

 protected:
  /**
   * Init method.
   */
  int DoInit( int argc, const char** argv );

  /**
   * Deinit method.
   */
  int DoDeinit();

  /**
   * Data source method.
   * @param [in] evtData       event data structure
   * @param [in] trigData	  trigger data structure
   * @param [in] outputPtr	  pointer to target buffer
   * @param [in,out] size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param [in] outputBlocks  list to receive output block descriptors
   * @return neg. error code if failed
   */
  int GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks );

  using AliHLTOfflineDataSource::GetEvent;

  virtual bool IsSelected(int /*equipmentId*/) const {return true;}

 protected:
  virtual int GetSpecificationFromEquipmentId(int id, AliHLTUInt32_t &specification) const;

 private:
  /** copy constructor prohibited */
  AliHLTRawReaderPublisherComponent(const AliHLTRawReaderPublisherComponent&);
  /** assignment operator prohibited */
  AliHLTRawReaderPublisherComponent& operator=(const AliHLTRawReaderPublisherComponent&);

  /** max output block size, estimated during DoInit */
  Int_t                   fMaxSize;                                //!transient

  /** detector string */
  TString                 fDetector;                               //!transient

  /** min equipment id */
  int                     fMinEquId;                               //!transient

  /** max equipment id */
  int                     fMaxEquId;                               //!transient

  /** verbosity */
  Int_t                   fVerbosity;                              //!transient

  /** data type */
  AliHLTComponentDataType fDataType;                               //!transient

  /** data specification */
  AliHLTUInt32_t          fSpecification;                          //!transient

  /** skip the generation of empty data blocks */
  Bool_t                  fSkipEmpty;                              //!transient

  ClassDef(AliHLTRawReaderPublisherComponent, 1);
};

#endif
