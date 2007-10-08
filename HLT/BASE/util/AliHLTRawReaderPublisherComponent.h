//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTRAWREADERPUBLISHERCOMPONENT_H
#define ALIHLTRAWREADERPUBLISHERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTRawReaderPublisherComponent.h
    @author Matthias Richter
    @date   
    @brief  A general data publisher component for the AliRawReader.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTOfflineDataSource.h"

/**
 * @class AliHLTRawReaderPublisherComponent
 * A general data publisher component for the AliRawReader.
 * The component publishs the data of a given detector and equipment ID.
 * Publication of several IDs, i.e. DDLs, requires derivation of
 * the data type and/or specification from the ID. This requires a child
 * class and implementation of @ref GetSpecificationFromEquipmentId.
 * 
 * Component ID: \b AliRawReaderPublisher <br>
 * Library: \b libAliHLTUtil.
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -detector     <i> detector name      </i>
 *      e.g. <tt> -detector TPC </tt>
 * \li -equipmentid  <i> id      </i>
 *      the equipmentid
 * \li -minid  <i> id      </i>
 *      the minimum equipmentid
 * \li -maxid  <i> id      </i>
 *      the maximum equipmentid
 * \li -verbose<br>
 *      print out some more info messages, mainly for the sake of tutorials
 * \li -datatype     <i> datatype   dataorigin </i> <br>
 *      data type ID and origin, e.g. <tt>-datatype DIGITS TPC </tt>
 * \li -dataspec     <i> specification </i> <br>
 *      data specification treated as decimal number or hex number if
 *      prepended by '0x'
 *
 * Optional arguments:<br>
 *
 *
 * @ingroup alihlt_system
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
		vector<AliHLTComponentBlockData>& outputBlocks );

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

  /** be verbose: info printouts */
  Bool_t                  fVerbose;                                //!transient

  /** data type */
  AliHLTComponentDataType fDataType;                               //!transient

  /** data specification */
  AliHLTUInt32_t          fSpecification;                          //!transient

  ClassDef(AliHLTRawReaderPublisherComponent, 0);
};

#endif
