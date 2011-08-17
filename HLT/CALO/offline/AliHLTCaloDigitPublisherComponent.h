//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTCALODIGITPUBLISHERCOMPONENT_H
#define ALIHLTCALODIGITPUBLISHERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTCaloDigitPublisherComponent.h
    @author Matthias Richter
    @date   
    @brief  TPC digit publisher component (input from offline).
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTOfflineDataSource.h"

class AliHLTCaloDigitHandler;
class AliHLTTPCFileHandler;

/**
 * @class AliHLTCaloDigitPublisherComponent
 * A digit publisher component for the TPC.
 * 
 * Component ID: \b CaloDigitPublisher <br>
 * Library: \b libAliHLTCalo.
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments:<br>
 *
 *
 * @ingroup alihlt_system
 */
class AliHLTCaloDigitPublisherComponent : public AliHLTOfflineDataSource {
 public:
  /** standard constructor */
  AliHLTCaloDigitPublisherComponent();
  /** destructor */
  virtual ~AliHLTCaloDigitPublisherComponent();

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
   * @param [out] constBase  Additive part, independent of the input data volume.
   * @param [out] inputMultiplier  Multiplication ratio.
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

  using AliHLTOfflineDataSource::GetEvent;

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

 private:
   
   /** Digit handler */
   AliHLTCaloDigitHandler *fDigitHandler;
   
   /** Data specification */
   AliHLTUInt32_t fSpecification;
   
   /** Data type */
   AliHLTComponentDataType fDataType;

   /** Module number */
   Int_t fModule;
   
  /** copy constructor prohibited */
  AliHLTCaloDigitPublisherComponent(const AliHLTCaloDigitPublisherComponent&);
  /** assignment operator prohibited */
  AliHLTCaloDigitPublisherComponent& operator=(const AliHLTCaloDigitPublisherComponent&);

  ClassDef(AliHLTCaloDigitPublisherComponent, 0);
};

#endif
