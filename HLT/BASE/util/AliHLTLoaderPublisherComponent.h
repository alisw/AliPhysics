//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTLOADERPUBLISHERCOMPONENT_H
#define ALIHLTLOADERPUBLISHERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTLoaderPublisherComponent.h
    @author Matthias Richter
    @date   
    @brief  A general tree publisher component for the AliLoader.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTOfflineDataSource.h"

class AliLoader;

/**
 * @class AliHLTLoaderPublisherComponent
 * A general tree publisher component for the AliLoader.
 * 
 * Component ID: \b AliLoaderPublisher <br>
 * Library: \b libAliHLTUtil.
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formating -->
 * \li -loader       <i> loader name      </i>
 *      e.g. <tt> -loader TPCLoader </tt>
 * \li -tree         <i> tree name </i> : digits (default), clusters     
 *      e.g. <tt> -tree digits </tt>
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
class AliHLTLoaderPublisherComponent : public AliHLTOfflineDataSource {
 public:
  /** standard constructor */
  AliHLTLoaderPublisherComponent();
  /** destructor */
  virtual ~AliHLTLoaderPublisherComponent();

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
   * @return
   */
  int GetEvent(const AliHLTComponentEventData& evtData,
	       AliHLTComponentTriggerData& trigData);

  using AliHLTOfflineDataSource::GetEvent;

 private:
  /** copy constructor prohibited */
  AliHLTLoaderPublisherComponent(const AliHLTLoaderPublisherComponent&);
  /** assignment operator prohibited */
  AliHLTLoaderPublisherComponent& operator=(const AliHLTLoaderPublisherComponent&);

  /**
   * Get tree of type specified in fTreeType from loader.
   */
  TTree* GetTree();

  /** max output block size, estimated during DoInit */
  Int_t                   fMaxSize;                                //!transient

  /** loader string */
  TString                 fLoaderType;                             //!transient

  /** tree string */
  TString                 fTreeType;                               //!transient

  /** be verbose: info printouts */
  Bool_t                  fVerbose;                                //!transient

  /** data type */
  AliHLTComponentDataType fDataType;                               //!transient

  /** data specification */
  AliHLTUInt32_t          fSpecification;                          //!transient

  /** instance of the AliLoader */
  AliLoader*              fpLoader;                                //!transient

  ClassDef(AliHLTLoaderPublisherComponent, 0);
};

#endif
