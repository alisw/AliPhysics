// -*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALIHLTROOTFILEPUBLISHERCOMPONENT_H
#define ALIHLTROOTFILEPUBLISHERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTRootFilePublisherComponent.h
    @author Matthias Richter
    @date   
    @brief  component for publishing of Root objects from a root file.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTFilePublisher.h"
#include <TList.h>

/**
 * @class AliHLTRootFilePublisherComponent
 * An HLT data source component which publishes data from one or a sequence
 * of files.<br>
 *
 * Component ID: \b RootFilePublisherComponent <br>
 * Library: \b libHLTBase (in order to use the component from the external
 * interface, it might be necessary to specify a dummy library with the
 * \em -componentlibrary argument).
 *
 * Mandatory arguments: <br>
 *
 * Optional arguments:<br>
 *
 * @see AliHLTFilePublisher for default arguments
 * @ingroup alihlt_component
 */
class AliHLTRootFilePublisherComponent : public AliHLTFilePublisher  {
 public:
  /** standard constructor */
  AliHLTRootFilePublisherComponent();
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTRootFilePublisherComponent(const AliHLTRootFilePublisherComponent&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTRootFilePublisherComponent& operator=(const AliHLTRootFilePublisherComponent&);
  /** destructor */
  virtual ~AliHLTRootFilePublisherComponent();

  const char* GetComponentID();
  AliHLTComponentDataType GetOutputDataType();
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

  /**
   * Open all files.
   * Opens all files from the file name list @ref fFileNames and adds TFile
   * opjects to the TFiles list.
   */
  int OpenFiles();

 protected:
  /**
   * Data processing method for the component.
   * The component uses the @ref alihltcomponent-high-level-interface
   * to put serialized Root object into the output stream. Despite of that it
   * implements the lox-level DumpEvent method in order to allow child classes
   * to use the low-level method.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @param outputPtr	  pointer to target buffer
   * @param size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param outputBlocks  list to receive output block descriptors
   * @return
   */
  int GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks );

  /**
   * Scan one argument and adjacent parameters.
   * Can be overloaded by child classes in order to add additional arguments
   * beyond the standard arguments of the file publisher. The method is called
   * whenever a non-standard argument is recognized.
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  virtual int ScanArgument(int argc, const char** argv);

 private:

  ClassDef(AliHLTRootFilePublisherComponent, 0)
};
#endif
