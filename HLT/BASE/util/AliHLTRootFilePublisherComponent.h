// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTROOTFILEPUBLISHERCOMPONENT_H
#define ALIHLTROOTFILEPUBLISHERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTRootFilePublisherComponent.h
    @author Matthias Richter, Jochen Thaeder
    @date   
    @brief  component for publishing of Root objects from a root file.
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTFilePublisher.h"
#include <TList.h>

/**
 * @class AliHLTRootFilePublisherComponent
 * An HLT data source component which publishes root objects from one 
 * or a sequence of root files. Be aware, one root file can contain 
 * several root objects. Either all objects or just one object can be selected.<br>
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ROOTFilePublisher <br>
 * Library: \b libAliHLTUtil.so     <br>
 * Input Data Types: none <br>
 * Output Data Types: according to arguments <br>
 *
 * <h2>Mandatory arguments:</h2>
 * @see AliHLTFilePublisher for mandatory defaultarguments
 *
 * <h2>Optional arguments:</h2>
 * @see AliHLTFilePublisher for optional default arguments
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -objectname   <i> objectname    </i>
 *      Name of the object in the root file to be fetched. This is set for 
 *      all events/files. If not given, all objects are fetched.
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
 * The component needs at least one argument \em -datafile or \em -datafilelist.
 * Both can occur multiple times. The \em -datatype and \em -dataspec
 * parameters are valid for all files until the next occurrence of
 * \em -datatype/spec.
 * All files are published within one event, unless the \em -nexevent specifies
 * where to break into multiple events. Be aware, one root file can contain 
 * several root objects. If \em -objectname is not used to select one, all 
 * objects are all published with the same datatype and specification.
 *
 * @ingroup alihlt_util_components
 */

class AliHLTRootFilePublisherComponent : public AliHLTFilePublisher  {
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
  
  /** standard constructor */
  AliHLTRootFilePublisherComponent();

  /** destructor */
  virtual ~AliHLTRootFilePublisherComponent();
  
  /*
   * ---------------------------------------------------------------------------------
   * Public functions to implement AliHLTComponent's interface.
   * These functions are required for the registration process
   * ---------------------------------------------------------------------------------
   */

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

 protected:
  /*
   * ---------------------------------------------------------------------------------
   * Protected functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  /**
   * Data processing method for the component.
   * The component uses the @ref alihltcomponent-high-level-interface
   * to put serialized Root object into the output stream. Despite of that it
   * implements the low-level DumpEvent method in order to allow child classes
   * to use the low-level method.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @param outputPtr     not used
   * @param size          not used
   * @param outputBlocks  not used
   * @return
   */
  Int_t GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks);

  using AliHLTFilePublisher::GetEvent;

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
  virtual Int_t ScanArgument(Int_t argc, const char** argv);

 private:

  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  /** copy constructor prohibited */
  AliHLTRootFilePublisherComponent(const AliHLTRootFilePublisherComponent&);

  /** assignment operator prohibited */
  AliHLTRootFilePublisherComponent& operator=(const AliHLTRootFilePublisherComponent&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** The current event */
  TObjLink *fpCurrentEvent;                  //! transient

  /** Name of the object which should be fetched 
   *  from the root file.
   */
  TString   fObjectName;                     //! objectname

  ClassDef(AliHLTRootFilePublisherComponent, 0)
};
#endif
