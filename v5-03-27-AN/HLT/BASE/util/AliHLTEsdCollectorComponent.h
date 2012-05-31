// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTESDCOLLECTORCOMPONENT_H
#define ALIHLTESDCOLLECTORCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTEsdCollectorComponent.h
    @author Matthias Richter
    @date   
    @brief  Collect ESDs of multiple events and write toi file
*/

#include "AliHLTFileWriter.h"
#include "TString.h"

class AliHLTEsdManager;

/**
 * @class AliHLTEsdCollectorComponent
 * The EsdCollector component merges ESDs from multiple events into one
 * ESD file per origin using the AliHLTEsdManager class.
 * \b Note: The component just merges ESDs of the same type/origin from
 * multiple events into one file. It does not implement merging of ESDs
 * from one event but several origins.
 *
 * The file name of the ESD file is derived from the origin of the ESD data
 * block.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b EsdCollector                                      <br>
 * Library: \b libAliHLTUtil.so					      <br>
 * Input Data Types: kAliHLTDataTypeESDObject, kAliHLTDataTypeESDTree <br>
 * Output Data Types: none					      <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * The only AliHLTFileWriter argument of relevance is the \em -directory
 * argument. See AliHLTFileWriter for full list of arguments. Note: The
 * file name of the ESD file is derieved from the origin of the ESD
 * data block.
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
class AliHLTEsdCollectorComponent : public AliHLTFileWriter
{
 public:
  /** standard constructor */
  AliHLTEsdCollectorComponent();
  /** destructor */
  virtual ~AliHLTEsdCollectorComponent();

  /**
   * The id of the component.
   * @return component id (string)
   */
  const char* GetComponentID() {return "EsdCollector";};

  /**
   * Spawn function.
   * @return new class instance
   */
  AliHLTComponent* Spawn() {return new AliHLTEsdCollectorComponent;}

 protected:
  // interface functions
  int InitWriter();
  int CloseWriter();
  int DumpEvent( const AliHLTComponentEventData& evtData,
		 const AliHLTComponentBlockData* blocks, 
		 AliHLTComponentTriggerData& trigData );
  
  using AliHLTFileWriter::DumpEvent;
  int ScanArgument(int argc, const char** argv);

private:
  /** copy constructor prohibited */
  AliHLTEsdCollectorComponent(const AliHLTEsdCollectorComponent&);
  /** assignment operator prohibited */
  AliHLTEsdCollectorComponent& operator=(const AliHLTEsdCollectorComponent&);

  /** the ESD manager instance writes the ESDs */
  AliHLTEsdManager* fpManager; //! transient
  /** name of the tree for ESD storage */
  TString fTreeName; //! transient

  ClassDef(AliHLTEsdCollectorComponent, 0)
};
#endif
