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

class AliHLTEsdManager;

/**
 * @class AliHLTEsdCollectorComponent
 * The EsdCollector component merges ESDs from multiple events into one
 * ESD file per origin using the AliHLTEsdManager class.
 * \b Note: The component just merges ESDs of the same type/origin from
 * multiple events into one file. It does not implement merging of ESDs
 * from one event but several origins.
 *
 * Component ID: \b EsdCollector <br>
 * Library: \b libAliHLTUtil.so
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * The only AliHLTFileWriter argument of relevance is the -directory argument.
 *
 * See AliHLTFileWriter for full list of arguments.
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

  ClassDef(AliHLTEsdCollectorComponent, 0)
};
#endif
