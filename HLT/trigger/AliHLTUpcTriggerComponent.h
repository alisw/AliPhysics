//-*- Mode: C++ -*-
// $Id$
//* This file is property of and copyright by the ALICE HLT Project       * 
//* ALICE Experiment at CERN, All rights reserved.                        *
//* See cxx source for full Copyright notice                              *
/// @file   AliHLTUpcTriggerComponent.h
/// @author Kyrre Skjerdal
/// @date   2010-04-16
/// @brief  HLT trigger component for Ultra-Peripheral Collisions (UPC)

#ifndef ALIHLTUPCTRIGGERCOMPONENT_H
#define ALIHLTUPCTRIGGERCOMPONENT_H

#include "AliHLTTrigger.h"
//#include "AliESDEvent.h"
class AliESDEvent;

/**
 * @class  AliHLTUpcTriggerComponent
 *
 * HLT trigger component for Ultra-Peripheral Collisions
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b UpcTrigger                             <br>
 * Library: \b libAliHLTTrigger.so                                        <br>
 * Input Data Types:  kAliHLTAllDataTypes                              <br>
 * Output Data Types: ::kAliHLTAnyDataType                                <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * By default, configuration is loaded from OCDB, can be overridden by
 * component arguments.
 *
 * <h2>Default CDB entries:</h2>
 *
 *
 * <h2>Performance:</h2>
 * 
 *
 * <h2>Memory consumption:</h2>
 * 
 *
 * <h2>Output size:</h2>
 * 
 *
 * \ingroup alihlt_trigger_components
 */

class AliHLTUpcTriggerComponent : public AliHLTTrigger
{
public:
  virtual ~AliHLTUpcTriggerComponent(){;}
  //Returning the trigger name (UpcTrigger)
  virtual const char* GetTriggerName() const;
  //Returning a new object of the class
  virtual AliHLTComponent* Spawn(); 
  //Gives the output data size
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

private:
  //The trigger code
  virtual int DoTrigger();
  //Check if the event has a reconstructed primary vertex
  Bool_t PrimaryVertexReconstructed(const AliESDEvent *event) const;
  ClassDef(AliHLTUpcTriggerComponent, 0)
};

#endif
