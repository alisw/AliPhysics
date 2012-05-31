//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTCOMPONENTCONFIGURATION_H
#define ALIHLTCOMPONENTCONFIGURATION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTComponentConfiguration.h
/// @author Matthias Richter
/// @date   2010-11-26
/// @brief  HLT configuration description for a single component.
/// @note   The class is used in Offline (AliRoot) context

#include "AliHLTConfiguration.h"

/**
 * @class AliHLTComponentConfiguration
 * @brief Description of an HLT component configuration.
 *
 * In addition to the base class AliHLTConfiguration the additional online
 * component parameters are available in this class
 * - node name
 * - shared memory settings
 *
 * This class is only used in the HLT offline environment, see @ref alihlt_system
 * for more details.
 *
 * @ingroup alihlt_system
 */
class AliHLTComponentConfiguration : public AliHLTConfiguration {
 public:
  /**
   * standard constructor. The configuration is automatically registered in the
   * global configuration manager
   */
  AliHLTComponentConfiguration();
  /**
   * constructor. The configuration is automatically registered in the
   * global configuration manager
   * @param id         unique id of the configuration
   * @param component  component id
   * @param sources    blank separated list of source configuration ids
   * @param arguments  argument string passed to the component at initialization
   */
  AliHLTComponentConfiguration(const char* id, const char* component,
		      const char* sources, const char* arguments);
  /** copy constructor */
  AliHLTComponentConfiguration(const AliHLTComponentConfiguration& src);
  /** assignment op */
  AliHLTComponentConfiguration& operator=(const AliHLTComponentConfiguration& src);
  /** destructor */
  virtual ~AliHLTComponentConfiguration();

  /**
   * Return the component library.
   */
  const char* GetComponentLibrary() const {return fLibrary.Data();}

  /**
   * Return the online command.
   */
  const char* GetOnlineCommand() const {return fOnlineCommand.Data();}

  /**
   * Return the online nodes.
   */
  const char* GetNodeSettings() const {return fNodeNames.Data();}
  
  void SetComponentLibrary(const char* library) {fLibrary=library;}

  void SetNodeNames(const char* nodes) {fNodeNames=nodes;}

  void AddNode(const char* node) {
    if (!node) return;
    if (!fNodeNames.IsNull()) fNodeNames+=" "; fNodeNames+=node;
  }

  /// set the online command string
  void SetOnlineCommand(const char* cmd);

  /**
   * overloaded from AliHLTConfiguration
   */
  virtual void PrintStatus() const;

  /**
   * overloaded from AliHLTConfiguration
   * options:
   *   status  - print status
   */
  virtual void Print(const char* option="") const;


 protected:
  
 private:
  /// component library
  TString fLibrary; // component library

  /// list of nodes of the component instances
  TString fNodeNames; // list of node names

  /// original command in the online configuration
  TString fOnlineCommand; // original online command

  ClassDef(AliHLTComponentConfiguration, 1);
};

#endif
