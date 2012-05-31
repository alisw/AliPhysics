// $Id$

#ifndef ALIHLTPLUGINBASE_H
#define ALIHLTPLUGINBASE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTPluginBase.h
    @author Matthias Richter
    @date   
    @brief  Base class for AliRoot HLT plugins.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TObject.h"

class AliHLTSystem;

/**
 * @class AliHLTPluginBase
 * Base class for AliRoot HLT plugins.
 *
 * In order to allow definitions of HLT chains in the same macro as
 * the simulation/reconstruction, a global instance of AliHLTSystem
 * is required to make the registration of configurations work.
 *
 * AliHLTPlugin, AliRawReaderHLT and AliHLTSimulation all use
 * the global AliHLTSystem instance hosted by this base class.
 */
class AliHLTPluginBase {
 public:
  AliHLTPluginBase();
  /** destructor */
  virtual ~AliHLTPluginBase();

  /**
   * Init the global AliHLTSystem instance.
   */
  static void InitInstance();

  /**
   * Get the global AliHLTSystem instance.
   */
  static AliHLTSystem* GetInstance();

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTPluginBase(const AliHLTPluginBase& src);
  /** assignment operator prohibited */
  AliHLTPluginBase& operator=(const AliHLTPluginBase& src);

  static AliHLTSystem* fpSystem; //! HLT steering object

  static int fNofInstances;

  ClassDef(AliHLTPluginBase, 0)   // base class for the HLT reconstruction
};

#endif
