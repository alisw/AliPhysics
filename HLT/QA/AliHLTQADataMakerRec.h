//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTQADATAMAKERREC_H
#define ALIHLTQADATAMAKERREC_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTQADataMakerRec.h
    @author Matthias Richter
    @date   2010-03-10
    @brief  Steering class for the HLT offline QA
*/

#include "AliHLTQADataMakerBase.h"
#include "TList.h"

/**
 * @class AliHLTQADataMakerRec
 * Steering class for HLT QA for reconstruction.
 *
 * HLT QA allows to define multiple HLT detector QA plugins. Each plugin
 * inherits through AliHLTQADataMakerBase from AliQADataMakerRec. Currently only
 * the AliQADataMakerRec interface is supported in the HLT QA. It seems that
 * AliQADataMakerSim is not relevant for HLT QA. However if so and at some
 * point it is required please inform the author of this class.
 *
 * AliHLTQADataMakerRec keeps a list of detector plugins and redirects the
 * different QA calls to all the plugins. At EndOfDetectorCycle all histograms
 * are collected from the plugins after EndOfDetectorCycle has been invoked
 * for every plugin.
 *
 * Detector plugins are added via the AliHLTModuleAgent. The optional function
 * AliHLTModuleAgent::GetQAPlugins() has to return a string of blank separated
 * class names.
 * 
 * HLT QA requires access to both the Esd and HLTEsd objects. Therefore the
 * Exec function is overloaded in AliHLTQADataMakerRec. A specific hnadling
 * in AliQAManager::RunOneEvent makes sure that an array of those objects
 * is passed, the call is then redirected to 
 * MakeESDs(AliESDEvent*, AliESDEvent*). Please note that the standard function
 * MakeESDs(AliESDEvent*) is usually not the place for HLT QA.
 */
class AliHLTQADataMakerRec: public AliHLTQADataMakerBase {

public:

  AliHLTQADataMakerRec();
  virtual ~AliHLTQADataMakerRec();

protected:
  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray** list);
  virtual void MakeRaws(AliRawReader * rawReader);
  virtual void MakeESDs(AliESDEvent * esd, AliESDEvent* hltesd);

  /// iterate over available agents and query class names of plugins
  int LoadAgents();

  /// load plugins from list of blank separated class names
  int LoadPlugins(const char* plugins=NULL);

  enum {
    kDigitsListInit    = 0x1,
    kESDsListInit      = 0x2,
    kRawsListInit      = 0x4,
    kRecPointsListInit = 0x8
  };

private:
  /** copy constructor prohibited */
  AliHLTQADataMakerRec(const AliHLTQADataMakerRec&);   
  /** assignment operator prohibited */
  AliHLTQADataMakerRec& operator = (const AliHLTQADataMakerRec&);

  TList fPlugins; //! list of HLT module QA plugins
  unsigned fFlags; //!

  ClassDef(AliHLTQADataMakerRec,0)  // HLT Quality Assurance Data Maker for reconstruction
};

#endif // ALIHLTQADATAMAKERREC_H
