//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTSAMPLEMONITORINGCOMPONENT_H
#define ALIHLTSAMPLEMONITORINGCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTSampleMonitoringComponent.h
    @author Matthias Richter
    @date   
    @brief  A sample monitoring component for the HLT.
*/

#include "AliHLTProcessor.h"

class TH1F;

/**
 * @class AliHLTSampleMonitoringComponent
 * An HLT sample component.
 * This component does not any data processing at all. It just
 * illustrates the existence of several components in ine library and
 * allows to set up a very simple chain with different components.
 *
 * The component generates two histograms and accumulates random data
 * as events are coming in. The two histograms are sent out via tree
 * different approaches, see output data types.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b Sample-MonitoringComponent <br>
 * Library: \b libAliHLTSample.so     <br>
 * Input Data Types: @ref kAliHLTAnyDataType <br>
 *    in fact, the component ignores all incoming data blocks
 *
 * Output Data Types:
 * \li {ROOT_TH1,EXPL}                                              <br>
 *     one histogram per data block, specification identifies the
 *     specific histogram.
 * \li {ROOTOBJA,EXPL}                                              <br>
 *     the two histograms are added to a TOBjArray which is pushed
 *     to the output stream
 * \li {ROOTTREE,EXPL}                                              <br>
 *     the two histograms are added to a TTree which is pushed
 *     to the output stream
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -push-histograms                              <br>
 *      push histograms idividually
 * \li -push-ttree (default)                         <br>
 *      push histograms embedded in TTree
 * \li -push-array                                   <br>
 *      push histograms embedded in TObjArray
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * The component has no default CDB entry.
 * It does not load any configuration from the global <tt>ConfigHLT</tt>
 * folder.
 * 
 * For re-configuration/steering the component expects a TObjString object
 * holding a string with the configuration parameters explained above.
 *
 * <h2>Performance:</h2>
 * The component does not any event data processing.
 *
 * <h2>Memory consumption:</h2>
 * The component does not any event data processing.
 *
 * <h2>Output size:</h2>
 * The component has no output data.
 *
 * The component implements the @ref alihltcomponent-low-level-interface.
 * for data processing.
 *
 * Using the latter case, a component can also be reconfigured. Special
 * events are propageted through the chain in order to trigger the re-
 * configuration. The component needs to implement the function
 * @ref Reconfigure(). The simplest version of a configuration object is
 * a string object (TObjString) containing configuration arguments.
 *
 * @ingroup alihlt_tutorial
 */
class AliHLTSampleMonitoringComponent : public AliHLTProcessor {
public:
  AliHLTSampleMonitoringComponent();
  virtual ~AliHLTSampleMonitoringComponent();

  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

 protected:
  // AliHLTComponent interface functions
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;

private:
  /** copy constructor prohibited */
  AliHLTSampleMonitoringComponent(const AliHLTSampleMonitoringComponent&);
  /** assignment operator prohibited */
  AliHLTSampleMonitoringComponent& operator=(const AliHLTSampleMonitoringComponent&);

  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   *
   * This function illustrates the scanning of an argument string. The string
   * was presumably fetched from the CDB.
   */
  int Configure(const char* arguments);

  /** send histograms individually as data blocks */
  bool fPushHistograms; // !transient

  /** send histograms TTree embedded */
  bool fPushTTree; // !transient

  /** send histograms Object Array embedded */
  bool fPushTObjArray; // !transient

  /** test histogram */
  TH1F* fHpx; //!transient

  /** test histogram */
  TH1F* fHpy; //!transient

  ClassDef(AliHLTSampleMonitoringComponent, 0)
};
#endif
