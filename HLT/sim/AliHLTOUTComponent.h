// @(#) $Id$

#ifndef ALIHLTOUTCOMPONENT_H
#define ALIHLTOUTCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTOUTComponent.h
    @author Matthias Richter
    @date   
    @brief  The HLTOUT data sink component similar to HLTOUT nodes
*/

// see class description below
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTOfflineDataSink.h"
#include <TString.h>

/**
 * @class AliHLTOUTComponent
 * The HLTOUT data sink component which models the behavior of the HLTOUT
 * nodes of the HLT cluster.
 * @ingroup alihlt_component
 */
class AliHLTOUTComponent : public AliHLTOfflineDataSink  {
 public:
  /** standard constructor */
  AliHLTOUTComponent();
  /** destructor */
  virtual ~AliHLTOUTComponent();

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponent* Spawn();

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
   * Data processing method for the component.
   * The function can be overloaded by other file writer components.
   * @param evtData       event data structure
   * @param blocks        input data block descriptors
   * @param trigData	  trigger data structure
   */
  int DumpEvent( const AliHLTComponentEventData& evtData,
		 const AliHLTComponentBlockData* blocks, 
		 AliHLTComponentTriggerData& trigData );

  /**
   * Fill ESD for one event.
   * @param eventNo       event No. \em Note: this is an internal enumeration of the
   *                      processed events.
   * @param runLoader     the AliRoot runloader
   * @param esd           an AliESDEvent instance
   * @return neg. error code if failed 
   */
  int FillESD(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd);

 private:
  /** copy constructor prohibited */
  AliHLTOUTComponent(const AliHLTOUTComponent&);
  /** assignment operator prohibited */
  AliHLTOUTComponent& operator=(const AliHLTOUTComponent&);

  ClassDef(AliHLTOUTComponent, 0)
};
#endif
