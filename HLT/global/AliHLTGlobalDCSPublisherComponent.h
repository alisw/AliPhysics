//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTGLOBALDCSPUBLISHERCOMPONENT_H
#define ALIHLTGLOBALDCSPUBLISHERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalDCSPublisherComponent.h
    @author Matthias Richter
    @date   20010-03-10
    @brief  DIM publisher component for global HLT data
*/
#include "AliHLTDataSink.h"

class AliHLTDimServer;

/**
 * @class AliHLTGlobalDCSPublisherComponent
 * DIM Publisher component for global HLT data.
 * It implements a DIM server which publishes global HLT data through the
 * following services:
 * - Vertex_X
 * - Vertex_Y
 * - Vertex_Z
 * - ResVertex_X
 * - ResVertex_Y
 * - ResVertex_Z
 * 
 */
class AliHLTGlobalDCSPublisherComponent : public AliHLTDataSink {
public:
  AliHLTGlobalDCSPublisherComponent();
  ~AliHLTGlobalDCSPublisherComponent();

  virtual const char* GetComponentID();
  virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  virtual AliHLTComponent* Spawn();

 protected:
  /// component initialization
  int DoInit( int argc, const char** argv );

  /// component cleanup
  int DoDeinit();

  /// Data processing method for the component.
  virtual int DumpEvent( const AliHLTComponentEventData& evtData,
			 AliHLTComponentTriggerData& trigData );

  using AliHLTDataSink::DumpEvent;

  int ScanConfigurationArgument(int argc, const char** argv);

private:
  /// copy constructor not permitted
  AliHLTGlobalDCSPublisherComponent(const AliHLTGlobalDCSPublisherComponent&);
  /// assignment operator not permitted
  AliHLTGlobalDCSPublisherComponent& operator=(const AliHLTGlobalDCSPublisherComponent&);

  AliHLTDimServer* fpServer; //! dim server instance

  ClassDef(AliHLTGlobalDCSPublisherComponent, 0)
};
#endif
