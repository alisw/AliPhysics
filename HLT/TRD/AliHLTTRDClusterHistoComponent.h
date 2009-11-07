//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDCLUSTERHISTOCOMPONENT_H
#define ALIHLTTRDCLUSTERHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *


#include "AliHLTProcessor.h"
#include "TH1D.h"
#include "TH2F.h"

/**
 * @class AliHLTTRDQHistoComponent
 * Component for ploting charge in clusters
 * 
 * Component ID: \b TRDQHisto <br>
 * Library: \b libAliHLTTRD.
 *
 * Mandatory arguments: <br>
 * 
 * 
 * Optional arguments: <br>
 * 
 *
 * @ingroup alihlt_tpc_components
 */
class TClonesArray;
class AliHLTTRDClusterHistoComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTTRDClusterHistoComponent();
  /** destructor */
  virtual ~AliHLTTRDClusterHistoComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

protected:

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 

  /** interface function, see AliHLTComponent for description */
  int DoInit( int argc, const char** argv );
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see AliHLTComponent for description */
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

  using AliHLTProcessor::DoEvent;

  int Configure(const char* arguments);
  
private:
  /** copy constructor prohibited */
  AliHLTTRDClusterHistoComponent(const AliHLTTRDClusterHistoComponent&);
  /** assignment operator prohibited */
  AliHLTTRDClusterHistoComponent& operator=(const AliHLTTRDClusterHistoComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 

  AliHLTUInt32_t fOutputSize;   // output size
  TClonesArray* fClusterArray;  // input array

  TH1D *fNClsDet;
  TH1D *fClsAmp;
  TH1D *fClsAmpDrift;
  TH1D *fClsTB;

  TH1D *fClsAmpDriftDet[540];
  TH1D *fClsAmpDist; 

  TH1D *fSClsDist;
  TH1D *fNScls;

  // kryptogramm
  TH2F fSlidingWindow[540];
  TH2F *fClusterDist;
  TH1D *fClusterCandCharge;

  ClassDef(AliHLTTRDClusterHistoComponent, 0);
};
#endif
