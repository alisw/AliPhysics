//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDTRACKHISTOCOMPONENT_H
#define ALIHLTTRDTRACKHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *


#include "AliHLTProcessor.h"

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

class TH1F;
class TClonesArray;
class AliHLTTRDTrackHistoComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTTRDTrackHistoComponent();
  /** destructor */
  virtual ~AliHLTTRDTrackHistoComponent();

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
  AliHLTTRDTrackHistoComponent(const AliHLTTRDTrackHistoComponent&);
  /** assignment operator prohibited */
  AliHLTTRDTrackHistoComponent& operator=(const AliHLTTRDTrackHistoComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 

  AliHLTUInt32_t fOutputSize;   // output size
  AliHLTUInt32_t fSpec;         // accumulated specification
  TClonesArray* fTracksArray;   // input array

  TH1F *fClPerTrkl;             // Number of clusters per tracklet
  TH1F *fTrklPerTrk;            // Number of tracklets per track
  TH1F *fEvSize;                // Event size in kbyte
  TH1F *fEtaDistrib;            // Eta distribution
  TH1F *fPhiDistrib;            // Phi distribution
  TH1F *fPtDistrib;             // Pt distribution

  ClassDef(AliHLTTRDTrackHistoComponent, 0);
};
#endif
