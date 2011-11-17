#ifndef ALIHLTMUONCLUSTERHISTOCOMPONENT_H
#define ALIHLTMUONCLUSTERHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTMUONClusterHistoComponent.h
    @author Arshad Ahmad <Arshad.Ahmad@cern.ch>
    @date  27 Nov 2009
    @brief  Component for ploting charge in clusters
*/

#include "TH1F.h"
#include "AliHLTMUONProcessor.h"
/**
 * @class AliHLTMUONClusterHistoComponent
 * Component for ploting charge in clusters
 * 
 * Component ID: \b MUONHisto <br>
 * Library: \b libAliHLTMUON.so <br>
 *
 * Mandatory arguments: <br>
 * 
 * 
 * Optional arguments: <br>
 * 
 *
 * @ingroup alihlt_dimuon_components
 */

class AliHLTMUONClusterHistoComponent : public AliHLTMUONProcessor
{
public:
  /** default constructor */
  AliHLTMUONClusterHistoComponent();
  /** destructor */
  virtual ~AliHLTMUONClusterHistoComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

private:

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component.

  /** interface function, see AliHLTComponent for description */
  int DoInit( int argc, const char** argv );
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see AliHLTComponent for description */
  /*   int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData ); */
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
              AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
              AliHLTUInt32_t& size, AliHLTComponentBlockDataList& outputBlocks);

  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;
  
  /** copy constructor prohibited */
  AliHLTMUONClusterHistoComponent(const AliHLTMUONClusterHistoComponent&);
  /** assignment operator prohibited */
  AliHLTMUONClusterHistoComponent& operator=(const AliHLTMUONClusterHistoComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  int Configure(const char* arguments);

  
  TH1F * fChargePerClusterBending;                                 //! transient
  TH1F * fChargePerClusterNonBending;                              //! transient
  TH1F * fNumberOfClusters;                                        //! transient
  
  Bool_t fPlotChargePerClusterBending;                             //! transient
  Bool_t fPlotChargePerClusterNonBending;                          //! transient
  Bool_t fPlotNClusters;                                           //! transient
 
  ClassDef(AliHLTMUONClusterHistoComponent, 0);

};
#endif //ALIHLTMUONCLUSTERHISTOCOMPONENT_H
