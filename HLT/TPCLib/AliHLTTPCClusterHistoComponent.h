// $Id$
#ifndef ALIHLTTPCCLUSTERHISTOCOMPONENT_H
#define ALIHLTTPCCLUSTERHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCQHistoComponent.h
    @author Gaute Ovrebekk
    @brief  Component for ploting charge in clusters
*/

#include "AliHLTProcessor.h"
#include "TH1F.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackSegmentData.h"

class AliHLTTPCConfMapper;

/**
 * @class AliHLTTPCQHistoComponent
 * Component for ploting charge in clusters
 * 
 * Component ID: \b TPCQHisto <br>
 * Library: \b libAliHLTTPC.
 *
 * Mandatory arguments: <br>
 * 
 * 
 * Optional arguments: <br>
 * 
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCClusterHistoComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTTPCClusterHistoComponent();
  /** destructor */
  virtual ~AliHLTTPCClusterHistoComponent();

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
  int DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData );

  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTTPCClusterHistoComponent(const AliHLTTPCClusterHistoComponent&);
  /** assignment operator prohibited */
  AliHLTTPCClusterHistoComponent& operator=(const AliHLTTPCClusterHistoComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  int Configure(const char* arguments);
  
  TH1F * fTotalClusterChargeOROCAll;                               //! transient
  TH1F * fTotalClusterChargeIROCAll;                               //! transient
  TH1F * fTotalClusterChargeROCSelection;                          //! transient
  TH1F * fTotalClusterChargePartitionSelection;                    //! transient
  TH1F * fQMaxPartitionAll;                                        //! transient
  TH1F * fQMaxROCAll;                                              //! transient
  TH1F * fNumberOfClusters;                                        //! transient
  
  Bool_t fPlotChargeOROCAll;                                       //! transient
  Bool_t fPlotChargeIROCAll;                                       //! transient
  Bool_t fPlotChargeROCSel;                                        //! transient
  Bool_t fPlotChargePartSel;                                       //! transient
  Bool_t fPlotQmaxPartAll;                                         //! transient
  Bool_t fPlotQmaxROCAll;                                          //! transient
  Bool_t fPlotNClusters;                                           //! transient
 
  ClassDef(AliHLTTPCClusterHistoComponent, 0);

};
#endif
