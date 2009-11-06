//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTITSCLUSTERHISTOCOMPONENT_H
#define ALIHLTITSCLUSTERHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTITSClusterHistoComponent.h
    @author Gaute Ovrebekk
    @brief  Component for ploting clusters
*/

#include "AliHLTProcessor.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliHLTITSSpacePointData.h"
#include "TClonesArray.h"
#include "AliITSRecPoint.h"

/**
 * @class AliHLTITSQHistoComponent
 * Component for ploting charge in clusters
 * 
 * Component ID: \b ITSClusterHisto <br>
 * Library: \b libAliHLTITS.so
 *
 * Mandatory arguments: <br>
 * 
 * 
 * Optional arguments: <br>
 * \li -plot-all <br>  
 *      will plot all histograms in the component.
 * \li -plot-xy <br>  
 *      will plot a historgam of the x and y projection of the clusters.
 * \li -plot-charge     <i> multiplicity  </i> <br> 
 *      will plot the charge of the clusters.
 * \li -plot-phieta <br>
 *      will plot the phi vs eta distrubution of the clusters.
 * @ingroup alihlt_tpc_components
 */
class AliHLTITSClusterHistoComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTITSClusterHistoComponent();
  /** destructor */
  virtual ~AliHLTITSClusterHistoComponent();

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
  AliHLTITSClusterHistoComponent(const AliHLTITSClusterHistoComponent&);
  /** assignment operator prohibited */
  AliHLTITSClusterHistoComponent& operator=(const AliHLTITSClusterHistoComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  int Configure(const char* arguments);
  
  TH2F * fXY;                              //! transient
  TH2F **fPhieta;                          //! transient
  TH1F * fCharge;                          //! transient
    
  Bool_t fPlotCharge;                      //! transient
  Bool_t fPlotXYPhiEta;                    //! transient

  /// maximum required size of the output buffer
  /// dynamically adjusted
  int fOutputSize;                         //! transient
   
  ClassDef(AliHLTITSClusterHistoComponent, 3);

};
#endif
