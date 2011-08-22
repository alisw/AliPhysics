//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCDEDXMONITORINGCOMPONENT_H
#define ALIHLTTPCDEDXMONITORINGCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCdEdxMonitoringComponent.h
/// @author Per-Ivar Lønne, Jochen Thaeder, Matthias Richter, Alexander Kalweit
/// @date   21.08.2011
/// @brief  Component for reading ESD from chain and produce a dEdx monitoring plot
///

/**
 * @class AliHLTTPCdEdxMonitoringComponent
 * A component meant to read ESD-files from the chain online
 * and produce a plot for dEdx monitoring.
 * 
 * description of what the component does in more detail
 * bla
 * bla
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCdEdxMonitoring <br>
 * Library: \b libAliHLTTPC.so <br>
 * Input Data Types: @ref kAliHLTDataTypeESDObject|kAliHLTDataOriginAny <br>
 * Output Data Types: @ref kAliHLTDataTypeHistogram|kAliHLTDataOriginHLT <br>
 *
 * <h2>Mandatory arguments:</h2>
 *
 * <h2>Optional arguments:</h2>
 * 
 * <h2>Configuration</h2>
 * \li -xbins <i> fxbins </i> <br>
 *     number of bins on x-axis
 *
 * <h2>Configuration</h2>
 * \li -xbins <i> fxmin </i> <br>
 *     minimum value of x-axis
 *
 * <h2>Configuration</h2>
 * \li -xbins <i> fxmax </i> <br>
 *     maximum value of x-axis
 *
 * <h2>Configuration</h2>
 * \li -ybins <i> fybins </i> <br>
 *     number of bins on y-axis
 *
 * <h2>Configuration</h2>
 * \li -ybins <i> fymin </i> <br>
 *     minimum value of y-axis
 *
 * <h2>Configuration</h2>
 * \li -ybins <i> fymax </i> <br>
 *     maximum value of y-axis
 *
 * <h2>Default CDB entries:</h2>
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory Consumption:</h2>
 *
 * <h2>Output size:</h2>
 * 4096
 *
 * @ingroup alihlt
 */

#include "AliHLTProcessor.h"

class AliESDtrackCuts; // For setting track cuts
class TH2F;

class AliHLTTPCdEdxMonitoringComponent : public AliHLTProcessor 
{
public:

  /*----------------------------------------------------*
   *          Constructor and destructor 
   *----------------------------------------------------*/
    
  AliHLTTPCdEdxMonitoringComponent();

  virtual ~AliHLTTPCdEdxMonitoringComponent();


  // AliHLTComponent interface functions
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  void GetOutputDataSize( unsigned long& constBase, Double_t& inputMultiplier );
  void GetOCDBObjectDescription( TMap* const targetMap);

  // Spawn function, return new class instance
  AliHLTComponent* Spawn();

 protected:
  // AliHLTComponent interface functions
  Int_t DoInit( Int_t argc, const char** argv );
  Int_t DoDeinit();
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  using AliHLTProcessor::DoEvent;
  Int_t ScanConfigurationArgument(Int_t argc, const char** argv);
  Int_t Reconfigure(const char* cdbEntry, const char* chainId);
  Int_t ReadPreprocessorValues(const char* modules);

private:
  /** copy constructor prohibited */
  AliHLTTPCdEdxMonitoringComponent(const AliHLTTPCdEdxMonitoringComponent&);
  /** assignment operator prohibited */
  AliHLTTPCdEdxMonitoringComponent& operator=(const AliHLTTPCdEdxMonitoringComponent&);

  // Sets standard trackcuts
  void SetDefaultConfiguration();

  // trackcuts
  AliESDtrackCuts *fESDTrackCuts; //! transistent

  // histogram
  TH2F *fHist; //! transistent
  Int_t fxbins; //! transistent
  Double_t fxmin; //! transistent
  Double_t fxmax; //! transistent
  Int_t fybins; //! transistent
  Double_t fymin; //! transistent
  Double_t fymax; //! transistent
  

  // sets some plotstyles
  void Plotstyle();

  ClassDef(AliHLTTPCdEdxMonitoringComponent, 0)
};
#endif
