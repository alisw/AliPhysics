// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCTRACKHISTOCOMPONENT_H
#define ALIHLTTPCTRACKHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCTrackHistoComponent.h
    @author Gaute Ovrebekk
    @date   
    @brief  Component for ploting charge in clusters
*/

#include "AliHLTProcessor.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackSegmentData.h"

class AliHLTTPCConfMapper;
class TH1F;

/**
 * @class AliHLTTPCTrackHistoComponent
 * Component for ploting proparties of Tracks. 
 * There has to be uesd one argument, or the component will not plot anything.
 * 
 * <h2>General properties:</h2> 
 *
 * Component ID: \b TPCTrackHisto <br>
 * Library: \b libAliHLTTPC.so <br>
 * Input Data Types: AliHLTTPCDefinitions::fgkClustersDataType,
 *                   AliHLTTPCDefinitions::fgkTrackSegmentsDataType or
 *                   AliHLTTPCDefinitions::fgkTracksDataType <br>
 * Output Data Types: @ref kAliHLTDataTypeHistogram <br> 
 *
 * <h2> Mandatory arguments: </h2>
 * \li One of the Optional arguments.
 * 
 * <h2> Optional arguments: </h2>
 * 
 * \li -plot-All <br>
 *      Plots all the Histograms (default kFALSE) 
 * \li -plot-nClusters <br>
 *      Plots Number of Clusters on Tracks (default kFALSE)
 * \li -plot-ChargeClusters <br>
 *      Plots Charge on all Clusters (default kFALSE)
 * \li -plot-ChargeUsedClusters <br>
 *      Plots Charge on Clusters used for Tracks (default kFALSE)
 * \li -plot-pT <br>
 *      Plots pT for Tracks (default kFALSE)
 * \li -plot-Residuals <br>
 *      Plots Residual of Tracks (default kFALSE)
 * \li -plot-Tgl <br>
 *      Plots Tgl for tracks (default kFALSE)
 * \li -plot-NClusters <br>
 *      Plots the number of clusters in the event (default kFALSE)
 * \li -plot-NUsedClusters <br>
 *      Plots th number of Used clusters in the event (default kFALSE)
 * \li -plot-NTracks <br>
 *      Plots the number of Tracks in the event (default kFALSE)
 * \li -plot-QMaxAll <br>
 *      Plots the Q Max for all clusters in the event (default kFALSE)
 * \li -plot-QMaxUsed <br>
 *      Plots the Q Max for clusters used on tracks (default kFALSE)
 * \li -reset-plots <br>
 *      Will reset the plots for every event (default kFALSE)
 *
 * <h2>Configuration:</h2>
 * 
 *
 * <h2>Default CDB entries:</h2>
 * The component has for now no CDB entry
 *
 * <h2>Performance:</h2>
 * Not Tested 
 *
 * <h2>Memory consumption:</h2>
 * Not Tested
 *
 * <h2>Output size:</h2>
 * The size of an histogram (588 bit) * the number of histograms you plot
 *
 * 
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCTrackHistoComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTTPCTrackHistoComponent();
  /** destructor */
  virtual ~AliHLTTPCTrackHistoComponent();

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

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTTPCTrackHistoComponent(const AliHLTTPCTrackHistoComponent&);
  /** assignment operator prohibited */
  AliHLTTPCTrackHistoComponent& operator=(const AliHLTTPCTrackHistoComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  int Configure(const char* arguments);
 
  void ReadTracks(const AliHLTComponentBlockData* iter,Int_t &tt);

  void PushHisto();
 
  TH1F * fHistoNClustersOnTracks;                                  //! transient
  TH1F * fHistoChargeAllClusters;                                  //! transient
  TH1F * fHistoChargeUsedClusters;                                 //! transient
  TH1F * fHistoPT;                                                 //! transient
  TH1F * fHistoResidual;                                           //! transient
  TH1F * fHistoTgl;                                                //! transient
  TH1F * fHistoNClusters;                                          //! transient
  TH1F * fHistoNUsedClusters;                                      //! transient
  TH1F * fHistoNTracks;                                            //! transient 
  TH1F * fHistoQMaxAllClusters;                                    //! transient
  TH1F * fHistoQMaxUsedClusters;                                   //! transient
  
  Bool_t fPlotAll;                                                 //! transient 
  Bool_t fPlotNClustersOnTracks;                                   //! transient 
  Bool_t fPlotChargeClusters;                                      //! transient 
  Bool_t fPlotChargeUsedClusters;                                  //! transient 
  Bool_t fPlotPT;                                                  //! transient 
  Bool_t fPlotResidual;                                            //! transient 
  Bool_t fPlotTgl;                                                 //! transient 
  Bool_t fPlotNClusters;                                           //! transient
  Bool_t fPlotNUsedClusters;                                       //! transient
  Bool_t fPlotNTracks;                                             //! transient
  Bool_t fPlotQMaxClusters;                                        //! transient
  Bool_t fPlotQMaxUsedClusters;                                    //! transient
  Bool_t fResetPlots;                                              //! transient

  vector<AliHLTTPCSpacePointData> fClusters;                       //! transient
  vector<AliHLTTPCTrackSegmentData> fTracks;                       //! transient
  
  vector<UInt_t> fTrackClusterID[36][6];                           //! transient

  ClassDef(AliHLTTPCTrackHistoComponent, 0);

};
#endif
