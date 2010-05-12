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
    @brief  Component for track histo
*/

#include "AliHLTProcessor.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackSegmentData.h"

class TNtuple;
class TH1F;
class TProfile;

/**
 * @class AliHLTTPCTrackHistoComponent
 * Component for plotting proparties of Tracks. 
 * The component gives out 2 NTuples. One for cluster and one for track properties
 * 
 * <h2>General properties:</h2> 
 *
 * Component ID: \b TPCTrackHisto <br>
 * Library: \b libAliHLTTPC.so <br>
 * Input Data Types: AliHLTTPCDefinitions::fgkClustersDataType,
 *                   AliHLTTPCDefinitions::fgkTracksDataType <br>
 * Output Data Types: ::kAliHLTDataTypeTNtuple <br> 
 *
 * <h2> Mandatory arguments: </h2>
 * \li No mandaroty arguments. 
 * 
 * <h2> Optional arguments: </h2>
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
 * Size variables in Ntuple
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
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();
  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);

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
  /** inherited from AliHLTComponent: handle re-configuration event */
  int Reconfigure(const char* cdbEntry, const char* chainId);
  /** inherited from AliHLTComponent, scan one argument and its parameters */
  int ScanConfigurationArgument(int argc, const char** argv);

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
  
  void ReadTracks(const AliHLTComponentBlockData* iter,Int_t &tt);

  void PushHisto();
 
  Int_t fNEvents;    //! transient
  Int_t fNtotTracks; //! transient

  Int_t fEvtMod;     //! number of events reached to reset the counter
  Int_t fBufferSize; //! size of circular buffer (number of entries) for the ntuples
  Bool_t fdEdx;      //! plot dEdx
  //  Bool_t fReset;  //! Reset track counter every certain events

  TH1F *fMeanMultiplicity; //! transient (mean multiplicity for every 20 evts vs. #evt by Z.Y.)
  TH1F *fMultiplicity;     //! transient (track multiplicity by Z.Y.)
  //TH1F *fdNdEta;           //! transient (dN/dEta)
     
  //TH2F *fNClusterVsXY;   //! transient (#Clusters vs. x, y positions, by ZY)
  //TH2F *fChargeVsXY;     //! transient (Charge distr. vs. x, y added by ZY)
  TProfile *fDeDxVsP;    //! transient (dEdX vs. p)

  TNtuple *fClusters;                             //! transient  
  TNtuple *fTracks;                               //! transient

  AliHLTTPCSpacePointData *fClustersArray[36][6]; //! transient
  UInt_t                   fNSpacePoints[36][6];  //! transient
  
  /** the default configuration entry for this component */
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliHLTTPCTrackHistoComponent, 5);

};
#endif
