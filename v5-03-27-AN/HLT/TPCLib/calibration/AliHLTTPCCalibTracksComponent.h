// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCALIBTRACKSCOMPONENT_H
#define ALIHLTTPCALIBTRACKSCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCCalibTracksComponent.h
    @author Kalliopi Kanaki
    @date   
    @brief  A tracks gain calibration component for the TPC.
*/

#include "AliHLTCalibrationProcessor.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "TObjArray.h"

class  AliTPCClusterParam;
class  AliTPCcalibTracksCuts;

class AliTPCcalibAlign;
class AliTPCcalibTracksGain;
class AliTPCcalibTracks;

class AliTPCseed;
class AliTPCclusterMI;
class AliHLTTPCOfflineCluster;

/**
 * @class AliHLTTPCCalibTracksComponent
 * 
 * This class is the calibration component for the AliTPCCalibPedestal class 
 * used for pedestal calibration of the TPC. 
 * 
 * It inherits from the AliHLTCalibrationProcessor and uses the high-level 
 * interface. The output is the class AliTPCCalibPedestal as a TObject.
 *
 * The component has the following component arguments:
 *   -enableanalysis : enable analysis before shipping data to FXS
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCCalibTracksComponent : public AliHLTCalibrationProcessor{

  public:
    /** constructor */
    AliHLTTPCCalibTracksComponent();
    /** destructor */
    virtual ~AliHLTTPCCalibTracksComponent();
    
    // Public functions to implement AliHLTComponent's interface.
    // These functions are required for the registration process

    const char* GetComponentID();
    void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
    AliHLTComponentDataType GetOutputDataType();  
    Int_t GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
    virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
    AliHLTComponent* Spawn();

  protected:

    using AliHLTCalibrationProcessor::ProcessCalibration;
    using AliHLTCalibrationProcessor::ShipDataToFXS;
    
    //int Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/);

    // Protected functions to implement AliHLTComponent's interface.
    // These functions provide initialization as well as the actual processing
    // capabilities of the component. 
    
    /** Initialize the calibration component. */
    Int_t InitCalibration();

    /** Scan commandline arguments of the calibration component. */
    Int_t ScanArgument( Int_t argc, const char** argv );

    /** DeInitialize the calibration component. */
    Int_t DeinitCalibration();

    /** Process the data in the calibration component. */
    Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

    /** Ship the data to the FXS at end of run or eventmodulo. */
    Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );
    
    Int_t Reconfigure(const char* cdbEntry, const char* chainId);

  private:
    /** copy constructor prohibited */
    AliHLTTPCCalibTracksComponent(const AliHLTTPCCalibTracksComponent&);
    /** assignment operator prohibited */
    AliHLTTPCCalibTracksComponent& operator=(const AliHLTTPCCalibTracksComponent&);
  
    Int_t Configure(const char* arguments);
    void ReadTracks(const AliHLTComponentBlockData* iter, Int_t &tt);
    
    AliTPCClusterParam    *fClustParam; //! TPC cluster parameters
    AliTPCcalibTracksCuts *fTrackCuts;  //! TPC track cuts 

    AliTPCcalibTracksGain *fCalibTracksGain; //!transient 
    AliTPCcalibAlign      *fCalibAlign;	     //!transient
    AliTPCcalibTracks     *fCalibTracks;     //!transient

    /** Minimum patch specifcation for this component */
    AliHLTUInt8_t fMinPatch;	 // see above

    /** Minimum patch specifcation for this component */
    AliHLTUInt8_t fMaxPatch;	 // see above

    /** The Specification for this component */
    AliHLTUInt32_t fSpecification; // see above

    /** Analyze calibration data before shipping to FXS */
    Bool_t fEnableAnalysis; //! transient
    
    /** seed argument for the Process(AliTPCseed*) function */
    AliTPCseed  *fSeed;							 // see above
    
    vector<AliHLTTPCSpacePointData> fClusters;                       //! transient
    vector<AliHLTTPCTrackSegmentData> fTracks;                       //! transient
    vector<UInt_t> fTrackClusterID[36][6];                           //! transient
    
    
    AliHLTTPCOfflineCluster *pConv;  //! transient
    TObjArray fOffArray;             //! transient
    
    Bool_t fReadMergedTracks; //! transient
    Bool_t fReadSliceTracks;  //! transient
   
    ClassDef(AliHLTTPCCalibTracksComponent, 0)

};
#endif
