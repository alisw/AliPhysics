#ifndef ALIMUONQAINDICES_H
#define ALIMUONQAINDICES_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// 
/// Definition of enums needed by the MUON QA data makers and checkers (both
/// tracker and trigger)
///
/// \author Laurent Aphecetche and Diego Stocco, Subatech

namespace AliMUONQAIndices
{
  /// Histogram indices common to raw and digits and/or recpoints.
  ///
  /// WARNING ! Those ones should not be reused anywhere else below.
  /// Numbers from 1 to 49 are thus reserved for ECommon enum !
  ///
  enum ECommon {
    kTrackerBusPatchOccupancy         =  1, ///< Bus patch occupancies
    kTrackerBusPatchParityErrors      =  2, ///< Parity errors during readout of the tracker
    kTrackerBusPatchTokenLostErrors   =  3, ///< Token lost errors during readout of the tracker
    kTrackerBusPatchPaddingErrors     =  4, ///< Padding errors during readout of the tracker
    kTrackerNofPhysicsEventsSeen      =  5, ///< Number of events seen 
    kTrackerNofGoodPhysicsEventsUsed  =  6, ///< Number of good physics events seen (and used)
    kTrackerBusPatchConfig            =  7, ///< Configuration of the tracker
    kTrackerDDLOccupancy              =  8, ///< DDL occupancy in percent
    kTrackerDDLNofEventsUsed          =  9, ///< nof of events per DDL (used) *WARNING* : same as above
    kTrackerDDLNofEventsSeen          = 10,  ///< nof of events per DDL (seen)
    kTrackerData                      = 11, ///< Accumulated data
    kTrackerIsThere                   = 12, ///< whether we're making QA of tracker or not
    kTriggerIsThere                   = 13  ///< whether we're making QA of trigger or not
  };
  
  /// Raw/digits histograms indices
  enum ERaw { 
    
    kTrackerReadoutStatusPerEvent     = 51, ///< as kTrackerReadoutStatus but normalized by the number of events
    kTrackerReadoutStatus             = 52, ///< Status of readout (errors, missing pads, etc...)
    kTrackerDDLEventSize             =  53, ///< event size per DDL
    kTrackerDDLEventSizePerEvent     =  54, ///< event size per DDL per event

    kTriggerScalersTime       = 60, ///< Trigger scalers acquisition time index
    kTriggerScalers           = 61, ///< Trigger scalers histogram per plane index
    kTriggerScalersDisplay    = 71, ///< Trigger scalers display histogram per plane index
    kTriggerCalibSummary      = 80, ///< Number of responding strips/boards and noisy strips 
    kTriggerCalibSummaryNorm  = 81, ///< Percentage of responding strips/boards and noisy strips
    kTriggerErrorLocalXPos = 82, ///< Local board: Number of XPos Error vs Local Board Id
    kTriggerErrorLocalYPos = 83, ///< Local board: Number of YPos Error vs Local Board Id
    kTriggerErrorLocalDev  = 84, ///< Local board: Number of Deviation Error vs Local Board
    kTriggerErrorLocalTriggerDec = 85, ///< Local board: Number of Trigger Decision (All Pt) Error vs Local Board Id
    kTriggerErrorLocalLPtLSB = 86, ///< Local board: Number of LSB Low Pt Error vs Local Board Id
    kTriggerErrorLocalLPtMSB = 87, ///< Local board: Number of MSB Low Pt Error vs Local Board Id
    kTriggerErrorLocalHPtLSB = 88, ///< Local board: Number of LSB High Pt Error vs Local Board Id
    kTriggerErrorLocalHPtMSB = 89, ///< Local board: Number of MSB High Pt Error vs Local Board Id
    kTriggerErrorLocalTrigY  = 90, ///< Local board: Number of TrigY Error vs Local Board Id
    kTriggerErrorLocal2RegionalLPtLSB  = 91, ///< Local to Regional: Number of LPt LSB error vs Local Board Id
    kTriggerErrorLocal2RegionalLPtMSB  = 92, ///< Local to Regional: Number of LPt MSB error vs Local Board Id
    kTriggerErrorLocal2RegionalHPtLSB  = 93, ///< Local to Regional: Number of HPt LSB error vs Local Board Id
    kTriggerErrorLocal2RegionalHPtMSB  = 94, ///< Local to Regional: Number of HPt MSB error vs Local Board Id
    kTriggerErrorOutGlobalFromInGlobal = 95, ///< Global board: Number of error vs output bit    with a re-reconstruction from Global inputs
    kTriggerErrorOutGlobalFromInLocal = 96, ///< Global board: Number of error vs output bit  with a re-reconstruction from Local inputs
    kTriggerErrorSummary      = 97,  ///< Number of errors for each trigger decision level (Local, Reg->Local, Reg, Reg->Glob, Global)
    kTriggerErrorSummaryNorm  = 98,  ///< Percentage of errors for each trigger decision level
    kTriggerErrorLocalYCopy     = 99, ///< Local board: Number of Y Copy Error vs Local Board Id
    kTriggerErrorLocalYCopyTest = 100, ///< Local Board: Number of Y copy error tests (for normalization)
    kTriggerErrorLocalYCopyNorm = 101, ///< Local Board: Number of Y Copy Error vs Local Board Id Normalized to the number of tests
    kTriggeredBoards          = 102,  ///< Triggered boards histogram index
    kTriggerBoardsDisplay     = 103,  ///< Triggered boards display histogram index
    kTriggerReadOutErrors     = 104,  ///< Number of read-out errors
    kTriggerReadOutErrorsNorm = 105,  ///< Percentage of read-out errors
    kTriggerGlobalOutput      = 106,  ///< Number of Global outputs and Global algo errors
    kTriggerGlobalOutputNorm  = 107,  ///< Percentage of Global outputs and Global algo errors
    kTriggerRawNAnalyzedEvents= 108,  ///< Number of analyzed events per event specie
    kTriggerLocalRatio4434           = 109,  ///< Ratio 44/34 vs Local Board Id
    kTriggerRatio4434AllEvents       = 110,  ///< Ratio 44/34 since the beginning of the run vs Event Number
    kTriggerRatio4434SinceLastUpdate = 111,  ///< Ratio 44/34 for the last kUpdateRatio4434 events vs Event Number
    kTriggerNumberOf34Dec            = 112,  ///< Number of Decision in coincidence 3/4 vs Local Board
    kTriggerNumberOf44Dec            = 113,  ///< Number of Decision in coincidence 4/4 vs Local Board
    kTriggerGlobalScalers            = 114,  ///< Number of L0 counts
    kTriggerGlobalScalersNorm        = 115   ///< L0 scaler rates
    
  };
  
  /// Rec points histograms indices
  enum ERecPoints { 
    kTrackerNumberOfClustersPerChamber    = 100, ///< Tracker: number of clusters per chamber
    kTrackerClusterMultiplicityPerChMean  = 101, ///< cluster size per Ch: mean
    kTrackerClusterMultiplicityPerChSigma = 102, ///< cluster size per Ch: dispersion
    kTrackerClusterChargePerChMean        = 103, ///< cluster charge per Ch: mean
    kTrackerClusterChargePerChSigma       = 104, ///< cluster charge per Ch: dispersion
    
    kTrackerRecPoints = 105, ///< Tracker : tracker data of clusters (all and mono-cathode ones)
    
    kTrackerClusterMultiplicityPerChamber = 200, ///< Tracker: cluster multiplicity per chamber
    kTrackerClusterChargePerChamber       = 300, ///< Tracker: cluster charge per chamber
    kTrackerClusterHitMapPerChamber       = 400, ///< Tracker: cluster position distribution per chamber
    
    kTrackerNumberOfClustersPerDE        = 500, ///< Tracker : number of clusters per DE		
    kTrackerClusterMultiplicityPerDEMean = 501, ///< cluster size per DE: mean
    kTrackerClusterChargePerDEMean       = 502, ///< cluster charge per DE: mean
    
    kTrackerClusterMultiplicityPerDE = 3000, ///< Tracker : cluster multiplicity per DE		
    kTrackerClusterChargePerDE       = 5000,  ///< Tracker : cluster charge per DE
    
    kTriggerNAnalyzedEvents           = 50, ///< Number of analyzed events per event specie
    kTriggerRPCtrips           = 51, ///< Trips in trigger chambers
    kTriggerRPChv              = 52  ///< Trigger chamber HV index
    
  };
  
  /// ESD histograms indices
  enum EESD { 
    kESDnTracks                 =  50,  ///< number of tracks
    kESDMatchTrig               =  51,  ///< number of tracks matched with trigger
    kESDMomentum                =  52,  ///< P distribution
    kESDPt                      =  53,  ///< Pt distribution
    kESDRapidity                =  54,  ///< rapidity distribution
    kESDChi2                    =  55,  ///< normalized chi2 distribution
    kESDProbChi2                =  56,  ///< distribution of probability of chi2
    
    kESDClusterHitMap           =  57,  ///< cluster position distribution in chamber i
    kESDnClustersPerTrack       =  67, ///< number of clusters per track
    kESDnClustersPerCh          =  68, ///< number of clusters per chamber per track
    kESDnClustersPerDE          =  69, ///< number of clusters per DE per track
    kESDClusterChargeInCh       =  70, ///< cluster charge distribution in chamber i
    kESDClusterChargePerChMean  =  80, ///< cluster charge per Ch: mean
    kESDClusterChargePerChSigma =  81, ///< cluster charge per Ch: dispersion
    kESDClusterChargePerDE      =  82, ///< cluster charge per DE: mean
    kESDClusterSizeInCh         =  83, ///< cluster size distribution in chamber i
    kESDClusterSizePerChMean    =  93, ///< cluster size per Ch: mean
    kESDClusterSizePerChSigma   =  94, ///< cluster size per Ch: dispersion
    kESDClusterSizePerDE        =  95, ///< cluster size per DE: mean
    
    kESDResidualXInCh           =  96, ///< cluster-track residual-X distribution in chamber i
    kESDResidualYInCh           = 106, ///< cluster-track residual-Y distribution in chamber i
    kESDResidualXPerChMean      = 116, ///< cluster-track residual-X per Ch: mean
    kESDResidualYPerChMean      = 117, ///< cluster-track residual-Y per Ch: mean
    kESDResidualXPerChSigma     = 118, ///< cluster-track residual-X per Ch: dispersion
    kESDResidualYPerChSigma     = 119, ///< cluster-track residual-Y per Ch: dispersion
    kESDResidualXPerDEMean      = 120, ///< cluster-track residual-X per DE: mean
    kESDResidualYPerDEMean      = 121, ///< cluster-track residual-Y per DE: mean
    kESDResidualXPerDESigma     = 122, ///< cluster-track residual-X per DE: dispersion
    kESDResidualYPerDESigma     = 123, ///< cluster-track residual-Y per DE: dispersion
    kESDLocalChi2XInCh          = 124, ///< local chi2-X distribution in chamber i
    kESDLocalChi2YInCh          = 134, ///< local chi2-Y distribution in chamber i
    kESDLocalChi2XPerChMean     = 144, ///< local chi2-X per Ch: mean
    kESDLocalChi2YPerChMean     = 145, ///< local chi2-Y per Ch: mean
    kESDLocalChi2XPerDEMean     = 146, ///< local chi2-X per DE: mean
    kESDLocalChi2YPerDEMean     = 147, ///< local chi2-Y per DE: mean
    kESDLocalChi2InCh           = 148, ///< local chi2-X distribution in chamber i
    kESDLocalChi2PerChMean      = 158, ///< local chi2 per Ch: mean
    kESDLocalChi2PerDEMean      = 159, ///< local chi2 per DE: mean
    
    kESDThetaX                  = 160, ///< thetaX distribution
    kESDThetaY                  = 161, ///< thetaY distribution
    
    kESDnTotClustersPerCh       = 1000, ///< total number of associated clusters per chamber
    kESDnTotClustersPerDE       = 1001, ///< total number of associated clusters per DE
    kESDnTotFullClustersPerDE   = 1002, ///< total number of associated clusters containing pad info per DE
    kESDSumClusterChargePerDE   = 1003, ///< sum of cluster charge per DE
    kESDSumClusterSizePerDE     = 1004, ///< sum of cluster size per DE
    kESDSumResidualXPerDE       = 1005, ///< sum of cluster-track residual-X per DE
    kESDSumResidualYPerDE       = 1006, ///< sum of cluster-track residual-Y per DE
    kESDSumResidualX2PerDE      = 1007, ///< sum of cluster-track residual-X**2 per DE
    kESDSumResidualY2PerDE      = 1008, ///< sum of cluster-track residual-Y**2 per DE
    kESDSumLocalChi2XPerDE      = 1009, ///< sum of local chi2-X per DE
    kESDSumLocalChi2YPerDE      = 1010, ///< sum of local chi2-Y per DE
    kESDSumLocalChi2PerDE       = 1011  ///< sum of local chi2 per DE
  };
  
  // Bins for summary histos
  enum {
    kTriggerRespStrips,    ///< Bin for % of responding trigger strips
    kTriggerRespLocal,     ///< Bin for % of responding trigger local boards
    kTriggerRespRegional,  ///< Bin for % of responding trigger regional boards
    kTriggerRespGlobal,    ///< Bin for % of responding trigger global boards
    kTriggerNoisyStrips,   ///< Bin for % of noisy trigger strips
    kNtrigCalibSummaryBins ///< Total number of bins for trigger calibration summary
  };
  
  // Bins for algorithm error histos
  enum {
    kAlgoLocalX,             ///< Bin for % of local board X pos errors
    kAlgoLocalY,             ///< Bin for % of local board Y pos errors
    kAlgoLocalLUT,           ///< Bin for % of local board deviation errors
    kAlgoLocalYCopy,         ///< Bin for % of local board Y copy errors
    kAlgoLocalToRegional,    ///< Bin for % of local to regional errors
    kAlgoRegional,           ///< Bin for % of regional board errors 
    kAlgoRegionalToGlobal,   ///< Bin for % of regional to global errors 
    kAlgoGlobalFromGlobal,   ///< Bin for % of global from global board errors 
    kAlgoGlobalFromLocal,    ///< Bin for % of global from local board errors 
    kAlgoGlobalFromRegional, ///< Bin for % of global from regional board errors 
    kNtrigAlgoErrorBins      ///< Total number of bins for trigger error summary
  };
  
  enum {
    kLocalStructError,    ///< Bin for % of errors in local struct
    kRegionalStructError, ///< Bin for % of errors in regional struct
    kGlobalStructError,   ///< Bin for % of errors in global struct
    kDarcStructError,     ///< Bin for % of errors in darc struct
    kNtrigStructErrorBins ///< Total number of bins for struct error summary
  };
    
  // Bins for tracker readout errors
  enum ETrackerReadoutStatus
  {
    kTrackerRawNofGlitchErrors = 0, ///< Bin for number of glitch errors
    kTrackerRawNofTokenLostErrors = 1, ///< Bin for number of token lost errors
    kTrackerRawNofParityErrors = 2, ///< Bin for number of parity errors
    kTrackerRawNofPaddingErrors = 3, ///< Bin for number of padding errors
    kTrackerRawNofEmptyEvents = 4, ///< Bin for number of empty events
    kTrackerRawNofMissingBusPatchesFromConfig = 5 , ///< Bin for number of missing bus patches (in config)
    kTrackerRawNofMissingBusPatchesFromDataStream = 6 ///< Bin for number of missing bus patches (in actual data)
  };
  
}

#endif
