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
  /// Raw histograms indices
  enum ERaw { 
    kTrackerData              = 3,  ///< Accumulated data
    kTrackerBusPatchOccupancy = 4, ///< Bus patch occupancies
    kTrackerBusPatchNofPads   = 5, ///< Number of pads per bus patch
    kTrackerBusPatchNofManus  = 6, ///< Number of manus per bus patch
    kTrackerBusPatchConfig    = 7, ///< Configuration of the tracker
    kTrackerBusPatchParityErrors    =  8, ///< Parity errors during readout of the tracker
    kTrackerBusPatchTokenLostErrors =  9, ///< Token lost errors during readout of the tracker
    kTrackerBusPatchPaddingErrors   = 10, ///< Padding errors during readout of the tracker
    kTrackerNofRawEventSeen         = 11, ///< Number of events seen (and used)
    kTrackerReadoutErrors           = 12,  ///< Integrated number of errors (and events for 1st bin)
    kTrackerDDLOccupancy            = 13, ///< DDL occupancy in percent
    kTrackerDDLNofEvents            = 14, ///< nof of events per DDL
    kTriggerScalersTime       = 22, ///< Trigger scalers acquisition time index
    kTriggerScalers           = 23, ///< Trigger scalers histogram per plane index
    kTriggerScalersDisplay    = 31, ///< Trigger scalers display histogram per plane index
    kTriggerCalibSummary      = 40, ///< Number of responding strips/boards and noisy strips 
    kTriggerCalibSummaryNorm  = 41, ///< Percentage of responding strips/boards and noisy strips
    kTriggerErrorLocalXPos = 50, ///< Local board: Number of XPos Error vs Local Board Id
    kTriggerErrorLocalYPos = 51, ///< Local board: Number of YPos Error vs Local Board Id
    kTriggerErrorLocalDev  = 52, ///< Local board: Number of Deviation Error vs Local Board
    kTriggerErrorLocalTriggerDec = 53, ///< Local board: Number of Trigger Decision (All Pt) Error vs Local Board Id
    kTriggerErrorLocalLPtLSB = 54, ///< Local board: Number of LSB Low Pt Error vs Local Board Id
    kTriggerErrorLocalLPtMSB = 55, ///< Local board: Number of MSB Low Pt Error vs Local Board Id
    kTriggerErrorLocalHPtLSB = 56, ///< Local board: Number of LSB High Pt Error vs Local Board Id
    kTriggerErrorLocalHPtMSB = 57, ///< Local board: Number of MSB High Pt Error vs Local Board Id
    kTriggerErrorLocalTrigY  = 58, ///< Local board: Number of TrigY Error vs Local Board Id
    kTriggerErrorLocal2RegionalLPtLSB  = 59, ///< Local to Regional: Number of LPt LSB error vs Local Board Id
    kTriggerErrorLocal2RegionalLPtMSB  = 60, ///< Local to Regional: Number of LPt MSB error vs Local Board Id
    kTriggerErrorLocal2RegionalHPtLSB  = 61, ///< Local to Regional: Number of HPt LSB error vs Local Board Id
    kTriggerErrorLocal2RegionalHPtMSB  = 62, ///< Local to Regional: Number of HPt MSB error vs Local Board Id
    kTriggerErrorOutGlobalFromInGlobal = 63, ///< Global board: Number of error vs output bit    with a re-reconstruction from Global inputs
    kTriggerErrorOutGlobalFromInLocal = 64, ///< Global board: Number of error vs output bit  with a re-reconstruction from Local inputs
    kTriggerErrorSummary      = 65,  ///< Number of errors for each trigger decision level (Local, Reg->Local, Reg, Reg->Glob, Global)
    kTriggerErrorSummaryNorm  = 66,  ///< Percentage of errors for each trigger decision level
    kTriggerErrorLocalYCopy     = 67, ///< Local board: Number of Y Copy Error vs Local Board Id
    kTriggerErrorLocalYCopyTest = 68, ///< Local Board: Number of Y copy error tests (for normalization)
    kTriggerErrorLocalYCopyNorm = 69, ///< Local Board: Number of Y Copy Error vs Local Board Id Normalized to the number of tests
    kTriggeredBoards          = 70,  ///< Triggered boards histogram index
    kTriggerBoardsDisplay     = 71,  ///< Triggered boards display histogram index
    kTriggerReadOutErrors     = 80,  ///< Number of read-out errors
    kTriggerReadOutErrorsNorm = 81,  ///< Percentage of read-out errors
    kTriggerGlobalOutput      = 90,  ///< Number of Global outputs and Global algo errors
    kTriggerGlobalOutputNorm  = 91,  ///< Percentage of Global outputs and Global algo errors
    kTriggerRawNAnalyzedEvents= 100,  ///< Number of analyzed events per event specie
    kTriggerLocalRatio4434           = 101,  ///< Ratio 44/34 vs Local Board Id
    kTriggerRatio4434AllEvents       = 102,  ///< Ratio 44/34 since the beginning of the run vs Event Number
    kTriggerRatio4434SinceLastUpdate = 103,  ///< Ratio 44/34 for the last kUpdateRatio4434 events vs Event Number
    kTriggerNumberOf34Dec            = 104,  ///< Number of Decision in coincidence 3/4 vs Local Board
    kTriggerNumberOf44Dec            = 105   ///< Number of Decision in coincidence 4/4 vs Local Board
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
    
    kTrackerNumberOfClustersPerDE        = 1000, ///< Tracker : number of clusters per DE		
    kTrackerClusterMultiplicityPerDEMean = 1001, ///< cluster size per DE: mean
    kTrackerClusterChargePerDEMean       = 1002, ///< cluster charge per DE: mean
    
    kTrackerClusterMultiplicityPerDE = 3000, ///< Tracker : cluster multiplicity per DE		
    kTrackerClusterChargePerDE       = 5000,  ///< Tracker : cluster charge per DE
    
    kTriggerNAnalyzedEvents           = 0, ///< Number of analyzed events per event specie
    kTriggerRPCtrips           = 1, ///< Trips in trigger chambers
    kTriggerRPChv              = 2  ///< Trigger chamber HV index
    
  };
  
  /// ESD histograms indices
  enum EESD { 
    kESDnTracks                 = 0,  ///< number of tracks
    kESDMatchTrig               = 1,  ///< number of tracks matched with trigger
    kESDMomentum                = 2,  ///< P distribution
    kESDPt                      = 3,  ///< Pt distribution
    kESDRapidity                = 4,  ///< rapidity distribution
    kESDChi2                    = 5,  ///< normalized chi2 distribution
    kESDProbChi2                = 6,  ///< distribution of probability of chi2
    
    kESDClusterHitMap           = 7,  ///< cluster position distribution in chamber i
    kESDnClustersPerTrack       = 17, ///< number of clusters per track
    kESDnClustersPerCh          = 18, ///< number of clusters per chamber per track
    kESDnClustersPerDE          = 19, ///< number of clusters per DE per track
    kESDClusterChargeInCh       = 20, ///< cluster charge distribution in chamber i
    kESDClusterChargePerChMean  = 30, ///< cluster charge per Ch: mean
    kESDClusterChargePerChSigma = 31, ///< cluster charge per Ch: dispersion
    kESDClusterChargePerDE      = 32, ///< cluster charge per DE: mean
    kESDClusterSizeInCh         = 33, ///< cluster size distribution in chamber i
    kESDClusterSizePerChMean    = 43, ///< cluster size per Ch: mean
    kESDClusterSizePerChSigma   = 44, ///< cluster size per Ch: dispersion
    kESDClusterSizePerDE        = 45, ///< cluster size per DE: mean
    
    kESDResidualXInCh           = 46, ///< cluster-track residual-X distribution in chamber i
    kESDResidualYInCh           = 56, ///< cluster-track residual-Y distribution in chamber i
    kESDResidualXPerChMean      = 66, ///< cluster-track residual-X per Ch: mean
    kESDResidualYPerChMean      = 67, ///< cluster-track residual-Y per Ch: mean
    kESDResidualXPerChSigma     = 68, ///< cluster-track residual-X per Ch: dispersion
    kESDResidualYPerChSigma     = 69, ///< cluster-track residual-Y per Ch: dispersion
    kESDResidualXPerDEMean      = 70, ///< cluster-track residual-X per DE: mean
    kESDResidualYPerDEMean      = 71, ///< cluster-track residual-Y per DE: mean
    kESDResidualXPerDESigma     = 72, ///< cluster-track residual-X per DE: dispersion
    kESDResidualYPerDESigma     = 73, ///< cluster-track residual-Y per DE: dispersion
    kESDLocalChi2XInCh          = 74, ///< local chi2-X distribution in chamber i
    kESDLocalChi2YInCh          = 84, ///< local chi2-Y distribution in chamber i
    kESDLocalChi2XPerChMean     = 94, ///< local chi2-X per Ch: mean
    kESDLocalChi2YPerChMean     = 95, ///< local chi2-Y per Ch: mean
    kESDLocalChi2XPerDEMean     = 96, ///< local chi2-X per DE: mean
    kESDLocalChi2YPerDEMean     = 97, ///< local chi2-Y per DE: mean
    kESDLocalChi2InCh           = 98, ///< local chi2-X distribution in chamber i
    kESDLocalChi2PerChMean      = 108, ///< local chi2 per Ch: mean
    kESDLocalChi2PerDEMean      = 109, ///< local chi2 per DE: mean
    
    kESDThetaX                  = 110, ///< thetaX distribution
    kESDThetaY                  = 111, ///< thetaY distribution
    
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
  enum ETrackerReadoutErrors
  {
    kTrackerRawNofGlitchErrors = 0, ///< Bin for number of glitch errors
    kTrackerRawNofTokenLostErrors = 1, ///< Bin for number of token lost errors
    kTrackerRawNofParityErrors = 2, ///< Bin for number of parity errors
    kTrackerRawNofPaddingErrors = 3 ///< Bin for number of padding errors
  };
  
}

#endif
