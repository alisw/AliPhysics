#ifndef ALIMUONQADATAMAKERREC_H
#define ALIMUONQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONQADataMakerRec
/// \brief MUON Quality assurance data maker
///

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"
#include "AliMUONRecoParam.h"

class AliMUONDigitMaker;
class AliMUONVClusterStore;
class AliMUONVDigitStore;
class AliMUONVStore;
class AliMUONVTrackerData;
class AliMUONVTrackerDataMaker;
class AliMUONCalibrationData;
class AliMUONTriggerElectronics;

class AliMUONQAMappingCheck;

class AliMUONQADataMakerRec: public AliQADataMakerRec {

public:
  AliMUONQADataMakerRec();         
  AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm);   
  AliMUONQADataMakerRec& operator=(const AliMUONQADataMakerRec& qadm);
  virtual ~AliMUONQADataMakerRec();
  
  AliMUONVTrackerData* GetTrackerData() const;

  virtual void InitRaws(); 
  virtual void InitRecPoints(); 

protected:
	
  virtual void StartOfDetectorCycle(); 

  virtual void InitDigits(); 
  virtual void InitESDs(); 
  
  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeDigits(TTree* dig); 
  virtual void MakeDigits()  {return;}
  virtual void MakeRecPoints(TTree* recpo); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  
  virtual void DefaultEndOfDetectorCycle(AliQAv1::TASKINDEX_t) {}

  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** list);

private:
  /// Raw histograms indices
  enum ERaw { 
    kTrackerData              = 3,  ///< Accumulated data
    kTrackerBusPatchOccupancy = 4, ///< Bus patch occupancies
    kTrackerBusPatchNofPads   = 5, ///< Number of pads per bus patch
    kTrackerBusPatchNofManus  = 6, ///< Number of manus per bus patch
    kTriggerScalers           = 22, ///< Trigger scalers histogram per plane index
    kTriggerScalersDisplay    = 30, ///< Trigger scalers display histogram per plane index
    kTriggerScalersTime       = 38, ///< Trigger scalers acquisition time index
    kTriggerRPCi              = 39, ///< Trigger chamber currents index
    kTriggerRPChv             = 43,  ///< Trigger chamber HV index
    kTriggerErrorLocalXPos = 50,  ///< Local board: Number of XPos Error vs Local Board Id
    kTriggerErrorLocalYPos = 51,  ///< Local board: Number of YPos Error vs Local Board Id
    kTriggerErrorLocalDev = 52,  ///< Local board: Number of Deviation Error vs Local Board
    kTriggerErrorLocalTriggerDec = 53,  ///< Local board: Number of Trigger Decision (All Pt) Error vs Local Board Id
    kTriggerErrorLocalLPtLSB = 54,  ///< Local board: Number of LSB Low Pt Error vs Local Board Id
    kTriggerErrorLocalLPtMSB = 55,  ///< Local board: Number of MSB Low Pt Error vs Local Board Id
    kTriggerErrorLocalHPtLSB = 56,  ///< Local board: Number of LSB High Pt Error vs Local Board Id
    kTriggerErrorLocalHPtMSB = 57,  ///< Local board: Number of MSB High Pt Error vs Local Board Id
    kTriggerErrorLocal2RegionalLPtLSB = 58,  ///< Local to Regional: Number of LPt LSB error vs Local Board Id
    kTriggerErrorLocal2RegionalLPtMSB = 59,  ///< Local to Regional: Number of LPt MSB error vs Local Board Id
    kTriggerErrorLocal2RegionalHPtLSB = 60,  ///< Local to Regional: Number of HPt LSB error vs Local Board Id
    kTriggerErrorLocal2RegionalHPtMSB = 61,  ///< Local to Regional: Number of HPt MSB error vs Local Board Id
    kTriggerErrorOutGlobalFromInGlobal = 62,  ///< Global board: Number of error vs output bit 
    kTriggerError = 63,  ///< percentage of error for each trigger decision level (Local, Reg->Local, Reg, Reg->Glob, Global)
    kTriggerErrorLocalTrigY = 64,  ///< Local board: Number of TrigY Error vs Local Board Id
    kTriggerErrorLocalYCopy = 65,  ///< Local board: Number of Y Copy Error vs Local Board Id

    kRawNAnalyzedEvents = 66 ///< Number of analyzed events per event specie
  };
         
  /// Rec points histograms indices
  enum ERecPoints { 
    kTriggerDigits             = 0,  ///< Trigger digits histogram per plane index
    kTriggeredBoards           = 8,  ///< Triggered boards histogram index
    kTriggerDigitsDisplay      = 9,  ///< Trigger digits display histogram per plane index
    kTriggerBoardsDisplay      = 17, ///< Triggered boards display histogram index

    kNAnalyzedEvents           = 90, ///< Number of analyzed events per event specie
    
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
    kTrackerClusterChargePerDE       = 5000  ///< Tracker : cluster charge per DE
    
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
  
private:
  void BeautifyTrackerBusPatchOccupancy(TH1& hbp);
  
  void EndOfDetectorCycleRaws(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleESDs();
 
  void InsertTrackerData(Int_t specie, TObjArray** list, TObject* object, 
                         Int_t indexNumber, Bool_t replace=kFALSE);

private:
	
  void Ctor();
  void DisplayTriggerInfo(AliQAv1::TASKINDEX_t task);
  Bool_t FillTriggerDCSHistos();
  void InitRecPointsTracker();
  void InitRecPointsTrigger();
  void MakeRawsTracker(AliRawReader* rawReader);
  void MakeRawsTrigger(AliRawReader* rawReader);
  void MakeRecPointsTracker(TTree* treeR);
  void MakeRecPointsTrigger(TTree* treeR);
  void RawTriggerInRegional2OutRegional();
  void RawTriggerInGlobal2OutGlobal();
  void RawTriggerMatchOutLocal();
  void RawTriggerMatchOutLocalInRegional();
  void RawTriggerMatchOutGlobalFromInGlobal();
	
  /// Return reco parameters
  const AliMUONRecoParam* GetMUONRecoParam() const { return dynamic_cast<const AliMUONRecoParam *>(fRecoParam); }
  
  Int_t fTriggerOutputLocalDataTriggerDec[235]; ///< Data Local Trigger decision for each active Local Board
  Int_t fTriggerOutputLocalDataTrigY[235]; ///< Data Local Trigger Y decision for each active Local Board
  Int_t fTriggerOutputLocalDataLPtDec[2][235]; ///< Data Local decision Low Pt for each active Local Board (2 Bits -> 0:LSB, 1:MSB)
  Int_t fTriggerOutputLocalDataHPtDec[2][235]; ///< Data Local decision High Pt for each active Local Board (2 Bits -> 0:LSB, 1:MSB)
  Int_t fTriggerOutputLocalDataXPos[235]; ///< Data Local XPos for each active Local Board
  Int_t fTriggerOutputLocalDataYPos[235]; ///< Data Local YPos for each active Local Board
  Int_t fTriggerOutputLocalDataDev[235]; ///< Data Local deviation for each active Local Board

  Int_t fTriggerOutputLocalRecTriggerDec[235]; ///< Reconstructed Local Trigger decision for each active Local Board
  Int_t fTriggerOutputLocalRecTrigY[235]; ///< Reconstructed Local Trigger Y decision for each active Local Board
  Int_t fTriggerOutputLocalRecLPtDec[2][235]; ///< Reconstructed Local decision Low Pt for each active Local Board (2 Bits -> 0:LSB, 1:MSB)
  Int_t fTriggerOutputLocalRecHPtDec[2][235]; ///< Reconstructed Local decision High Pt for each active Local Board (2 Bits -> 0:LSB, 1:MSB)
  Int_t fTriggerOutputLocalRecXPos[235]; ///< Reconstructed Local XPos for each active Local Board
  Int_t fTriggerOutputLocalRecYPos[235]; ///< Reconstructed Local YPos for each active Local Board
  Int_t fTriggerOutputLocalRecDev[235]; ///< Reconstructed Local deviation for each active Local Board

  Int_t fTriggerInputRegionalDataLPt[2][235]; ///< Data Regional Input LPt for each Local board
  Int_t fTriggerInputRegionalDataHPt[2][235]; ///< Data Regional Input HPt for each Local board
  Int_t fTriggerOutputRegionalData[16]; ///< Data Regional Trigger decision for each Regional Board (1R:0, 2R:1, ... , 1L:8, ...) -> 4 bits LPt, 4 bits HPt
  Int_t fTriggerInputRegionalRecLPt[2][16][16]; ///< Reconstructed Regional Input LPt for each Regional Board ([bit][reg][loc]) (reg -> 1R:0, 2R:1, ... , 1L:8, ...)
  Int_t fTriggerInputRegionalRecHPt[2][16][16]; ///< Reconstructed Regional Input HPt for each Regional Board ([bit][reg][loc]) (reg -> 1R:0, 2R:1, ... , 1L:8, ...)
  Int_t fTriggerOutputRegionalRec[16]; ///< Reconstructed Regional Trigger decision for each Regional Board (8 Bits)

  Int_t fTriggerInputGlobalDataLPt[16][4]; ///< Data Global inputs LPt (1R:0, 2R:1, ... , 1L:8, ...)
  Int_t fTriggerInputGlobalDataHPt[16][4]; ///< Data Global inputs HPt (1R:0, 2R:1, ... , 1L:8, ...)
  Int_t fTriggerOutputGlobalData[6]; ///< Data Global outputs
  Int_t fTriggerOutputGlobalRecFromGlobalInput[6]; //< Reconstructed Global outputs from Global inputs
  Int_t fTriggerOutputGlobalRecFromLocalInput[6]; //< Reconstructed Global outputs from Local inputs
  Int_t fTriggerOutputGlobalRecFromLocalOutput[6]; //< Reconstructed Global outputs from Local outputs
  Int_t fgitmp[4]; //< Tempory used to store Global inputs
  Int_t fgotmp[6]; //< Tempory used to store Global outputs

  Int_t fTriggerPatternX1[243][16]; ///< Local pattern X1
  Int_t fTriggerPatternX2[243][16]; ///< Local pattern X2
  Int_t fTriggerPatternX3[243][16]; ///< Local pattern X3
  Int_t fTriggerPatternX4[243][16]; ///< Local pattern X4
  Int_t fTriggerPatternY1[243][16]; ///< Local pattern Y1
  Int_t fTriggerPatternY2[243][16]; ///< Local pattern Y2
  Int_t fTriggerPatternY3[243][16]; ///< Local pattern Y3
  Int_t fTriggerPatternY4[243][16]; ///< Local pattern Y4

  Bool_t fTriggerErrorLocalYCopy[235]; ///< True if Y copy error for Local Board i
  
  AliMUONVDigitStore*   fDigitStore; //!< pointer to digits store
  AliMUONDigitMaker*    fDigitMaker;  //!< pointer to digit maker
  AliMUONVClusterStore* fClusterStore; //!< pointer to cluster store
	
  AliMUONVTrackerDataMaker* fTrackerDataMaker; //!< tracker data accumulation (Raw)
  
  AliMUONQAMappingCheck* fMappingCheckRecPoints; //!< mapping cross-checker (RecPoints)
  
  AliMUONCalibrationData* fCalibrationData; //!< Used to load Local, Regional and Global masks

  AliMUONTriggerElectronics* fTriggerProcessor; //!< trigger processore to re-compute response
  
  ClassDef(AliMUONQADataMakerRec,9)  // MUON Quality assurance data maker

};
#endif
