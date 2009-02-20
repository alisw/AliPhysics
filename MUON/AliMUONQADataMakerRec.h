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
class AliMUONVDigitStore;
class AliMUONVTriggerStore;
class AliMUONVClusterStore;
class AliMUONVTrackerData;
class AliMUONVTrackerDataMaker;

class AliMUONQADataMakerRec: public AliQADataMakerRec {

public:
  AliMUONQADataMakerRec();         
  AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm);   
  AliMUONQADataMakerRec& operator=(const AliMUONQADataMakerRec& qadm);
  virtual ~AliMUONQADataMakerRec();
  
  AliMUONVTrackerData* GetTrackerData() const;
  
protected:
	
  virtual void StartOfDetectorCycle(); 

  virtual void InitRaws(); 
  virtual void InitRecPoints(); 
  virtual void InitESDs(); 
  
  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeRecPoints(TTree* recpo); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  
  virtual void DefaultEndOfDetectorCycle(AliQA::TASKINDEX_t) {}

  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray** list);

private:
  /// Raw histograms indices
  enum ERaw { 
    kTrackerData           = 3,  ///< Accumulated data
    kTrackerBusPatchOccupancy = 4, ///< Bus patch occupancies
    kTriggerScalersBP      = 22, ///< Trigger scalers on BP histogram per chamber index
    kTriggerScalersNBP     = 23, ///< Trigger scalers on NBP histogram per chamber index
    kTriggerScalersDisplay = 24, ///< Trigger scalers display histogram per chamber index
    kTriggerScalersTime    = 32  ///< Trigger scalers acquisition time index
  };
         
  /// Rec points histograms indices
  enum ERecPoints { 
    kTriggerDigitsBendPlane    = 0,  ///< Trigger digits on BP histogram index
    kTriggerDigitsNonBendPlane = 1,  ///< Trigger digits on NBP histogram index
    kTriggeredBoards           = 2,  ///< Triggered boards histogram index
    kTriggerDigitsDisplay      = 3,  ///< Trigger digits display histogram per plane index
    kTriggerBoardsDisplay      = 11, ///< Triggered boards display histogram index
    kTriggerRPCi               = 12, ///< Trigger chamber currents index
    kTriggerRPChv              = 13, ///< Trigger chamber HV index
    kTriggerIDisplay           = 14, ///< Trigger chamber currents display histogram index
    kTriggerHVDisplay          = 18, ///< Trigger chamber HV display histogram index
    
    kTrackerNumberOfClustersPerChamber    = 100, ///< Tracker: # of clusters per chamber
    kTrackerClusterMultiplicityPerChamber = 200, ///< Tracker: cluster multiplicity per chamber
    kTrackerClusterChargePerChamber       = 300, ///< Tracker: cluster charge per chamber
				
    kTrackerNumberOfClustersPerDE    = 1000, ///< Tracker : number of clusters per DE		
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
    kESDSumLocalChi2YPerDE      = 1010  ///< sum of local chi2-X per DE
  };
  
private:
	
  void Ctor();
  void DisplayTriggerInfo(AliQA::TASKINDEX_t task);
  Bool_t FillTriggerDCSHistos();
  void InitRecPointsTracker();
  void InitRecPointsTrigger();
  void MakeRawsTracker(AliRawReader* rawReader);
  void MakeRawsTrigger(AliRawReader* rawReader);
  void MakeRecPointsTracker(TTree* treeR);
  void MakeRecPointsTrigger(TTree* treeR);
	
  /// Return reco parameters
  const AliMUONRecoParam* GetRecoParam() const { return dynamic_cast<const AliMUONRecoParam *>(fRecoParam); }
  
  Bool_t  fIsInitRaws;       //!<  info if InitRaws() went ok
  Bool_t  fIsInitRecPointsTracker;  //!<  info if InitRecPoints() went ok
  Bool_t  fIsInitRecPointsTrigger;  //!<  info if InitRecPoints() went ok
  Bool_t  fIsInitESDs;       //!<  info if InitESDs() went ok
  
  AliMUONVDigitStore*   fDigitStore; //!< pointer to digits store
  AliMUONVTriggerStore* fTriggerStore; //!< pointer to trigger store
  AliMUONDigitMaker*    fDigitMaker;  //!< pointer to digit maker
  AliMUONVClusterStore* fClusterStore; //!< pointer to cluster store
	
  AliMUONVTrackerDataMaker* fTrackerDataMaker; //!< tracker data accumulation
  
  ClassDef(AliMUONQADataMakerRec,5)  // MUON Quality assurance data maker

};
#endif
