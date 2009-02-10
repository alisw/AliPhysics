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
class AliMUONVTrackerDataMaker;

class AliMUONQADataMakerRec: public AliQADataMakerRec {

public:
  AliMUONQADataMakerRec();         
  AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm);   
  AliMUONQADataMakerRec& operator=(const AliMUONQADataMakerRec& qadm);
  virtual ~AliMUONQADataMakerRec();
  
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
    kTriggerScalersDisplay = 24  ///< Trigger scalers display histogram per chamber index
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
    kESDnTracks             = 0,  ///< number of tracks
    kESDMatchTrig           = 1,  ///< number of tracks matched with trigger
    kESDMomentum            = 2,  ///< P distribution
    kESDPt                  = 3,  ///< Pt distribution
    kESDRapidity            = 4,  ///< rapidity distribution
    kESDChi2                = 5,  ///< normalized chi2 distribution
    
    kESDClusterHitMap       = 6,  ///< cluster position distribution in chamber i
    kESDnClustersPerTrack   = 16, ///< number of clusters per track
    kESDnClustersPerCh      = 17, ///< number of clusters per chamber
    kESDnClustersPerDE      = 18, ///< number of clusters per DE
    kESDClusterCharge       = 19, ///< cluster charge distribution
    kESDClusterChargeInCh   = 20, ///< cluster charge distribution in chamber i
    kESDClusterChargePerDE  = 30, ///< cluster charge per DE: mean +- dispersion
    kESDClusterMult         = 31, ///< cluster multiplicity distribution
    kESDClusterMultInCh     = 32, ///< cluster multiplicity distribution in chamber i
    kESDClusterMultPerDE    = 42, ///< cluster multiplicity per DE: mean +- dispersion
    
    kESDResidualX           = 43, ///< cluster-track residual-X distribution
    kESDResidualY           = 44, ///< cluster-track residual-Y distribution
    kESDResidualXInCh       = 45, ///< cluster-track residual-X distribution in chamber i
    kESDResidualYInCh       = 55, ///< cluster-track residual-Y distribution in chamber i
    kESDResidualXPerDEMean  = 65, ///< cluster-track residual-X per DE: mean
    kESDResidualYPerDEMean  = 66, ///< cluster-track residual-Y per DE: mean
    kESDResidualXPerDESigma = 67, ///< cluster-track residual-X per DE: sigma
    kESDResidualYPerDESigma = 68  ///< cluster-track residual-Y per DE: sigma
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
	
  ClassDef(AliMUONQADataMakerRec,4)  // MUON Quality assurance data maker

};
#endif
