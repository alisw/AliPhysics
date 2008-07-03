#ifndef AliMUONQADataMakerRec_H
#define AliMUONQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONQADataMakerRec
/// \brief MUON Quality assurance data maker
///

// --- AliRoot header files ---
class AliMUONDigitMaker;
class AliMUONVDigitStore;
class AliMUONVTriggerStore;
class AliMUONVClusterStore;

#include "AliQADataMakerRec.h"

class AliMUONVTrackerDataMaker;

class AliMUONQADataMakerRec: public AliQADataMakerRec {

public:
  AliMUONQADataMakerRec();         
  AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm);   
  AliMUONQADataMakerRec& operator=(const AliMUONQADataMakerRec& qadm);
  virtual ~AliMUONQADataMakerRec();
  
private:
  /// Raw histograms indices
  enum ERaw { 
		kTrackerData           = 3,   ///< Accumulated data
    kTriggerScalersBP      = 22,  ///< Trigger scalers on BP histogram per chamber index
    kTriggerScalersNBP     = 23,  ///< Trigger scalers on NBP histogram per chamber index
    kTriggerScalersDisplay = 24   ///< Trigger scalers display histogram per chamber index
  };
         
  /// Rec points histograms indices
  enum ERecPoints { 
    
		kTriggerDigitsBendPlane    = 0,  ///< Trigger digits on BP histogram index
    kTriggerDigitsNonBendPlane = 1,  ///< Trigger digits on NBP histogram index
    kTriggeredBoards           = 2,  ///< Triggered boards histogram index
    kTriggerDigitsDisplay      = 3,  ///< Trigger digits display histogram per plane index
    kTriggerBoardsDisplay      = 11,  ///< Triggered boards display histogram index

		kTrackerNumberOfClustersPerChamber = 100, ///< Tracker: # of clusters per chamber
		kTrackerClusterMultiplicityPerChamber = 200, ///< Tracker: cluster multiplicity per chamber
		kTrackerClusterChargePerChamber = 300, ///< Tracker: cluster charge per chamber
				
		kTrackerNumberOfClustersPerDE    = 1000, ///< Tracker : number of clusters per DE		
		kTrackerClusterMultiplicityPerDE = 3000, ///< Tracker : cluster multiplicity per DE		
		kTrackerClusterChargePerDE       = 5000 ///< Tracker : cluster charge per DE
		
  };
          
  /// ESD histograms indices
  enum EESD { 
    kESDnTracks       = 0,  ///< ESD nTrack histogram index
    kESDMomentum      = 1,  ///< ESD momentum histogram index
    kESDPt            = 2,  ///< ESD Pt histogram index
    kESDRapidity      = 3,  ///< ESD Rapidity histogram index
    kESDClusterHitMap = 4   ///< ESD Cluster hit map histogram index
  };

protected:
	
  virtual void StartOfDetectorCycle(); 

  virtual void InitRaws(); 
  virtual void InitRecPoints(); 
  virtual void InitESDs(); 
  
  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeRecPoints(TTree* recpo); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  
  virtual void DefaultEndOfDetectorCycle(AliQA::TASKINDEX_t) {}

	virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray* list);

private:
	
	void Ctor();
  void DisplayTriggerInfo(AliQA::TASKINDEX_t task);
	void InitRecPointsTracker();
	void InitRecPointsTrigger();
	void MakeRawsTracker(AliRawReader* rawReader);
	void MakeRawsTrigger(AliRawReader* rawReader);
  void MakeRecPointsTracker(TTree* treeR);
  void MakeRecPointsTrigger(TTree* treeR);
	
  Bool_t  fIsInitRaws;       //!<  info if InitRaws() went ok
  Bool_t  fIsInitRecPointsTracker;  //!<  info if InitRecPoints() went ok
  Bool_t  fIsInitRecPointsTrigger;  //!<  info if InitRecPoints() went ok
  Bool_t  fIsInitESDs;       //!<  info if InitESDs() went ok
  
  AliMUONVDigitStore*   fDigitStore; //!< pointer to digits store
  AliMUONVTriggerStore* fTriggerStore; //!< pointer to trigger store
  AliMUONDigitMaker*    fDigitMaker;  //!< pointer to digit maker
  AliMUONVClusterStore* fClusterStore; //!< pointer to cluster store
	
	AliMUONVTrackerDataMaker* fTrackerDataMaker; //!< tracker data accumulation
	
  ClassDef(AliMUONQADataMakerRec,3)  // MUON Quality assurance data maker

};
#endif
