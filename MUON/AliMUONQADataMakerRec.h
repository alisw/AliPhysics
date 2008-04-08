#ifndef AliMUONQADataMakerRec_H
#define AliMUONQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup rec
/// \class AliMUONQADataMakerRec
/// \brief MUON Quality assurance data maker
///
//  Author Christian Finck

// dummy function for simulation part
// to avoid circular dependencie

// --- ROOT system ---
class TObjArray;
class TArrayF;

// --- AliRoot header files ---
class AliMUONVTrackStore;
class AliMUONDigitMaker;
class AliMUONVDigitStore;
class AliMUONVTriggerStore;

#include "AliQADataMakerRec.h"

class AliMUONQADataMakerRec: public AliQADataMakerRec {

public:
  AliMUONQADataMakerRec();         
  AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm);   
  AliMUONQADataMakerRec& operator=(const AliMUONQADataMakerRec& qadm);
  virtual ~AliMUONQADataMakerRec();
  
private:
  /// Raw histograms indices
  enum ERaw { 
    kRawBusPatch           = 0,   ///< Raw bus patch histogram index
    kRawCharge             = 1,   ///< Raw charge histogram index
    kRawBuspatchDDL        = 2,   ///< Raw buspatch hit map histogram per DDL index
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
    kTriggerBoardsDisplay      = 11  ///< Triggered boards display histogram index
  };
          
  /// ESD histograms indices
  enum EESD { 
    kESDnTracks       = 0,  ///< ESD nTrack histogram index
    kESDMomentum      = 1,  ///< ESD momentum histogram index
    kESDPt            = 2,  ///< ESD Pt histogram index
    kESDRapidity      = 3,  ///< ESD Rapidity histogram index
    kESDClusterHitMap = 4   ///< ESD Cluster hit map histogram index
  };

  virtual void   StartOfDetectorCycle(); 

  virtual void   InitRaws(); 
  virtual void   InitRecPoints(); 
  virtual void   InitESDs(); 
  
  virtual void   MakeRaws(AliRawReader* rawReader); 
  virtual void   MakeRecPoints(TTree* recpo); 
  virtual void   MakeESDs(AliESDEvent* esd) ;
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray* list);

  void DisplayTriggerInfo(AliQA::TASKINDEX_t task);
  void InitDisplayHistos(AliQA::TASKINDEX_t task);
  Bool_t AddSortedPoint(Float_t currVal, TArrayF& position, const Float_t kResetValue);
  
  AliMUONVDigitStore* fDigitStore; //!< pointer to digits store
  AliMUONVTriggerStore* fTriggerStore; //!< pointer to trigger store
  AliMUONDigitMaker* fDigitMaker;  //!< pointer to digit maker

  ClassDef(AliMUONQADataMakerRec,2)  // MUON Quality assurance data maker

};
#endif
