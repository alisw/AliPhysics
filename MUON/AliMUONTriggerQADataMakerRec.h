#ifndef ALIMUONTRIGGERQADATAMAKERREC_H
#define ALIMUONTRIGGERQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONTriggerQADataMakerRec
/// \brief MUON Quality assurance data maker for MTR
///

// --- AliRoot header files ---
#include "AliMUONVQADataMakerRec.h"

class AliMUONCalibrationData;
class AliMUONDigitMaker;
class AliMUONVClusterStore;
class AliMUONTriggerElectronics;
class AliMUONVDigitStore;
class TObjectArray;
class TMap;
class AliMpDCSNamer;
class AliMUONVTriggerStore;
class AliMUONGlobalTrigger;

class AliMUONTriggerQADataMakerRec: public AliMUONVQADataMakerRec {

public:
  AliMUONTriggerQADataMakerRec(AliQADataMakerRec* master);         
  virtual ~AliMUONTriggerQADataMakerRec();
  
  virtual void InitRaws(); 
  virtual void InitRecPoints(); 
  virtual void InitDigits(); 
  virtual void InitESDs(); 
  
  void EndOfDetectorCycleRaws(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleESDs(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleDigits(Int_t , TObjArray** ) {}
    
  virtual void MakeRaws(AliRawReader* rawReader); 
  
  virtual void MakeDigits(TTree* dig); 
  virtual void MakeRecPoints(TTree* recpo); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  
public:

  /// Raw histograms indices
  
  enum ERaw {
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
    kTriggerErrorOutGlobalFromInGlobal = 63, ///< Global board: Number of error vs output bit 
    kTriggerErrorSummary      = 64,  ///< Number of errors for each trigger decision level (Local, Reg->Local, Reg, Reg->Glob, Global)
    kTriggerErrorSummaryNorm  = 65,  ///< Percentage of errors for each trigger decision level
    kTriggerErrorLocalYCopy     = 67, ///< Local board: Number of Y Copy Error vs Local Board Id
    kTriggerErrorLocalYCopyTest = 68, ///< Local Board: Number of Y copy error tests (for normalization)
    kTriggerErrorLocalYCopyNorm = 69, ///< Local Board: Number of Y Copy Error vs Local Board Id Normalized to the number of tests
    kTriggeredBoards          = 70,  ///< Triggered boards histogram index
    kTriggerBoardsDisplay     = 71,  ///< Triggered boards display histogram index
    kTriggerReadOutErrors     = 80,  ///< Number of read-out errors
    kTriggerReadOutErrorsNorm = 81,  ///< Percentage of read-out errors
    kTriggerGlobalOutput      = 90,  ///< Number of Global outputs and Global algo errors
    kTriggerGlobalOutputNorm  = 91,  ///< Percentage of Global outputs and Global algo errors
    kRawNAnalyzedEvents       = 100  ///< Number of analyzed events per event specie
  };
         
  /// Rec points histograms indices
  enum ERecPoints { 
    kNAnalyzedEvents           = 0, ///< Number of analyzed events per event specie
    kTriggerRPCtrips           = 1, ///< Trips in trigger chambers
    kTriggerRPChv              = 2  ///< Trigger chamber HV index
  };
  
//  /// ESD histograms indices
//  enum EESD { 
//  };

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
  
private:

  AliMUONTriggerQADataMakerRec(const AliMUONTriggerQADataMakerRec& qadm);   
  AliMUONTriggerQADataMakerRec& operator=(const AliMUONTriggerQADataMakerRec& qadm);

  void DisplayTriggerInfo();
  Bool_t FillTriggerDCSHistos();
  TObjArray* GetDCSValues(Int_t iMeas, Int_t detElemId,
			  TMap* triggerDcsMap, AliMpDCSNamer& triggerDcsNamer);
  void RawTriggerInRegional2OutRegional();
  UChar_t RawTriggerInGlobal2OutGlobal(UInt_t globalInput[4]);
  void RawTriggerMatchOutLocal(AliMUONVTriggerStore& inputTriggerStore, AliMUONVTriggerStore& recoTriggerStore);
  //void RawTriggerMatchOutLocalInRegional();
  void RawTriggerMatchOutGlobalFromInGlobal(AliMUONGlobalTrigger& inputLocalTrigger,
					    AliMUONGlobalTrigger& recoGlobalTrigger);

  //Int_t fTriggerOutputRegionalData[16]; ///< Data Regional Trigger decision for each Regional Board (1R:0, 2R:1, ... , 1L:8, ...) -> 4 bits LPt, 4 bits HPt
  //Int_t fTriggerInputRegionalRecLPt[2][16][16]; ///< Reconstructed Regional Input LPt for each Regional Board ([bit][reg][loc]) (reg -> 1R:0, 2R:1, ... , 1L:8, ...)
  //Int_t fTriggerInputRegionalRecHPt[2][16][16]; ///< Reconstructed Regional Input HPt for each Regional Board ([bit][reg][loc]) (reg -> 1R:0, 2R:1, ... , 1L:8, ...)
  //Int_t fTriggerOutputRegionalRec[16]; ///< Reconstructed Regional Trigger decision for each Regional Board (8 Bits)

  //Int_t fTriggerInputGlobalDataLPt[16][4]; ///< Data Global inputs LPt (1R:0, 2R:1, ... , 1L:8, ...)
  //Int_t fTriggerInputGlobalDataHPt[16][4]; ///< Data Global inputs HPt (1R:0, 2R:1, ... , 1L:8, ...)
  //Int_t fTriggerOutputGlobalRecFromLocalInput[6]; //< Reconstructed Global outputs from Local inputs
  //Int_t fTriggerOutputGlobalRecFromLocalOutput[6]; //< Reconstructed Global outputs from Local outputs
  
  AliMUONDigitMaker* fDigitMaker; //!< pointer to digit maker
  AliMUONCalibrationData* fCalibrationData; //!< Used to load Local, Regional and Global masks
  AliMUONTriggerElectronics* fTriggerProcessor; //!< trigger processore to re-compute response
  AliMUONVDigitStore* fDigitStore; //!< pointer to digits store
  
  ClassDef(AliMUONTriggerQADataMakerRec,1)  // MUON Quality assurance data maker

};

#endif
