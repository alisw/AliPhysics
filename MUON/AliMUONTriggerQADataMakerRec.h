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
    kTriggerScalers           = 22, ///< Trigger scalers histogram per plane index
    kTriggerScalersDisplay    = 30, ///< Trigger scalers display histogram per plane index
    kTriggerScalersTime       = 38, ///< Trigger scalers acquisition time index
    kTriggerCalibSummary      = 39, ///< Summary of responding strips/boards and noisy strips 
    kTriggeredBoards          = 40,  ///< Triggered boards histogram index
    kTriggerBoardsDisplay     = 41, ///< Triggered boards display histogram index
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

    kRawNAnalyzedEvents = 66, ///< Number of analyzed events per event specie
    kTriggerReadOutErrors = 67, ///< Number of read-out errors
    kTriggerGlobalOutput = 68 //< Histo including Global outputs and Global algo errors
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
  void RawTriggerInGlobal2OutGlobal();
  void RawTriggerMatchOutLocal();
  void RawTriggerMatchOutLocalInRegional();
  void RawTriggerMatchOutGlobalFromInGlobal();
	
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
  
  AliMUONDigitMaker* fDigitMaker; //!< pointer to digit maker
  AliMUONCalibrationData* fCalibrationData; //!< Used to load Local, Regional and Global masks
  AliMUONTriggerElectronics* fTriggerProcessor; //!< trigger processore to re-compute response
  AliMUONVDigitStore* fDigitStore; //!< pointer to digits store
  
  ClassDef(AliMUONTriggerQADataMakerRec,1)  // MUON Quality assurance data maker

};

#endif
