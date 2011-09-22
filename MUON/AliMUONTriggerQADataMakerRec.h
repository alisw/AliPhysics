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
  
  /// Empty implementation
  void EndOfDetectorCycleDigits(Int_t , TObjArray** ) {}
    
  virtual void MakeRaws(AliRawReader* rawReader); 
  
  virtual void MakeDigits(TTree* dig); 
  virtual void MakeRecPoints(TTree* recpo); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  
  void ResetDetectorRaws(TObjArray* list);
  
private:
  /// Not implemented
  AliMUONTriggerQADataMakerRec(const AliMUONTriggerQADataMakerRec& qadm);   
  /// Not implemented
  AliMUONTriggerQADataMakerRec& operator=(const AliMUONTriggerQADataMakerRec& qadm);

  void DisplayTriggerInfo();
  void FillRatio4434Histos(Int_t evtInterval);
  Bool_t FillTriggerDCSHistos();
  TObjArray* GetDCSValues(Int_t iMeas, Int_t detElemId,
			  TMap* triggerDcsMap, AliMpDCSNamer& triggerDcsNamer);
  UChar_t RawTriggerInGlobal2OutGlobal(UInt_t globalInput[4]);
  void RawTriggerMatchOutLocal(const AliMUONVTriggerStore& inputTriggerStore, const AliMUONVTriggerStore& recoTriggerStore);
  //void RawTriggerMatchOutLocalInRegional();
  void RawTriggerMatchOutGlobal(AliMUONGlobalTrigger& inputLocalTrigger,
				AliMUONGlobalTrigger& recoGlobalTrigger,
				Char_t histo);
  AliMUONTriggerElectronics* TriggerElectronics();
  AliMUONCalibrationData* CalibrationData();

  //Int_t fTriggerOutputRegionalData[16]; ///< Data Regional Trigger decision for each Regional Board (1R:0, 2R:1, ... , 1L:8, ...) -> 4 bits LPt, 4 bits HPt
  //Int_t fTriggerInputRegionalRecLPt[2][16][16]; ///< Reconstructed Regional Input LPt for each Regional Board ([bit][reg][loc]) (reg -> 1R:0, 2R:1, ... , 1L:8, ...)
  //Int_t fTriggerInputRegionalRecHPt[2][16][16]; ///< Reconstructed Regional Input HPt for each Regional Board ([bit][reg][loc]) (reg -> 1R:0, 2R:1, ... , 1L:8, ...)
  //Int_t fTriggerOutputRegionalRec[16]; ///< Reconstructed Regional Trigger decision for each Regional Board (8 Bits)

  //Int_t fTriggerInputGlobalDataLPt[16][4]; ///< Data Global inputs LPt (1R:0, 2R:1, ... , 1L:8, ...)
  //Int_t fTriggerInputGlobalDataHPt[16][4]; ///< Data Global inputs HPt (1R:0, 2R:1, ... , 1L:8, ...)
  //Int_t fTriggerOutputGlobalRecFromLocalInput[6]; //< Reconstructed Global outputs from Local inputs
  //Int_t fTriggerOutputGlobalRecFromLocalOutput[6]; //< Reconstructed Global outputs from Local outputs

  static const Int_t fgkUpdateRatio4434=50; ///< Event interval between 2 update of the Ratio4434 histos
  
  AliMUONDigitMaker* fDigitMaker; //!< pointer to digit maker
  AliMUONCalibrationData* fCalibrationData; //!< Used to load Local, Regional and Global masks
  AliMUONTriggerElectronics* fTriggerProcessor; //!< trigger processore to re-compute response
  AliMUONVDigitStore* fDigitStore; //!< pointer to digits store
  
  ClassDef(AliMUONTriggerQADataMakerRec,1)  // MUON Quality assurance data maker

};

#endif
