#ifndef ALIPHOSTRIGGERRAWDIGIPRODUCER_H
#define ALIPHOSTRIGGERRAWDIGIPRODUCER_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

class AliPHOSTriggerRawReader;
class AliPHOSTriggerParameters;
class AliRawReader;
class AliCaloRawStreamV3;

#include "TString.h"
#include "TClonesArray.h"

#include <vector>

class AliPHOSTriggerRawDigiProducer
{
 public:
  
  AliPHOSTriggerRawDigiProducer();
  AliPHOSTriggerRawDigiProducer(AliRawReader *rawReader);
  
  virtual ~AliPHOSTriggerRawDigiProducer();
  
  void ProcessEvent(TClonesArray* tdigits);

  void ProcessL0(TClonesArray* tdigits);
  void ProcessL1(TClonesArray* tdigits);

  void GetGammaPatchXY(Int_t itru, Int_t ieta, Int_t iphi, Int_t& x, Int_t& y);
  void GetL1GammaPatchModuleXZ(Int_t itru, Int_t xglob, Int_t yglob, Int_t& module, Int_t& x, Int_t& z);

  void SetTriggerParameters(AliPHOSTriggerParameters* parameters) {fParameters = parameters;}
  void SetAnalyseModule(int mod, bool analyse = true) {fModules[mod] = analyse;}
    
  static int Get2x2Max(AliPHOSTriggerRawReader*, AliPHOSTriggerParameters*, int mod, int xIdx, int zIdx);
  static int Get2x2Signal(AliPHOSTriggerRawReader*, AliPHOSTriggerParameters*, int mod, int xIdx, int zIdx, int timeBin);
  static int Get4x4Max(AliPHOSTriggerRawReader*, AliPHOSTriggerParameters*, int mod, int TRURow, int branch, int xIdx, int zIdx);
  static int Get4x4Signal(AliPHOSTriggerRawReader*, AliPHOSTriggerParameters*, int mod, int TRURow, int branch, int xIdx, int zIdx, int timeBin);
  
  static bool Is2x2Active(AliPHOSTriggerRawReader*, int mod, int xIdx, int zIdx);
  static bool Is2x2Active(AliPHOSTriggerRawReader*, int mod, int xIdx, int zIdx, int timeBin);
    
  const static int kNMods = 5; // number of PHOS modules
  const static int kNTRURows = 4; // number of TRU rows
  const static int kNBranches = 2; // number of branches
  const static int kN2x2X = 32; // (=64/2) Number of 2x2 in X direction
  const static int kN2x2Z = 28; // (=56/2) Number of 2x2 in Z direction
  const static Int_t kN2x2XPrTRURow = 8; // (=64/2/4) Number of 2x2 pr. row
  const static Int_t kN2x2ZPrBranch = 14; // (=56/2/2) Number of 2x2 pr. branch
  const static Int_t kN4x4XPrTRURow = 7; // (=64/2/4 -1) Number of 4x4 pr. row
  const static Int_t kN4x4ZPrBranch = 13; // (=56/2/2) -1 Number of 4x4 pr. branch  
  const static int kNTRUTimeBins = 128; // number of TRU time bins
  const static int kNDefaultNEMCTimeBins = 62;

 private:  
  AliPHOSTriggerRawDigiProducer(const AliPHOSTriggerRawDigiProducer &tdp); // not implemented
  AliPHOSTriggerRawDigiProducer& operator= (const AliPHOSTriggerRawDigiProducer &tdp); // not implemented

 protected:
  std::vector<bool> fModules; // , per module: should analyser analyse module
  UShort_t fSaturationThreshold;
  AliPHOSTriggerParameters* fParameters;
  
private:
  AliRawReader            * fRawReader;       //! Raw data reader
  AliCaloRawStreamV3      * fRawStream;       //! Calorimeter decoder of ALTRO format
  AliPHOSTriggerRawReader * fTriggerReader;   //! TriggerRawReader
  static const Int_t fgkSTUDDL = 20; //! DDL ID of the PHOS STU

  ClassDef(AliPHOSTriggerRawDigiProducer,3)
};

#endif
