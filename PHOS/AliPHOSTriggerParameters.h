#ifndef ALIPHOSTRIGGERPARAMETERS_H
#define ALIPHOSTRIGGERPARAMETERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Henrik Qvigstad <henrik.qvigstad@cern.ch> 17/10-2011
/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for PHOS Trigger Parameters                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TNamed.h"

/*  class for TRU Calib Data: Pedestals, etc...
 */
class AliPHOSTriggerParameters : public TNamed
{

public:
  AliPHOSTriggerParameters();
  AliPHOSTriggerParameters(const char* name);
  AliPHOSTriggerParameters(const AliPHOSTriggerParameters& );
  AliPHOSTriggerParameters& operator= (const AliPHOSTriggerParameters& );
  virtual ~AliPHOSTriggerParameters();
  
  // Getters
  UShort_t GetTRUPedestal(Int_t mod, Int_t TRURow, Int_t branch, Int_t xIdx, Int_t zIdx) const;
  Bool_t GetTRUReadoutOn(Int_t mod, Int_t TRURow, Int_t branch) const;
  Bool_t GetTRUSignalReadoutOn(Int_t mod, Int_t TRURow, Int_t branch) const;
  UShort_t GetTRUSignalTimeBinFrom(Int_t mod, Int_t TRURow, Int_t branch) const;
  UShort_t GetTRUSignalTimeBinTo(Int_t mod, Int_t TRURow, Int_t branch) const;
  UShort_t GetTRUThreshold(Int_t mod, Int_t TRURow, Int_t branch) const;
  UShort_t GetTRUMaskChannel(Int_t mod, Int_t TRURow, Int_t branch) const;
  const UShort_t* GetTORMaskArray(Int_t mod, Int_t tor) const;
  const UShort_t* GetTORReadoutMask(Int_t mod, Int_t tor) const;
  
  
  // Setters
  void SetTRUPedestal(UShort_t pedestal, Int_t mod, Int_t TRURow, Int_t branch, Int_t xIdx, Int_t zIdx);
  void SetTRUReadoutOn(Bool_t isOn, Int_t mod, Int_t TRURow, Int_t branch);
  void SetTRUSignalReadoutOn(Bool_t isOn, Int_t mod, Int_t TRURow, Int_t branch);
  void SetTRUSignalTimeBinFrom(UShort_t fromBin, Int_t mod, Int_t TRURow, Int_t branch);
  void SetTRUSignalTimeBinTo(UShort_t toBin, Int_t mod, Int_t TRURow, Int_t branch);
  void SetTRUThreshold(UShort_t threshold, Int_t mod, Int_t TRURow, Int_t branch);
  void SetTRUMaskChannel(UShort_t mask, Int_t mod, Int_t TRURow, Int_t branch);
  void SetTORMaskArray(const UShort_t ma[3], Int_t mod, Int_t tor);
  void SetTORReadoutMask(const UShort_t rm[2], Int_t mod, Int_t tor);
  
  // Misc
  virtual void Print(Option_t *option = "") const; 
  void Reset();

  // Constants
  static const Int_t kNMods     = 5; // Number of PHOS Modules
  static const Int_t kNTORs     = 2; // Number of TORs per Module
  static const Int_t kNTRURows  = 4; // Number of TRU rows
  static const Int_t kNBranches = 2; // Number of Branches
  static const Int_t kNTRUX     = 8; // Number of 2x2 per TRU in x
  static const Int_t kNTRUZ     = 14; // Number of 2x2 per TRU in z
  static const UShort_t kIdealTRUPedestal    = 512; // Ideal TRU Pedestal
  static const Int_t    kDefaultNTRUTimeBins = 128; // Number of timebins
  
protected:
  // TRU Parameters:
  UShort_t fTRUPedestals          [kNMods][kNTRURows][kNBranches][kNTRUX][kNTRUZ]; // TRU Pedestals
  Bool_t   fTRUTriggerBitReadoutOn[kNMods][kNTRURows][kNBranches]; // TRU TriggerBit Readout is on
  Bool_t   fTRUSignalReadoutOn    [kNMods][kNTRURows][kNBranches]; // TRU Signal Readout is on
  UChar_t  fTRUSignalTimeBinFrom  [kNMods][kNTRURows][kNBranches]; // TRU from (including) timebin
  UChar_t  fTRUSignalTimeBinTo    [kNMods][kNTRURows][kNBranches]; // TRU to (including) timebin
  UShort_t fTRUThreshold          [kNMods][kNTRURows][kNBranches]; // TRU Threshold
  UShort_t fTRUMaskChannel        [kNMods][kNTRURows][kNBranches]; // TRU Mask Channel
  
  // TOR Parameters:
  UShort_t fTORMaskArray[kNMods][kNTORs][3]; // TOR Mask Array
  UShort_t fTORReadoutMask[kNMods][kNTORs][2]; // TOR Readout Mask
  
  ClassDef(AliPHOSTriggerParameters, 0) // PHOS Trigger Parameters
};

#endif
