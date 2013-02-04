/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Henrik Qvigstad <henrik.qvigstad@cern.ch>
/* $Id$ */


#ifndef ALITRUPEDESTALOUTPUT_H
#define ALITRUPEDESTALOUTPUT_H

#include <TObject.h>
#include <Rtypes.h>

class TH1I;
class TH1F;
class TH2F;

class AliTRUPedestalOutput : public TObject
{

public:
  AliTRUPedestalOutput();
  virtual ~AliTRUPedestalOutput();

  void SetRun(Int_t run);
  void EventAdded();

  UInt_t GetEventsAdded() { return fNEvents; }

  // Histograms Getters:
  TH1F* GetPedestals();
  TH1F* GetPedestalRMS();
  TH1I* GetPedestalSamples();
  TH2F* GetPedestals2d(UInt_t mod);
  TH2F* GetPedestalRMS2d(UInt_t mod);
  TH1F* GetPedestalsId();
  TH1F* GetPedestals_branch(UInt_t mod, UInt_t row, UInt_t branch);

  TH1I* GetTRUSignals(UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z);

  // Other Getters:
  Double_t GetPedestal(UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z);
  Double_t GetPedestalError(UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z);
  Double_t GetRMS(UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z);
  Double_t GetSamples(UInt_t mod, UInt_t row, UInt_t branch, UInt_t x, UInt_t z);
  
  // Constants
  const static UInt_t kNMods         = 5;
  const static UInt_t kNTRURows      = 4;
  const static UInt_t kNBranches     = 2;
  const static UInt_t kN2x2X         = 64/2;
  const static UInt_t kN2x2Z         = 56/2;
  const static UInt_t kN2x2XPrTRURow = kN2x2X / kNTRURows;
  const static UInt_t kN2x2ZPrBranch = kN2x2Z / kNBranches;
  const static UInt_t kNTRUTimeBins  = 128;
  const static UInt_t kNEMCTimeBins  = 62;

private:
  AliTRUPedestalOutput ( const AliTRUPedestalOutput& other ); // not impl.
  AliTRUPedestalOutput& operator= ( const AliTRUPedestalOutput& other ); // not impl.
  
  Int_t fRun;
  UInt_t fNEvents;
  
  // Event Global Histograms:
  TH1F* fPedestals;   //! Pedestals
  TH1F* fPedestalRMS; //!
  TH1I* fPedestalSamples; //!
  TH2F* fPedestals2d[kNMods]; //!
  TH2F* fPedestalRMS2d[kNMods]; //!
  TH1F* fPedestalsId; //! Pedestals v Id
  TH1F* fPedestals_branch[kNMods][kNTRURows][kNBranches]; //! Pedestals, pr mod
  // TH2F* fPedestalsId_branch[kNMods][kNTRURows][kNBranches]; //! Pedestals v Id, pr mod

  // Regular Histograms
  TH1I* fTRUSignals[kNMods][kNTRURows][kNBranches][kN2x2XPrTRURow][kN2x2ZPrBranch]; //->

  
  ClassDef(AliTRUPedestalOutput, 0)
};

#endif // ALITRUPEDESTALOUTPUT_H
