/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Henrik Qvigstad <henrik.qvigstad@cern.ch>
/* $Id$ */


#ifndef ALITRUPEDESTALANALYSIS_H
#define ALITRUPEDESTALANALYSIS_H

#include <Rtypes.h>

class AliCaloRawStreamV3;
class AliTRUPedestalOutput;
class AliPHOSTriggerRawReader;

class AliTRUPedestalAnalysis
{
public:
  AliTRUPedestalAnalysis();
  virtual ~AliTRUPedestalAnalysis();

  void ProcessEvent(AliCaloRawStreamV3* );
  AliTRUPedestalOutput* GetOutput() { return fOutput; }

  static UInt_t Get2x2Max(AliPHOSTriggerRawReader* reader, int mod, int row, int branch, int x, int z);

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
  AliTRUPedestalAnalysis ( const AliTRUPedestalAnalysis& other ); // not impl.
  AliTRUPedestalAnalysis& operator= ( const AliTRUPedestalAnalysis& other ); // not impl.

  AliTRUPedestalOutput* fOutput;
  AliPHOSTriggerRawReader* fTriggerReader;
};

#endif // ALITRUPEDESTALANALYSIS_H
