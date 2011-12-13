#ifndef ALIPHOSTRIGGERRAWREADER_H
#define ALIPHOSTRIGGERRAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

#include "TObject.h"

class AliCaloRawStreamV3;
class AliPHOSTRURawReader;

/*
 *  Class for reading the Trigger Data Stream from Raw.
 *  Author: Henrik Qvigstad (henrik.qvigstad@cern.ch)
 */
class AliPHOSTriggerRawReader : public TObject
{
 public:
  AliPHOSTriggerRawReader();
 ~AliPHOSTriggerRawReader();
  
  AliPHOSTRURawReader* GetTRU(Int_t mod, Int_t truRow, Int_t branch);
  
  void ReadFromStream(AliCaloRawStreamV3* );
  void Reset();
  
 private:
  AliPHOSTriggerRawReader(const AliPHOSTriggerRawReader&); // not implemented
  AliPHOSTriggerRawReader& operator= (const AliPHOSTriggerRawReader&); // not implemented

 private:
  // constants
  const static Int_t kNMods = 5; // n. mods
  const static Int_t kNTRURows = 4; // n. tru rows
  const static Int_t kNBranches = 2; // n. branches

  AliPHOSTRURawReader* fTRUs[kNMods][kNTRURows][kNBranches]; // TRU raw readers [mod][truRow][branch]


  ClassDef(AliPHOSTriggerRawReader, 0)
};

#endif 
