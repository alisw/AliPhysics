#ifndef ALITOFFEEREADER_H
#define ALITOFFEEREADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides the TOF FEE reader.                 //
//                                                           //
///////////////////////////////////////////////////////////////

/* 
 * author: Roberto Preghenella (R+), Roberto.Preghenella@bo.infn.it
 */

#include "TObject.h"
//#include "AliTOFGeometry.h"

class AliTOFFEEConfig;

class AliTOFFEEReader :
public TObject 
{

 public:

  AliTOFFEEReader(); // default constructor
  AliTOFFEEReader(const AliTOFFEEReader &source); // copy constructor
  AliTOFFEEReader &operator=(const AliTOFFEEReader &source); // operator =
  virtual ~AliTOFFEEReader(); // default destructor

  /* getters */
  static Int_t GetNumberOfDDLs() {return fgkNumberOfDDLs;}; // get number of DDLs
  static Int_t GetNumberOfTRMs() {return fgkNumberOfTRMs;}; // get number of TRMs
  static Int_t GetNumberOfChains() {return fgkNumberOfChains;}; // get number of chains
  static Int_t GetNumberOfTDCs() {return fgkNumberOfTDCs;}; // get number of TDCs
  static Int_t GetNumberOfChannels() {return fgkNumberOfChannels;}; // get number of channels
  static Int_t GetNumberOfIndexes() {return fgkNumberOfIndexes;}; // get number of indexes
  AliTOFFEEConfig *GetFEEConfig() const {return fFEEConfig;}; // get FEE config
  Bool_t GetChannelEnabled(Int_t iIndex) const {return iIndex < GetNumberOfIndexes() ? fChannelEnabled[iIndex] : kFALSE;}; // get channel enabled
  Int_t GetMatchingWindow(Int_t iIndex) const {return iIndex < GetNumberOfIndexes() ? fMatchingWindow[iIndex] : 0;}; // get matching window
  
  /* setters */
  
  /* methods */
  void LoadFEEConfig(const Char_t *FileName) const; // load FEE config
  void DumpFEEConfig(); // dump FEE config
  Int_t ParseFEEConfig(); // parse FEE config
  void ResetChannelEnabledArray(); // reset channel enabled array
  void Reset(); // reset channel enabled array
  Bool_t IsChannelEnabled(Int_t iDDL, Int_t iTRM, Int_t iChain, Int_t iTDC, Int_t iChannel) const; // is channel enabled
  Bool_t IsChannelEnabled(Int_t iIndex) const {return GetChannelEnabled(iIndex);}; // is channel enabled
  Int_t GetMatchingWindow(Int_t iDDL, Int_t iTRM, Int_t iChain, Int_t iTDC, Int_t iChannel) const; // get matching window
  
 private:

  static const Int_t fgkNumberOfDDLs = 72; // number of DDLs
  static const Int_t fgkNumberOfTRMs = 10; // number of TRMs
  static const Int_t fgkNumberOfChains = 2; // number of chains
  static const Int_t fgkNumberOfTDCs = 15; // number of TDCs
  static const Int_t fgkNumberOfChannels = 8; // number of channels
  static const Int_t fgkNumberOfIndexes = 157248; // number of indexes

  AliTOFFEEConfig *fFEEConfig; // FEE config
  Bool_t fChannelEnabled[fgkNumberOfIndexes]; // channel enabled
  Int_t fMatchingWindow[fgkNumberOfIndexes]; // matching window

  ClassDef(AliTOFFEEReader, 1);

};

#endif /* ALITOFFEEREADER_H */
