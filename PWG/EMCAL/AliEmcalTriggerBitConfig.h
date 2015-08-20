#ifndef ALIEMCALTRIGGERBITCONFIG_H
#define ALIEMCALTRIGGERBITCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliLog.h"
#include <TNamed.h>

class AliEmcalTriggerBitConfig  : public TNamed {
public:
  AliEmcalTriggerBitConfig();
  AliEmcalTriggerBitConfig(Int_t l0bit, Int_t j1bit, Int_t j2bit, Int_t g1bit, Int_t g2bit, Int_t mcoffset);
  virtual ~AliEmcalTriggerBitConfig() {}

  void Initialise(const AliEmcalTriggerBitConfig &ref);

  Int_t GetLevel0Bit() const { if(fL0Bit < 0) AliFatal("Invalid trigger configuration: Level0 bit < 0"); return fL0Bit; }
  Int_t GetJetHighBit() const { if(fJHighBit < 0) AliFatal("Invalid trigger configuration: Jet high bit < 0"); return fJHighBit; }
  Int_t GetJetLowBit() const { if(fJLowBit < 0) AliFatal("Invalid trigger configuration: Jet low bit < 0"); return fJLowBit; }
  Int_t GetGammaHighBit() const { if(fGHighBit < 0) AliFatal("Invalid trigger configuration: Gamma high bit < 0"); return fGHighBit; }
  Int_t GetGammaLowBit() const { if(fGLowBit < 0) AliFatal("Invalid trigger configuration: Gamma low bit < 0"); return fGLowBit; }
  Int_t GetTriggerTypesEnd() const {if(fTriggerTypesEnd < 0) AliFatal("Invalid trigger configuration: MC Offset bit < 0"); return fTriggerTypesEnd; }

protected:
 Int_t fL0Bit;      // Level0 bit
 Int_t fJHighBit;   // Jet High bit
 Int_t fJLowBit;    // Jet Low bit
 Int_t fGHighBit;   // Gamma High bit
 Int_t fGLowBit;    // Gamma Low bit
 Int_t fTriggerTypesEnd;   // Monte-Carlo offset

 ClassDef(AliEmcalTriggerBitConfig, 1);
};

class AliEmcalTriggerBitConfigOld : public AliEmcalTriggerBitConfig{
public:
  AliEmcalTriggerBitConfigOld();
  virtual ~AliEmcalTriggerBitConfigOld() {}

  ClassDef(AliEmcalTriggerBitConfigOld, 1);
};

class AliEmcalTriggerBitConfigNew : public AliEmcalTriggerBitConfig{
public:
  AliEmcalTriggerBitConfigNew();
  virtual ~AliEmcalTriggerBitConfigNew() {}

  ClassDef(AliEmcalTriggerBitConfigNew, 1);
};

#endif /* ALIEMCALTRIGGERBITCONFIG_H */
