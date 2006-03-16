#ifndef AliALTROMAPPING_H
#define AliALTROMAPPING_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
// Class used to setup the mapping of hardware adresses //
// in ALTRO to pad-rows and pad indeces.                //
// The mapping is defined in an external mapping files  //
// separately for TPC,PHOS and FMD.                     //
//////////////////////////////////////////////////////////

#include <TObject.h>

class AliAltroMapping: public TObject {
 public:
  AliAltroMapping(const char *mappingFile);
  virtual ~AliAltroMapping();

  AliAltroMapping(const AliAltroMapping& mapping);
  AliAltroMapping& operator = (const AliAltroMapping& mapping);

  Int_t GetHWAdress(Int_t padrow, Int_t pad) const;
  Int_t GetPadRow(Int_t hwAdress) const;
  Int_t GetPad(Int_t hwAdress) const;

 private:
  Bool_t ReadMapping(const char *mappingFile);

  Int_t     fNumberOfChannels; // Number of ALTRO channels
  Int_t     fMaxHWAdress;      // Maximum HW adress
  Int_t     fMinPadRow;        // Minimum Index of pad-row
  Int_t     fMaxPadRow;        // Maximum Index of pad-row
  Int_t     fMaxPad;           // Maximum Index of pad inside row
  Short_t **fMapping;          // Array which connects hardware adresses to pad and pad-row indeces
  Short_t **fInvMapping;       // Inverse of fMapping

  ClassDef(AliAltroMapping,0)  // Altro mapping handler class
};

#endif
