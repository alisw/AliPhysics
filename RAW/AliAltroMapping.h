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

  virtual Int_t GetHWAdress(Int_t padrow, Int_t pad, Int_t sector) const = 0;
  virtual Int_t GetPadRow(Int_t hwAdress) const = 0;
  virtual Int_t GetPad(Int_t hwAdress) const = 0;
  virtual Int_t GetSector(Int_t hwAdress) const = 0;

 protected:
  Bool_t OpenMappingFile(const char *mappingFile);
  Bool_t CloseMappingFile();
  virtual Bool_t ReadMapping() = 0;
  virtual void   DeleteMappingArrays() = 0;

  ifstream *fIn;               // External mapping file
  Int_t     fNumberOfChannels; // Number of ALTRO channels
  Int_t     fMaxHWAdress;      // Maximum HW adress

  ClassDef(AliAltroMapping,0)  // Altro mapping handler class
};

#endif
