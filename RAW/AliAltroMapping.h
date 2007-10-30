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
  AliAltroMapping();
  AliAltroMapping(const char *mappingFile);
  virtual ~AliAltroMapping();

  virtual Int_t GetHWAddress(Int_t padrow, Int_t pad, Int_t sector) const = 0;
  virtual Int_t GetPadRow(Int_t hwAddress) const = 0;
  virtual Int_t GetPad(Int_t hwAddress) const = 0;
  virtual Int_t GetSector(Int_t hwAddress) const = 0;

 protected:
  Bool_t OpenMappingFile(const char *mappingFile);
  Bool_t CloseMappingFile();
  virtual Bool_t ReadMapping() = 0;
  virtual void   DeleteMappingArrays() = 0;

  ifstream *fIn;               //! External mapping file
  Int_t     fNumberOfChannels; // Number of ALTRO channels
  Int_t     fMaxHWAddress;     // Maximum HW adress
  Int_t     fMappingSize;      // Maximum size of the mapping array, used by the streamer of derived classes
  Int_t     fInvMappingSize;   // Maximum size of the inverse mapping arrays, used by the streamer of derived classes

 private:
  AliAltroMapping(const AliAltroMapping& mapping);
  AliAltroMapping& operator = (const AliAltroMapping& mapping);

  ClassDef(AliAltroMapping,3)  // Altro mapping handler class
};

#endif
