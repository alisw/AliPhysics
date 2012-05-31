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

  virtual Int_t GetHWAddress(Int_t padrow, Int_t pad, Int_t sector) = 0;
  virtual Int_t GetPadRow(Int_t hwAddress) const = 0;
  virtual Int_t GetPad(Int_t hwAddress) const = 0;
  virtual Int_t GetSector(Int_t hwAddress) const = 0;

 protected:
  void           CloseMappingFile();
  virtual Bool_t ReadMapping() = 0;
  virtual Bool_t CreateInvMapping() = 0;

  ifstream *fIn;               //! External mapping file
  Int_t     fNumberOfChannels; // Number of ALTRO channels
  Int_t     fMaxHWAddress;     // Maximum HW adress
  Int_t     fMappingSize;      // Size of the mapping array, used by the streamer of derived classes
  Short_t  *fMapping;          //[fMappingSize] Array which connects hardware adresses to detector element indices

 private:
  Bool_t    OpenMappingFile(const char *mappingFile);

  AliAltroMapping(const AliAltroMapping& mapping);
  AliAltroMapping& operator = (const AliAltroMapping& mapping);

  ClassDef(AliAltroMapping,4)  // Altro mapping handler class
};

#endif
