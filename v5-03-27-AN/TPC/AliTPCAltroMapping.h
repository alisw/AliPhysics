#ifndef ALITPCALTROMAPPING_H
#define ALITPCALTROMAPPING_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
// Class used to setup the mapping of hardware adresses //
// in ALTRO to pad-rows and pad indeces.                //
// The mapping is defined in an external mapping files  //
// separately. The class derives from the base altro    //
// mapping class defined in the RAW package.            //
//////////////////////////////////////////////////////////

#include "AliAltroMapping.h"

class AliTPCAltroMapping: public AliAltroMapping {
 public:
  AliTPCAltroMapping();
  AliTPCAltroMapping(const char *mappingFile);
  virtual ~AliTPCAltroMapping();

  virtual Int_t GetHWAddress(Int_t padrow, Int_t pad, Int_t sector);
  virtual Int_t GetPadRow(Int_t hwAddress) const;
  virtual Int_t GetPad(Int_t hwAddress) const;
  virtual Int_t GetSector(Int_t hwAddress) const;

 protected:
  virtual Bool_t ReadMapping();
  virtual Bool_t CreateInvMapping();

  Int_t     fMinPadRow;        // Minimum Index of pad-row
  Int_t     fMaxPadRow;        // Maximum Index of pad-row
  Int_t     fMaxPad;           // Maximum Index of pad inside row
  Short_t  *fInvMapping;       //! Inverse of fMapping

 private:

  AliTPCAltroMapping(const AliTPCAltroMapping& mapping);
  AliTPCAltroMapping& operator = (const AliTPCAltroMapping& mapping);

  ClassDef(AliTPCAltroMapping,3)  // Altro mapping handler class
};

#endif
