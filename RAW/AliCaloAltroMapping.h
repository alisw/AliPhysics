#ifndef ALICALOALTROMAPPING_H
#define ALICALOALTROMAPPING_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////
// Class used to setup the mapping of hardware adresses //
// in ALTRO to Calo cells (column and row index +       //
// low/high gain.                                       //
// The mapping is defined in an external mapping files  //
// separately. The class derives from the base altro    //
// mapping class defined in the RAW package.            //
//////////////////////////////////////////////////////////

/// Exported from PHOS to be used also by EMCAL
/// November 2006 Gustavo Conesa Balbastre

#include "AliAltroMapping.h"

class AliCaloAltroMapping: public AliAltroMapping {
 public:
  AliCaloAltroMapping();
  AliCaloAltroMapping(const char *mappingFile);
  virtual ~AliCaloAltroMapping();

  // In case of PHOS/EMCAL the relevant segmentation is row-column-gain
  // or eta-phi-gain
  virtual Int_t GetHWAddress(Int_t row, Int_t column, Int_t gain);
  // Get Row (not pad-row as in the base class)
  virtual Int_t GetPadRow(Int_t hwAddress) const;
  // Get Column (not pad as in the base class)
  virtual Int_t GetPad(Int_t hwAddress) const;
  // Get Gain (0/1) (not sector as in the base class)
  virtual Int_t GetSector(Int_t hwAddress) const;

 protected:
  virtual Bool_t ReadMapping();
  virtual Bool_t CreateInvMapping();

  Int_t     fMinRow;        // Minimum row index
  Int_t     fMaxRow;        // Maximum row index
  Int_t     fMinCol;        // Minimum column index
  Int_t     fMaxCol;        // Maximum column index
  Short_t  *fInvMappingLow; //! Inverse of fMapping (Low gain)
  Short_t  *fInvMappingHigh;//! Inverse of fMapping (High gain)

 private:

  AliCaloAltroMapping(const AliCaloAltroMapping& mapping);
  AliCaloAltroMapping& operator = (const AliCaloAltroMapping& mapping);

  ClassDef(AliCaloAltroMapping,3)  // Altro mapping handler class
};

#endif
