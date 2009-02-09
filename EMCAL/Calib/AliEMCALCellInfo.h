#ifndef ALIEMCALCELLINFO_H
#define ALIEMCALCELLINFO_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: AliEMCALCellInfo.h 24500 2008-03-13 23:39:38Z jklay $ */

//_________________________________________________________________________
// This is a table of emcal indexes.
// Simplify a work with cell indexes 
// Initial version was created with TTable staff
// 
//*-- Authors: Aleksei Pavlinov (WSU)
#include <TNamed.h>

// Aug 1, 2007; Corr. Sep 05
// cellInfo -> AliEMCALCellIndexes - Oct 15, 2007

class AliEMCALCellIndexes : public TObject { 
  // See AliEMCALGeometry
  // Indexes information
  friend class AliEMCALCellInfo;
  friend class AliEMCALPi0Calibration;
 public:
  virtual const char* GetName() const {return Form("Ind%5.5i",fAbsId);}
  AliEMCALCellIndexes();
  virtual ~AliEMCALCellIndexes() {};

 protected:
  Int_t fAbsId;   // abs id of cell as in Geant
  // Geant numbering tree - see AliEMCALv2
  Int_t fNSupMod; // index of super module (SM)
  Int_t fNModule; // index of module in SM
  Int_t fNIphi;   // phi index of tower(cell) in module
  Int_t fNIeta;   // eta index of tower(cell) in module
  // Inside SM - used in cluster finder
  Int_t fIPhi;    // phi index of tower(cell) in SM 
  Int_t fIEta;    // eta index of tower(cell) in SM
  Int_t fIPhim;   // phi index of module in SM
  Int_t fIEtam;   // eta index of module in SM
  // Coordinate information should be include too ??
  ClassDef(AliEMCALCellIndexes,2) // Cell indexes information
};

class AliEMCALGeometry;
class TObjArray;

class AliEMCALCellInfo : public TNamed {
 public:
  AliEMCALCellInfo(); // default constractor
  AliEMCALCellInfo(const AliEMCALCellInfo& info); //copy constructor
  AliEMCALCellInfo(const char* name, const Int_t nrow);
  virtual ~AliEMCALCellInfo();

  AliEMCALCellInfo & operator = (const AliEMCALCellInfo  & /*rvalue*/) {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented");
    return *this;
  };
  // 
  void AddAt(AliEMCALCellIndexes* r);
  AliEMCALCellIndexes* GetTable(Int_t i) const;

  void PrintTable(int ind1=-1, int ind2=-1) const;  //*MENU*

  static AliEMCALCellInfo *GetTableForGeometry(const char* geoName);
  static AliEMCALCellInfo *GetTableForGeometry(AliEMCALGeometry *g);
  //
 protected:
  TObjArray *fTable; // Array of  AliEMCALCellIndexes
  Int_t fCurrentInd; // Current index

  ClassDef(AliEMCALCellInfo,2) // Table of emcal indexes  
};

#endif // ALIEMCALCELLINFO_H
