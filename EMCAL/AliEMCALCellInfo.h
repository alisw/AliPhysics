#ifndef ALIEMCALCELLINFO_H
#define ALIEMCALCELLINFO_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
// This is a table of emcal indexes.
// Simplify a work with cell indexes 
// Initial version was created with TTable staff
// 

//*-- Authors: Aleksei Pavlinov (WSU)
#include <TNamed.h>

// Aug 1, 2007; Corr. Sep 05
class cellInfo : public TObject { 
  // See AliEMCALGeometry
  // Indexes information
 public:
  virtual const char* GetName() const {return Form("Ind%5.5i",absId);}
  cellInfo();
  virtual ~cellInfo() {};

  Int_t absId;   // abs id of cell as in Geant
  // Geant numbering tree - see AliEMCALv2
  Int_t nSupMod; // index of super module (SM)
  Int_t nModule; // index of module in SM
  Int_t nIphi;   // phi index of tower(cell) in module
  Int_t nIeta;   // eta index of tower(cell) in module
  // Inside SM - ised in cluster finder
  Int_t iPhi;    // phi index of tower(cell) in SM 
  Int_t iEta;    // eta index of tower(cell) in SM
  Int_t iPhim;   // phi index of module in SM
  Int_t iEtam;   // eta index of module in SM
  // Coordinate information should be include too ??
  ClassDef(cellInfo,1) // Cell indexes information
};

class AliEMCALGeometry;
class TObjArray;

class AliEMCALCellInfo : public TNamed {
 public:
  AliEMCALCellInfo(); // default constractor
  AliEMCALCellInfo(const char* name, const Int_t nrow);
  virtual ~AliEMCALCellInfo();

  AliEMCALCellInfo & operator = (const AliEMCALCellInfo  & /*rvalue*/) {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented");
    return *this;
  };
  // 
  void AddAt(cellInfo* r);
  cellInfo* GetTable(Int_t i) const;

  void PrintTable(int ind1=-1, int ind2=-1) const;  //*MENU*

  static AliEMCALCellInfo *GetTableForGeometry(const char* geoName);
  static AliEMCALCellInfo *GetTableForGeometry(AliEMCALGeometry *g);
  //
 protected:
  TObjArray *fTable;
  Int_t fCurrentInd;

  ClassDef(AliEMCALCellInfo,2) // Table of emcal indexes  
};

#endif // ALIEMCALCELLINFO_H
