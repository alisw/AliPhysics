#ifndef ALIITS_FOEFFICIENCYSPDCOLUMN_H
#define ALIITS_FOEFFICIENCYSPDCOLUMN_H

/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to store Fast-OR efficiency values in OCDB.  //
// One value per pixel chip column in this daughter class.         //
// The values are the probability that a pixel hit will generate a //
// fast-OR signal.                                                 //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliITSFOEfficiencySPD.h"

class AliITSFOEfficiencySPDColumn : public AliITSFOEfficiencySPD {

 public:
  AliITSFOEfficiencySPDColumn();
  AliITSFOEfficiencySPDColumn(const AliITSFOEfficiencySPDColumn& foEff);
  virtual ~AliITSFOEfficiencySPDColumn();
  AliITSFOEfficiencySPDColumn& operator=(const AliITSFOEfficiencySPDColumn& foEff);

  virtual void    ResetValues();
  virtual void    SetColumnEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, Float_t value);

  virtual Float_t GetEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t /*row*/) const
    {return GetColumnEfficiency(eq,hs,chip,col);}
  virtual Float_t GetColumnEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col) const;

 protected:
  Float_t fColumnEfficiency[20][6][10][32]; // efficiency values per chip column

  ClassDef(AliITSFOEfficiencySPDColumn,1)   // FO Efficiency Per Column
};

#endif
