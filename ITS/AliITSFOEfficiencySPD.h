#ifndef ALIITS_FOEFFICIENCYSPD_H
#define ALIITS_FOEFFICIENCYSPD_H

/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to store Fast-OR efficiency values in OCDB.  //
// One value per pixel chip in this base class (if per column      //
// accuracy is needed, use AliITSFOEfficiencySPDColumn class).     //
// The values are the probability that a pixel hit will generate a //
// fast-OR signal.                                                 //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TError.h>

class AliITSFOEfficiencySPD : public TObject {

 public:
  AliITSFOEfficiencySPD();
  AliITSFOEfficiencySPD(const AliITSFOEfficiencySPD& foEff);
  virtual ~AliITSFOEfficiencySPD();
  AliITSFOEfficiencySPD& operator=(const AliITSFOEfficiencySPD& foEff);

  virtual void    ResetValues();
  virtual void    SetChipEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, Float_t value);
  virtual void    SetColumnEfficiency(UInt_t /*eq*/, UInt_t /*hs*/, UInt_t /*chip*/, UInt_t /*col*/, Float_t /*value*/) 
    {Error("AliITSFOEfficiencySPD::SetColumnEfficiency","You need daughter class to set column efficiencies!");}

  virtual Float_t GetChipEfficiency(UInt_t eq, UInt_t hs, UInt_t chip) const;
  virtual Float_t GetColumnEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t /*col*/) const
    {return GetChipEfficiency(eq,hs,chip);}
  virtual Float_t GetEfficiency(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t /*col*/, UInt_t /*row*/) const
    {return GetChipEfficiency(eq,hs,chip);}

 protected:
  Float_t fChipEfficiency[20][6][10]; // efficiency values per chip

  ClassDef(AliITSFOEfficiencySPD,1) // FO Efficiency Base
};

#endif
