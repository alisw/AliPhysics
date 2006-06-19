/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup base
/// \class AliMUONTriggerEfficiencyCells
/// \brief Store and give access to the trigger chamber efficiency.
///
//  Author: Diego Stocco; INFN Torino

#ifndef ALIMUONTRIGGEREFFICIENCYCELLS_H
#define ALIMUONTRIGGEREFFICIENCYCELLS_H

#include "TObject.h"
#include "TArrayI.h"
#include "TVector2.h"
#include "TString.h"

class AliMUONTriggerEfficiencyCells : public TObject
{
public:
  AliMUONTriggerEfficiencyCells();
  AliMUONTriggerEfficiencyCells(const char* filename);

  virtual ~AliMUONTriggerEfficiencyCells();

  Float_t GetCellEfficiency(Int_t detElemId, Int_t cathode, Float_t x, Float_t y);
  void GetCellEfficiency(Int_t detElemId, Float_t x, Float_t y, Float_t &eff1, Float_t &eff2);
    
  Bool_t IsTriggered(Int_t detElemId, Int_t cathode, Float_t x, Float_t y);
  void IsTriggered(Int_t detElemId, Float_t x, Float_t y, Bool_t &trig1, Bool_t &trig2);

  TVector2 ChangeReferenceFrame(Float_t x, Float_t y, Float_t x0, Float_t y0);

  void Reset();
    
protected:
    TArrayI CellByCoord(Int_t detElemId, Float_t x, Float_t y);
    void ReadFile(const char* filename="$ALICE_ROOT/MUON/data/efficiencyCells.dat");

private:
    Int_t FindChamberIndex(Int_t detElemId);
    Int_t FindSlatIndex(Int_t detElemId);
    
    static const Int_t fgkNofCells=80; ///< Number of cells
    
    /// Cell content [trig. chambers][RPCs][cathode][cellsX][cellsY]
    Float_t fCellContent[4][18][2][fgkNofCells][fgkNofCells]; //[trig. chambers][RPCs][cathode][cellsX][cellsY]

    Float_t fCellSize[4][18][2]; ///< the size of the cells
    Int_t fCellNumber[4][18][2]; ///< id of the cells
    
    ClassDef(AliMUONTriggerEfficiencyCells,1) // Trigger efficiency store
};
#endif
