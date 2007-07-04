/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup calib
/// \class AliMUONTriggerEfficiencyCells
/// \brief Store and give access to the trigger chamber efficiency.
///
//  Author: Diego Stocco; INFN Torino

#ifndef ALIMUONTRIGGEREFFICIENCYCELLS_H
#define ALIMUONTRIGGEREFFICIENCYCELLS_H

#include "TObject.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "TVector2.h"
#include "TString.h"
#include "TMatrix.h"

class AliMUONTriggerEfficiencyCells : public TObject
{
public:
  AliMUONTriggerEfficiencyCells();
  AliMUONTriggerEfficiencyCells(const char* filename);

  virtual ~AliMUONTriggerEfficiencyCells();

  void GetCellEfficiency(Int_t detElemId, Float_t x, Float_t y, Float_t &eff1, Float_t &eff2) const;
  void GetCellEfficiency(Int_t detElemId, Int_t localBoard, Float_t &eff1, Float_t &eff2) const;
    
  void IsTriggered(Int_t detElemId, Float_t x, Float_t y, Bool_t &trig1, Bool_t &trig2) const;
  void IsTriggered(Int_t detElemId, Int_t localBoard, Bool_t &trig1, Bool_t &trig2) const;

  TVector2 ChangeReferenceFrame(Float_t x, Float_t y, Float_t x0, Float_t y0);

  void Reset();
    
protected:
    TArrayI CellByCoord(Int_t detElemId, Float_t x, Float_t y) const;
    void ReadFile(const char* filename="$ALICE_ROOT/MUON/data/efficiencyCells.dat");

private:
    void CheckConstants() const;
    Int_t FindChamberIndex(Int_t detElemId) const;
    Int_t FindSlatIndex(Int_t detElemId) const;
    void ReadFileXY(ifstream &file);
    void ReadFileBoards(ifstream &file);
    void ReadHistoBoards(const char* filename="MUON.TriggerEfficiencyMap.root");
    
    static const Int_t fgkNcells=80;   ///< Number of cells
    static const Int_t fgkNcathodes=2; ///<Number of cathodes
    static const Int_t fgkNchambers=4; ///<Number of chambers
    static const Int_t fgkNplanes=8;   ///<Number of planes
    static const Int_t fgkNslats=18;   ///<Number of slats

    
    TMatrixF fCellContent[fgkNplanes][fgkNslats]; ///< the cells content
    TArrayF  fCellSize[fgkNplanes];    ///< the size of the cells
    TArrayI  fCellNumber[fgkNplanes];  ///< id of the cells
    TArrayF  fBoardContent[fgkNplanes];///< the boards content

    ClassDef(AliMUONTriggerEfficiencyCells,3) // Trigger efficiency store
};
#endif
