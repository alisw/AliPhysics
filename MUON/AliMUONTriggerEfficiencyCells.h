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
#include "TMatrix.h"
#include "TH1F.h"
#include "TList.h"

class AliMUONTriggerEfficiencyCells : public TObject
{
public:
  AliMUONTriggerEfficiencyCells();
  AliMUONTriggerEfficiencyCells(const Char_t* filename);
  AliMUONTriggerEfficiencyCells(TList *countHistoList, TList *noCountHistoList);

  AliMUONTriggerEfficiencyCells(const AliMUONTriggerEfficiencyCells& other); // copy constructor
  AliMUONTriggerEfficiencyCells& operator=(const AliMUONTriggerEfficiencyCells& other); // assignment operator

  virtual ~AliMUONTriggerEfficiencyCells();

  void GetCellEfficiency(Int_t detElemId, Float_t x, Float_t y, Float_t &eff1, Float_t &eff2) const;
  void GetCellEfficiency(Int_t detElemId, Int_t localBoard, Float_t &eff1, Float_t &eff2) const;
    
  void IsTriggered(Int_t detElemId, Float_t x, Float_t y, Bool_t &trig1, Bool_t &trig2) const;
  void IsTriggered(Int_t detElemId, Int_t localBoard, Bool_t &trig1, Bool_t &trig2) const;

  void DisplayEfficiency(Bool_t perSlat=kFALSE, const Char_t* geoFilename="geometry.root");
  Bool_t SumRunEfficiency(const AliMUONTriggerEfficiencyCells &other);

  void SetFiredStrips(TList *firedStrips){fFiredStrips = firedStrips;}
  void CheckFiredStrips(const Char_t *geoFilename="geometry.root");
                                // Check for strips with lower counts than others:
                                // syntomatic of possible read-out problems in boards
  void Reset();
    
protected:
    TArrayI CellByCoord(Int_t detElemId, Float_t x, Float_t y) const;
    TVector2 ChangeReferenceFrame(Float_t x, Float_t y, Float_t x0, Float_t y0);
    void ReadFile(const Char_t* filename="$ALICE_ROOT/MUON/data/efficiencyCells.dat");
    void CalculateEfficiency(Int_t trigger44, Int_t trigger34,
			     Float_t &efficiency, Float_t &error,
			     Bool_t failuresAsInput);


private:
    void CheckConstants() const;
    Int_t FindChamberIndex(Int_t detElemId) const;
    Int_t FindSlatIndex(Int_t detElemId) const;
    void ReadFileXY(ifstream &file);
    void ReadFileBoards(ifstream &file);
    void ReadHistoBoards(const Char_t* filename="MUON.TriggerEfficiencyMap.root");
    void InitHistos();
    void FillHistosFromList();
    Bool_t GetListsForCheck(const Char_t* geoFilename="geometry.root");
    
    static const Int_t fgkNcells=80;   ///< Number of cells
    static const Int_t fgkNcathodes=2; ///<Number of cathodes
    static const Int_t fgkNchambers=4; ///<Number of chambers
    static const Int_t fgkNplanes=8;   ///<Number of planes
    static const Int_t fgkNslats=18;   ///<Number of slats

    
    TMatrixF fCellContent[fgkNplanes][fgkNslats]; ///< the cells content
    TArrayF  fCellSize[fgkNplanes];    ///< the size of the cells
    TArrayI  fCellNumber[fgkNplanes];  ///< id of the cells
    
    TH1F *fBoardEfficiency[fgkNplanes];///< the boards content
    TH1F *fSlatEfficiency[fgkNplanes];///< the slats content

    TList *fCountHistoList; ///<list of efficiency numerators
    TList *fNoCountHistoList; ///<list of efficiency denominators
    TList *fFiredStrips; ///<list of fired strips for efficiency check

    TList *fDisplayHistoList; //! list of efficiency histograms for display
    TList *fBoardLabelList; //! list of board labels for display
    TList *fFiredFitHistoList; //! list of fired strips for checks
    TList *fFiredDisplayHistoList; //! list of fired strips for display

    ClassDef(AliMUONTriggerEfficiencyCells,5) // Trigger efficiency store
};
#endif
