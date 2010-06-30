/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONTriggerChamberEfficiency
/// \brief Calculate, apply and possibly draw trigger chamber efficiency.
///
//  Author: Diego Stocco; Subatech, Nantes

#ifndef ALIMUONTRIGGERCHAMBEREFFICIENCY_H
#define ALIMUONTRIGGERCHAMBEREFFICIENCY_H

#include "TObject.h"
class TH1;
class TList;
class TObjArray;
class TGraphAsymmErrors;
class AliMUONTriggerEfficiencyCells;

class AliMUONTriggerChamberEfficiency : public TObject
{
public:
  AliMUONTriggerChamberEfficiency(AliMUONTriggerEfficiencyCells* effMap);
  AliMUONTriggerChamberEfficiency(const Char_t* filename,
				  const Char_t* listname="triggerChamberEff");

  AliMUONTriggerChamberEfficiency(const AliMUONTriggerChamberEfficiency& other); // copy constructor
  AliMUONTriggerChamberEfficiency& operator=(const AliMUONTriggerChamberEfficiency& other); // assignment operator

  virtual ~AliMUONTriggerChamberEfficiency();

  Float_t GetCellEfficiency(Int_t detElemId, Int_t localBoard, Int_t hType) const;
    
  void IsTriggered(Int_t detElemId, Int_t localBoard, Bool_t &trigBend, Bool_t &trigNonBend) const;

  Bool_t LowStatisticsSettings(Bool_t useMeanValues=kTRUE);

  // Methods for display
  void DisplayEfficiency(Bool_t perSlat=kFALSE, Bool_t show2Dhisto = kTRUE);

  enum{
    kHboardEff,     ///< Efficiency per board index
    kHslatEff       ///< Efficiency per slat index
  };

private:
    Int_t FindChamberIndex(Int_t detElemId) const;
    void FillFromList(Bool_t useMeanValues = kFALSE);

    Int_t GetIndex(Int_t histoType, Int_t countType, 
		   Int_t chamber = -1) const;

    TGraphAsymmErrors* GetEfficiencyGraph(TH1* histoNum, TH1* histoDen);

    Bool_t fIsOwner; ///< Owner of efficiency map
    AliMUONTriggerEfficiencyCells* fEfficiencyMap; ///< Efficiency map

    TObjArray* fEfficiencyObjects; ///< Collect all efficiency 
    TList* fDisplayList; //!< List of objects for display

    ClassDef(AliMUONTriggerChamberEfficiency,0) // Trigger efficiency store
};
#endif
