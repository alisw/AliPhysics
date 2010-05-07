/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONTriggerEfficiencyCells
/// \brief Store and give access to the trigger chamber efficiency.
///
//  Author: Diego Stocco; INFN Torino

#ifndef ALIMUONTRIGGEREFFICIENCYCELLS_H
#define ALIMUONTRIGGEREFFICIENCYCELLS_H

#include "TObject.h"
class TH1F;
class TList;

class AliMUONTriggerEfficiencyCells : public TObject
{
public:
  AliMUONTriggerEfficiencyCells();
  AliMUONTriggerEfficiencyCells(const Char_t* filename, const Char_t* listname="triggerChamberEff");
  AliMUONTriggerEfficiencyCells(TList *countHistoList);

  AliMUONTriggerEfficiencyCells(const AliMUONTriggerEfficiencyCells& other); // copy constructor
  AliMUONTriggerEfficiencyCells& operator=(const AliMUONTriggerEfficiencyCells& other); // assignment operator

  virtual ~AliMUONTriggerEfficiencyCells();

  enum {
    kBendingEff,     ///< Bending plane fired
    kNonBendingEff,  ///< Non-bending plane fired
    kBothPlanesEff,  ///< Both planes fired
    kAllTracks,      ///< tracks used for calculation
    kNcounts         ///< Number of count type
  };

  enum {
    kHboardCount,   ///< Counts per board index 
    kHslatCount,    ///< Counts per slat index
    kHchamberCount  ///< Counts per chamber index
  };

  const Char_t* GetHistoName(Int_t histoType, Int_t countType, 
			     Int_t chamber = -1);

  /// Get list of histograms
  TList* GetHistoList() { return fCountHistoList; }

  TH1F* GetOldEffHisto(Int_t hType, Int_t ich, Int_t icath); // obsolete

protected:
    void ResetHistos(Bool_t deleteObjects = kFALSE);

    void ReadFile(const Char_t* filename, 
		  const Char_t* listname);

private:
    void CheckConstants() const;

    static const Int_t fgkNcathodes=2; ///<Number of cathodes
    static const Int_t fgkNchambers=4; ///<Number of chambers
    static const Int_t fgkNplanes=8;   ///<Number of planes
    
    TH1F *fBoardEfficiency[fgkNplanes];///< the boards content (obsolete)
    TH1F *fSlatEfficiency[fgkNplanes];///< the slats content (obsolete)

    TList *fCountHistoList; ///< list of histograms for efficiency calculation
    TList *fNoCountHistoList; ///<list of efficiency denominators (obsolete)
    TList *fFiredStrips; ///<list of fired strips for efficiency check (obsolete)

    ClassDef(AliMUONTriggerEfficiencyCells,6) // Trigger efficiency store
};
#endif
