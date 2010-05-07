#ifndef ALIMUONTRIGGERDISPLAY_H
#define ALIMUONTRIGGERDISPLAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONTriggerDisplay
/// \brief Converts histograms as a function of strip/board/slat number in human readable histograms
///
//  Author Diego Stocco

// --- ROOT system ---

#include "TH2.h"

class TArrayD;
class TString;
class TGraph;

class AliMUONTriggerDisplay: public TObject {

public:
  AliMUONTriggerDisplay();
  virtual ~AliMUONTriggerDisplay();
  
  /// Display element inidices (strip,board,slat)
  enum EDisplayType {
    kDisplayStrips,    ///< Draw strips
    kDisplayBoards,    ///< Draw boards
    kDisplaySlats      ///< Draw slats
  };

  /// Display options inidices
  enum EDisplayOption {
    kDefaultDisplay,   ///< Default display
    kNumbered,         ///< Histogram filled with board numbers
    kShowZeroes,       ///< Displays strip/board/slat content even if it is 0
    kNormalizeToArea   ///< Draw input histo divided by element area
  };

  TH2* GetEmptyDisplayHisto(TString displayHistoName, EDisplayType displayType,
			    Int_t cathode, Int_t chamber=11,
			    TString displayHistoTitle="");

  TH2* GetBoardNumberHisto(TString displayHistoName,
			   Int_t chamber=11,TString displayHistoTitle="");
  
  TH2* GetDisplayHistogram(TH1* inputHisto, TString displayHistoName,
			   EDisplayType displayType, Int_t cathode,
			   Int_t chamber=11, TString displayHistoTitle="",
			   EDisplayOption displayOpt=kDefaultDisplay);

  TH2* GetDisplayHistogram(TGraph* inputGraph, TString displayHistoName,
			   EDisplayType displayType, Int_t cathode,
			   Int_t chamber=11, TString displayHistoTitle="",
			   EDisplayOption displayOpt=kDefaultDisplay);
  
  Bool_t FillDisplayHistogram(TH1* inputHisto, TH2* displayHisto,
			      EDisplayType displayType, Int_t cathode,
			      Int_t chamber=11,EDisplayOption displayOpt=kDefaultDisplay);
  
  Bool_t FillDisplayHistogram(TGraph* inputGraph, TH2* displayHisto,
			      EDisplayType displayType, Int_t cathode,
			      Int_t chamber=11,EDisplayOption displayOpt=kDefaultDisplay);
  
private:
  Bool_t AddSortedPoint(Float_t currVal, TArrayD& position, const Float_t kResetValue);
  /// Return index
  Int_t GetIndex(Int_t chamber, Int_t cathode) { return 2*chamber + cathode;}

  Bool_t InitOrDisplayTriggerInfo(TObject* inputHisto, TH2* displayHisto,
				  EDisplayType displayType,
				  Int_t cathode, Int_t chamber,
				  TString displayHistoName, TString displayHistoTitle,
				  EDisplayOption displayOpt=kDefaultDisplay);

  void FillBins(TObject* inputHisto, TH2* displayHisto,
		Int_t iElement1, Int_t iElement2,
		Float_t x1, Float_t x2, Float_t y1, Float_t y2,
		const Float_t kShiftX, const Float_t kShiftY,
		EDisplayOption displayOpt);

  ClassDef(AliMUONTriggerDisplay,0)  // Class for converting trigger histograms

};
#endif
