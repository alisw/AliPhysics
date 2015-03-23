#ifndef ALIMUONTRACKERDATACOMPAREDIALOG_H
#define ALIMUONTRACKERDATACOMPAREDIALOG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONTrackerDataCompareDialog
/// \brief
/// 
/// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include "TGFrame.h"
#endif

#include "AliMUONTrackerData.h"


class TGButtonGroup;
class TGComboBox;
class TGCompositeFrame;
class TGTextEntry;

class AliMUONTrackerDataCompareDialog : public TGTransientFrame
{
public:
  AliMUONTrackerDataCompareDialog(const TGWindow* p = 0, const TGWindow* main = 0, UInt_t w = 1, UInt_t h = 1);
  virtual ~AliMUONTrackerDataCompareDialog();

  void DoOK();
  void DoCancel();
  
private:
  
    /// not defined
    AliMUONTrackerDataCompareDialog(const AliMUONTrackerDataCompareDialog& rhs);
  /// not defined
  AliMUONTrackerDataCompareDialog& operator=(const AliMUONTrackerDataCompareDialog& rhs);

  void CompareData(const char* d1name, const char* d2name, AliMUONTrackerData::EDiffType difftype) const;
  
private:

    TGCompositeFrame* fF1; ///< frame for data source 1 selection
  TGComboBox* fData1; ///< to select first data
  TGCompositeFrame* fF2; ///< frame for data source 2 selection
  TGComboBox* fData2; ///< to select second data
  TGCompositeFrame* fF3; ///< frame for difference type selection
  TGComboBox* fDiffType; ///< to select the kind of difference to make
  TGCompositeFrame* fF4; ///< frame for output basename selection
  TGTextEntry* fBasename; ///< basename of resulting (diff-ed) data
  TGCompositeFrame* fButtonFrame; ///< to hold OK and Cancel buttons
  TGTextButton* fOK; ///< ok button
  TGTextButton* fCancel; ///< cancel button
  
  ClassDef(AliMUONTrackerDataCompareDialog,2) // Dialog to select two data sources to compare
};

#endif
