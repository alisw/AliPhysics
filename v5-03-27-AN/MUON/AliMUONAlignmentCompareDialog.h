#ifndef ALIMUONALIGNMENTCOMPAREDIALOG_H
#define ALIMUONALIGNMENTCOMPAREDIALOG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONAlignmentCompareDialog
/// \brief
/// 
/// Authors Philippe Pillot, Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include "TGFrame.h"
#endif

class AliMUONVTrackerData;
class TGButtonGroup;
class TGNumberEntry;
class TGCompositeFrame;
class TGTextEntry;

class AliMUONAlignmentCompareDialog : public TGTransientFrame
{
public:
  AliMUONAlignmentCompareDialog(const TGWindow* p = 0, const TGWindow* main = 0, UInt_t w = 1, UInt_t h = 1);
  virtual ~AliMUONAlignmentCompareDialog();

  void DoOK();
  void DoCancel();
  
private:
  
    /// not defined
    AliMUONAlignmentCompareDialog(const AliMUONAlignmentCompareDialog& rhs);
    /// not defined
    AliMUONAlignmentCompareDialog& operator=(const AliMUONAlignmentCompareDialog& rhs);

    AliMUONVTrackerData* CompareAlignment(const char* ocdbPathForAlign1, Int_t run1,
                                          const char* ocdbPathForAlign2, Int_t run2);
  
    void AddInput(TGCompositeFrame* frame, const char* msg, TGTextEntry*& text, TGNumberEntry*&     run);

private:

  TGCompositeFrame* fF1; ///< frame for align 1 selection
  TGTextEntry* fOCDBPath1; ///< to select first alignment path
  TGNumberEntry* fRun1; ///< to select first run
  TGCompositeFrame* fF2; ///< frame for align 2 selection
  TGTextEntry* fOCDBPath2; ///< to select second alignment path
  TGNumberEntry* fRun2; ///< to select second run
  TGCompositeFrame* fF3; ///< frame for difference type selection
  TGTextEntry* fBasename; ///< basename of resulting (diff-ed) data
  TGCompositeFrame* fButtonFrame; ///< to hold OK and Cancel buttons
  TGTextButton* fOK; ///< ok button
  TGTextButton* fCancel; ///< cancel button
    
  ClassDef(AliMUONAlignmentCompareDialog,1) // Dialog to select two data sources to compare
};

#endif
