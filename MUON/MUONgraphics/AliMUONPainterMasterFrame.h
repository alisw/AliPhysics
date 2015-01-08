#ifndef ALIMUONPAINTERMASTERFRAME_H
#define ALIMUONPAINTERMASTERFRAME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterMasterFrame
/// \brief The main window for the offline "a la mood" display
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include <TGFrame.h>
#endif
#ifndef ROOT_TArrayI
#  include <TArrayI.h>
#endif

class AliMUONAttPainter;
class AliMUONPainterMatrix;
class AliMUONPainterMatrixFrame;
class AliMUONAttPainterSelectorFrame;
class AliMUONVPainter;
class TGButton;
class TGButtonGroup;
class TGComboBox;
class TGLabel;
class TObjArray;

class AliMUONPainterMasterFrame : public TGCompositeFrame
{
public:
  AliMUONPainterMasterFrame(const TGWindow* p, UInt_t w, UInt_t h, AliMUONPainterMatrix* matrix);
  virtual ~AliMUONPainterMasterFrame();

  void Backward();

  void ChangeTitle(const char* newTitle);
  
  void Clicked(AliMUONVPainter* painter, Double_t* values);

  void ShiftClicked(AliMUONVPainter* painter, Double_t* values);

  void Forward();
  
  void PainterMatrixWantToShow(AliMUONPainterMatrix* group);

  void Update();

  void AttributesChanged(const AliMUONAttPainter* newValues);
  
  void SaveAs(const char* filename = "", Option_t* option = "") const;
  
  void PrintAs() const;

  void PrintMe() const;
  
private:
  /// not implemented
  AliMUONPainterMasterFrame(const AliMUONPainterMasterFrame& rhs);
  /// not implemented
  AliMUONPainterMasterFrame& operator=(const AliMUONPainterMasterFrame& rhs);
  
  void AddPainterMatrix(AliMUONPainterMatrix* group);
  void MakeTopPainterMatrix(UInt_t w, UInt_t h, AliMUONPainterMatrix* matrix);
  void SetNavigation(Int_t i);
  void ShowPainterMatrix(AliMUONPainterMatrix* group);  
  void UpdateNavigation();
  void UpdateAttributes(const AliMUONPainterMatrix& painterMatrix);
  
private:
  TGHorizontalFrame* fNavigationFrame; ///< top frame for navigation
  AliMUONPainterMatrixFrame* fPainterMatrixFrame; ///< main frame with painters
  
  TGButton* fBackButton; ///< navigation back 
  TGButton* fForwardButton; ///< navigation forward
  TGLabel* fGroupTitle; ///< top title
  TGButton* fPrintMeButton; ///< print button
  TGButton* fPrintAsButton; ///< print... button
  
  TArrayI fNavigation; ///< navigation "history"
    
  Int_t fCurrentNavigationPosition; ///< current position in navigation history

  AliMUONAttPainterSelectorFrame* fAttPainterSelectorFrame; ///< view type selection frame
  
  static const Int_t fgkBorderSize; ///< border sizes to use when placing frames

  ClassDef(AliMUONPainterMasterFrame,0) // Main window of display
};

#endif
