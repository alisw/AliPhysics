#ifndef ALIMUONPAINTERMATRIXFRAME_H
#define ALIMUONPAINTERMATRIXFRAME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterMatrixFrame
/// \brief Widget to plot a matrix of painters
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include "TGFrame.h"
#endif
#include <float.h>

class AliMUONPainterColorSlider;
class AliMUONPainterPlotSelector;
class AliMUONPainterMatrix;
class AliMUONPainterHighlighter;
class AliMUONVPainter;
class AliMUONVTrackerData;
class TGButtonGroup;
class TGToolTip;
class TList;
class TObjArray;
class TRootEmbeddedCanvas;

class AliMUONPainterMatrixFrame : public TGCompositeFrame
{
public:
  AliMUONPainterMatrixFrame(const TGWindow* window, UInt_t w, UInt_t h);
  virtual ~AliMUONPainterMatrixFrame();
  
  void Clear(Option_t* opt="");

  void DataRangeAutoRequested();
  
  void DataRangeWasChanged(Double_t* range);
    
  void DataSourceWasChanged(const char* name, AliMUONVTrackerData* data, Int_t dataIndex);

  void EventInfo(Int_t event, Int_t px, Int_t py, TObject* selected);

  void MouseEnter(AliMUONVPainter* painter); // *SIGNAL*

  void MouseMotion(AliMUONVPainter* painter, Double_t* position); // *SIGNAL*

  void MouseLeave(const AliMUONVPainter* painter); // *SIGNAL*
  
  void ResponderButtonWasClicked(Int_t id); 

  void OutlineButtonWasClicked(Int_t id); 
  
  void Use(AliMUONPainterMatrix* group);
  
  void TitleHasChanged(const char* newTitle); // *SIGNAL*
  
  void Update(); 

  /// Get the matrix pointer
  AliMUONPainterMatrix* Matrix() const { return fPainterMatrix; }
  
  void SaveAs(const char* filename="", Option_t* option="") const;

  void UpdateInterface(Bool_t fromScratch);
  
private:
  /// not implemented
  AliMUONPainterMatrixFrame(const AliMUONPainterMatrixFrame& rhs);
  /// not implemented
  AliMUONPainterMatrixFrame& operator=(const AliMUONPainterMatrixFrame& rhs);

  void ChangeTitle(const TString& title);
  
  void ChangeTitle(AliMUONVPainter* painter, const char* basename=0x0,
                   Double_t x=FLT_MAX, Double_t y=FLT_MAX);
  
  void CreateButtons();
  
  void UpdateDataRange();
  
  void ViewModified();
  
private:
  AliMUONPainterMatrix* fPainterMatrix; ///< the matrix we plot (not owner)
  TRootEmbeddedCanvas* fView; ///< the canvas used to plot
  TGHorizontalFrame* fInterface;  ///< the interface frame
  TGButtonGroup* fResponderButtons; ///< the responder buttons
  TGButtonGroup* fOutlineButtons; ///< the outline buttons
  
  AliMUONPainterPlotSelector* fPlotSelector; ///< the data source selection
    
  AliMUONPainterHighlighter* fPainterHighlighter; ///< the highlighter
  
  UInt_t fCanvasWidth; ///< canvas width
  UInt_t fCanvasHeight; ///< canvas height
  
  TGCompositeFrame* fMainFrame; ///< our main frame
  
  AliMUONPainterColorSlider* fColorSlider; ///< color slider (for data)
  
  ClassDef(AliMUONPainterMatrixFrame,1) // Widget for drawing painter matrix
};

#endif
