/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONPainterMatrixFrame.h"

#include "AliMUONPainterColorSlider.h"
#include "AliMUONPainterMatrix.h"
#include "AliMUONPainterGroup.h"
#include "AliMUONPainterHighlighter.h"
#include "AliMUONPainterInterfaceHelper.h"
#include "AliMUONPainterPlotSelector.h"
#include "AliMUONPainterRegistry.h"
#include "AliMUONVTrackerData.h"
#include "AliMUONVPainter.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGFrame.h>
#include <TGListBox.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRootEmbeddedCanvas.h>
#include <TString.h>
#include <cassert>
#include <float.h>

/// \class AliMUONPainterMatrixFrame
///
/// A widget to draw a painter matrix, and the corresponding interface
/// to select what to outline or paint, and which part of the painter
/// is responding to mouse events
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterMatrixFrame)
///\endcond

//_____________________________________________________________________________
AliMUONPainterMatrixFrame::AliMUONPainterMatrixFrame(const TGWindow* window, 
                                                   UInt_t w, UInt_t h)
: TGCompositeFrame(window,w,h,kVerticalFrame|kDoubleBorder),
  fPainterMatrix(0x0),
  fView(0x0),
  fInterface(0x0),
  fResponderButtons(0x0),
  fOutlineButtons(0x0),
  fPlotSelector(0x0),
  fPainterHighlighter(new AliMUONPainterHighlighter),
  fCanvasWidth(0),
  fCanvasHeight(0),
  fMainFrame(0x0),
  fColorSlider(0x0)
{
    /// ctor
    const Int_t kBorderSize = 10;

    UInt_t wi = w - kBorderSize*0;
    UInt_t hi = h - kBorderSize*1;
    
    fCanvasWidth = wi;
    fCanvasHeight = (UInt_t)(hi*0.75);
    
    fMainFrame = new TGHorizontalFrame(this,fCanvasWidth,hi);
    
    const Int_t kColorWidth = 100;
    
    fColorSlider = new AliMUONPainterColorSlider(fMainFrame,kColorWidth,fCanvasHeight);
    
    fView = new TRootEmbeddedCanvas("ec",fMainFrame,fCanvasWidth-kColorWidth,fCanvasHeight,kChildFrame);
    
    fInterface = new TGHorizontalFrame(this,fCanvasWidth);
    
    fMainFrame->AddFrame(fView, new TGLayoutHints(kLHintsExpandX));
    fMainFrame->AddFrame(fColorSlider,new TGLayoutHints(kLHintsTop|kLHintsRight|kLHintsCenterY,kBorderSize/2));

    AliMUONPainterInterfaceHelper::SetBackgroundColor("MatrixFrame.ColorSlider",*fColorSlider);
    
    fResponderButtons = new TGButtonGroup(fInterface,"Responder");
    
    fOutlineButtons = new TGButtonGroup(fInterface,"Outline");
    
    fInterface->AddFrame(fResponderButtons);
    fInterface->AddFrame(fOutlineButtons);

    fPlotSelector = 
      new AliMUONPainterPlotSelector(fInterface);//,wi,interfaceHeight);
    
    fInterface->AddFrame(fPlotSelector);//,new TGLayoutHints(kLHintsRight|kLHintsExpandX));

    fOutlineButtons->Show();
    fResponderButtons->Show();
        
    AddFrame(fMainFrame,new TGLayoutHints(kLHintsExpandX|kLHintsTop,
                                         0,0,0,0));

    AddFrame(fInterface,new TGLayoutHints(kLHintsExpandX|kLHintsBottom,
             0,0,kBorderSize,0));


    // Set the connections
    
    fPlotSelector->Connect("DataSourceWasChanged(const char*,AliMUONVTrackerData*,Int_t)",                                
                                 "AliMUONPainterMatrixFrame",
                                 this,
                                 "DataSourceWasChanged(const char*,AliMUONVTrackerData*,Int_t)");
    
    fColorSlider->Connect("DataRangeWasChanged(Double_t*)",
                          "AliMUONPainterMatrixFrame",
                         this,
                         "DataRangeWasChanged(Double_t*)");

    fColorSlider->Connect("DataRangeAutoRequested()",
                          "AliMUONPainterMatrixFrame",
                         this,
                         "DataRangeAutoRequested()");
    
    // Set the colors (mainly for debugging frame layout)

    AliMUONPainterInterfaceHelper::SetBackgroundColor("MatrixFrame.Main",*this);
    
    fMainFrame->HideFrame(fColorSlider);
    
    fMainFrame->Resize();
}

//_____________________________________________________________________________
AliMUONPainterMatrixFrame::~AliMUONPainterMatrixFrame()
{
  /// dtor
  delete fPainterHighlighter;
  AliError("Please write a decent dtor for this class !");
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::ChangeTitle(const TString& title)
{
  /// Change title
  
  TitleHasChanged(title.Data());
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::ChangeTitle(AliMUONVPainter* painter, 
                                       const char* basename,
                                       Double_t x, Double_t y)
{
  /// Change the title according to painter
  
  TString name;
  
  if (painter) 
  {

    if ( basename ) name = basename;
    else name = painter->PathName();
    
    AliMUONVPainter* master = painter->Master();
    
    AliMUONPainterGroup* group = master->PlotterGroup();

    AliDebug(1,Form("Painter is %s plotterGroup is %p %s",
                    painter->PathName().Data(),
                    group,
                    ( group ? group->Type() : "")));
    
    
    if ( group && group->Data() ) 
    {
      name += "\n";
      name += painter->Describe(*(group->Data()),group->DataIndex(),x,y);
    }
  }
  else
  {
    name = fPainterMatrix->Name();
  }
  
  TitleHasChanged(name.Data());
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::Clear(Option_t*)
{
  /// Clear the view(s)
  
  fPainterMatrix = 0x0;

  AliMUONPainterInterfaceHelper::ClearButtons(*fOutlineButtons);
  AliMUONPainterInterfaceHelper::ClearButtons(*fResponderButtons);
  
  fView->GetCanvas()->SetEditable(kTRUE);
  fView->GetCanvas()->Clear();
  fView->GetCanvas()->Modified();
  fView->GetCanvas()->Update();
  fView->GetCanvas()->SetEditable(kFALSE);

  Layout();  
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrixFrame::CreateButtons()
{
  /// Create the interface buttons
  
  AliDebug(1,"");
  
//  AliMUONVPainter* painter = fPainterMatrix->Painter(0);
  
  /// create buttons    
  TObjArray types;
  
  fPainterMatrix->GetTypes(types);
  
  TIter nextType(&types);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(nextType()) ) )
  {
    AliMUONPainterInterfaceHelper::AddRadioButton(*fResponderButtons,str->String());
    AliMUONPainterInterfaceHelper::AddCheckButton(*fOutlineButtons,str->String());
  }
  
  fOutlineButtons->Connect("Clicked(Int_t)","AliMUONPainterMatrixFrame",
                        this,"OutlineButtonWasClicked(Int_t)");
  
  fResponderButtons->Connect("Clicked(Int_t)","AliMUONPainterMatrixFrame",
                             this,"ResponderButtonWasClicked(Int_t)");        
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrixFrame::DataRangeAutoRequested()
{
  /// Get there when the "Auto" button below the color slider is clicked,
  /// to compute the data range actually painted.
  
  Double_t dataMin, dataMax;

  AliDebug(1,"");
  
  fPainterMatrix->ComputeDataRange();

  fPainterMatrix->GetDataRange(dataMin,dataMax);
  
  AliDebug(1,Form("dataMin,Max for SetRange=%e,%e",dataMin,dataMax));
  
  Bool_t emit(kTRUE);
  
  fColorSlider->SetRange(dataMin,dataMax,emit);
  
  Update();
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrixFrame::DataRangeWasChanged(Double_t* range)
{
  /// Get there when the data range is changed
  
  AliDebug(1,Form("range=%e,%e",range[0],range[1]));

  fPainterMatrix->SetDataRange(range[0],range[1]);
  
  if ( !fColorSlider->IsLocked() )
  {
    Update();
  }
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::DataSourceWasChanged(const char* type,
                                               AliMUONVTrackerData* data,
                                               Int_t indexInData)
{
  /// Update what to plot
 
  TString pattern(type);

  AliDebug(1,Form("type=%s data=%s index=%d",type,
                  (data ? data->GetName() : "null"),indexInData));
    
  AliMUONVTrackerData* d = data;
  
  if ( !d || !data || indexInData < 0 || pattern == "" )
  {
    pattern = "*";
    d = 0;
    indexInData = -1;
  }
  
  fPainterMatrix->SetData(pattern,d,indexInData);
    
  Update();
  
  ChangeTitle(fPainterMatrix->Painter(0));
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::EventInfo(Int_t event, Int_t px ,Int_t py, TObject* object)
{
  /// Used to detect entering/leaving a given painter 
  
  if (!gPad || !object) return;
  
//  cout << "EventInfo : event " << event << " px " << px << " py " << py
//    << " object " << object << " " << object->GetName() << endl;
// 
  if ( event == 7 ) 
  {
    if ( object->InheritsFrom("AliMUONVPainter") )
    {
      AliMUONVPainter* p = static_cast<AliMUONVPainter*>(object);
      p->ExecuteEvent(7,px,py);
      return;
    }      
  }
  
  if ( event == kMouseLeave )
  {
    if ( object->InheritsFrom("AliMUONVPainter") )
    {
      AliMUONVPainter* p = static_cast<AliMUONVPainter*>(object);
      MouseLeave(p);
      fPainterHighlighter->SetPainter(0x0);
      gPad->Modified();
      gPad->Update();
    }
  }

  if ( event == kMouseEnter )
  {
    if ( object->InheritsFrom("AliMUONVPainter") )
    {
      AliMUONVPainter* painter = static_cast<AliMUONVPainter*>(object);
      if ( painter->IsResponder() && !painter->HandleMouseMotion() )
      {
        MouseEnter(static_cast<AliMUONVPainter*>(object));
        fPainterHighlighter->SetPainter(painter);
        gPad->Modified();
        gPad->Update();
      }
      else if ( !painter->HandleMouseMotion() )
      {
        MouseEnter(static_cast<AliMUONVPainter*>(object)); 
      }
    }
  }
  
   if ( event == kMouseMotion ) 
  {
    if ( object->InheritsFrom("AliMUONVPainter") )
    {
      AliMUONVPainter* painter = static_cast<AliMUONVPainter*>(object);

      if ( painter->HandleMouseMotion() && painter->IsResponder() )
      {
        Double_t pos[2];
        TVirtualPad* padsave = gPad;
        painter->Pad()->cd();
        painter->PixelToPad(px,py,pos[0],pos[1]);
        MouseMotion(static_cast<AliMUONVPainter*>(object),pos);
        fPainterHighlighter->SetPainter(painter,pos[0],pos[1]);
        gPad->Modified();
        gPad->Update();
        gPad = padsave;
      }
    }    
  }
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::MouseEnter(AliMUONVPainter* painter)
{
  /// Emit a signal to notify that mouse pointer is entering a given painter

  AliDebug(1,Form("painter=%p %s",painter,painter->PathName().Data()));
  
  ChangeTitle(painter);

  Long_t params[] = { (Long_t)painter };
  
  Emit("MouseEnter(AliMUONVPainter*)",params);  
}


//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::MouseLeave(const AliMUONVPainter* painter)
{
  /// Emit a signal to notify that mouse pointer is leaving a given painter
  
  ChangeTitle(fPainterMatrix->Name());

  Long_t params[] = { (Long_t)painter };
  
  Emit("MouseLeave(AliMUONVPainter*)",params);
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::MouseMotion(AliMUONVPainter* painter, Double_t* position)
{
  /// Emit a signal to notify that mouse pointer is moving within a given painter
  
  ChangeTitle(painter,painter->NameAtPosition(position[0],position[1]),
              position[0],position[1]);
  
  Long_t params[] = { (Long_t)painter, (Long_t)position };
  
  Emit("MouseMotion(AliMUONVPainter*,Double_t*)",params);
}


//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::ResponderButtonWasClicked(Int_t id)
{
  /// One responder button was clicked
  
  TGTextButton* button = static_cast<TGTextButton*>(fResponderButtons->GetButton(id));
  TString pattern = button->GetString();
  
//  AliInfo(Form("id=%d button=%d %s",id,button->IsOn(),pattern.Data()));

  assert(button->IsOn()==1);
  
  fPainterMatrix->SetResponder(pattern.Data());
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::OutlineButtonWasClicked(Int_t id)
{
  /// One outline button was clicked
  
  TGTextButton* button = static_cast<TGTextButton*>(fOutlineButtons->GetButton(id));
  TString pattern = button->GetString();
  
  fPainterMatrix->SetOutlined(pattern.Data(),button->IsOn());
  
  ViewModified();
  fView->GetCanvas()->Update();
    
  // Update the interface (e.g. list of possible responders can have 
  // changed due to visibility change)
  UpdateInterface(kFALSE);
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrixFrame::SaveAs(const char* filename, Option_t* option) const
{
  /// Save painter matrix (in the sense of "print") in filename
  
  TCanvas* d = fPainterMatrix->CreateCanvas();
  
  d->SaveAs(filename,option);
  
  delete d;
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::TitleHasChanged(const char* title)
{
  /// Emit the TitleHasChanged signal
  
  Long_t params[] = { (Long_t)title };
  Emit("TitleHasChanged(const char*)",params);
}


//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::Update()
{
  /// Force update of all canvases

  UpdateDataRange();
  
  fView->GetCanvas()->SetEditable(kTRUE);
  
  Bool_t colorSlider = ( fPainterMatrix->Data() != 0x0 );
  
  ViewModified();

  fView->GetCanvas()->SetEditable(kFALSE);
  
  AliDebug(1,Form("colorSlider=%d",colorSlider));
  
  if ( colorSlider )
  {
    fMainFrame->ShowFrame(fColorSlider);
  }
  else
  {
    fMainFrame->HideFrame(fColorSlider);
  }
  
  fMainFrame->Layout();
  Layout();
}

//_____________________________________________________________________________
void
AliMUONPainterMatrixFrame::UpdateDataRange()
{
  /// Update the data range

  if ( fColorSlider->IsLocked() ) 
  {
    fColorSlider->SetRange(0,0,kTRUE);
    return;
  }
  
  Double_t min, max;

  fPainterMatrix->GetDataRange(min,max);

  AliDebug(1,Form("min %e max %e",min,max));

  if ( min > max ) 
  {
    fPainterMatrix->ComputeDataRange();
    fPainterMatrix->GetDataRange(min,max);
  }

  fColorSlider->SetRange(min,max,kFALSE);
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrixFrame::UpdateInterface(Bool_t fromScratch)
{
  /// Update the full interface
  
  if ( fromScratch || fOutlineButtons->GetCount() == 0 ) 
  {
    CreateButtons();
  }
  
  AliMUONPainterInterfaceHelper::Unselect(*fResponderButtons,"*");
  AliMUONPainterInterfaceHelper::Unselect(*fOutlineButtons,"*");
  
  AliMUONVPainter* painter = fPainterMatrix->Painter(0);
  
  TObjArray types;
  types.SetOwner(kTRUE);
  
  fPainterMatrix->GetTypes(types);
  
  // update button states
  TIter next(&types);
  TObjString* otype;
  
  TString theResponder;
  
  while ( ( otype = static_cast<TObjString*>(next()) ) )
  {    
    AliMUONPainterGroup* group = painter->Group(otype->String());
    
    if ( group && group->IsOutlined() ) 
    {
      AliMUONPainterInterfaceHelper::Select(*fOutlineButtons,otype->String().Data());
    }    
  }
  
  if ( painter ) 
  {
    AliMUONPainterGroup* responderGroup = painter->ResponderGroup();
  
    if (responderGroup)
    {
      AliMUONPainterInterfaceHelper::Select(*fResponderButtons,responderGroup->Type());
    }
  }
  
  // update data source view
  
  fPlotSelector->Update(*fPainterMatrix);
    
  fResponderButtons->Show();
  fOutlineButtons->Show();

  Layout();
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrixFrame::Use(AliMUONPainterMatrix* group)
{
  /// Change the matrix used
  
  Clear();
  
  fPainterMatrix = group;
  
  fView->GetCanvas()->SetEditable(kTRUE);
  
  fView->GetCanvas()->Divide(fPainterMatrix->Nx(),fPainterMatrix->Ny());
  
  for ( Int_t i = 0; i < fPainterMatrix->Size(); ++i ) 
  {
    AliMUONVPainter* painter = fPainterMatrix->Painter(i);
    fView->GetCanvas()->cd(i+1);
    painter->Draw("R");
    fPainterHighlighter->SetPainter(0x0);
    fPainterHighlighter->Draw();
  }  

  Update();
  
  UpdateInterface(kTRUE);
  
  ChangeTitle(fPainterMatrix->Name());
  
  fView->GetCanvas()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
                              "AliMUONPainterMatrixFrame",this, 
                              "EventInfo(Int_t,Int_t,Int_t,TObject*)"); 
  
}

//_____________________________________________________________________________
void 
AliMUONPainterMatrixFrame::ViewModified()
{
  /// Update our canvas
  
  for ( Int_t i = 0; i < fPainterMatrix->Size(); ++i ) 
  {
    fView->GetCanvas()->GetPad(i+1)->Modified();
  }
  fView->GetCanvas()->Modified();
  fView->GetCanvas()->Update();
}


