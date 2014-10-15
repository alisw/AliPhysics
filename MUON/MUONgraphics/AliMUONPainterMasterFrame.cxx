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

#include "AliMUONPainterMasterFrame.h"

#include "AliMUONChamberPainter.h"
#include "AliMUONPainterGroup.h"
#include "AliMUONPainterMatrix.h"
#include "AliMUONPainterMatrixFrame.h"
#include "AliMUONPainterInterfaceHelper.h"
#include "AliMUONPainterRegistry.h"
#include "AliMUONAttPainterSelectorFrame.h"
#include "AliMUONVPainter.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TGComboBox.h>
#include <TGFileDialog.h>
#include <TGLabel.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGButtonGroup.h>
#include <TGMsgBox.h>
#include <TSystem.h>

/// \class AliMUONPainterMasterFrame
///
/// Main window of the 2D display
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterMasterFrame)
///\endcond

const Int_t AliMUONPainterMasterFrame::fgkBorderSize = 10;

//_____________________________________________________________________________
AliMUONPainterMasterFrame::AliMUONPainterMasterFrame(const TGWindow* p, 
                                                     UInt_t w, UInt_t h, AliMUONPainterMatrix* matrix)
: TGCompositeFrame(p,w,h,kVerticalFrame),
fNavigationFrame(0x0),
fPainterMatrixFrame(0x0),
fBackButton(0x0),
fForwardButton(0x0),
fGroupTitle(0x0),
fPrintMeButton(0x0),
fPrintAsButton(0x0),
fNavigation(),
fCurrentNavigationPosition(-1),
fAttPainterSelectorFrame(0x0)
{  
  /// ctor
    
  UInt_t wi = w - fgkBorderSize*2;
  UInt_t hi = h - fgkBorderSize*3;
  
  fNavigationFrame = new TGHorizontalFrame(this,wi);
  
  AddFrame(fNavigationFrame,new TGLayoutHints(kLHintsExpandX|kLHintsTop,
                                              fgkBorderSize,fgkBorderSize,
                                              fgkBorderSize,fgkBorderSize));
  
  fBackButton = new TGPictureButton(fNavigationFrame,
                                    gClient->GetPicture("tb_back.xpm"));
  
  fForwardButton = new TGPictureButton(fNavigationFrame,
                                       gClient->GetPicture("tb_forw.xpm"));    

  fPrintMeButton = new TGTextButton(fNavigationFrame,"Print");
  fPrintAsButton = new TGTextButton(fNavigationFrame,"Print As...");

  fAttPainterSelectorFrame = new AliMUONAttPainterSelectorFrame(fNavigationFrame,w/2,20);
  
  fGroupTitle = new TGLabel(fNavigationFrame,"");
  
  fNavigationFrame->AddFrame(fBackButton,new TGLayoutHints(kLHintsCenterY));
  fNavigationFrame->AddFrame(fForwardButton,new TGLayoutHints(kLHintsCenterY));
  
  fNavigationFrame->AddFrame(fAttPainterSelectorFrame,new TGLayoutHints(kLHintsCenterY,10));

  fNavigationFrame->AddFrame(fPrintMeButton,new TGLayoutHints(kLHintsCenterY,10));
  fNavigationFrame->AddFrame(fPrintAsButton,new TGLayoutHints(kLHintsCenterY,10));

  fAttPainterSelectorFrame->Connect("Clicked(AliMUONAttPainter*)",
                                    "AliMUONPainterMasterFrame",
                                    this,
                                    "AttributesChanged(AliMUONAttPainter*)");
  
  fNavigationFrame->AddFrame(fGroupTitle,new TGLayoutHints(kLHintsExpandX|kLHintsCenterX|kLHintsCenterY,10));
  
  fForwardButton->Connect("Clicked()","AliMUONPainterMasterFrame",
                          this,
                          "Forward()");

  fBackButton->Connect("Clicked()","AliMUONPainterMasterFrame",
                          this,
                          "Backward()");
    
  fPrintMeButton->Connect("Clicked()","AliMUONPainterMasterFrame",
                        this,
                        "PrintMe()");

  fPrintAsButton->Connect("Clicked()","AliMUONPainterMasterFrame",
                        this,
                        "PrintAs()");
  
  UInt_t w1 = wi;
  //  UInt_t h1 = hi - fNavigationFrame->GetHeight() - 3*fgkBorderSize;
  UInt_t h1 = hi - 7*12;
  
  MakeTopPainterMatrix(w1,h1,matrix);

  AddFrame(fPainterMatrixFrame,new TGLayoutHints(kLHintsExpandX,
                                                fgkBorderSize,fgkBorderSize,
                                                0,fgkBorderSize));
  
  AliMUONPainterInterfaceHelper::SetBackgroundColor("MasterFrame.Navigation",*fNavigationFrame);
  AliMUONPainterInterfaceHelper::SetBackgroundColor("MasterFrame.Main",*this);
  
  AliDebug(1,Form("fNavigation=%p",&fNavigation));
  
  AliMUONPainterRegistry::Instance()->Connect("PainterMatrixWantToShow(AliMUONPainterMatrix*)",
                                              "AliMUONPainterMasterFrame",
                                              this,
                                              "PainterMatrixWantToShow(AliMUONPainterMatrix*)");
  
  fPainterMatrixFrame->DataSourceWasChanged(matrix->DataPattern().Data(),matrix->Data(),matrix->DataIndex());
}

//_____________________________________________________________________________
AliMUONPainterMasterFrame::~AliMUONPainterMasterFrame()
{
  /// dtor
  Cleanup();
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::AddPainterMatrix(AliMUONPainterMatrix* painterMatrix)
{
  /// array is adopted (by the registry)

  AliDebug(1,Form("matrix=%p %s",painterMatrix,painterMatrix->GetName()));
  
  Int_t i = AliMUONPainterRegistry::Instance()->Register(painterMatrix);

  SetNavigation(i);
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::PainterMatrixWantToShow(AliMUONPainterMatrix* group)
{
  /// FIXME: should check whether we are the active window before
  /// responding to this message ?

  AliDebug(1,Form("group=%p %s",group,group->GetName()));
  
  Int_t i = AliMUONPainterRegistry::Instance()->FindIndexOf(group);

  Int_t alreadyThere(-1);
  
  for ( Int_t j = 0; j < fNavigation.GetSize(); ++j )
  {
    if ( fNavigation[j] == i ) alreadyThere = j;
  }
  
  if (alreadyThere<0) 
  {
    SetNavigation(i);
  }
  else
  {
    fCurrentNavigationPosition = alreadyThere;
  }
  
  ShowPainterMatrix(group);  
}
     
//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::PrintAs() const
{
  /// Handle the PrintAs button
  
  TGFileInfo fileInfo;
  
  new TGFileDialog(gClient->GetRoot(),gClient->GetRoot(),
                   kFDSave,&fileInfo);
  
  if ( fileInfo.fFilename ) 
  {
    SaveAs(gSystem->ExpandPathName(Form("%s",fileInfo.fFilename)),"RECREATE");
  }
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::PrintMe() const
{
  /// Handle the PrintMe button
  
  SaveAs(gSystem->ExpandPathName(Form("%s.png",fPainterMatrixFrame->Matrix()->GetName())),"RECREATE");
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::SetNavigation(Int_t i)
{
  /// Change navigation position
  
  ++fCurrentNavigationPosition;
  fNavigation.Set(fCurrentNavigationPosition+1);
  fNavigation[fCurrentNavigationPosition] = i;
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::ShowPainterMatrix(AliMUONPainterMatrix* painterMatrix)
{
  /// Change the painter matrix we show
  
  fPainterMatrixFrame->Use(painterMatrix);
  
  painterMatrix->Connect("Clicked(AliMUONVPainter*,Double_t*)",
                         "AliMUONPainterMasterFrame",this,
                         "Clicked(AliMUONVPainter*,Double_t*)");

  painterMatrix->Connect("ShiftClicked(AliMUONVPainter*,Double_t*)",
                         "AliMUONPainterMasterFrame",this,
                         "ShiftClicked(AliMUONVPainter*,Double_t*)");
  
  fPainterMatrixFrame->Connect("TitleHasChanged(const char*)",
                              "AliMUONPainterMasterFrame",this,
                              "ChangeTitle(const char*)");
  UpdateNavigation();
  
  UpdateAttributes(*(fPainterMatrixFrame->Matrix()));
  
  AliMUONPainterRegistry::Instance()->AddToHistory(painterMatrix);
  
  Layout();
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::ChangeTitle(const char* newTitle)
{
  /// Change the top title
  
  fGroupTitle->SetText(newTitle);
  fGroupTitle->Resize();
  Layout();
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::Backward()
{
  /// Move back one step in the history
  --fCurrentNavigationPosition;
  
  AliMUONPainterMatrix* group = 
    AliMUONPainterRegistry::Instance()->PainterMatrix(fNavigation[fCurrentNavigationPosition]);
  
  ShowPainterMatrix(group);
  
  UpdateNavigation();
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::Forward()
{
  /// Move forward one step in history
  
  ++fCurrentNavigationPosition;
  
  AliMUONPainterMatrix* group = 
    AliMUONPainterRegistry::Instance()->PainterMatrix(fNavigation[fCurrentNavigationPosition]);
  
  ShowPainterMatrix(group);
  
  UpdateNavigation();
}

//_____________________________________________________________________________
void 
AliMUONPainterMasterFrame::Clicked(AliMUONVPainter* painter, Double_t* values)
{
  /// A given painter was (singly) clicked

  if ( painter->CanBeDetached() )
  {
    fPainterMatrixFrame->MouseLeave(painter);
  
    AliMUONPainterMatrix* matrix = new AliMUONPainterMatrix(painter->Name().Data());

    AliMUONVPainter* p = painter->Detach();

    p->SetResponder(1);

    matrix->Adopt(p);
  
    AddPainterMatrix(matrix);
    ShowPainterMatrix(matrix);
  }
  else
  {
    painter->DrawHistogram(values);
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterMasterFrame::ShiftClicked(AliMUONVPainter* painter, Double_t*)
{
  /// A given painter was shift-clicked
  
  if ( !painter->CanBeDetached() ) return;
  
  AliMUONPainterMatrix* currentMatrix = fPainterMatrixFrame->Matrix();
  
  AliMUONAttPainter a = painter->Attributes();
  
  TString basename(Form("%s-DUAL",painter->GetName()));
  
  TString newName = AliMUONPainterMatrix::NameIt(currentMatrix->Whatname(),basename.Data(),a);
  
  AliMUONPainterMatrix* matrix = AliMUONPainterRegistry::Instance()->PainterMatrix(newName.Data());
  
  if (!matrix)
  {
    // No. So we must make a new matrix painter from the existing one,
    // and add to this new matrix the painters of the other one, but
    // using the new attributes...
    
    // create "opposite" attributes
    AliMUONAttPainter a1(a);
    AliMUONAttPainter a2(a);
  
    a2.Invert();
    
    a1.SetCathodeAndPlaneDisabled(kTRUE);
    a2.SetCathodeAndPlaneDisabled(kTRUE);
    
    AliMUONVPainter* p1 = AliMUONVPainter::CreatePainter(painter->ClassName(),
                                                         a1,
                                                         painter->ID0(),
                                                         painter->ID1());
    
    AliMUONVPainter* p2 = AliMUONVPainter::CreatePainter(painter->ClassName(),
                                                         a2,
                                                         painter->ID0(),
                                                         painter->ID1());
    
    if (!p1 || !p2)
    {
      Int_t ret;
      new TGMsgBox(gClient->GetRoot(), this,
                   "Invalid combination", "Cannot create 2 views from this painter",
                   kMBIconExclamation, kMBOk, &ret);
      PainterMatrixWantToShow(currentMatrix);
      delete p1;
      delete p2;
      return;
    }
    
    p1->UpdateGroupsFrom(*(painter->Master()));
    p2->UpdateGroupsFrom(*(painter->Master()));
    
    p1->SetResponder(1);
    p2->SetResponder(1);
    
    Int_t nx(2);
    Int_t ny(1);
    
    AliMpArea area(painter->Area());
    
    if ( area.GetDimensionX() > 1.2*area.GetDimensionY() ) 
    {
      nx = 1;
      ny = 2;
    }
    
    matrix = new AliMUONPainterMatrix(basename.Data(),nx,ny);
    
    matrix->Adopt(p1);
    matrix->Adopt(p2);
    
    AddPainterMatrix(matrix);
  }
  
  matrix->SetData(currentMatrix->DataPattern(),
                  currentMatrix->Data(),
                  currentMatrix->DataIndex());
  
  fPainterMatrixFrame->MouseLeave(painter);
  
  PainterMatrixWantToShow(matrix);
}

//_____________________________________________________________________________
void 
AliMUONPainterMasterFrame::SaveAs(const char* filename, Option_t* option) const
{
  /// Save painter matrix (in the sense of "print") in filename
  fPainterMatrixFrame->SaveAs(filename,option);
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::Update()
{
  /// Update ourselves
  
  fPainterMatrixFrame->Update();
  fPainterMatrixFrame->UpdateInterface(kFALSE);
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::UpdateAttributes(const AliMUONPainterMatrix& painterMatrix)
{
  /// Update the view buttons from the matrix we actually plot
  
  fAttPainterSelectorFrame->Update(painterMatrix.Attributes());
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::MakeTopPainterMatrix(UInt_t w, UInt_t h, AliMUONPainterMatrix* matrix)
{
  /// Create the first painter matrix that appears when we are create
  /// FIXME: how to make this more flexible ?
  
  fPainterMatrixFrame = new AliMUONPainterMatrixFrame(this,w,h);

  if (matrix)
  {
    PainterMatrixWantToShow(matrix);
  }
  else
  {
    AliError("Cannot work without a painterMatrix");
  }
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::UpdateNavigation()
{
  /// Update navigation frame

  fBackButton->SetEnabled(kTRUE);
  fForwardButton->SetEnabled(kTRUE);

  if ( fCurrentNavigationPosition == 0 ) 
  {
    fBackButton->SetEnabled(kFALSE);
  }
  if ( fCurrentNavigationPosition == fNavigation.GetSize()-1 ) 
  {
    fForwardButton->SetEnabled(kFALSE);
  }
}

//_____________________________________________________________________________
void
AliMUONPainterMasterFrame::AttributesChanged(const AliMUONAttPainter* newValues)
{
  /// Attributes changed (e.g. from cath0 to cath1 or bending to nonbending, etc...)
  
  AliMUONPainterMatrix* currentMatrix = fPainterMatrixFrame->Matrix();
  
  AliMUONAttPainter a = currentMatrix->Validate(*newValues);
  
  if (!a.IsValid())
  {
    Int_t ret;
    new TGMsgBox(gClient->GetRoot(), this,
                 "Invalid combination", "Change of attributes not possible for this object",
                 kMBIconExclamation, kMBOk, &ret);
    PainterMatrixWantToShow(currentMatrix);
    return;
  }
  
  // First check if we already have this matrix available
  
  TString newName = AliMUONPainterMatrix::NameIt(currentMatrix->Whatname(),currentMatrix->Basename(),a);
  
  AliMUONPainterMatrix* matrix = AliMUONPainterRegistry::Instance()->PainterMatrix(newName.Data());

  if (!matrix)
  {
    // No. So we must make a new matrix painter from the existing one,
    // and add to this new matrix the painters of the other one, but
    // using the new attributes...
    
    matrix = currentMatrix->Clone(a);
  
    AddPainterMatrix(matrix);
  }
  
  matrix->SetData(currentMatrix->DataPattern(),
                  currentMatrix->Data(),
                  currentMatrix->DataIndex());
  
  PainterMatrixWantToShow(matrix);
}
