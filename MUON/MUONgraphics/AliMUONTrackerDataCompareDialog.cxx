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

#include "AliMUONTrackerDataCompareDialog.h"

/// \class AliMUONTrackerDataCompareDialog
///
/// Widget to select 2 VTrackerData objects (D1,D2) to be compared
///
/// The type of differences that can be used are : 
///
/// - Difference = plain difference D1-D2
/// - Absolute difference = absolute value of the preceeding = |D1-D2|
/// - Relative difference = relative difference = (D1-D2)/D1
/// - Absolute relative difference = absolute value of preceeding = |(D1-D2)/D1|
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONTrackerData.h"
#include "AliMUONTrackerDataWrapper.h"
#include "AliMUONVTrackerData.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include <TGComboBox.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TTimer.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrackerDataCompareDialog)
/// \endcond



//_____________________________________________________________________________
AliMUONTrackerDataCompareDialog::AliMUONTrackerDataCompareDialog(const TGWindow* p, const TGWindow* main, UInt_t w, UInt_t h)
: TGTransientFrame(p,main,w,h),
fF1(new TGHorizontalFrame(this)),
fData1(new TGComboBox(fF1)),
fF2(new TGHorizontalFrame(this)),
fData2(new TGComboBox(fF2)),
fF3(new TGHorizontalFrame(this)),
fDiffType(new TGComboBox(fF3)),
fF4(new TGHorizontalFrame(this)),
fBasename(new TGTextEntry(fF4)),
fButtonFrame(new TGHorizontalFrame(this)),
fOK(new TGTextButton(fButtonFrame,"OK")),
fCancel(new TGTextButton(fButtonFrame,"Cancel"))
{
  /// ctor
  
  SetCleanup(kDeepCleanup);
  
  AliMUONPainterDataRegistry* reg = AliMUONPainterDataRegistry::Instance();
  
  for ( Int_t i = 0; i < reg->NumberOfDataSources(); ++i ) 
  {
    AliMUONVTrackerData* data = reg->DataSource(i);
    fData1->AddEntry(data->GetName(),i);
    fData2->AddEntry(data->GetName(),i);
  }
  
  fDiffType->AddEntry("Difference",AliMUONTrackerData::kDifference);
  fDiffType->AddEntry("Absolute difference",AliMUONTrackerData::kAbsoluteDifference);
  fDiffType->AddEntry("Relative difference",AliMUONTrackerData::kRelativeDifference);
  fDiffType->AddEntry("Absolute relative difference",AliMUONTrackerData::kAbsoluteRelativeDifference);
  fDiffType->AddEntry("All four",AliMUONTrackerData::kAll);

  fData1->Select(0);
  fData2->Select(0);
  fDiffType->Select(4);
  
  fF1->AddFrame(new TGLabel(fF1,"First data"),new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
  fF1->AddFrame(fData1,new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsTop,5,5,5,5));

  fF2->AddFrame(new TGLabel(fF2,"Second data"),new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
  fF2->AddFrame(fData2,new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsTop,5,5,5,5));

  fF3->AddFrame(new TGLabel(fF3,"Difference type"),new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
  fF3->AddFrame(fDiffType,new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsTop,5,5,5,5));

  fF4->AddFrame(new TGLabel(fF4,"Output basename"),new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
  fF4->AddFrame(fBasename,new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsTop,5,5,5,5));

  AddFrame(fF1,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  AddFrame(fF2,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  AddFrame(fF3,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  AddFrame(fF4,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  
  fButtonFrame->AddFrame(fOK,new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
  fButtonFrame->AddFrame(fCancel,new TGLayoutHints(kLHintsRight|kLHintsTop,5,5,5,5));

  AddFrame(fButtonFrame,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  
  fData1->Resize(200,20);
  fData2->Resize(200,20);
  fDiffType->Resize(200,20);
  
  fOK->Connect("Clicked()", "AliMUONTrackerDataCompareDialog",this,"DoOK()");
  fCancel->Connect("Clicked()","AliMUONTrackerDataCompareDialog",this,"DoCancel()");
}

//_____________________________________________________________________________
AliMUONTrackerDataCompareDialog::~AliMUONTrackerDataCompareDialog()
{
  /// dtor
}

//______________________________________________________________________________
void
AliMUONTrackerDataCompareDialog::DoOK()
{
  /// Do the job.
  
  TGTextLBEntry* t1 = static_cast<TGTextLBEntry*>(fData1->GetSelectedEntry());
  TString s1 = t1->GetText()->GetString();
  TGTextLBEntry* t2 = static_cast<TGTextLBEntry*>(fData2->GetSelectedEntry());
  TString s2 = t2->GetText()->GetString();
  
  AliMUONTrackerData::EDiffType nd = static_cast<AliMUONTrackerData::EDiffType>(fDiffType->GetSelected());
  
  if ( nd == AliMUONTrackerData::kAll )
  {
    CompareData(s1.Data(),s2.Data(),AliMUONTrackerData::kDifference);
    CompareData(s1.Data(),s2.Data(),AliMUONTrackerData::kRelativeDifference);
    CompareData(s1.Data(),s2.Data(),AliMUONTrackerData::kAbsoluteDifference);
    CompareData(s1.Data(),s2.Data(),AliMUONTrackerData::kAbsoluteRelativeDifference);
  }
  else
  {
    CompareData(s1.Data(),s2.Data(),nd);
  }
  
  TTimer::SingleShot(150,"AliMUONTrackerDataCompareDialog",this,"CloseWindow()");
}

//______________________________________________________________________________
void
AliMUONTrackerDataCompareDialog::DoCancel()
{
  /// Kills the dialog
  TTimer::SingleShot(150,"AliMUONTrackerDataCompareDialog",this,"CloseWindow()");
}

//______________________________________________________________________________
void
AliMUONTrackerDataCompareDialog::CompareData(const char* d1name,
                                             const char* d2name,
                                             AliMUONTrackerData::EDiffType difftype) const
{
  /// Compare two data sources
  
  AliMUONPainterDataRegistry* reg = AliMUONPainterDataRegistry::Instance();
  
  AliMUONVTrackerData* d1 = reg->DataSource(d1name);
  if (!d1)
  {
    AliError(Form("Cannot find data source %s",d1name));
    return;
  }
  
  AliMUONVTrackerData* d2 = reg->DataSource(d2name);
  if (!d2)
  {
    AliError(Form("Cannot find data source %s",d2name));
    return;
  }
  
  TString basename = fBasename->GetText();
  
  AliMUONVTrackerData* d = AliMUONTrackerData::CompareData(*d1,*d2,basename.Data(),difftype);
  
  if (d)
  {
    AliMUONVTrackerDataMaker* dw = new AliMUONTrackerDataWrapper(d);
    
    AliMUONPainterDataRegistry::Instance()->Register(dw);
  }
}
