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

const Int_t AliMUONTrackerDataCompareDialog::fgkDifference(1);
const Int_t AliMUONTrackerDataCompareDialog::fgkAbsoluteDifference(2);
const Int_t AliMUONTrackerDataCompareDialog::fgkRelativeDifference(3);
const Int_t AliMUONTrackerDataCompareDialog::fgkAbsoluteRelativeDifference(4);
const Int_t AliMUONTrackerDataCompareDialog::fgkAll(5);

namespace
{
  
#define PRECISION 1E-12
  
  Double_t Difference(Double_t v1, Double_t v2)
  {
    Double_t d = v1-v2;
    return TMath::Abs(d) < PRECISION ? 0.0 : d;
  }
    
  Double_t AbsoluteDifference(Double_t v1, Double_t v2)
  {
    return TMath::Abs(Difference(v1,v2));
  }
  
  
  Double_t RelativeDifference(Double_t v1, Double_t v2)
  {
    if ( TMath::Abs(v1) < PRECISION ) return 0.0;
    return (v1-v2)/v1;
  }
  
  Double_t AbsoluteRelativeDifference(Double_t v1, Double_t v2)
  {
    return TMath::Abs(RelativeDifference(v1,v2));
  }
}


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
  
  fDiffType->AddEntry("Difference",fgkDifference);
  fDiffType->AddEntry("Absolute difference",fgkAbsoluteDifference);
  fDiffType->AddEntry("Relative difference",fgkRelativeDifference);
  fDiffType->AddEntry("Absolute relative difference",fgkAbsoluteRelativeDifference);
  fDiffType->AddEntry("All four",fgkAll);

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
  
  Int_t nd = fDiffType->GetSelected();
  
  if ( nd == fgkAll ) 
  {
    CompareData(s1.Data(),s2.Data(),fgkDifference);
    CompareData(s1.Data(),s2.Data(),fgkRelativeDifference);
    CompareData(s1.Data(),s2.Data(),fgkAbsoluteDifference);
    CompareData(s1.Data(),s2.Data(),fgkAbsoluteRelativeDifference);
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
                                             Int_t difftype) const
{
  /// Compare two data sources
  
  AliMUONPainterDataRegistry* reg = AliMUONPainterDataRegistry::Instance();
  
  AliMUONVTrackerData* d1 = reg->DataSource(d1name);
  if (!d1)
  {
    AliError(Form("Cannot find data source %s",d1name));
  }
  
  AliMUONVTrackerData* d2 = reg->DataSource(d2name);
  if (!d2)
  {
    AliError(Form("Cannot find data source %s",d2name));
  }
  
  Double_t (*difffunction)(Double_t,Double_t)=0x0;
  TString suffix("unknown");
  
  if ( difftype == fgkDifference ) 
  {
    difffunction = Difference;
    suffix = "D";
  }
  if ( difftype == fgkAbsoluteDifference ) 
  {
    difffunction = AbsoluteDifference;
    suffix = "AD";
  }
  if ( difftype == fgkRelativeDifference ) 
  {
    difffunction = RelativeDifference;
    suffix = "RD";
  }
  if ( difftype == fgkAbsoluteRelativeDifference ) 
  {
    difffunction = AbsoluteRelativeDifference;
    suffix = "ARD";
  }
  
  TString basename = fBasename->GetText(); 
  
  AliMUONVTrackerData* d = CompareData(*d1,*d2,Form("%s:%s",basename.Data(),suffix.Data()),difffunction);
  
  AliMUONVTrackerDataMaker* dw = new AliMUONTrackerDataWrapper(d);
  
  AliMUONPainterDataRegistry::Instance()->Register(dw);
}

//______________________________________________________________________________
AliMUONVTrackerData*
AliMUONTrackerDataCompareDialog::CompareData(const AliMUONVTrackerData& d1,
                                             const AliMUONVTrackerData& d2,
                                             const char* outname,
                                             Double_t(*diff)(Double_t,Double_t)) const
{
  /// Compare two data objects, using the diff method
  
  if ( d1.NumberOfDimensions() != d2.NumberOfDimensions() ) 
  {
    AliError("Cannot compare data of incompatible dimensions");
    return 0x0;
  }
  
  AliMpManuIterator it;
  Int_t detElemId, manuId;
  
  AliMUONVStore* store = new AliMUON2DMap(kTRUE);
  
  while ( it.Next(detElemId,manuId) )
  {
    if ( d1.HasDetectionElement(detElemId) && d2.HasDetectionElement(detElemId) &&
         d1.HasManu(detElemId,manuId) && d2.HasManu(detElemId,manuId) )
    {
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
      
      AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(store->FindObject(detElemId,manuId));
      
      if (!param)
	    {
	      param = new AliMUONCalibParamND(d1.ExternalDimension(),64,detElemId,manuId,
                                        AliMUONVCalibParam::InvalidFloatValue());
	      store->Add(param);
	    }
      
      for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i ) 
	    {
	      if ( de->IsConnectedChannel(manuId,i) )
        {
          for ( Int_t k = 0; k < d1.ExternalDimension(); ++k ) 
          {
            
            Double_t d = diff(d1.Channel(detElemId,manuId,i,k),
                                  d2.Channel(detElemId,manuId,i,k));
          
            param->SetValueAsDouble(i,k,d);
          }
        }
	    }
    }
  }
  
  AliMUONVTrackerData* d = new AliMUONTrackerData(outname,outname,d1.ExternalDimension(),kTRUE);
  for ( Int_t k = 0; k < d1.ExternalDimension(); ++k ) 
  {
    d->SetDimensionName(k,Form("D:%s",d1.ExternalDimensionName(k).Data()));
  }
  d->Add(*store);
  
  return d;
}

