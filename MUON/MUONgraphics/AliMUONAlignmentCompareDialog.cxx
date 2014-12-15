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

#include "AliMUONAlignmentCompareDialog.h"

/// \class AliMUONAlignmentCompareDialog
///
/// Widget to select 2 alignments objects from the OCDB (A1,A2) to be compared
///
/// \author Philippe Pillot, Laurent Aphecetche
///

#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONTrackerData.h"
#include "AliMUONTrackerDataWrapper.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpVSegmentation.h"
#include <TGComboBox.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGTextEntry.h>
#include <TGeoManager.h>
#include <TTimer.h>

/// \cond CLASSIMP
ClassImp(AliMUONAlignmentCompareDialog)
/// \endcond

namespace
{
  
#define PRECISION 1E-12
  
  Double_t Difference(Double_t v1, Double_t v2)
  {
    Double_t d = v1-v2;
    return TMath::Abs(d) < PRECISION ? 0.0 : d;
  }  
}

//_____________________________________________________________________________
AliMUONAlignmentCompareDialog::AliMUONAlignmentCompareDialog(const TGWindow* p, const TGWindow* main, UInt_t w, UInt_t h)
: TGTransientFrame(p,main,w,h),
fF1(new TGVerticalFrame(this)),
fOCDBPath1(0x0),
fRun1(0x0),
fF2(new TGVerticalFrame(this)),
fOCDBPath2(0x0),
fRun2(0x0),
fF3(new TGHorizontalFrame(this)),
fBasename(new TGTextEntry(fF3)),
fButtonFrame(new TGHorizontalFrame(this)),
fOK(new TGTextButton(fButtonFrame,"OK")),
fCancel(new TGTextButton(fButtonFrame,"Cancel"))
{
  /// ctor
  
  SetCleanup(kDeepCleanup);
    
  AddInput(fF1,"First alignment",fOCDBPath1,fRun1);
  AddInput(fF2,"Second alignment",fOCDBPath2,fRun2);
    
  fF3->AddFrame(new TGLabel(fF3,"Output basename"),new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
  fF3->AddFrame(fBasename,new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsTop,5,5,5,5));

  AddFrame(fF1,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  AddFrame(fF2,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  AddFrame(fF3,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  
  fButtonFrame->AddFrame(fOK,new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
  fButtonFrame->AddFrame(fCancel,new TGLayoutHints(kLHintsRight|kLHintsTop,5,5,5,5));

  AddFrame(fButtonFrame,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
  
  fOK->Connect("Clicked()", "AliMUONAlignmentCompareDialog",this,"DoOK()");
  fCancel->Connect("Clicked()","AliMUONAlignmentCompareDialog",this,"DoCancel()");
}

//_____________________________________________________________________________
AliMUONAlignmentCompareDialog::~AliMUONAlignmentCompareDialog()
{
  /// dtor
}

//_____________________________________________________________________________
void AliMUONAlignmentCompareDialog::AddInput(TGCompositeFrame* frame, const char* msg,
                                             TGTextEntry*& text, TGNumberEntry*& run)
{
    
    TGHorizontalFrame* hf1 = new TGHorizontalFrame(frame);
    
    hf1->AddFrame(new TGLabel(hf1,TString(msg) + " ocdb path"),new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
    
    text = new TGTextEntry(hf1,"alien://folder=/alice/data/2012/OCDB");
    
    hf1->AddFrame(text,new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsTop,5,5,5,5));
    
    TGHorizontalFrame* hf2 = new TGHorizontalFrame(frame);
    
    hf2->AddFrame(new TGLabel(hf2,TString(msg) + " run number"),new TGLayoutHints(kLHintsLeft|kLHintsTop,5,5,5,5));
    
    run = new TGNumberEntry(hf2);
    
    hf2->AddFrame(run,new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsTop,5,5,5,5));
    
    frame->AddFrame(hf1,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
    frame->AddFrame(hf2,new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsTop,5,5,5,5));
}

//______________________________________________________________________________
void
AliMUONAlignmentCompareDialog::DoOK()
{
  /// Do the job.
  
  AliMUONVTrackerData* d = CompareAlignment(fOCDBPath1->GetText(),fRun1->GetNumber(),
                                            fOCDBPath2->GetText(),fRun2->GetNumber());
                   
  if (!d) return;
  
  TString basename = fBasename->GetText(); 
  
  AliMUONVTrackerDataMaker* dw = new AliMUONTrackerDataWrapper(d);
  
  AliMUONPainterDataRegistry::Instance()->Register(dw);
  
  TTimer::SingleShot(150,"AliMUONAlignmentCompareDialog",this,"CloseWindow()");
}

//______________________________________________________________________________
void
AliMUONAlignmentCompareDialog::DoCancel()
{
  /// Kills the dialog
  TTimer::SingleShot(150,"AliMUONAlignmentCompareDialog",this,"CloseWindow()");
}

//______________________________________________________________________________
AliMUONVTrackerData*
AliMUONAlignmentCompareDialog::CompareAlignment(const char* ocdbPathForAlign1, Int_t run1,
                                                const char* ocdbPathForAlign2, Int_t run2)
{
  // ocdb access
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage("raw://");
  
  // get geometry transformers
  AliMUONGeometryTransformer geoTransformer[2];

  const char* align[2] = {
    ocdbPathForAlign1,
    ocdbPathForAlign2
  };
  
  Int_t runs[] = {
     run1,run2
  };
  
  for (Int_t i = 0; i < 2; i++) 
  {
    cdbm->UnloadFromCache("GRP/Geometry/Data");
    cdbm->UnloadFromCache("MUON/Align/Data");
    AliGeomManager::GetGeometry()->UnlockGeometry();
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return 0x0;
    cdbm->SetSpecificStorage("MUON/Align/Data",align[i]);
    cdbm->SetRun(runs[i]);
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    geoTransformer[i].LoadGeometryData();
  }
  
  // store for cluster shifts
  AliMUON2DMap shiftStore(kTRUE);
  
  // loop over chamber
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    
    // loop over DEs
    AliMpDEIterator nextDE;
    nextDE.First(iCh);
    while (!nextDE.IsDone()) {
      
      Int_t deId = nextDE.CurrentDE()->GetId();
      
      // loop over cathods
      for (Int_t icath = 0; icath < 2; icath++) {
        
        const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::GetCathodType(icath));
        
        // loop over pads
        AliMpVPadIterator *nextPad = seg->CreateIterator();
        nextPad->First();
        while (!nextPad->IsDone()) {
          
          AliMpPad pad = nextPad->CurrentItem();
          Int_t manuId = pad.GetManuId();
          Int_t manuChannel = pad.GetManuChannel();
          
          // local position
          Double_t xl = pad.GetPositionX();
          Double_t yl = pad.GetPositionY();
          Double_t zl = 0.;
          
          // position with first alignment
          Double_t x1, y1, z1;
          geoTransformer[0].Local2Global(deId,xl,yl,zl,x1,y1,z1);
          
          // position with second alignment
          Double_t x2, y2, z2;
          geoTransformer[1].Local2Global(deId,xl,yl,zl,x2,y2,z2);
          
          // pad shift
          Double_t dx = ::Difference(x2,x1);
          Double_t dy = ::Difference(y2,y1);
          Double_t dz = ::Difference(z2,z1);
          
          // store pad shifts
          AliMUONVCalibParam* p = static_cast<AliMUONVCalibParam*>(shiftStore.FindObject(deId,manuId));
          if (!p) {
            p = new AliMUONCalibParamND(3,AliMpConstants::ManuNofChannels(),deId,manuId,0.);
            shiftStore.Add(p);
          }
          p->SetValueAsDouble(manuChannel,0,dx);
          p->SetValueAsDouble(manuChannel,1,dy);
          p->SetValueAsDouble(manuChannel,2,dz);
          
          nextPad->Next();
        }
        
        delete nextPad;
        
      }
      
      nextDE.Next();
    }
    
  }
  
  // create tracker data
  AliMUONTrackerData* data = new AliMUONTrackerData(fBasename->GetText(),fBasename->GetText(),3,kTRUE);
  data->SetDimensionName(0,"dx"); // max shift in x
  data->SetDimensionName(1,"dy"); // max shift in y
  data->SetDimensionName(2,"dz"); // max shift in z
  data->Add(shiftStore);
  
  return data;
}