// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file muon_digits.C
/// \brief  Macro to visualise digits from MUON spectrometer 
/// (both tracker and trigger).
///
/// Use muon_digits() in order to run it
///
/// Needs that alieve_init() is already called
///
/// \author P. Pillot, L. Aphecetche; Subatech

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TTree.h>
#include <TStyle.h>
#include <TEveManager.h>
#include <TEveQuadSet.h>

#include <AliLog.h>
#include <AliMUONGeometryTransformer.h>
#include <AliMUONVDigit.h>
#include <AliMUONVDigitStore.h>
#include <AliMpPad.h>
#include <AliMpSegmentation.h>
#include <AliMpVSegmentation.h>
#include <AliMpCDB.h>
#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

//______________________________________________________________________________
void add_muon_digits(TIter* next, TEveQuadSet* bending, TEveQuadSet* nonBending, Bool_t fromRaw)
{
  // load mapping
  AliMpCDB::LoadAll(kFALSE);
  
  // load geometry
  static AliMUONGeometryTransformer* gMUONGeometryTransformer = 0x0;
  if (!gMUONGeometryTransformer) 
  {
    AliEveEventManager::Instance()->AssertGeometry();
    gMUONGeometryTransformer = new AliMUONGeometryTransformer();
    gMUONGeometryTransformer->LoadGeometryData();
  }
  
  // loop over digits and produce corresponding graphic objects
  AliMUONVDigit* digit;
  while ( ( digit = static_cast<AliMUONVDigit*>((*next)() ) ) )
  {
    if (!digit->IsTrigger() && !fromRaw && digit->Charge() < 1.e-3) continue;
    
    Int_t detElemId = digit->DetElemId();
    Int_t manuId = digit->ManuId();
    
    const AliMpVSegmentation* vseg =
      AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(digit->Cathode()));
    if (!vseg) 
    {
      AliInfoGeneral("muon_digits.C", Form("Could not get segmentation for DE %4d MANU %4d",detElemId,manuId));
      continue; // should not happen, unless we got a readout error and thus a bad de,manu pair
    }
    
    AliMpPad pad = vseg->PadByLocation(manuId,digit->ManuChannel());
    
    Double_t local[] = { pad.GetPositionX(), pad.GetPositionY(), 0.0 };
    Double_t global[] = { 0.0, 0.0, 0.0 };
    
    gMUONGeometryTransformer->Local2Global(detElemId,
                                           local[0], local[1], local[2],
                                           global[0], global[1], global[2]);
    
    TEveQuadSet* pads = bending;
    if (vseg->PlaneType()==AliMp::kNonBendingPlane) pads = nonBending;
    
    pads->AddQuad(global[0]-pad.GetDimensionX(),global[1]-pad.GetDimensionY(),global[2],
		  2.*pad.GetDimensionX(),2.*pad.GetDimensionY());
    
    if (fromRaw && !digit->IsTrigger()) pads->QuadValue(digit->ADC());
    else pads->QuadValue((Int_t) digit->Charge());
  }
  
}

//______________________________________________________________________________
void muon_digits()
{
    printf("*** DIGITS MUON ***");
    
  // load digits
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("MUON");
  TTree* dt = rl->GetTreeD("MUON", kFALSE);
  if (!dt) return;
  AliMUONVDigitStore *digitStore = AliMUONVDigitStore::Create(*dt);
  digitStore->Clear();
  digitStore->Connect(*dt,0);
  dt->GetEvent(0);
  rl->UnloadDigits("MUON");
  
  if (digitStore->GetSize() == 0 && !gEve->GetKeepEmptyCont()) {
    delete digitStore;
    return;
  }
  
  // container for graphic representation of digits
  TEveElementList* cont = new TEveElementList("MUON Digits");
  
  TEveQuadSet* bending = new TEveQuadSet(TEveQuadSet::kQT_RectangleXY, kFALSE, 32);
  bending->SetName("Bending");
  bending->SetRenderMode(TEveDigitSet::kRM_Fill);
  bending->SetPickable(kFALSE);
  cont->AddElement(bending);
  
  TEveQuadSet* nonBending = new TEveQuadSet(TEveQuadSet::kQT_RectangleXY, kFALSE, 32);
  nonBending->SetName("Non bending");
  nonBending->SetRenderMode(TEveDigitSet::kRM_Line);
  nonBending->SetPickable(kFALSE);
  cont->AddElement(nonBending);
  
  // add digits to the containers
  TIter next(digitStore->CreateIterator());
  add_muon_digits(&next, bending, nonBending, kFALSE);
  delete digitStore;
  
  // set containers' title
  Int_t nDigitB = bending->GetPlex()->Size();
  Int_t nDigitNB = nonBending->GetPlex()->Size();
  cont->SetTitle(Form("N=%d",nDigitB+nDigitNB));
  bending->SetTitle(Form("N=%d",nDigitB));
  nonBending->SetTitle(Form("N=%d",nDigitNB));
  
  // automatic scaling
  gStyle->SetPalette(1);
  bending->AssertPalette();
  nonBending->AssertPalette();
  
  // add graphic containers
  gEve->DisableRedraw();
  gEve->AddElement(cont);
  gEve->EnableRedraw();
  gEve->Redraw3D();
}
