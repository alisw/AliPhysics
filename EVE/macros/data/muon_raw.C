// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file muon_raw.C
/// \brief Macro to visualise rootified raw-data from MUON spectrometer 
/// (both tracker and trigger).
///
/// Use muon_raw() in order to run it
///
/// Needs that alieve_init() is already called
///
/// \author P. Pillot, L. Aphecetche; Subatech

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <Riostream.h>
#include <TROOT.h>
#include <TEveManager.h>
#include <TEveUtil.h>
#include <TEveQuadSet.h>

#include <AliLog.h>
#include <AliMUONDigitMaker.h>
#include <AliMUONDigitStoreV2R.h>
#include <AliMpCDB.h>
#include <AliRawReader.h>
#include <AliEveEventManager.h>
#endif

void muon_raw()
{
    printf("*** RAW MUON ***");
    
  // load mapping
  AliMpCDB::LoadAll(kFALSE);

  // load raw data
  AliRawReader* reader = AliEveEventManager::AssertRawReader();
/*  if ( reader->GetEventHeader() ) 
    AliInfoGeneral("muon_raw.C", Form("RUN %d EVENT %d", reader->GetRunNumber(),reader->GetEventIndex()) );
  else
    AliInfoGeneral("muon_raw.C", "NO EVENT HEADER ?");
*/  
  // convert raw to digits
  AliMUONDigitMaker digitMaker;
  digitMaker.SetMakeTriggerDigits(kTRUE);
  AliMUONDigitStoreV2R digitStore;
  digitMaker.Raw2Digits(reader,&digitStore);
  if (digitStore.GetSize() == 0 && !gEve->GetKeepEmptyCont()) return;
  
  // container for graphic representation of digits
  TEveElementList* cont = new TEveElementList("MUON Raw digits");
  cont->SetTitle(Form("N=%d",digitStore.GetSize()));
  
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
  TEveUtil::LoadMacro("muon_digits.C+");
  TIter next(digitStore.CreateIterator());
  gROOT->ProcessLine(Form("add_muon_digits((TIter*)%p, (TEveQuadSet*)%p, (TEveQuadSet*)%p, kTRUE);",
			  &next, bending, nonBending));
  
  // set containers' title
  bending->SetTitle(Form("N=%d",bending->GetPlex()->Size()));
  nonBending->SetTitle(Form("N=%d",nonBending->GetPlex()->Size()));
  
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
