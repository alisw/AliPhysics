// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TEveUtil.h>

#include <AliRawReader.h>
#include <AliEveEventManager.h>
#include <AliEveITSDigitsInfo.h>

#include "its_common_foos.C"
#endif

// Load ITS raw-data.
// Argument mode is a bitwise or determining which layers to import:
//    1,  2 : SPD
//    4,  8 : SDD
//   16, 32 : SSD
// By default import all layers.

void its_raw(Int_t mode            = 63,
             Bool_t check_empty    = kTRUE,
             Bool_t scaled_modules = kFALSE)
{
    printf("*** RAW ITS ***");
    
  AliRawReader *rawReader = AliEveEventManager::AssertRawReader();

  TEveUtil::LoadMacro("its_common_foos.C");

  AliEveITSDigitsInfo* di = new AliEveITSDigitsInfo();
  di->ReadRaw(rawReader,mode);
  // di->Dump();

  gStyle->SetPalette(1, 0);

  its_display_raw_digits(di, mode, check_empty, scaled_modules);
}
