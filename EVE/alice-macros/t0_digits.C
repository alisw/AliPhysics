// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TStyle.h>

#include <AliRunLoader.h>
#include <AliT0digit.h>
#include <AliEveEventManager.h>
#include <AliEveITSDigitsInfo.h>
#include <AliEveT0Module.h>
#endif

void t0_digits()
{
  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();

  rl->LoadDigits("T0");
  TTree* dt = rl->GetTreeD("T0", false);

  AliT0digit *digits = 0;
  dt->SetBranchAddress("T0", &digits);
  dt->GetEntry(0);

  gStyle->SetPalette(1, 0);

  AliEveT0Module::MakeModules(digits);
}

