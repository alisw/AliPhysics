// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include <TEveManager.h>
#include <TEveElement.h>

#include <AliEveITSModule.h>
#include <AliEveITSDigitsInfo.h>
#include <AliEveITSScaledModule.h>
#endif


void its_common_foos()
{}

AliEveITSModule* its_make_module(Int_t i, TEveElement* parent,
                                 AliEveITSDigitsInfo*  di,
                                 AliEveDigitScaleInfo* si,
                                 Bool_t check_empty,
				 Bool_t scaled_modules)
{
  AliEveITSModule* m = 0;

  Int_t det_id = 0;
  if (i > 239 && i < 500) det_id = 1;
  else if (i >= 500)      det_id = 2;

  if (!check_empty || di->HasData(i, det_id) || di->IsDead(i, det_id))
  {
    if (scaled_modules)
      m = new AliEveITSScaledModule(i, di, si);
    else
      m = new AliEveITSModule(i, di);

    // Before 5.26 ROOT did not draw frames of empty quad-sets.
    // Bypass until we move there.
    if (!di->HasData(i, det_id))
    {
      m->AddQuad(0,0,0,0);
      m->RefitPlex();
    }

    if (parent)
      parent->AddElement(m);
  }

  return m;
}

void its_display_raw_digits(AliEveITSDigitsInfo* di, Int_t mode,
                            Bool_t check_empty    = kTRUE,
                            Bool_t scaled_modules = kFALSE)
{
  const TString bsSector = "Sector";
  const TString bsStave  = "Stave";
  const TString bsLadder = "Ladder";
  TString sSector;
  TString sStave;
  TString sLadder;

  Int_t  i=0;
  Long_t nsec, nstave, nlad, nMod;

  gEve->DisableRedraw();

  AliEveDigitScaleInfo* si = 0;
  if (scaled_modules)
  {
    si = new AliEveDigitScaleInfo;
  }

  if (mode & 1)
  {
    TEveElementList* l = new TEveElementList("SPD0");
    l->SetTitle("SPDs' first layer");
    l->SetMainColor(kRed);
    gEve->AddElement(l);

    for (nsec=0; nsec<10; nsec++)
    {
      TEveElementList* relSector = new TEveElementList(bsSector + nsec);
      relSector->SetMainColor(kRed);
      l->AddElement(relSector);

      for (nstave=0; nstave<2; nstave++)
      {
	TEveElementList* relStave = new TEveElementList(bsStave + nstave);
	relStave->SetMainColor(kRed);
	relSector->AddElement(relStave);

	for (nMod=0; nMod<4; ++nMod, ++i)
	{
          its_make_module(i, relStave, di, si, check_empty, scaled_modules);
	}
      }
    }
  }
  else
  {
    i += 10*2*4;
  }

  if (mode & 2)
  {
    TEveElementList* l = new TEveElementList("SPD1");
    l->SetTitle("SPDs' second layer");
    l->SetMainColor(kRed);
    gEve->AddElement(l);

    for (nsec=0; nsec<10; nsec++)
    {
      TEveElementList* relSector = new TEveElementList(bsSector + nsec);
      relSector->SetMainColor(kRed);
      l->AddElement(relSector);

      for (nstave=0; nstave<4; nstave++)
      {
	TEveElementList* relStave = new TEveElementList(bsStave + nstave);
	relStave->SetMainColor(kRed);
	relSector->AddElement(relStave);

	for (nMod=0; nMod<4;  ++nMod, ++i)
	{
          its_make_module(i, relStave, di, si, check_empty, scaled_modules);
	}
      }
    }
  }
  else
  {
    i += 10*4*4;
  }

  if (mode & 4)
  {
    TEveElementList* l = new TEveElementList("SDD2");
    l->SetTitle("SDDs' first layer");
    l->SetMainColor(kBlue);
    gEve->AddElement(l);

    for (nlad=0; nlad<14; nlad++)
    {
      TEveElementList* relLadder = new TEveElementList(bsLadder + nlad);
      relLadder->SetMainColor(kBlue);
      l->AddElement(relLadder);
      for (nMod=0; nMod<6; ++nMod, ++i)
      {
        its_make_module(i, relLadder, di, si, check_empty, scaled_modules);
      }
    }
  }
  else
  {
    i += 14*6;
  }

  if (mode & 8)
  {
    TEveElementList* l = new TEveElementList("SDD3");
    l->SetTitle("SDDs' second layer");
    l->SetMainColor(kBlue);
    gEve->AddElement(l);

    for (nlad=0; nlad<22; nlad++)
    {
      TEveElementList* relLadder = new TEveElementList(bsLadder + nlad);
      relLadder->SetMainColor(kBlue);
      l->AddElement(relLadder);
      for (nMod=0; nMod<8;  ++nMod, ++i)
      {
        its_make_module(i, relLadder, di, si, check_empty, scaled_modules);
      }
    }
  }
  else
  {
    i += 22*8;
  }

  if (mode & 16)
  {
    TEveElementList* l = new TEveElementList("SSD4");
    l->SetTitle("SSDs' first layer");
    l->SetMainColor(kGreen);
    gEve->AddElement(l);

    for (nlad=0; nlad<34; nlad++)
    {
      TEveElementList* relLadder = new TEveElementList(bsLadder + nlad);
      relLadder->SetMainColor(kGreen);
      l->AddElement(relLadder);
      for (nMod=0; nMod<22; ++nMod, ++i)
      {
        its_make_module(i, relLadder, di, si, check_empty, scaled_modules);
      }
    }
  }
  else
  {
    i += 34*22;
  }

  if (mode & 32)
  {
    TEveElementList* l = new TEveElementList("SSD5");
    l->SetTitle("SSDs' second layer");
    l->SetMainColor(kGreen);
    gEve->AddElement(l);

    for (nlad=0; nlad<38; nlad++)
    {
      TEveElementList* relLadder = new TEveElementList(bsLadder + nlad);
      relLadder->SetMainColor(kGreen);
      l->AddElement(relLadder);
      for (nMod=0; nMod<25; ++nMod, ++i)
      {
        its_make_module(i, relLadder, di, si, check_empty, scaled_modules);
      }
    }
  }
  else
  {
    i += 38*25;
  }

  gEve->EnableRedraw();
}
