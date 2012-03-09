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
#include <TString.h>
#include <TEveElement.h>
#include <TEveFrameBox.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>
#include <TEveTrans.h>

#include <AliHMPIDDigit.h>
#include <AliHMPIDv3.h>
#include <AliCluster3D.h>
#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#include <AliEvePMDModule.h>
#endif

void pmd_digits(Int_t mode = 0)
{
  Bool_t drawBackPolygons = kFALSE;

  gStyle->SetPalette(1, 0);


  TString spl;

  TString sddl;
  TString bsddl="DDL";


  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("PMD");
  TTree* pmdt = rl->GetTreeD("PMD", false);

  //  cout << pmdt->GetEntries() << endl;

  gEve->DisableRedraw();

  TEveElementList* l = new TEveElementList("PMD");
  // l->SetTitle("tooltip");
  // l->SetMainColor(3);
  gEve->AddElement(l);

  TEveRGBAPalette* pal = new TEveRGBAPalette(20, 1000);
  pal->SetLimits(0, 1024);

  Int_t NSM         = 0;
  Int_t istartDDL   = 0;
  Int_t iendDDL     = 0;
  Int_t modnumber   = 0;
  Int_t istartPlane = 0;
  Int_t iendPlane   = 0;
  Float_t zpos      = 0;

  switch(mode)
    {
    case 0:
      istartPlane = 0;
      iendPlane   = 1;
      printf("--- Visualization is set for PREshower Plane ---\n");
      break;

    case 1:
      istartPlane = 1;
      iendPlane   = 2;
      printf("--- Visualization is set for CPV Plane ---\n");
      break;

    case 2:
      istartPlane = 0;
      iendPlane   = 2;
      printf("--- Visualization is set for both Planes ---\n");
      break;

    default:
      printf("--- Not set for any Plane ---\n");
    }

  for (Int_t ipl = istartPlane; ipl < iendPlane; ipl++)
    {

      if (ipl == 0)
	{
	  spl       = "PRE";
	  istartDDL = 0;
	  iendDDL   = 4;
	  zpos      = 365.;
	}
      if (ipl == 1)
	{
	  spl       = "CPV";
	  istartDDL = 4;
	  iendDDL   = 6;
	  zpos      = 360.;
	}

      TEveElementList* lplane = new TEveElementList(spl.Data());
      //  l->SetMainColor(3);
      gEve->AddElement(lplane, l);

      for (Int_t iddl = istartDDL; iddl < iendDDL; iddl++)
	{
	  sddl = bsddl;
	  sddl += iddl;
	  TEveElementList* lddl = new TEveElementList(sddl.Data());
	  //  l->SetMainColor(3);
	  gEve->AddElement(lddl, lplane);

	  modnumber = iddl*6;

	  if (iddl < 4)
	    {
	      NSM = 6;
	    }
	  else if (iddl >=4 && iddl < 6)
	    {
	      NSM = 12;
	    }

	  for (Int_t ism = 0; ism < NSM; ism++)
	    {

	      AliEvePMDModule *lmodule = new AliEvePMDModule();
	      lmodule->SetPosition(0.,0.,zpos);
	      lmodule->DisplayDigitsData(modnumber, pmdt);
              lmodule->SetPalette(pal);
              if (drawBackPolygons)
              {
                TEveFrameBox* b = lmodule->GetFrame();
                b->SetFrameWidth(1.5);
                b->SetFrameColor(1);
                b->SetBackColor ((Color_t) (kTeal - 9));
                b->SetFrameFill (kFALSE);
                b->SetDrawBack  (kTRUE);
              }
	      gEve->AddElement(lmodule, lddl);
	      modnumber++;
	      if (iddl == 4 && modnumber == 30) modnumber = 42;
	    }

	}

    }

  gEve->EnableRedraw();
}

// ---------------------------------------------------------------------- //
