// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
// ---------------------------------------------------------------------- //

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TObjArray.h>
#include <TString.h>
#include <TStyle.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEvePointSet.h>

#include <AliPMDRawStream.h>
#include <AliRawReaderRoot.h>
#include <AliCluster.h>
#include <AliEveEventManager.h>
#include <AliEvePMDModule.h>
#endif

void pmd_raw(Int_t mode = 0)
{
  gStyle->SetPalette(1, 0);


  TObjArray *pmdddlcont = new TObjArray();

  TString spl;

  TString sddl;
  TString bsddl="DDL";

  // Use this to get data consistent with current event:
  // AliRawReader *reader = AliEveEventManager::AssertRawReader();

  Int_t ievt = 159;
  AliRawReaderRoot reader("raw.root",ievt);
  AliPMDRawStream stream(&reader);

  gEve->DisableRedraw();

  TEveElementList* l = new TEveElementList("PMD");
  //  l->SetTitle("PMD");
  //  l->SetMainColor(3);
  gEve->AddElement(l);

  Int_t NSM       = 0;
  Int_t istartDDL = 0;
  Int_t iendDDL   = 0;
  Int_t modnumber = 0;
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
      printf("--- Visualization is set for both the Plane ---\n");
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
	  spl = "CPV";
	  istartDDL = 4;
	  iendDDL   = 6;
	  zpos      = 360.;
	}

      TEveElementList* lplane = new TEveElementList(spl.Data());
      //  l->SetMainColor(3);
      gEve->AddElement(lplane, l);

      Int_t iddl = -1;

      while ((iddl = stream.DdlData(pmdddlcont)) >=0) {
	  if (iddl >= istartDDL && iddl < iendDDL){
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
		  lmodule->DisplayRawData(modnumber,pmdddlcont);
		  gEve->AddElement(lmodule, lddl);
		  modnumber++;
		  if (iddl == 4 && modnumber == 30) modnumber = 42;
	      }
	      
	      pmdddlcont->Delete();
	  }
      }
      
      gEve->EnableRedraw();
    }
  
}
// ---------------------------------------------------------------------- //
  
