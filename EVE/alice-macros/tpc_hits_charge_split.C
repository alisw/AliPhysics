// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveTreeTools.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif

void tpc_hits_charge_split(const char *varexp    =
			"TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ"
			":log(TPC2.fArray.fCharge)",
			const char *selection = "TPC2.fArray.fR>80")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("TPC");

  TTree* ht = rl->GetTreeH("TPC", false);

  TEvePointSetArray* l = new TEvePointSetArray("TPC hits - Log-Charge Slices", "");
  l->SetSourceCS(TEvePointSelectorConsumer::kTVT_RPhiZ);
  l->SetMarkerColor(3);
  l->SetMarkerStyle(20); // full circle
  l->SetMarkerSize(.5);

  gEve->AddElement(l);
  l->InitBins("Log Charge", 20, 0, 5);

  TEvePointSelector ps(ht, l, varexp, selection);
  ps.Select();

  l->CloseBins();

  gEve->Redraw3D();
}
