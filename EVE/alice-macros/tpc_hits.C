// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEvePointSet*
tpc_hits(const char  *varexp    = "TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",
	 const char  *selection = "TPC2.fArray.fR>80",
         TEveElement *cont      = 0)
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadHits("TPC");

  TTree* ht = rl->GetTreeH("TPC", false);

  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  TEvePointSet* points = new TEvePointSet(Form("TPC Hits '%s'", selection));
  char form[1000];
  sprintf(form,"TPC Hits '%s'", selection);
  TEvePointSet* points = new TEvePointSet(form);
  points->SetSourceCS(TEvePointSelectorConsumer::kTVT_RPhiZ);

  TEvePointSelector ps(ht, points, varexp, selection);
  ps.Select();

  if (points->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
    Warning("tpc_hits", Form("No hits match '%s'", selection));
    delete points;
    return 0;
  }

  //PH  points->SetTitle(Form("N=%d", points->Size()));
  sprintf(form,"N=%d", points->Size());
  points->SetTitle(form);
  points->SetMarkerSize(.5);
  points->SetMarkerColor((Color_t)3);

  gEve->AddElement(points, cont);
  gEve->Redraw3D();

  return points;
}
