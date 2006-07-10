// $Id$

//void
//muon_hits(const char *varexp    = "fX:fY:fZ",
//	 const char *selection = "",
//	 Option_t *option      = "goff")
void muon_hits(Int_t iShowCha = 1)
{
  char name[128];
  char title[128];

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("MUON");

  TTree* ht = rl->GetTreeH("MUON", false);
  //ht->Draw(varexp, selection, option);
  //ht->Print();

  for (Int_t iSta = 1; iSta <= 7; iSta++) {

  for (Int_t iCha = 1; iCha <= 2; iCha++) {

  Int_t iChamber = (iSta-1) * 2 + iCha; 

  if (iChamber != iShowCha) continue;

  sprintf(name,"M-ST%1dCH%1d/Hits",iSta,iCha);
  Reve::RenderElementList* l = new Reve::RenderElementList(name);
  if (iSta <= 5) {
    sprintf(title,"Station %1d chamber %1d (tracking)",iSta,iCha);
  } else {
    sprintf(title,"Station %1d chamber %1d (trigger)",iSta,iCha);
  }
  
  l->SetTitle(title);
  l->SetMainColor((Color_t)4);
  gReve->AddRenderElement(l);
      
  AliMpDEIterator ite;
  for ( ite.First(iChamber-1); ! ite.IsDone(); ite.Next() ) {
	
  Int_t detElemId = ite.CurrentDE();

  Int_t nTracks = ht->GetEntries();
  printf("Found %d tracks. \n",nTracks);
  for (Int_t it = 0; it < nTracks; it++) {

    TClonesArray *hits = 0;
    ht->SetBranchAddress("MUONHits",&hits);

    ht->GetEvent(it);
    
    Int_t nHits = hits->GetEntriesFast();
    //printf("Found %d hits in track %d. \n",nHits,it);

    Int_t nstep = 10;
    Reve::PointSet* points = new Reve::PointSet(nHits*nstep);
    sprintf(name,"Track %d: hits",it);
    points->SetName(name);
    sprintf(title,"TreeH entry %d",it);
    points->SetTitle(title);
    Int_t np = 0;    

    AliMUONHit *hit;
    Float_t x, y, z;
    for (Int_t ih = 0; ih < nHits; ih++) {
      
      hit = (AliMUONHit*)hits->UncheckedAt(ih);
      if (detElemId != hit->DetElemId()) continue;
      
      x = hit->X();
      y = hit->Y();
      z = hit->Z();
      
      Float_t dstep = 2.0*TMath::Pi() / (Float_t)nstep;
      Float_t d;
      Float_t r = 1.0;
      for (Int_t istep = 1; istep < (nstep+1); istep++) {
	
	d = istep * dstep;
	Float_t xc = x + r * TMath::Cos(d);
	Float_t yc = y + r * TMath::Sin(d);

	points->SetPoint(np,z,yc,xc);
	np++;
      
      }

    }
    
    if (np > 0) {
      gReve->AddRenderElement(l, points);      
    } else {
      delete points;
    }

  }

  }  // end detElemId loop

  }  // end chamber loop

  }  // end station loop

  gReve->Redraw3D();

  return;
}
