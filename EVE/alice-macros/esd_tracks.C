// $Id$

Reve::Track* esd_make_track(Reve::TrackRnrStyle* rnrStyle,
			    AliESDtrack* at,
			    AliExternalTrackParam* tp=0)
{
  // Helper function
  Double_t        pbuf[3], vbuf[3];
  Reve::RecTrack  rt;

  if(tp == 0) tp = at;

  rt.label  =         at->GetLabel();
  rt.status = (Int_t) at->GetStatus();
  rt.sign   = tp->GetSign();
  tp->GetXYZ(vbuf);
  rt.V.Set(vbuf);
  tp->GetPxPyPz(pbuf);
  rt.P.Set(pbuf);
  Double_t ep = at->GetP(), mc = at->GetMass();
  rt.beta = ep/TMath::Sqrt(ep*ep + mc*mc);
 
  Reve::Track* track = new Reve::Track(&rt, rnrStyle);
  track->SetName(Form("ESDTrack %d", rt.label));
  track->SetTitle(Form("pT=%.3f, pZ=%.3f; V=(%.3f, %.3f, %.3f)",
		       rt.sign*TMath::Hypot(rt.P.x, rt.P.y), rt.P.z,
		       rt.V.x, rt.V.y, rt.V.z));
  return track;
}

Bool_t gkFixFailedITSExtr = kFALSE;

Reve::TrackList* esd_tracks(Double_t min_pt=0.1, Double_t max_pt=100)
{
  AliESD* esd = Alieve::Event::AssertESD();

  Double_t minptsq = min_pt*min_pt;
  Double_t maxptsq = max_pt*max_pt;
  Double_t ptsq;

  Reve::TrackList* cont = new Reve::TrackList("ESD Tracks"); 
  cont->SetMainColor(Color_t(6));
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->SetMagField( esd->GetMagneticField() );

  gReve->AddRenderElement(cont);

  Int_t    count = 0;
  Double_t pbuf[3];
  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++) {
    AliESDtrack* at = esd->GetTrack(n);

    // Here would be sweet to have TObjectFormula.
    at->GetPxPyPz(pbuf);
    ptsq = pbuf[0]*pbuf[0] + pbuf[1]*pbuf[1];
    if(ptsq < minptsq || ptsq > maxptsq)
      continue;

    ++count;

    // If gkFixFailedITSExtr is TRUE (FALSE by default) and
    // if ITS refit failed, take track parameters at inner TPC radius.
    AliExternalTrackParam* tp = at;
    if (gkFixFailedITSExtr && !at->IsOn(AliESDtrack::kITSrefit)) {
       tp = at->GetInnerParam();
    }

    Reve::Track* track = esd_make_track(rnrStyle, at, tp);
    gReve->AddRenderElement(cont, track);
  }

  const Text_t* tooltip = Form("pT ~ (%.2lf, %.2lf), N=%d", min_pt, max_pt, count);
  cont->SetTitle(tooltip); // Not broadcasted automatically ...
  cont->UpdateItems();

  cont->MakeTracks();
  cont->MakeMarkers();
  gReve->Redraw3D();

  return cont;
}

/**************************************************************************/
// esd_tracks_from_array()
/**************************************************************************/

Reve::TrackList* esd_tracks_from_array(TCollection* col, AliESD* esd=0)
{
  // Retrieves AliESDTrack's from collection.
  // See example usage with AliAnalysisTrackCuts in the next function.

  if(esd == 0) esd = Alieve::Event::AssertESD();

  Reve::TrackList* cont = new Reve::TrackList("ESD Tracks"); 
  cont->SetMainColor(Color_t(6));
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->SetMagField( esd->GetMagneticField() );

  gReve->AddRenderElement(cont);

  Int_t    count = 0;
  TIter    next(col);
  TObject *obj;
  while((obj = next()) != 0) {
    if(obj->IsA()->InheritsFrom("AliESDtrack") == kFALSE) {
      Warning("Object '%s', '%s' is not an AliESDtrack.",
	      obj->GetName(), obj->GetTitle());
      continue;
    }

    ++count;
    AliESDtrack* at = (AliESDtrack*) obj;

    Reve::Track* track = esd_make_track(rnrStyle, at);
    gReve->AddRenderElement(cont, track);
  }

  const Text_t* tooltip = Form("N=%d", count);
  cont->SetTitle(tooltip); // Not broadcasted automatically ...
  cont->UpdateItems();

  cont->MakeTracks();
  cont->MakeMarkers();
  gReve->Redraw3D();

  return cont;
}

void esd_tracks_alianalcuts_demo()
{
  AliESD* esd = Alieve::Event::AssertESD();
  gSystem->Load("libANALYSIS");

  AliAnalysisTrackCuts atc;
  atc.SetPtRange(0.1, 5);
  atc.SetRapRange(-1, 1);

  esd_tracks_from_array(atc.GetAcceptedParticles(esd), esd);
}

/**************************************************************************/
// esd_tracks_vertex_cut
/**************************************************************************/

Float_t get_sigma_to_vertex(AliESDtrack* esdTrack)
{
  // Taken from: PWG0/esdTrackCuts/AliESDtrackCuts.cxx
  // Float_t AliESDtrackCuts::GetSigmaToVertex(AliESDtrack* esdTrack)

  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // -----------------------------------
  // How to get to a n-sigma cut?
  //
  // The accumulated statistics from 0 to d is
  //
  // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
  // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
  //
  // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-x**2)/2)
  // Can this be expressed in a different way?

  if (bRes[0] == 0 || bRes[1] ==0)
    return -1;

  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  // stupid rounding problem screws up everything:
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  if (TMath::Exp(-d * d / 2) < 1e-10)
    return 1000;

  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  return d;
}

Reve::RenderElementList* esd_tracks_vertex_cut()
{
  // Import ESD tracks, separate them into three containers according to
  // primary-vertex cut.

  AliESD* esd = Alieve::Event::AssertESD();

  Reve::RenderElementList* cont = new Reve::RenderElementList("ESD Tracks", 0, kTRUE);
  gReve->AddRenderElement(cont);
  Reve::TrackList *tl[4];
  Int_t            tc[4];
  Int_t            count = 0;

  tl[0] = new Reve::TrackList("Sigma < 3");
  printf("%p %s\n", tl[0], tl[0]->IsA()->GetName());
  tc[0] = 0;
  tl[0]->GetRnrStyle()->SetMagField( esd->GetMagneticField() );
  tl[0]->SetMainColor(Color_t(3));
  gReve->AddRenderElement(cont, tl[0]);

  tl[1] = new Reve::TrackList("3 < Sigma < 5");
  tc[1] = 0;
  tl[1]->GetRnrStyle()->SetMagField( esd->GetMagneticField() );
  tl[1]->SetMainColor(Color_t(7));
  gReve->AddRenderElement(cont, tl[1]);

  tl[2] = new Reve::TrackList("5 < Sigma");
  tc[2] = 0;
  tl[2]->GetRnrStyle()->SetMagField( esd->GetMagneticField() );
  tl[2]->SetMainColor(Color_t(46));
  gReve->AddRenderElement(cont, tl[2]);

  tl[3] = new Reve::TrackList("ITS refit failed");
  tc[3] = 0;
  tl[3]->GetRnrStyle()->SetMagField( esd->GetMagneticField() );
  tl[3]->SetMainColor(Color_t(41));
  gReve->AddRenderElement(cont, tl[3]);

  for (Int_t n=0; n<esd->GetNumberOfTracks(); n++) {
    AliESDtrack* at = esd->GetTrack(n);

    Float_t s  = get_sigma_to_vertex(at);
    Int_t   ti;
    if      (s <  3) ti = 0;
    else if (s <= 5) ti = 1;
    else             ti = 2;

    AliExternalTrackParam* tp = at;
    // If ITS refit failed, take track parameters at inner TPC radius and
    // put track in a special container.
    // This ignores state of gkFixFailedITSExtr (used in esd_tracks()).
    // Use BOTH functions to compare results.
    if (!at->IsOn(AliESDtrack::kITSrefit)) {
      tp = at->GetInnerParam();
      ti = 3;
    }

    Reve::TrackList* tlist = tl[ti];
    ++tc[ti];
    ++count;

    Reve::Track* track = esd_make_track(tlist->GetRnrStyle(), at, tp);
    track->SetName(Form("track %d, sigma=%5.3f", at->GetLabel(), s));
    gReve->AddRenderElement(tlist, track);
  }

  for (Int_t ti=0; ti<4; ++ti) {
    Reve::TrackList* tlist = tl[ti];
    const Text_t* tooltip = Form("N tracks=%d", tc[ti]);
    tlist->SetTitle(tooltip); // Not broadcasted automatically ...
    tlist->UpdateItems();

    tlist->MakeTracks();
    tlist->MakeMarkers();
  }
  cont->SetTitle(Form("N all tracks = %d", count));
  cont->UpdateItems();
  gReve->Redraw3D();

  return cont;
}
