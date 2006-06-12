// $Id$

Reve::Track* esd_make_track(Reve::TrackRnrStyle* rnrStyle,
			    AliESDtrack* at,
			    AliExternalTrackParam* tp=0)
{
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

Bool_t gkFixFailedITSExtr = kTRUE;

Reve::TrackList* esd_tracks(Double_t min_pt=0.1, Double_t max_pt=100)
{
  AliESD* esd = Alieve::Event::AssertESD();
  AliPID  pid;

  Double_t minptsq = min_pt*min_pt;
  Double_t maxptsq = max_pt*max_pt;
  Double_t ptsq;

  TGListTree*     l_tree = gReve->GetListTree();
  TGListTreeItem* parent = gReve->GetEventTreeItem();

  Reve::TrackList* cont = new Reve::TrackList("ESD Tracks"); 
  cont->SetMainColor(Color_t(6));
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->SetMagField( esd->GetMagneticField() );

  TGListTreeItem *holder = gReve->AddRenderElement(cont);

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

    // If ITS refit failed, take track parameters at inner TPC radius.
    AliExternalTrackParam* tp = at;
    if (gkFixFailedITSExtr && !at->IsOn(AliESDtrack::kITSrefit)) {
       tp = at->GetInnerParam();
    }

    Reve::Track* track = esd_make_track(rnrStyle, at, tp);
    cont->AddElement(track);
    gReve->AddRenderElement(holder, track);
  }

  const Text_t* tooltip = Form("pT ~ (%.2lf, %.2lf), N=%d", min_pt, max_pt, count);
  cont->SetTitle(tooltip); // Not broadcasted automatically ...
  holder->SetTipText(tooltip);

  cont->MakeTracks();
  cont->MakeMarkers();
  gReve->DrawRenderElement(cont);

  return cont;
}

Reve::TrackList* esd_tracks_from_array(TCollection* col, AliESD* esd=0)
{
  // Retrieves AliESDTrack's from collection.
  // See example usage with AliAnalysisTrackCuts in the next function.

  if(esd == 0) esd = Alieve::Event::AssertESD();
  AliPID  pid;

  TGListTree*     l_tree = gReve->GetListTree();
  TGListTreeItem* parent = gReve->GetEventTreeItem();

  Reve::TrackList* cont = new Reve::TrackList("ESD Tracks"); 
  cont->SetMainColor(Color_t(6));
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->SetMagField( esd->GetMagneticField() );

  TGListTreeItem *holder = gReve->AddRenderElement(cont);

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
    cont->AddElement(track);
    gReve->AddRenderElement(holder, track);
  }

  const Text_t* tooltip = Form("N=%d", count);
  cont->SetTitle(tooltip); // Not broadcasted automatically ...
  holder->SetTipText(tooltip);

  cont->MakeTracks();
  cont->MakeMarkers();
  gReve->DrawRenderElement(cont);

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
