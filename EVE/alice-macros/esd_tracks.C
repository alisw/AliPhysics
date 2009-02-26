// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Use inner-tpc track params when its refit failed.
Bool_t g_esd_tracks_use_ip_on_failed_its_refit = kFALSE;

TString esd_track_title(AliESDtrack* t)
{
  TString s;

  Int_t label = t->GetLabel(), index = t->GetID();
  TString idx(index == kMinInt ? "<undef>" : Form("%d", index));
  TString lbl(label == kMinInt ? "<undef>" : Form("%d", label));

  Double_t p[3], v[3];
  t->GetXYZ(v);
  t->GetPxPyPz(p);

  s = Form("Index=%s, Label=%s\nChg=%d, Pdg=%d\n"
	   "pT=%.3f, pZ=%.3f\nV=(%.3f, %.3f, %.3f)\n",
	   idx.Data(), lbl.Data(), t->Charge(), 0,
	   TMath::Sqrt(p[0]*p[0] + p[1]*p[1]), p[2], v[0], v[1], v[2]);

  Int_t   o;
  s += "Det(in,out,refit,pid):\n";
  o  = AliESDtrack::kITSin;
  s += Form("ITS (%d,%d,%d,%d)  ",  t->IsOn(o), t->IsOn(o<<1), t->IsOn(o<<2), t->IsOn(o<<3));
  o  = AliESDtrack::kTPCin;
  s += Form("TPC(%d,%d,%d,%d)\n",   t->IsOn(o), t->IsOn(o<<1), t->IsOn(o<<2), t->IsOn(o<<3));
  o  = AliESDtrack::kTRDin;
  s += Form("TRD(%d,%d,%d,%d) ",    t->IsOn(o), t->IsOn(o<<1), t->IsOn(o<<2), t->IsOn(o<<3));
  o  = AliESDtrack::kTOFin;
  s += Form("TOF(%d,%d,%d,%d)\n",   t->IsOn(o), t->IsOn(o<<1), t->IsOn(o<<2), t->IsOn(o<<3));
  o  = AliESDtrack::kHMPIDout;
  s += Form("HMPID(out=%d,pid=%d)\n", t->IsOn(o), t->IsOn(o<<1));
  s += Form("ESD pid=%d", t->IsOn(AliESDtrack::kESDpid));

  if (t->IsOn(AliESDtrack::kESDpid))
  {
    Double_t pid[5];
    t->GetESDpid(pid);
    s += Form("\n[%.2f %.2f %.2f %.2f %.2f]", pid[0], pid[1], pid[2], pid[3], pid[4]);
  }

  return s;
}

void esd_track_add_param(AliEveTrack* track, AliExternalTrackParam* tp)
{
  // Add additional track parameters as a path-mark to track.

  if (tp == 0)
    return;

  Double_t pbuf[3], vbuf[3];
  tp->GetXYZ(vbuf);
  tp->GetPxPyPz(pbuf);

  TEvePathMark pm(TEvePathMark::kReference);
  pm.fV.Set(vbuf);
  pm.fP.Set(pbuf);
  track->AddPathMark(pm);
}

TEveTrack* esd_make_track(AliESDtrack *at, TEveTrackList* cont)
{
  // Make a standard track representation and put it into given container.

  // Choose which parameters to use a track's starting point.
  // If gkFixFailedITSExtr is TRUE (FALSE by default) and
  // if ITS refit failed, take track parameters at inner TPC radius.
  AliExternalTrackParam* tp = at;

  Bool_t innerTaken = kFALSE;
  if ( ! at->IsOn(AliESDtrack::kITSrefit) && g_esd_tracks_use_ip_on_failed_its_refit)
  {
    tp = at->GetInnerParam();
    innerTaken = kTRUE;
  }

  Double_t pbuf[3], vbuf[3];
  tp->GetXYZ(vbuf);
  tp->GetPxPyPz(pbuf);

  TEveRecTrack rt;
  rt.fLabel  = at->GetLabel();
  rt.fIndex  = at->GetID();
  rt.fStatus = (Int_t) at->GetStatus();
  rt.fSign   = at->GetSign();
  rt.fV.Set(vbuf);
  rt.fP.Set(pbuf);
  Double_t ep = at->GetP(), mc = at->GetMass();
  rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);

  AliEveTrack* track = new AliEveTrack(&rt, cont->GetPropagator());
  track->SetAttLineAttMarker(cont);
  track->SetName(Form("AliEveTrack %d", at->GetID()));
  track->SetElementTitle(esd_track_title(at));
  track->SetSourceObject(at);

  // Add inner/outer track parameters as path-marks.
  if (at->IsOn(AliESDtrack::kTPCrefit))
  {
    if ( ! innerTaken)
    {
      esd_track_add_param(track, at->GetInnerParam());
    }
    esd_track_add_param(track, at->GetOuterParam());
  }

  return track;
}

TEveTrackList* esd_tracks()
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();

  TEveTrackList* cont = new TEveTrackList("ESD Tracks");
  cont->SetMainColor(6);
  TEveTrackPropagator* trkProp = cont->GetPropagator();
  trkProp->SetMagField(0.1*esd->GetMagneticField());
  trkProp->SetMaxR    (520);

  gEve->AddElement(cont);

  Int_t count = 0;
  for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
  {
    ++count;
    AliEveTrack* track = esd_make_track(esd->GetTrack(n), cont);

    cont->AddElement(track);
  }
  cont->SetTitle(Form("N=%d", count));
  cont->MakeTracks();

  gEve->Redraw3D();

  return cont;
}

TEveElementList* esd_tracks_MI()
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();

  TEveElementList* cont = new TEveElementList("ESD Tracks MI");
  gEve->AddElement(cont);

  Int_t count = 0;
  for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
  {
    ++count;
    AliESDtrack* at = esd->GetTrack(n);
    TEveLine* l = new TEveLine; 
    l->SetLineColor(5);
    at->FillPolymarker(l, esd->GetMagneticField(), 0, 250, 5);
    l->SetName(Form("ESDTrackMI %d", at->GetID()));
    l->SetElementTitle(esd_track_title(at));
    l->SetSourceObject(at);

    cont->AddElement(l);
  }
  cont->SetTitle(Form("N=%d", count));

  gEve->Redraw3D();

  return cont;
}

/******************************************************************************/
// esd_tracks_from_array()
/******************************************************************************/

TEveTrackList* esd_tracks_from_array(TCollection* col, AliESDEvent* esd=0)
{
  // Retrieves AliESDTrack's from collection.
  // See example usage with AliAnalysisTrackCuts in the next function.

  if (esd == 0) esd = AliEveEventManager::AssertESD();

  TEveTrackList* cont = new TEveTrackList("ESD Tracks");
  cont->SetMainColor(6);
  TEveTrackPropagator* trkProp = cont->GetPropagator();
  trkProp->SetMagField(0.1*esd->GetMagneticField());
  trkProp->SetMaxR    (520);

  gEve->AddElement(cont);

  Int_t    count = 0;
  TIter    next(col);
  TObject *obj;
  while ((obj = next()) != 0)
  {
    if (obj->IsA()->InheritsFrom("AliESDtrack") == kFALSE)
    {
      Warning("esd_tracks_from_array", "Object '%s', '%s' is not an AliESDtrack.",
	      obj->GetName(), obj->GetTitle());
      continue;
    }

    ++count;
    AliESDtrack* at = (AliESDtrack*) obj;

    AliEveTrack* track = esd_make_track(at, cont);
    cont->AddElement(track);
  }
  cont->SetTitle(Form("N=%d", count));
  cont->MakeTracks();

  gEve->Redraw3D();

  return cont;
}

void esd_tracks_alianalcuts_demo()
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();
  gSystem->Load("libANALYSIS");

  AliAnalysisTrackCuts atc;
  atc.SetPtRange(0.1, 5);
  atc.SetRapRange(-1, 1);

  esd_tracks_from_array(atc.GetAcceptedParticles(esd), esd);
}

/******************************************************************************/
// esd_tracks_by_category
/******************************************************************************/

Float_t get_sigma_to_vertex(AliESDtrack* esdTrack)
{
  // Taken from: PWG0/esdTrackCuts/AliESDtrackCuts.cxx
  // Float_t AliESDtrackCuts::GetSigmaToVertex(AliESDtrack* esdTrack)

  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);
  if (bCov[0] <= 0 || bCov[2] <= 0)
  {
    printf("Estimated b resolution lower or equal zero!\n");
    bCov[0] = bCov[2] = 0;
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

  if (bRes[0] == 0 || bRes[1] == 0)
    return -1;

  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  // stupid rounding problem screws up everything:
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  if (TMath::Exp(-d * d / 2) < 1e-10)
    return 1000;

  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  return d;
}

TEveElementList* g_esd_tracks_by_category_container = 0;

TEveElementList* esd_tracks_by_category()
{
  // Import ESD tracks, separate them into several containers
  // according to primary-vertex cut and ITS&TPC refit status.

  AliESDEvent* esd = AliEveEventManager::AssertESD();

  TEveElementList* cont = new TEveElementList("ESD Tracks by category");
  gEve->AddElement(cont);

  g_esd_tracks_by_category_container = cont;

  const Int_t   nCont = 7;
  const Float_t maxR  = 520;
  const Float_t magF  = 0.1*esd->GetMagneticField();

  TEveTrackList *tl[nCont];
  Int_t          tc[nCont];
  Int_t          count = 0;

  tl[0] = new TEveTrackList("Sigma < 3");
  tc[0] = 0;
  tl[0]->GetPropagator()->SetMagField(magF);
  tl[0]->GetPropagator()->SetMaxR    (maxR);
  tl[0]->SetMainColor(3);
  cont->AddElement(tl[0]);

  tl[1] = new TEveTrackList("3 < Sigma < 5");
  tc[1] = 0;
  tl[1]->GetPropagator()->SetMagField(magF);
  tl[1]->GetPropagator()->SetMaxR    (maxR);
  tl[1]->SetMainColor(7);
  cont->AddElement(tl[1]);

  tl[2] = new TEveTrackList("5 < Sigma");
  tc[2] = 0;
  tl[2]->GetPropagator()->SetMagField(magF);
  tl[2]->GetPropagator()->SetMaxR    (maxR);
  tl[2]->SetMainColor(46);
  cont->AddElement(tl[2]);

  tl[3] = new TEveTrackList("no ITS refit; Sigma < 5");
  tc[3] = 0;
  tl[3]->GetPropagator()->SetMagField(magF);
  tl[3]->GetPropagator()->SetMaxR    (maxR);
  tl[3]->SetMainColor(41);
  cont->AddElement(tl[3]);

  tl[4] = new TEveTrackList("no ITS refit; Sigma > 5");
  tc[4] = 0;
  tl[4]->GetPropagator()->SetMagField(magF);
  tl[4]->GetPropagator()->SetMaxR    (maxR);
  tl[4]->SetMainColor(48);
  cont->AddElement(tl[4]);

  tl[5] = new TEveTrackList("no TPC refit");
  tc[5] = 0;
  tl[5]->GetPropagator()->SetMagField(magF);
  tl[5]->GetPropagator()->SetMaxR    (maxR);
  tl[5]->SetMainColor(kRed);
  cont->AddElement(tl[5]);

  tl[6] = new TEveTrackList("ITS stand-alone");
  tc[6] = 0;
  tl[6]->GetPropagator()->SetMagField(magF);
  tl[6]->GetPropagator()->SetMaxR    (maxR);
  tl[6]->SetMainColor(kMagenta - 9);
  cont->AddElement(tl[6]);

  for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
  {
    AliESDtrack* at = esd->GetTrack(n);

    Float_t s  = get_sigma_to_vertex(at);
    Int_t   ti;
    if      (s <  3) ti = 0;
    else if (s <= 5) ti = 1;
    else             ti = 2;

    if (at->IsOn(AliESDtrack::kITSin) && ! at->IsOn(AliESDtrack::kTPCin))
    {
      ti = 6;
    }
    else if (at->IsOn(AliESDtrack::kTPCin) && ! at->IsOn(AliESDtrack::kTPCrefit))
    {
      ti = 5;
    }
    else if ( ! at->IsOn(AliESDtrack::kITSrefit))
    {
      ti = (ti == 2) ? 4 : 3;
    }

    TEveTrackList* tlist = tl[ti];
    ++tc[ti];
    ++count;

    AliEveTrack* track = esd_make_track(at, tlist);
    track->SetName(Form("AliEveTrack idx=%d, sigma=%5.3f", at->GetID(), s));
    tlist->AddElement(track);
  }

  for (Int_t ti = 0; ti < nCont; ++ti)
  {
    TEveTrackList* tlist = tl[ti];
    tlist->SetName(Form("%s [%d]", tlist->GetName(), tlist->NumChildren()));
    tlist->SetTitle(Form("N tracks=%d", tc[ti]));

    tlist->MakeTracks();
  }
  cont->SetTitle(Form("N all tracks = %d", count));
  // ??? The following does not always work:
  cont->FindListTreeItem(gEve->GetListTree())->SetOpen(kTRUE);

  gEve->Redraw3D();

  return cont;
}
