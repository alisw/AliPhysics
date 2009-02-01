// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Import tracks from kinematics-tree / particle-stack.
// Preliminary/minimal solution.

#include "TParticlePDG.h"

TEveTrackList*
kine_tracks(Double_t min_pt  = 0,     Double_t min_p   = 0,
	    Bool_t   pdg_col = kTRUE, Bool_t   recurse = kTRUE,
	    Bool_t   use_track_refs = kTRUE)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();
  if (!stack) {
    Error("kine_tracks.C", "can not get kinematics.");
    return 0;
  }

  gEve->DisableRedraw();

  TEveTrackList* cont = new TEveTrackList("Kine Tracks");
  cont->SetMainColor(3);
  TEveTrackPropagator* rnrStyle = cont->GetPropagator();
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  rnrStyle->SetMagField(fld ? -0.1*fld->SolenoidField() : 0);

  gEve->AddElement(cont);
  Int_t count = 0;
  Int_t N = stack->GetNtrack();
  Int_t Np = stack->GetNprimary();
  for (Int_t i = 0; i < Np; ++i)
  {
    TParticle* p = stack->Particle(i);
    if (p->GetStatusCode() <= 1)
    {
      if (p->Pt() < min_pt && p->P() < min_p) continue;

      ++count;
      TEveTrack* track = new TEveTrack(p, i, rnrStyle);

      //PH The line below is replaced waiting for a fix in Root
      //PH which permits to use variable siza arguments in CINT
      //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
      //PH    track->SetName(Form("%s [%d]", p->GetName(), i));
      char form[1000];
      sprintf(form,"%s [%d]", p->GetName(), i);
      track->SetName(form);
      track->SetStdTitle();
      Int_t ml = p->GetMother(0);
      if (ml != -1)
      {
        track->SetTitle(Form("%s\nMother label=%d\nMother Pdg=%d",
                             track->GetElementTitle(),
                             ml, stack->Particle(ml)->GetPdgCode()));
      }
      set_track_color(track, pdg_col);

      gEve->AddElement(track, cont);

      if (recurse)
	kine_daughters(track, stack, min_pt, min_p, pdg_col, recurse);
    }
  }

  // set path marks
  AliEveKineTools kt;
  kt.SetDaughterPathMarks(cont, stack, recurse);
  if (use_track_refs && rl->LoadTrackRefs() == 0)
  {
    kt.SetTrackReferences(cont, rl->TreeTR(), recurse);
    rnrStyle->SetEditPathMarks(kTRUE);
  }
  kt.SortPathMarks(cont, recurse);

  //PH  const Text_t* tooltip = Form("min pT=%.2lf, min P=%.2lf), N=%d", min_pt, min_p, count);
  char tooltip[1000];
  sprintf(tooltip,"min pT=%.2lf, min P=%.2lf), N=%d", min_pt, min_p, count);
  cont->SetTitle(tooltip); // Not broadcasted automatically ...

  cont->MakeTracks(recurse);
  gEve->EnableRedraw();
  gEve->Redraw3D();

  return cont;
}

void kine_daughters(TEveTrack* parent,  AliStack* stack,
		    Double_t     min_pt,  Double_t  min_p,
		    Bool_t       pdg_col, Bool_t    recurse)
{
  TParticle *p = stack->Particle(parent->GetLabel());
  if (p->GetNDaughters() > 0)
  {
    TEveTrackPropagator* rs = parent->GetPropagator();
    for (int d=p->GetFirstDaughter(); d>0 && d<=p->GetLastDaughter(); ++d)
    {
      TParticle* dp = stack->Particle(d);
      if (dp->Pt() < min_pt && dp->P() < min_p) continue;

      TEveTrack* dtrack = new TEveTrack(dp, d, rs);
      char form[1000];
      sprintf(form,"%s [%d]", dp->GetName(), d);
      dtrack->SetName(form);
      dtrack->SetStdTitle();
      set_track_color(dtrack, pdg_col);

      gEve->AddElement(dtrack, parent);

      if (recurse)
	kine_daughters(dtrack, stack, min_pt, min_p, pdg_col, recurse);
    }
  }
}

Color_t set_track_color(TEveTrack* t, Bool_t pdg_col)
{
  if (pdg_col)
    t->SetMainColor(get_pdg_color(t->GetPdg()));
  else
    t->SetMainColor(30);
}

Color_t get_pdg_color(Int_t pdg)
{
  // PDG color indices
  static const Color_t DefCol   = 30;
  static const Color_t ECol     = 5;
  static const Color_t MuCol    = 6;
  static const Color_t GammaCol = 7;
  static const Color_t MesCol1  = 3;
  static const Color_t MesCol2  = 38;
  static const Color_t BarCol   = 10;

  Int_t pdga  = TMath::Abs(pdg);
  Color_t col = DefCol;

  // elementary  particles
  if (pdga < 100) {
    switch (pdga) {
      case 11:
	col = ECol; break;
      case 12:
	col = MuCol; break;
      case 22:
	col = GammaCol; break;
    }
  }
  // mesons and barions
  else if (pdga < 100000) {
    Int_t i  = pdga;
    Int_t i0 = i%10; i /= 10;
    Int_t i1 = i%10; i /= 10;
    Int_t i2 = i%10; i /= 10;
    Int_t i3 = i%10; i /= 10;
    Int_t i4 = i%10;
    //printf("pdg(%d) quark indices (%d,%d,%d,%d,%d) \n",pdg, i4,i3,i2, i1, i0);
    // meson
    if ((i3 == 0) && ( i4 < 2)){
      col = MesCol1; // quarks: i1,i2 (spin = i0)
      if(i1 == 3 || i2 == 3)
	col = MesCol2;
    } // barion
    else if ( i2 >= i1 && i3 >= i2 ) {
      col = BarCol; // quarks: i1,i2, i3 (spin = i0))
    }
  }

  return col;
}

/******************************************************************************/

TEveElement*
kine_track(Int_t  label,
	   Bool_t import_mother    = kTRUE, Bool_t import_daughters = kTRUE,
	   Bool_t pdg_col          = kTRUE, Bool_t recurse          = kTRUE,
           TEveElement* cont = 0)

{
  // Create mother and daughters tracks with given label.
  // mother     -> particle with label
  // daughters  -> daughters of label

  if (label < 0) {
    Warning("kine_track", "label not set.");
    return 0;
  }

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();
  if (label >= stack->GetNtrack())
  {
    Warning("kine_track", "label out of range.");
    return 0;
  }

  TParticle* p = stack->Particle(label);

  if (import_mother || (import_daughters && p->GetNDaughters()))
  {
    TEveTrack* toptrack = 0;
    TEveTrackList* tracklist = 0;
    TEveTrackPropagator* rs = 0;

    if (cont == 0)
    {
      TEveTrackList* tlist = new TEveTrackList
	(Form("Kinematics of %d", label, p->GetNDaughters()));
      cont = tlist;

      TEveTrackPropagator* rnrStyle = tlist->GetPropagator();
      // !!! Watch the '-', apparently different sign convention then for ESD.
      AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
      rnrStyle->SetMagField( fld ? -0.1*fld->SolenoidField() : 0 );
      char tooltip[1000];
      sprintf(tooltip,"Ndaughters=%d", p->GetNDaughters());
      tlist->SetTitle(tooltip);
      rnrStyle->fMaxOrbs = 2;
      rnrStyle->SetEditPathMarks(kTRUE);

      gEve->AddElement(cont);
      rs = tlist->GetPropagator();
    }
    else
    {
      // check if container is TEveTrackList or TEveTrack (has rnr-style)
      TEveTrack* t = dynamic_cast<TEveTrack*>(cont);
      if (t) {
	rs = t->GetPropagator();
      } else {
        TEveTrackList* l = dynamic_cast<TEveTrackList*>(cont);
        if (l)
	  rs = l->GetPropagator();
        else
	  Error("kine_tracks.C", "TrackRenderStyle not set.");
      }
    }

    if (import_mother)
    {
      TEveTrack* track = new TEveTrack(p, label, rs);
      char form[1000];
      sprintf(form,"%s [%d]", p->GetName(), label);
      track->SetName(form);
      track->SetStdTitle();
      set_track_color(track, pdg_col);

      track->MakeTrack();
      gEve->AddElement(track, cont);
      cont = track;
    }

    if (import_daughters && p->GetNDaughters())
    {
      for (int d=p->GetFirstDaughter(); d>0 && d<=p->GetLastDaughter(); ++d)
      {
	TParticle* dp = stack->Particle(d);
	TEveTrack* track = new TEveTrack(dp, d, rs);
	char form[1000];
	sprintf(form,"%s [%d]", dp->GetName(), d);
	track->SetName(form);
	track->SetStdTitle();
	set_track_color(track, pdg_col);

        track->MakeTrack();
	gEve->AddElement(track, cont);

	if (recurse)
	  kine_daughters(track, stack, 0, 0, pdg_col, recurse);
      }
    }
  }

  gEve->Redraw3D();
  return cont;
}

