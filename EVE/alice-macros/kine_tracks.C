// Import tracks from kinematics-tree / particle-stack.
// Preliminary/minimal solution.
#include "TParticlePDG.h"

// PDG color indices
static Color_t DefCol   = 30;
static Color_t ECol     = 5;
static Color_t MuCol    = 6;
static Color_t GammaCol = 7; 
static Color_t MesCol1  = 3;
static Color_t MesCol2  = 38;
static Color_t BarCol   = 10;

Reve::TrackList*
kine_tracks(Double_t min_pt  = 0.1,   Double_t min_p   = 0.2,
	    Bool_t   pdg_col = kTRUE, Bool_t   recurse = kTRUE,
	    Bool_t   use_track_refs = kTRUE)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();
  if (!stack) {
    Error("kine_tracks.C", "can not get kinematics.");
    return 0;
  }

  gReve->DisableRedraw();
 
  Reve::TrackList* cont = new Reve::TrackList("Kine Tracks"); 
  cont->SetMainColor(Color_t(3));
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  // !!! Watch the '-', apparently different sign convention then for ESD.
  rnrStyle->SetMagField( - gAlice->Field()->SolenoidField() );

  gReve->AddRenderElement(cont);
  Int_t count = 0;
  Int_t N = stack->GetNtrack();
  for (Int_t i=0; i<N; ++i) 
  {
    if(stack->IsPhysicalPrimary(i)) 
    {
      TParticle* p = stack->Particle(i);
      if (p->Pt() < min_pt && p->P() < min_p) continue;

      ++count;
      Reve::Track* track = new Reve::Track(p, i, rnrStyle);
  
      //PH The line below is replaced waiting for a fix in Root
      //PH which permits to use variable siza arguments in CINT
      //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
      //PH    track->SetName(Form("%s [%d]", p->GetName(), i));
      char form[1000];
      sprintf(form,"%s [%d]", p->GetName(), i);
      track->SetName(form);
      if (pdg_col && p->GetPDG())
	track->SetMainColor(get_pdg_color(p->GetPDG()->PdgCode()));
      gReve->AddRenderElement(track, cont);

      if (recurse)
	kine_daughters(track, stack, min_pt, min_p, pdg_col, recurse);
    }
  }

  // set path marks
  Alieve::KineTools kt; 
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
  cont->UpdateItems();

  cont->MakeTracks(recurse);
  gReve->EnableRedraw();
  gReve->Redraw3D();

  return cont;
}

void kine_daughters(Reve::Track* parent,  AliStack* stack,
		    Double_t     min_pt,  Double_t  min_p,
		    Bool_t       pdg_col, Bool_t    recurse)
{
  TParticle *p = stack->Particle(parent->GetLabel());
  if (p->GetNDaughters() > 0) 
  {
    Reve::TrackRnrStyle* rs = parent->GetRnrStyle();
    for (int d=p->GetFirstDaughter(); d>0 && d<=p->GetLastDaughter(); ++d) 
    {	
      TParticle* dp = stack->Particle(d);
      if (dp->Pt() < min_pt && dp->P() < min_p) continue;

      Reve::Track* dtrack = new Reve::Track(dp, d, rs);  
      char form[1000];
      sprintf(form,"%s [%d]", dp->GetName(), d);
      dtrack->SetName(form);
      if (pdg_col && dp->GetPDG())
	dtrack->SetMainColor(get_pdg_color(dp->GetPDG()->PdgCode()));
      gReve->AddRenderElement(dtrack, parent);

      if (recurse)
	kine_daughters(dtrack, stack, min_pt, min_p, pdg_col, recurse);
    }
  }
}

Color_t get_pdg_color(Int_t pdg)
{
  Int_t pdga = TMath::Abs(pdg);
  Color_t col = Reve::DefCol;

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
  else if (pdga < 100000){ 
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


// Create mother and daughters tracks with given label.

Reve::RenderElement*
kine_track(Int_t  label,
	   Bool_t import_mother    = kTRUE,
           Bool_t import_daughters = kTRUE,
           Reve::RenderElement* cont = 0)

{
  if (label < 0) {
    Warning("kine_track", "label not set.");
    return 0;
  }
 
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();
  TParticle* p = stack->Particle(label);

  if (import_mother || (import_daughters && p->GetNDaughters()))
  {
    Track* toptrack = 0;
    TrackList* tracklist = 0;  
    TrackRnrStyle* rs = 0;

    if (cont == 0)
    {
      Reve::TrackList* tlist = new Reve::TrackList
	(Form("Kinematics of %d", label, p->GetNDaughters()));
      cont = tlist;

      Reve::TrackRnrStyle* rnrStyle = tlist->GetRnrStyle();
      // !!! Watch the '-', apparently different sign convention then for ESD.
      rnrStyle->SetMagField( - gAlice->Field()->SolenoidField() );
      char tooltip[1000];
      sprintf(tooltip,"Ndaughters=%d", p->GetNDaughters());
      tlist->SetTitle(tooltip);
      tlist->SelectByPt(0.2, 100);
      rnrStyle->fMaxOrbs = 8;
      rnrStyle->SetEditPathMarks(kTRUE);
      gReve->AddRenderElement(cont);
      rs = tlist->GetRnrStyle();
    }
    else
    {
      // check if argument is TrackList
      Reve::Track* t = dynamic_cast<Reve::Track*>(cont);
      if(t) 
      {
	rs = t->GetRnrStyle();
      }
      else
      {
        Reve::TrackList* l = dynamic_cast<Reve::TrackList*>(cont);
        if(l)
	{
	  rs = l->GetRnrStyle();
	}
        else
	{
	  Error("kine_tracks.C", "TrackRenderStyle not set.");
	}
      }
    }

    if (import_mother)
    {
      Track* track = new Reve::Track(p, label, rs);  
      char form[1000];
      sprintf(form,"%s [%d]", p->GetName(), label);
      track->SetName(form);
      gReve->AddRenderElement(track, cont);

    }

    if (import_daughters && p->GetNDaughters()) 
    {
      for (int d=p->GetFirstDaughter(); d>0 && d<=p->GetLastDaughter(); ++d) 
      {	
	TParticle* dp = stack->Particle(d);
	Track* track = new Reve::Track(dp, d, rs);  
	char form[1000];
	sprintf(form,"%s [%d]", dp->GetName(), d);
	track->SetName(form);
        track->MakeTrack();
	gReve->AddRenderElement(track, cont);
      }
    }
  }

  cont->UpdateItems();
  gReve->Redraw3D();
  return cont;
}

