// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void trackcount_init()
{
  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init(".", -1);

  TEveUtil::LoadMacro("primary_vertex.C");
  TEveUtil::LoadMacro("esd_tracks.C");
  TEveUtil::LoadMacro("its_clusters.C+");
  TEveUtil::LoadMacro("tpc_clusters.C+");

  TEveTrackCounter* g_trkcnt = new TEveTrackCounter("Primary Counter");
  gEve->AddGlobalElement(g_trkcnt);

  gAliEveEvent->AddNewEventCommand("on_new_event();");
  gAliEveEvent->GotoEvent(0);

  gEve->Redraw3D(kTRUE);
}

void on_new_event()
{
  TEvePointSet* itsc = its_clusters();
  itsc->SetMarkerColor(5);

  TEvePointSet* tpcc = tpc_clusters();
  tpcc->SetMarkerColor(4);

  primary_vertex(1, 1);

  TEveElementList* cont = esd_tracks_vertex_cut();
  TGListTree* lt = gEve->GetListTree();
  TGListTreeItem* ti = cont->FindListTreeItem(lt);
  lt->OpenItem(ti);

  // Here we expect five TEveTrackList containers.
  // First two have reasonable primaries (sigma-to-prim-vertex < 5).
  // Other three are almost certainly secondaries.
  Int_t count = 1;
  TEveTrackCounter* g_trkcnt = TEveTrackCounter::fgInstance;
  g_trkcnt->Reset();
  g_trkcnt->SetEventId(gAliEveEvent->GetEventId());
  TEveElement::List_i i = cont->BeginChildren();
  while (i != cont->EndChildren()) {
    TEveTrackList* l = dynamic_cast<TEveTrackList*>(*i);
    if (l != 0) {
      l->SetLineWidth(2);
      g_trkcnt->RegisterTracks(l, (count <= 2));
      ++count;
    }
    ++i;
  }
}

TParticle* id(Int_t label=0, Bool_t showParents=kTRUE)
{
  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  printf("Number primaries %d, all particles %d, label %d\n",
	 stack->GetNprimary(), stack->GetNtrack(), label);
  if (label < 0 || label >= stack->GetNtrack()) {
    printf("  Label exceeds available range.\n");
    return 0;
  }

  TParticle* part = stack->Particle(label);
  if(part != 0) {
    part->Print();
    if(showParents) {
      while(part->GetMother(0) >= 0) {
	part = stack->Particle(part->GetMother(0));
	part->Print();
      }
    }
  }

  return stack->Particle(label);
}
