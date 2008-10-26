void embed_init()
{
  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init("Signal", -1);

  AliEveEventManager::GetMaster()->AddNewEventCommand("main_event()");

  // ------------------------------------------------------------------------

  Info("embed_init", "Opening background event ...");
  // Need strings as CINT can't properly convert "raw-string" to const TString&.
  TString name("Background Event");
  TString path("Background");
  AliEveEventManager* bkg = AliEveEventManager::AddDependentManager(name, path);
  bkg->IncDenyDestroy();
  bkg->AddNewEventCommand("background_event()");
  gEve->AddToListTree(bkg, kTRUE);

  TEveScene* bs = gEve->SpawnNewScene("Background");
  bs->AddElement(bkg);

  gEve->GetDefViewer()->AddScene(bs);

  // ------------------------------------------------------------------------

  TEveUtil::LoadMacro("SplitGLView.C+"); // Needed for dependencies.
  TEveUtil::LoadMacro("its_clusters.C+");
  TEveUtil::LoadMacro("tpc_clusters.C+");

  // ------------------------------------------------------------------------

  TEveBrowser* browser = gEve->GetBrowser();

  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
  browser->StopEmbedding("EventCtrl");

  // ------------------------------------------------------------------------

  AliEveEventManager::GetMaster()->GotoEvent(0);
  gEve->Redraw3D(kTRUE);
}

void main_event()
{
  printf("Main Event - post load\n");

  its_clusters();
  tpc_clusters();
}

void background_event()
{
  printf("Background Event - post load\n");

  TEvePointSet* c;

  c = its_clusters();
  c->SetMarkerColor(kOrange);

  c = tpc_clusters();
  c->SetMarkerColor(kCyan);
}
