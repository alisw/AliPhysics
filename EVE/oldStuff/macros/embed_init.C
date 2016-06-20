TEveViewer *gSignalView     = 0;
TEveViewer *gBackgroundView = 0;

void embed_init()
{
  TEveUtil::LoadMacro("alieve_init.C");
  alieve_init("Signal", -1);

  AliEveEventManager::Instance()->AddNewEventCommand("main_event()");

  // ------------------------------------------------------------------------

  Info("embed_init", "Opening background event ...");
  AliEveEventManager* bkg =
    AliEveEventManager::AddDependentManager("Background Event", "Background");
  bkg->IncDenyDestroy();
  bkg->AddNewEventCommand("background_event()");
  gEve->AddToListTree(bkg, kTRUE);

  TEveScene* bs = gEve->SpawnNewScene("Background");
  bs->AddElement(bkg);

  gEve->GetDefaultViewer()->AddScene(bs);

  // ------------------------------------------------------------------------

  TEveBrowser* browser = gEve->GetBrowser();

  TEveWindowSlot *slot = 0;
  TEveWindowPack *pack = 0;

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  pack = slot->MakePack();
  pack->SetElementName("Parallel View");
  pack->SetHorizontal();
  pack->SetShowTitleBar(kFALSE);

  pack->NewSlot()->MakeCurrent();
  gSignalView = gEve->SpawnNewViewer("Signal View", "");
  gSignalView->AddScene(gEve->GetEventScene());

  pack->NewSlot()->MakeCurrent();
  gBackgroundView = gEve->SpawnNewViewer("Background View", "");
  gBackgroundView->AddScene(bs);

  // ------------------------------------------------------------------------

  TEveUtil::LoadMacro("its_clusters.C+");
  TEveUtil::LoadMacro("tpc_clusters.C+");

  // ------------------------------------------------------------------------

  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::Instance());
  browser->StopEmbedding("EventCtrl");

  // ------------------------------------------------------------------------

  AliEveEventManager::Instance()->GotoEvent(0);
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
