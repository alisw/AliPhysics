// Stand-Alone reve initialization.
// If class Reve::RGTopFrame is not know it attempts to load
// "libReve".

void reve_sa_init(int mode=1)
{
  TClass* c = gROOT->GetClass("Reve::RGTopFrame");
  if (!c) {
    Info("reve_sa_init", "tring to load libReve.");
    //gSystem->Load("libTreePlayer");
    //gSystem->Load("libGeomPainter");
    //gSystem->Load("libGed");
    gSystem->Load("libRGL");
    gSystem->Load("libEG");
    if (gSystem->Load("libReve") == -1)
      Warning("reve_sa_init", "loading of libReve failed.");
  }

  Reve::SetupEnvironment();
  Reve::RGTopFrame::SpawnGui(mode);
}
