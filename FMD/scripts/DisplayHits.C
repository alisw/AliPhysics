//____________________________________________________________________
//
// $Id$
//
// Draw hits in the specialised FMD event display 
//
void
DisplayHits()
{
  gSystem->Load("libFMDutil.so");
  AliFMDDisplay* d = new AliFMDDisplay;
  d->AddLoad(AliFMDInput::kHits);
  d->AddLoad(AliFMDInput::kKinematics);
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
