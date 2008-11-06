void
DrawGeometry(const char* vol="ALIC")
{
  AliLog::SetModuleDebugLevel("FMD", 1);
  gAlice->InitMC("$ALICE_ROOT/FMD/Config.C");
  gGeoManager->GetVolume(vol)->Draw("ogl");
  new TBrowser("geoBrowser", gGeoManager);
}

