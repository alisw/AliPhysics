void
TestAlignable(const char* file="geom.root", Int_t lvl=35)
{
  TGeoManager::Import(file);
  if (!gGeoManager) return;
  
  AliFMDGeometry* g = AliFMDGeometry::Instance();
  g->Init();
  
  AliLog::SetModuleDebugLevel("FMD", lvl);
  

  g->SetAlignableVolumes();
}

