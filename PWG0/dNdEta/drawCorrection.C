void drawCorrection()
{
  gSystem->Load("libdNdEta.so");
  
  dNdEtaCorrection* dNdEtaMap = new dNdEtaCorrection();
  dNdEtaMap->LoadCorrection("correction_map.root");
  
  dNdEtaMap->DrawHistograms();
}
