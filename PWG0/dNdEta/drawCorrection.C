/* $Id$ */

void drawCorrection()
{
  gSystem->Load("libPWG0base");

  dNdEtaCorrection* dNdEtaMap = new dNdEtaCorrection();
  dNdEtaMap->LoadCorrection("correction_map.root");
  
  dNdEtaMap->DrawHistograms();
}
