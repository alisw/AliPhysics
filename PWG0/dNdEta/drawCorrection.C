/* $Id$ */

void drawCorrection()
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* dNdEtaMap = new AlidNdEtaCorrection();
  dNdEtaMap->LoadCorrection("correction_map.root");

  dNdEtaMap->DrawHistograms();

  dNdEtaMap->GetMeasuredFraction(0.3, -1, kTRUE);
}
