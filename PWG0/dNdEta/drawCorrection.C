/* $Id$ */

void drawCorrection()
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* dNdEtaMap = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
  dNdEtaMap->LoadCorrection("correction_map.root");

  dNdEtaMap->DrawHistograms();

  dNdEtaMap->GetMeasuredFraction(0.3, -1, kTRUE);
}
