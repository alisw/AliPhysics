// this macro combines the correction and the analysis and draws them

void CreatedNdEta(Bool_t correct = kTRUE, const Char_t* filename = "analysis_esd.root")
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* dNdEtaCorrection = 0;
  if (correct)
  {
    dNdEtaCorrection = new AlidNdEtaCorrection();
    dNdEtaCorrection->LoadHistograms("correction_map.root","dndeta_correction");
    //dNdEtaCorrection->RemoveEdges(2, 0, 2);
  }

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");

  TFile* file = TFile::Open(filename);
  if (!file)
  {
    cout << "Error. File out.root not found" << endl;
    return;
  }
  fdNdEtaAnalysis->LoadHistograms();

  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3);

  fdNdEtaAnalysis->DrawHistograms();
}


