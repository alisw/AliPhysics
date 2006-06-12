// this macro combines the correction and the analysis and draws them

void CreatedNdEta(Bool_t correct = kTRUE)
{
  gSystem->Load("libPWG0base");

  dNdEtaCorrection* dNdEtaCorrection = 0;
  if (correct)
  {
    dNdEtaCorrection = new dNdEtaCorrection();
    dNdEtaCorrection->LoadHistograms("correction_map.root","dndeta_correction");
    dNdEtaCorrection->RemoveEdges(2, 0, 2);
  }

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");

  TFile* file = TFile::Open("out.root");
  if (!file)
  {
    cout << "Error. File out.root not found" << endl;
    return;
  }
  fdNdEtaAnalysis->LoadHistograms();

  fdNdEtaAnalysis->Finish(dNdEtaCorrection);

  fdNdEtaAnalysis->DrawHistograms();
}


