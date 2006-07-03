// this macro combines the correction and the analysis and draws them

void CreatedNdEta(Bool_t correct = kTRUE, const Char_t* filename = "analysis_esd.root", const char* object = "dndeta")
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* dNdEtaCorrection = 0;
  if (correct)
  {
    dNdEtaCorrection = new AlidNdEtaCorrection("dndeta_correction", "dndeta_correction");
    dNdEtaCorrection->LoadHistograms("correction_map.root","dndeta_correction");
    //dNdEtaCorrection->RemoveEdges(2, 0, 2);
  }

  fdNdEtaAnalysis = new dNdEtaAnalysis(object, object);

  TFile* file = TFile::Open(filename);
  if (!file)
  {
    cout << "Error. File out.root not found" << endl;
    return;
  }
  fdNdEtaAnalysis->LoadHistograms();

  fdNdEtaAnalysis->Finish(dNdEtaCorrection, (correct) ? 0.3 : -1);

  fdNdEtaAnalysis->DrawHistograms();
}


