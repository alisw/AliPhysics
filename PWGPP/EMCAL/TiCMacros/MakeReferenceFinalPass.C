void MakeReferenceFinalPass(TString output="Reference_LHC15i_final.root") {

  char tmpin[256];
  char tmpout[256];

    sprintf(tmpin,"AnalysisResults.root");
    //sprintf(tmpout,"Reference_LHC15i_calib_final.root");
    //AliAnalysisTaskEMCALTimeCalib::ProduceCalibConsts(tmpin,tmpout,kTRUE);
    AliAnalysisTaskEMCALTimeCalib::ProduceCalibConsts(tmpin,output.Data(),kTRUE);

}

