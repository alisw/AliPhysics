void MakeReferencePass1() {

  char tmpin[256];
  char tmpout[256];
  Int_t runno[]={
     259867, 259703, 259091, 260187, 259954, 259471, 259470, 259379, 259095

  };
  Int_t nruns=9;//172;

  for (Int_t i=0;i<nruns;i++){
    sprintf(tmpin,"%d/AnalysisResults.root",runno[i]);
    sprintf(tmpout,"Reference_%d.pass1.root",runno[i]);
    AliAnalysisTaskEMCALTimeCalib::ProduceCalibConsts(tmpin,tmpout,kFALSE);
  }



}

