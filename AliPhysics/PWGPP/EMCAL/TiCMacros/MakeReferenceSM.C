void MakeReferenceSM(TString outname="ReferenceSM_LHC15o_mcp1_v1.root") {

  char tmpin[256];
  char tmpout[256];
  Int_t runno[]={
   238097, 238164, 238170, 238176, 238179, 238184, 238432, 238451, 238454, 238455, 238456, 238457, 238458, 238459, 238460, 238468, 238469, 238470, 238474
  };
  Int_t nruns=19;

  for (Int_t i=0;i<nruns;i++){
    sprintf(tmpin,"Reference_%d.pass1.root",runno[i]);
    TFile *file=new TFile(tmpin);
    if(!file) {

      file->ls();
      delete file;

      continue;
    }
    //sprintf(tmpout,"ReferenceSM_LHC15o_muon_calo_pass1_v1.root");
    sprintf(tmpout,outname.Data());
    AliAnalysisTaskEMCALTimeCalib::ProduceOffsetForSMsV2(runno[i],tmpin,tmpout,kFALSE,kTRUE);
  }



}

