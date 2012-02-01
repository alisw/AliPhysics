{
  //load libraries
  //gSystem->SetBuildDir("/tmp");
  gROOT->LoadMacro("AliCentralityByFunction.cxx+");

  const char *finname ="/home/alberica/analysis/zdc/out5/output_f8TOT_correlations.root"; // name input file
  const char *foutname="test_AliCentralityByFunction.root"; // name output file

  float percentXsec=100.0;

  AliCentralityByFunction *mPM = new AliCentralityByFunction();
  
  mPM.AddHisto("hmultV0vsmultFMD_all"); 
  mPM.SetFitFunction("hmultV0vsmultFMD_all","fitf_pol3", 0,30000);
  
  mPM.AddHisto("hNtrackletsvsmultV0_all");  
  mPM.SetFitFunction("hNtrackletsvsmultV0_all","fitf_pol2",0,10000);
  
  mPM.AddHisto("hEzemvsEzdc_all");
  mPM.SetFitFunction("hEzemvsEzdc_all","fitf_pol6",0,3.1);

  // mPM.AddHisto("hNtracksvsEzdc_all");
  // mPM.AddHisto("hNtrackletsvsEzdc_all");
  // mPM.AddHisto("hNclusters0vsEzdc_all");
  // mPM.AddHisto("hmultV0vsEzdc_all");
  // mPM.AddHisto("hmultFMDvsEzdc_all");
  // mPM.AddHisto("hNtracksvsEzem_all");
  // mPM.AddHisto("hNtrackletsvsEzem_all");
  // mPM.AddHisto("hNclusters0vsEzem_all");
  // mPM.AddHisto("hmultV0vsEzem_all");
  // mPM.AddHisto("hmultFMDvsEzem_all");
  // mPM.AddHisto("hNtracksvsmultV0_all");    
  // mPM.AddHisto("hNclusters0vsmultV0_all");
  // mPM.AddHisto("hNtracksvsmultFMD_all");
  // mPM.AddHisto("hNtrackletsvsmultFMD_all");
  // mPM.AddHisto("hNclusters0vsmultFMD_all");
  
  mPM.SetPercentileCrossSection(percentXsec);
  mPM.SetPercentileFile(foutname);

  mPM.MakePercentiles(finname);


}
