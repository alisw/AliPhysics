{
  //load libraries
  //gSystem->SetBuildDir("/tmp");
  gROOT->LoadMacro("AliCentralityBy1D.cxx+");

  const char *finname ="/home/alberica/analysis/zdc/out5/output_f8TOT.root"; // name input file
  const char *foutname="test_AliCentralityBy1D_95.root"; // name output file

  float percentXsec=95.0;

  AliCentralityBy1D *mPM = new AliCentralityBy1D();
  mPM.AddHisto("hNtracklets");
  mPM.AddHisto("hNclusters0");
  mPM.AddHisto("hNtracks");
  mPM.AddHisto("hmultV0");
  mPM.AddHisto("hmultFMD");

  mPM.SetPercentileCrossSection(percentXsec);
  mPM.SetPercentileFile(foutname);

  mPM.MakePercentiles(finname);


}
