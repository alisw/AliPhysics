{
  //load libraries
  //gSystem->SetBuildDir("/tmp");
  gROOT->LoadMacro("AliCentralityGlauberFit.cxx+");

  const char *finname ="/home/alberica/analysis/zdc/out5/output_f8TOT.root"; // name input file
  const char *foutname="test_AliCentralityGlauberFit.root"; // name output file

  float percentXsec=95.0;

  AliCentralityGlauberFit *mPM = new AliCentralityGlauberFit();
  mPM.AddHisto("hmultV0");
  mPM.SetGlauberParam(5,60.0,80.0,2,0.5,1.5,7,0.9,1.2); 
  // mPM.AddHisto("hNclusters0");
  // mPM.AddHisto("hNtracks");
  // mPM.AddHisto("hmultV0");
  // mPM.AddHisto("hmultFMD");

  //  mPM.SetPercentileCrossSection(percentXsec);
  mPM.SetOutputFile(foutname);

  mPM.MakeFits(finname);


}
