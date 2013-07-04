void CreateLookup(const char* spaceChargeFile, const char* outputFile, Int_t phiBin) {
  TFile fn(spaceChargeFile);
  gROOT->cd();
  fSpaceCharge=(AliTPCSpaceCharge3D*)fn.Get("map");
  fn.Close();

  AliTPCCorrectionLookupTable l;
  l.SetupDefaultLimits();
  l.CreateLookupTableSinglePhi(*fSpaceCharge, phiBin, 1);

//   TString outFile(gSystem->BaseName(spaceChargeFile));
//   outFile.ReplaceAll(".root","");
//   outFile.Append(Form("_%03d",phiBin));
//   outFile.Append(".root");
//   outFile.Prepend("/");
//   outFile.Prepend(outputDir);

  TFile fout(outputFile,"recreate");
  l.Write("map");
  fout.Close();
}
