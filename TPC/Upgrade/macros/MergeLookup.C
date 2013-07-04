void MergeLookup(const char* files, const char* outfile) {
  AliTPCCorrectionLookupTable l;
  l.MergePhiTables(files);

  TFile f(outfile,"recreate");
  l.Write("map");
  f.Close();
}
