void date2root(const char* prefix="2006run2211")
{
  // Converts raw DATE file to ROOT file.

  TString dfile(prefix);
  dfile += ".raw";

  TString rfile(prefix);
  rfile += ".root";

  AliSimulation sim;
  sim.ConvertDateToRoot(dfile.Data(),rfile.Data());

}
