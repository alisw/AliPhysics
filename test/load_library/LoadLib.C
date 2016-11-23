Int_t LoadLib(const TString &libName)
{
  if (libName.Contains("AliPythia") ||
      libName.Contains("TEPEMGEN") ||
      libName.Contains("TPHIC"))
    gSystem->Load("libpythia6");

  if (libName.Contains("MONITOR") ||
      libName.Contains("Eve") ||
      libName.Contains("HLT") ||
      libName.Contains("ConverterCollada")) {
    // Might load system-installed libraries (e.g. libdim.so)
    gSystem->AddDynamicPath("/lib:/lib64:/usr/lib:/usr/lib64:/usr/local/lib:/usr/local/lib64:/usr/lib/x86_64-linux-gnu");
  }

  return gSystem->Load(libName);
}
