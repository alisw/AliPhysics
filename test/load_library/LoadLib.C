Int_t LoadLib(const TString &libName)
{
  if (libName.Contains("AliPythia") ||
      libName.Contains("TEPEMGEN") ||
      libName.Contains("TPHIC"))
    gSystem->Load("libpythia6");

  return gSystem->Load(libName);
}
