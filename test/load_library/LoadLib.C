Int_t LoadLib(const TString &libName)
{
  if (libName.Contains("PWGPP") ||
      libName.Contains("PWGUDdiffractive")) {
    // Might need system-installed libraries (e.g. libzmq.so)
    gSystem->AddDynamicPath("/lib:/lib64:/usr/lib:/usr/lib64:/usr/local/lib:/usr/local/lib64:/usr/lib/x86_64-linux-gnu");
  }

  return gSystem->Load(libName);
}
