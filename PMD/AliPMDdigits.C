void AliPMDdigits()
{
  // This macro converts sdigits to digits
  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  manager->SetInputStream(0,"galice.root");
  AliPMDDigitizer *dpmd  = new AliPMDDigitizer(manager);
  manager->Exec("");
}
