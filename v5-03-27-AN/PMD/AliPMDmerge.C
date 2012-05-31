void AliPMDmerge()
{
  // This macro is for event merging
  AliRunDigitizer * manager = new AliRunDigitizer(2,1);
  manager->SetInputStream(0,"run1/galice.root");
  manager->SetInputStream(1,"run2/bg.root");
  manager->SetOutputFile("digits.root");
  AliPMDDigitizer *dpmd  = new AliPMDDigitizer(manager);
  manager->SetNrOfEventsToWrite(1);
  manager->Exec("");
}
