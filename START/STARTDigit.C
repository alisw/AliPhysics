void STARTDigit () 
{
// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  manager->SetInputStream(0,"galice.root");
  //  manager->SetOutputFile("digits.root");
  AliSTARTDigitizer *START = new AliSTARTDigitizer(manager);
  manager->Exec("");
}














