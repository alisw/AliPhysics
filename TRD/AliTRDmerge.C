void AliTRDmerge()
{
  /////////////////////////////////////////////////////////////////////////
  //
  // Merges different file with summable digits
  //
  /////////////////////////////////////////////////////////////////////////

  Char_t *fileSignal     = "galice_signal.root";
  Char_t *fileBackground = "galice_background.root";

  if (gAlice) {
    printf("<AliTRDmerge> Delete the AliRun object from memory.\n");
    delete gAlice;
    gAlice = 0;
  }

  AliRunDigitizer *manager = new AliRunDigitizer(2,1);
  manager->SetInputStream(0,fileSignal);
  manager->SetInputStream(1,fileBackground);

  AliTRDdigitizer *digitizer = new AliTRDdigitizer(manager
                                                  ,"TRDdigitizer"
                                                  ,"TRD digitizer class");

  // Define the parameter object
  // If no external parameter object is defined, 
  // default parameter will be used
  AliTRDparameter *parameter = new AliTRDparameter("TRDparameter"
						  ,"TRD parameter class");
  digitizer->SetParameter(parameter);

  // Do the merging
  manager->Exec("deb");

}
