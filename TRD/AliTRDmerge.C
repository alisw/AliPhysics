void AliTRDmerge()
{
  /////////////////////////////////////////////////////////////////////////
  //
  // Merges different file with summable digits
  //
  /////////////////////////////////////////////////////////////////////////

  AliRunDigitizer *manager = new AliRunDigitizer(2,1);
  manager->SetInputStream(0,"galice_signal.root");
  manager->SetInputStream(1,"galice_background.root");

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
