void merge()
{
  /////////////////////////////////////////////////////////////////////////
  //
  // Merges different file with summable digits
  //
  /////////////////////////////////////////////////////////////////////////

  AliRunDigitizer *manager = new AliRunDigitizer(2,1);
  manager->SetInputStream(0,"galice_1.root");
  manager->SetInputStream(1,"galice_2.root");

  AliTRDdigitizer *digitizer = new AliTRDdigitizer(manager
                                                  ,"TRDdigitizer"
                                                  ,"TRD digitizer class");

  manager->Exec("deb");

}
