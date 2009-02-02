void AliTOFhits2sdigits(Int_t evNumber1=-1, Int_t numberOfEvents=0)
{

  /////////////////////////////////////////////////////////////////////////
  //
  // Creates TOF summable digits from the hits for all event in the header file
  //
  // Use case:
  // start root
  // // load the macro
  // root[0] .L AliTOFhits2sdigits.C
  // root[1] AliTOFhits2sdigits()
  //
  // By default, it creates sdigits for all the events in the header file.
  //
  // If you want create sdigits only the 3th event (existing in the header file)
  // you can use the following line:
  //
  // root[0] .L AliTOFhits2sdigits.C
  // root[1] AliTOFhits2sdigits(3,1)
  //
  // Created by: F. Pierella
  // Updated to the new I/O: C. Zampolli
  //
  // Report problems to decaro@sa.infn.it
  //
  /////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs

  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  } 
  
  if (gAlice)
    {
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice = 0x0;
    }
  
  // Create the TOF sdigitzer and sdigitize all events by default
  AliTOFSDigitizer *sdigitizer = new AliTOFSDigitizer("galice.root",evNumber1,numberOfEvents);

  // Activate this line if you want to print the parameters
  // used in sdigitization
  // sdigitizer->PrintParameters();

  // e.g. Activate this line if you want to sdigitize only hits 
  // with the plate number 3 and the sector number 15
  // pay attention that sector must be in the range [0,17]
  //                and plate  must be in the range [0,4]
  // by default we sdigitize hits of all plates in all sectors
  // sdigitizer->SelectSectorAndPlate(15,3);

  // performs sdigitization with "all" verbose option
  // "tim" option is also available for benchmarking only

  sdigitizer->Exec("all");

  sdigitizer = 0x0;
  delete sdigitizer;

  if (gAlice)
    {
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice = 0x0;
    }

}
