void AliTRDdigits2cluster() 
{

///////////////////////////////////////////////////////////////////////// 
//
// Creates cluster from the digit information. 
//
///////////////////////////////////////////////////////////////////////// 

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  }

  // Input and output file names
  Char_t *infile  = "galice.root";
  Char_t *outfile = "AliTRDclusters.root";

  // Create the clusterizer
  AliTRDclusterizerV1 *clusterizer = 
    new AliTRDclusterizerV1("clusterizer","Clusterizer class"); 

  // Define the parameter object
  // If no external parameter object is defined, 
  // default parameter will be used
  AliTRDparameter *parameter = new AliTRDparameter("TRDparameter"
						  ,"TRD parameter class");
  clusterizer->SetParameter(parameter);

  // Open the AliRoot file 
  clusterizer->Open(infile,outfile,0);

  // Load the digits
  clusterizer->ReadDigits();
 
  // Find the cluster
  clusterizer->MakeClusters();

  // Write the cluster tree into file AliTRDclusters.root
  clusterizer->WriteClusters(-1);

  // Save the clusterizer class in the AliROOT file
  // clusterizer->Write();

}
