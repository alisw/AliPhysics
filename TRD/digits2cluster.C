void digits2cluster() 
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
  Char_t *outfile = "TRDclusters.root";

  // Create the clusterizer
  AliTRDclusterizerV1 *clusterizer = 
    new AliTRDclusterizerV1("clusterizer","Clusterizer class"); 

  // Set the parameter
  clusterizer->SetClusMaxThresh(0);
  clusterizer->SetClusSigThresh(0);
  //clusterizer->SetVerbose(1);
  clusterizer->Dump();

  // Open the AliRoot file 
  clusterizer->Open(infile,0);
  //clusterizer->Open(infile,outfile,0);

  // Load the digits
  clusterizer->ReadDigits();
 
  // Find the cluster
  clusterizer->MakeClusters();

  // Write the cluster tree into file AliTRDclusters.root
  clusterizer->WriteClusters(-1);

  // Save the clusterizer class in the AliROOT file
  // clusterizer->Write();

}
