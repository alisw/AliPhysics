void slowClusterCreate() {

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
  AliTRDclusterizerV1 *Clusterizer = 
    new AliTRDclusterizerV1("clusterizer","slow clusterizer class"); 

  // Define output file name
  Clusterizer->Init(outfile);

  // Open the AliRoot file 
  Clusterizer->Open(infile);

  // Load the digits
  Clusterizer->ReadDigits();
 
  // Find the cluster
  Clusterizer->MakeClusters();

  // Write the cluster tree into file AliTRDclusters.root
  Clusterizer->WriteClusters(-1);

  // Save the clusterizer class in the AliROOT file
  // Clusterizer->Write();

}
