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

  // Input (and output) file name
  Char_t *alifile = "galice_c_v1.root";

  // Create the clusterizer
  AliTRDclusterizerV1 *Clusterizer = 
    new AliTRDclusterizerV1("clusterizer","slow clusterizer class"); 

  // Open the AliRoot file 
  Clusterizer->Open(alifile);

  // Load the digits
  Clusterizer->ReadDigits();

  // Find the cluster
  Clusterizer->MakeCluster();

  // Write the cluster into the input file
  Clusterizer->WriteCluster();

  // Save the clusterizer class in the AliROOT file
  Clusterizer->Write();

}
