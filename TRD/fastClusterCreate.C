void fastClusterCreate() {

///////////////////////////////////////////////////////////////////////// 
//
// Creates cluster from the hit information (fast simulator). 
// An additional hit-tree is added to the input file.
//
///////////////////////////////////////////////////////////////////////// 

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  // Input (and output) file name
  Char_t *alifile = "galice_r_v0.root";

  // Create the clusterizer
  AliTRDclusterizerV0 *Clusterizer = 
    new AliTRDclusterizerV0("clusterizer","fast clusterizer class"); 

  // Open the AliRoot file 
  Clusterizer->Open(alifile);

  // Find the cluster
  Clusterizer->MakeCluster();

  // Write the cluster into the input file
  Clusterizer->WriteCluster();

  // Save the clusterizer class in the AliROOT file
  Clusterizer->Write();

}
