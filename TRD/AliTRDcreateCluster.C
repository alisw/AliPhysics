Int_t AliTRDcreateCluster()
{
  //
  // Creates the cluster from the digits
  //

  Int_t rc = 0;

  // Create the clusterizer
  AliTRDclusterizerV1 *clusterizer =
    new AliTRDclusterizerV1("clusterizer","Clusterizer class");

  // Set the parameter
  clusterizer->SetClusMaxThresh(0);
  clusterizer->SetClusSigThresh(0);
  clusterizer->Dump();
 
  // Open the file
  if (!(clusterizer->Open("TRD_test.root",0))) {
    rc = 1;
    return rc;
  }    

  // Load the digits
  if (!(clusterizer->ReadDigits())) {
    rc = 2;
    return rc;
  }    

  // Find the cluster
  if (!(clusterizer->MakeClusters())) {
    rc = 3;
    return rc;
  }

  // Write the cluster tree into the file 
  if (!(clusterizer->WriteClusters(-1))) {
    rc = 4;
    return rc;
  }

  // Save the clusterizer class in the file
  if (!(clusterizer->Write())) {
    rc = 5;
    return rc;
  }

  return rc;

}
