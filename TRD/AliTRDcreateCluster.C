Int_t AliTRDcreateCluster()
{
  //
  // Creates the cluster from the digits
  //

  Int_t rc = 0;

  // Create the clusterizer
  AliTRDclusterizerV1 *Clusterizer =
    new AliTRDclusterizerV1("clusterizer","Clusterizer class");

  // Set the parameter
  Clusterizer->SetClusMaxThresh(0.0);
  Clusterizer->SetClusSigThresh(0.0);
  Clusterizer->Dump();
 
  // Open the file
  if (!(Clusterizer->Open("TRD_test.root",0))) {
    rc = 1;
    return rc;
  }    

  // Load the digits
  if (!(Clusterizer->ReadDigits())) {
    rc = 2;
    return rc;
  }    

  // Find the cluster
  if (!(Clusterizer->MakeClusters())) {
    rc = 3;
    return rc;
  }

  // Write the cluster tree into the file 
  if (!(Clusterizer->WriteClusters(-1))) {
    rc = 4;
    return rc;
  }

  // Save the clusterizer class in the file
  if (!(Clusterizer->Write())) {
    rc = 5;
    return rc;
  }

  return rc;

}
