void AliTRDdigits2clusterMI() 
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
  AliTRDclusterizerMI *clusterizer = 
    new AliTRDclusterizerMI("clusterizer","Clusterizer class"); 

  // Read the parameter
  TFile *parfile = TFile::Open(infile);
  AliTRDparameter *par = (AliTRDparameter *) parfile->Get("TRDparameter"); 
  par->ReInit();
  clusterizer->SetParameter(par);

  // Set the parameter
  clusterizer->SetVerbose(1);

  // Open the AliRoot file 
  clusterizer->Open(infile,0);
  //clusterizer->Open(infile,outfile,0);


  // Load the digits
  clusterizer->ReadDigits();
  clusterizer->Dump();
 
  // Find the cluster
  clusterizer->MakeClusters();

  // Write the cluster tree into file AliTRDclusters.root
  clusterizer->WriteClusters(-1);

  // Save the clusterizer class in the AliROOT file
  // clusterizer->Write();

}
