void AliTRDDigitsDDL2Clusters() 
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

  // Create the clusterizer
  AliTRDclusterizerV1 *clusterizer = 
    new AliTRDclusterizerV1("clusterizer","Clusterizer class"); 

  // Read the parameter
  TFile *parfile = TFile::Open(infile);
  AliTRDparameter *par = (AliTRDparameter *) parfile->Get("TRDparameter"); 
  par->ReInit();
  clusterizer->SetParameter(par);

  // Set the parameter
  clusterizer->SetVerbose(1);

  //Number of events
  TTree * te = (TTree*)parfile->Get("TE");
  Int_t nev = (Int_t)te->GetEntries();

  for(Int_t iev=0;iev<nev;iev++) {

    AliRawReaderFile rawReader(iev);

    // Open the AliRoot file 
    clusterizer->Open(infile,iev);

    // Load the digits
    clusterizer->ReadDigits(&rawReader);
    clusterizer->Dump();
 
    // Find the cluster
    clusterizer->MakeClusters();

    // Write the cluster tree into file AliTRDclusters.root
    clusterizer->WriteClusters(-1);

  }

}
