void readVertex(Int_t evNumber=1) 
{
 
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  
  // Connect the Root Galice file containing Geometry, Kine and Hits

  TFile *file =  (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
  if (!file) file = new TFile("production.root","UPDATE");
  
  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }
  char nameTD[8],nameTR[8];

  TH1F *hVertex = new TH1F("hVertex","Z position of vertex",100,-100,100);
  TH1F *hRealVertex = new TH1F("hRealVertex","Z position of vertex",
			       100,-100,100);
  
  digits = new AliT0digit();
  vertex = new AliT0vertex();
 // Event ------------------------- LOOP  
  for (Int_t j=0; j<evNumber; j++){
    
    sprintf(nameTR,"T0_R_%d",j);
    printf("%s\n",nameTR);
    TObject *tr = (TObject*)gDirectory->Get(nameTR);
    cout<<" tr "<<tr<<endl;
    vertex->Read(nameTR);
    
    hVertex->Fill(vertex->GetVertex());
    printf(" Z position %f\n",vertex->GetVertex());
    gAlice->GetEvent(j);
    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* genHeader = header->GenEventHeader();
    TArrayF *o = new TArrayF(3); 
    genHeader->PrimaryVertex(*o);
    cout<<" o "<<o<<endl;
    Float_t zRealVertex=o->At(2);
    cout<<" zRealVertex "<<zRealVertex<<endl;
    hRealVertex->Fill(zRealVertex);
       
  }
    Hfile = new TFile("figs.root","RECREATE","Histograms for T0 Vertex");
   printf("Writting histograms to root file \n");
   Hfile->cd();

 //Create a canvas, set the view range, show histograms
   
  gStyle->SetOptStat(111111);
  hVertex->SetXTitle("vertex position, mm");
  hVertex->SetYTitle("number of events");
  hVertex->Write();
  hRealVertex->SetXTitle("vertex position, mm");
  hRealVertex->SetYTitle("number of events");
  hRealVertex->Write();

} // end of macro




