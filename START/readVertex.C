void readVertex(Int_t evNumber=1) 
{
 
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  
  // Connect the Root Galice file containing Geometry, Kine and Hits

  TFile *file =  (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
  if (!file) file = new TFile("galice.root","UPDATE");
  
  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }
  char nameTD[8],nameTR[8];

  TH1F *hVertex = new TH1F("hVertex","Z position of vertex",100,-350,350);
  
  digits = new AliSTARTdigit();
  vertex = new AliSTARTvertex();
  TBranch *bd;
  TBranch *bRec;

 // Event ------------------------- LOOP  
  for (j=0; j<evNumber; j++){
    sprintf(nameTD,"TreeD%d",j);
    printf("%s\n",nameTD);
    TTree *TD = (TTree*)gDirectory->Get(nameTD);
    bd = TD->GetBranch("START");
    bd->SetAddress(&digits);
    bd->GetEvent(0); 
    printf(" Digits: %d \n ",digits->GetTime()); 
    
    
    sprintf(nameTR,"TreeR%d",j);
    printf("%s\n",nameTR);
    TTree *TR = (TTree*)gDirectory->Get(nameTR);
    bRec = TR->GetBranch("START");
    bRec->SetAddress(&vertex);
    bRec->GetEvent(0);
    if(digits->GetTime()!=999999){
      hVertex->Fill(vertex->GetVertex());
      // printf(" Z position %f\n",vertex->GetVertex());
    }
  }
    Hfile = new TFile("figs1.root","RECREATE","Histograms for START Vertex");
   printf("Writting histograms to root file \n");
   Hfile->cd();

 //Create a canvas, set the view range, show histograms
   
  gStyle->SetOptStat(111111);
  TCanvas *c1 = new TCanvas("c1","Alice START Time (vertex)",400,10,600,600);

  hVertex->SetXTitle("vertex position, mm");
  hVertex->SetYTitle("number of events");
  hVertex->Write();

} // end of macro




