void readDigits(Int_t evNumber=1) 
{
 
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  
  // Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file =  (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
  //TFile *file =  (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
  if (!file) file = new TFile("galice.root","UPDATE");
  
  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }
  char nameTD[8],nameTR[8];

  TH1F *hTimediff = new TH1F("hTimediff","Time difference",100,0,256);
  TH1F *hTimePs = new TH1F("hTimePs","Time in Ps",100,-2000,2000);
  
  digits = new AliSTARTdigit();
  TBranch *bd;


 // Event ------------------------- LOOP  
  for (j=0; j<evNumber; j++){
    sprintf(nameTD,"TreeD%d",j);
    printf("%s\n",nameTD);
    TTree *TD = (TTree*)gDirectory->Get(nameTD);
    bd = TD->GetBranch("START");
    bd->SetAddress(&digits);
    bd->GetEvent(0); 
    //    printf(" Digits: %d \n ",digits->GetTime()); 
    
    printf("time %d\n",digits->GetTime());
    
    if(digits->GetTime()!=999999){
      Int_t timediff = digits->GetTime();
      Double_t timePs=(timediff-128)*10.; // time in Ps channel_width =10ps
      //      cout<<"timediff "<<timediff<<" timePs "<<timePs<<endl;
      hTimediff->Fill(timediff);
      hTimePs->Fill(timePs);
    }
  }

    Hfile = new TFile("figs.root","UPDATE","Histograms for STASRT digits");
   printf("Writting histograms to root file \n");
   Hfile->cd();
//Create a canvas, set the view range, show histograms
  gStyle->SetOptStat(111111);
  TCanvas *c1 = new TCanvas("c1","Alice START Time ",400,10,600,600);
  hTimePs->SetXTitle("arriving time, ps");
  hTimePs->SetYTitle("number of events");
  hTimePs->Write();

    
 
} // end of macro




