// This macro displays the hits belonging to a track for selected detectors
// Input: in the tracks contains the interesting tracks
//        ntracks is the number of interesing tracks
//        The default values correspond to "Show everything"
// Note: For the moment it works only with HIJING events, the PYTHIA is 
//       still not supported 
void display (const char *filename="galice.root",Int_t nevent=0, Int_t * tracks=0, Int_t ntracks=0) {
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice = 0;
   }
      
// Connect the Root Galice file containing Geometry, Kine and Hits
   AliRunLoader *rl = 0x0;
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if(file){
     cout<<"galice.root is already open \n";
   }

   
   rl = AliRunLoader::Open(filename,"DISPLAYED EVENT");
   
   if (rl == 0x0)
    {
      cerr<<"Error <display.C()>: can not get Run Loader. Exiting"<<endl;
      return;
    }
// Get AliRun object from file or create it if not on file

   rl->LoadgAlice();
 
   gAlice = rl->GetAliRun();
   if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    return;
  }
    
// Load data
   rl->GetEvent(nevent);
   rl->LoadKinematics();
   rl->LoadHeader();
   rl->LoadHits();

// Create Event Display object
   AliDisplay *edisplay = new AliDisplay(750);
   if (ntracks>0) edisplay->SetTracksToDisplay(tracks, ntracks);
   
// Display the requested event
   edisplay->ShowNextEvent(0);

// Define the buttons to switch on/off the existing modules
   Float_t nodet=0.;
   TObjArray *moduli = gAlice->Modules();
   Int_t nomod=moduli->GetEntriesFast();
   AliModule *modu;
   for (Int_t j=0; j<nomod; j++){
     modu=(AliModule*)moduli->At(j);
     char *avoid=strstr("BODY MAG ABSO DIPO HALL FRAME SHIL PIPE",modu->GetName());
     if(avoid)continue;
     nodet++;
   }
   TDialogCanvas *dialog = new TDialogCanvas("Modules"," ",150,30*nodet);
   Float_t yval1=1./nodet*0.9*0.05;
   Float_t yval2=1./nodet*0.9*0.95;
   char action[50];
   char title[30];
   char bname[30];
   TButton *butto1;
   for (Int_t j=0; j<nomod; j++){
     modu=(AliModule*)moduli->At(j);
     char *avoid=strstr(" BODY MAG ABSO DIPO HALL FRAME SHIL PIPE",modu->GetName());
     if(avoid)continue;
     sprintf(action,"swioff(\"%s\")",modu->GetName());
     sprintf(title,"%s is on",modu->GetName());
     sprintf(bname,"but%s",modu->GetName());
     butto1 = new TButton(title,action,.05,yval1,.95,yval2);
     butto1->SetName(bname);
     butto1->SetFillColor(3);
     butto1->Draw();
     yval1+=1./nodet;
     yval2+=1./nodet;
   }
} 

void swioff(const char *dete){
  gAlice->Display()->DisableDetector(dete);
  gAlice->Display()->Pad()->Modified();
  gAlice->Display()->Pad()->Update();
  char bname[30];
  char action[50];
  char title[30];
  sprintf(bname,"but%s",dete);
  TDialogCanvas *dia = (TDialogCanvas *)gROOT->FindObject("Modules");
  TButton *bt = (TButton *)dia->FindObject(bname);
  bt->SetFillColor(2);
  sprintf(action,"swion(\"%s\")",dete);
  bt->SetMethod(action);
  sprintf(title,"%s is off",dete);
  bt->SetTitle(title);
  dia->Draw();
}

void swion(const char *dete){
  gAlice->Display()->EnableDetector(dete);
  gAlice->Display()->Pad()->Modified();
  gAlice->Display()->Pad()->Update();
  TDialogCanvas *dia = (TDialogCanvas *)gROOT->FindObject("Modules");
  char bname[30];
  char action[50];
  char title[30];
  sprintf(bname,"but%s",dete);
  TButton *bt = (TButton *)dia->FindObject(bname);
  bt->SetFillColor(3);
  sprintf(action,"swioff(\"%s\")",dete);
  bt->SetMethod(action);
  sprintf(title,"%s is on",dete);
  bt->SetTitle(title);
  dia->Draw();
}

