void display (const char *filename="galice.root",Int_t nevent=0) {
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
      delete gAlice;
      gAlice = 0;
   }
      
// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if(file){
     cout<<"galice.root is already open \n";
   }
   else {
     if(!file)file=TFile::Open(filename);
   }

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   
// Create Event Display object
   AliDisplay *edisplay = new AliDisplay(750);

// Display the requested event
   gAlice->GetEvent(nevent);
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

