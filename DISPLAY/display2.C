// This macro displays the hits belonging to a track for selected detectors
// Input: in the tracks contains the interesting tracks
//        ntracks is the number of interesing tracks
//        The default values correspond to "Show everything"
// Note: For the moment it works only with HIJING events, the PYTHIA is 
//       still not supported 
//#include <ClassTable.h>

void display2 (const char *filename="galice.root",Int_t nevent=0, Int_t * tracks=0, Int_t ntracks=0) {
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
	  delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice = 0;
   }
    //gROOT->LoadMacro("AliDisplay2.C");  
   gSystem->Load("libAliL3Src");
	gSystem->Load("libDISPLAY");
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
   
// Create Event Display object
   AliDisplay2 *edisplay = new AliDisplay2(gClient->GetRoot(), 900, 700);
//   if (ntracks>0) edisplay->SetTracksToDisplay(tracks, ntracks);

// Display the requested event
//   gAlice->GetEvent(nevent);
  rl->GetEvent(nevent);
   rl->LoadKinematics();
   rl->LoadHeader();
   rl->LoadHits();
   //	edisplay->FindModules();
   edisplay->ShowNextEvent(0);  
} 



