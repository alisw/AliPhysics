#include <iostream.h>

void MUONrecodisplay(Int_t evNumber=0)
{
////////////////////////////////////////////////////////////////////
// The display pops up a context menu when the right-button is    //
//	clicked on the main pad. The following functions are available //
//	* SetDrawHits()	- switches on or off Geant hits ;				//
//	* CutMomentum()	- displays only tracks within Pmin - Pmax 	//
//	* ListTracks()	- prints ID and momentum info. for all				//
//	tracks within momentum range Pmin,Pmax;								//
//	* Highlight()	- shows only one selected reco. track				//
//	and its best matching Geant track;										//
//	*UnHighlight()	- self explaining;										//
//	*SetDrawHits()	- switch on or off Geant hits							//
////////////////////////////////////////////////////////////////////
   if (evNumber<0) return;


//Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
	loadlibs();
   }
    
   
// Connect the Root Galice file containing ...
    
   TFile *galice_file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!galice_file) galice_file = new TFile("galice.root");
   if (!galice_file) {
      cout << "File galice.root not found\n";
      return;
   }
     
// Get AliRun object from file or create it if not on file
   if (!gAlice) {
	gAlice = (AliRun*)galice_file->Get("gAlice");
	if (gAlice) {
         cout << "AliRun object found on file\n";
      } else {
         cout << "AliRun object not found on file !!!\n";
         return;
      }
   }
    
   AliMUONRecoDisplay *display = new AliMUONRecoDisplay(evNumber);   

   display->DisableDetector("ITS");
   display->DisableDetector("TPC");
   display->DisableDetector("TOF");
   display->DisableDetector("RICH");
   display->DisableDetector("ZDC");
   display->DisableDetector("CASTOR");
   display->DisableDetector("TRD");
   display->DisableDetector("FMD");
   display->DisableDetector("PHOS");
   display->DisableDetector("PMD");

   display->ShowNextEvent(0);

}
