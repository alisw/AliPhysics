Int_t giNevents;//should be global for menu->AddButton(...)

void rMenu(Int_t iNevents=1)// How many events to generate, digitise etc.
{ 
   giNevents=iNevents;
   cout<<giNevents<<" event(s) are requested\n";
   
   gROOT->ProcessLine(".x loadlibs.C"); // load ALICE libraries if that are not yet loaded. see rich/loadlibs.C
     
   TControlBar *menu = new TControlBar("vertical","RICH");
   
   menu->AddButton("Help for RICH",      "gSystem->Exec(\"less RICHHelp.txt\");", "Explains how to use RICH menus");
   menu->AddButton("Run",                "gAlice->Run(giNevents)","Process an Alice event - WARNING: Overwrites previous data file!");
   menu->AddButton("Run Lego",           ".x RICHRunLego.C","Special runs to generate the radl/absl lego plots");
   menu->AddButton("Merge and Digitise", "MergeAndDigitise()","Merge with background file and digitise");
   menu->AddButton("Clusterize Event",   ".x RICHrawclusters.C(0,events-1)","Reconstruct clusters");
   menu->AddButton("3D Hough Pat. Rec.", ".x RICHdetect.C(0,events-1)","Lisbon");
   menu->AddButton("1D Hough Pat. Rec.", ".x RICHpatrec.C(0,events-1)","Bari");
   menu->AddButton("Diagnostics",        "Diagnostics()",             "Miscellaneous diagnostics");
   menu->AddButton("Display",            ".x RICHdisplay.C","Display run");
   menu->AddButton("Geometry Browser",   "Gui()","Browse the GEANT geometry - WARNING: Overwrites previous data file!");
   menu->AddButton("File Browser",       "TBrowser new;","Browse data files");
   menu->AddButton("Quit AliRoot",       ".q","Close session");

   menu->Show();
}//void rmenu(

void RICHInit(Int_t events)
{
  gAlice->Init("Config.C");
  ((TGeant3*)gMC)->InitHIGZ();
  gAlice->Run(events);
}

void Gui()
{
  gAlice->Init("Config.C");
  ((TGeant3*)gMC)->InitHIGZ();
  gROOT->ProcessLine(".x TGeant3GUI.C");
}

void MergeAndDigitise(Int_t iFirstEventN=0,Int_t iLastEventN=giNevents-1, Int_t merging=0)
{
   TFile *pFile = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if(!pFile){
      cout<<"MergeAndDigitise> Opening file galice.root\n"; 
      pFile = new TFile("galice.root","UPDATE");
   }
	 
   if(!gAlice) gAlice=(AliRun*)pFile->Get("gAlice");
   cout<<"MergerAndDigitise> "<<gAlice->GetEvNumber()<<" events found\n";

   AliRICHMerger* pRICHMerger = new AliRICHMerger();
      
   pRICHMerger->SetMode(merging);
   pRICHMerger->SetSignalEventNumber(0);
   pRICHMerger->SetBackgroundEventNumber(0);
   pRICHMerger->SetBackgroundFileName("bg.root");
   
   AliRICH *pRICH= (AliRICH*) gAlice->GetDetector("RICH");
   pRICH->Print();   
        
}//void MergeAndDigitise(Int_t evNumber1=0,Int_t evNumber2=0, Int_t merging=0)


void Diagnostics()
{
   TControlBar *menu = new TControlBar("vertical","RICH diag");
   
   menu->AddButton("Single Ring Hits",".x RICHpadtest.C(1,0,events-1","Hits in central region");
   menu->AddButton("Single Ring Spectra",".x RICHpadtest.C(2,0,events-1","Photon spectra");
   menu->AddButton("Single Ring Statistics",".x RICHpadtest.C(3,0,events-1","Production and clusters");
   menu->AddButton("Single Ring Reconstruction",".x RICHpadtest.C(4,0,events-1","Generated and reconstructed values");
   menu->AddButton("Single Ring Gain Variation",".x RICHgainvar.C(0,events-1","Gain variation along wires");
   menu->AddButton("Full Event Hits",".x RICHpadtest.C(5,0,events-1","Hits in seven modules");
   menu->AddButton("Full Event Spectra",".x RICHpadtestC.C","Individual particles' spectra and fluxes");
   menu->AddButton("Full Event Occupancy",".x RICHoccupancy.C","Mean and per chamber occupancy");
   
   menu->Show();  
}//void Diagnostics()
