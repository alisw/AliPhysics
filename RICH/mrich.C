Int_t events=1;

mrich()
{ 

  //gAlice=0;
  //char *confi="Config.C";
  //Int_t events=10;
  
   TControlBar *menu = new TControlBar("vertical","RICH menu");
   menu->AddButton("      Help for RICH      ","RICHHelp()", "Explains how to use RICH menus");
   menu->AddButton("Configure",            "gSystem->Exec(\"rconfig\"); gSystem->Exit(0);","Configure the simulation");
   menu->AddButton("Run",               "RICHInit(events)","Process an Alice event - WARNING: Overwrites previous data file!");
   //menu->AddButton("Run Lego",          "gAlice->RunLego()","Special run to generate the radl/absl lego plots");
   menu->AddButton("Digitise Event",             ".x RICHdigit.C(0,events-1)","Digitise event");
   menu->AddButton("Clusterize Event",      ".x RICHrawclusters.C(0,events-1)","Reconstruct clusters");
   // TODO: add the diaglevel here before the script
   menu->AddButton("3D Hough Pat. Rec.",      ".x RICHdetect.C(0,events-1)","Lisbon");
   menu->AddButton("1D Hough Pat. Rec.",      ".x RICHpatrec.C(0,events-1)","Bari");
   menu->AddButton("Diagnostics",       ".x RICHpadtest.C(3,0,events-1)","Miscellaneous diagnostics");
   menu->AddButton("Display",           ".x RICHdisplay.C","Display run");
   menu->AddButton("Geometry Browser",           "Gui()","Browse the GEANT geometry - WARNING: Overwrites previous data file!");
   menu->AddButton("File Browser",           "TBrowser new;","Browse data files");
//   menu->AddButton("Test",              ".x RICHtest.C","bla bla");
//   menu->AddButton("Edit Configuration",".x RICHConfig.C","Interactive Configuration");
//   menu->AddButton("Draw",              ".x DrawRICH.C","bla bla");
//   menu->AddButton("View",              ".x ViewRICH.C","does nothing???");
   menu->AddButton("Quit",             ".q","Close session");
//   menu->AddButton("Reset",             "RICHReset()","Close and Restart AliRoot");
   
   //gROOT->SaveContext();

   //gAlice->Init(config);
   //((TGeant3*)gMC)->InitHIGZ();

   menu->Show();
}

void RICHHelp()
{
   gSystem->Exec("xemacs RICHHelp.C &");
}

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
  //printf("Doesn't work yet\n");
}


void RICHReset()
{
//   gSystem->Exec("aliroot mrich.C &");

   // This doesn't work for Win (eheheh) and aliroot must be in the path
   gSystem->Exec("xterm +ls -e aliroot mrich.C &");
   gSystem->Exit(0);
}


