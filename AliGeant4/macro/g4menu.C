#include <iostream.h>

void g4menu()
{
  // load Geant4 libraries
  if (!gInterpreter->IsLoaded("g4libs.C")) gROOT->LoadMacro("g4libs.C");
  gInterpreter->ProcessLine("g4libs()");

  // load AliRoot core libraries
  gInterpreter->ProcessLine("steerlibs()");

  // menu
  TControlBar* menu = new TControlBar("vertical","Alice Geant4 menu");
  
  menu->AddButton("Geant4",   "CreateGeant4()", "Create Geant4 Monte Carlo");
  menu->AddButton("Geant4UI", "StartGeant4UI()","Go to Geant4 Interactive session");
  menu->AddButton("Init",     "gAlice->Init()", "Initialize AliRun");
  menu->AddButton("Run",      "gAlice->Run()",  "Process Alice run");
  menu->AddButton("Run Lego", "gAlice->RunLego()", "Process special lego run");
  gROOT->SaveContext();
  menu->Show();
}

void CreateGeant4()
{
  if (!gMC) {
    // AliRunConfiguration for Geant4
    AliRunConfiguration* runConfiguration 
      = new AliRunConfiguration();
  
    // TGeant4
    TGeant4* geant4
      = new TGeant4("TGeant4", "The Geant4 Monte Carlo", runConfiguration);

    cout << "Geant4 has been created." << endl;
   }      
  else {
    cout << "Monte Carlo has been already created." << endl;
  }       
}    

void StartGeant4UI()
{
  if (gMC) {
    // release Root terminal control

    // go into non-raw term mode
    Getlinem(kCleanUp, 0);
    
    // add test if gMC is TGeant4
    TGeant4* g4 = (TGeant4*)gMC;
  
    g4->StartGeantUI();

    // new Root prompt
    Getlinem(kInit, ((TRint*)gROOT->GetApplication())->GetPrompt());  
  }
  else {  
    cout << "Monte Carlo has not been yet created." << endl;
  }       
}  

