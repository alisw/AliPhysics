// $Id$
//
// Root macro that opens a mini GUI for running aliroot with Geant4.
//      
//  To run aliroot with Geant4 using the g4menu.C:
//  aliroot
//  root [0] .x g4menu.C
//  --> Select "Init" and then "Run" button
//	    
// The the bar enables to start Geant4 interactive session:
//  --> Select "Geant4UI" button and use Geant4 interactive commands;
//      In case TGeant4 has not yet been created you need first
//      select "Geant4" button before selecting  "Geant4UI". 	    
// To go back to Root UI, type exit.


#include <iostream>

void g4menu()
{

  // Load Geant4 libraries 
  if (!gInterpreter->IsLoaded("$ALICE/geant4_vmc/examples/macro/g4libs.C"))
    gROOT->LoadMacro("$ALICE/geant4_vmc/examples/macro/g4libs.C");
  gInterpreter->ProcessLine("g4libs()");

  // Menu
  TControlBar* menu = new TControlBar("vertical","Alice Geant4 menu");
  
  menu->AddButton("Init",     "gAlice->Init(\"g4Config.C\")", "Initialize \" AliRun \"");
  menu->AddButton("Run",      "gAlice->Run()",  "Process Alice run");
  menu->AddButton("Geant4",   "CreateGeant4()", "Create Geant4 only (without initializing AliRun)");
  menu->AddButton("Geant4UI", "StartGeant4UI()","Go to Geant4 Interactive session");
  menu->AddButton("XML",      "GenerateXML()","Generate XML (AGDD) file with geometry description");
  menu->AddButton("Quit",     "Quit()", "Quit aliroot");
  gROOT->SaveContext();
  menu->Show();
}

void CreateGeant4()
{  
  if (!gMC) {
  
     // TG4RunConfiguration for Geant4
     TG4RunConfiguration* runConfiguration 
      = new TG4RunConfiguration(true);
  
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

void GenerateXML()
{
  if (gMC) {
    // release Root terminal control

    // go into non-raw term mode
    //Getlinem(kCleanUp, 0);
    
    // add test if gMC is TGeant4
    TGeant4* g4 = (TGeant4*)gMC;
    
    g4->ProcessGeantCommand("/xml/generateAGDD");

    // new Root prompt
    //Getlinem(kInit, ((TRint*)gROOT->GetApplication())->GetPrompt());  
  }
  else {  
    cout << "Monte Carlo has not been yet created." << endl;
  } 
}        

void Quit()
{
  delete gAlice->GetRunLoader();
  delete gAlice;
  
  exit(0);
}  
