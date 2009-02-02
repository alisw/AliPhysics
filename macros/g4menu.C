// $Id$
//
// Root macro that opens a mini GUI for running aliroot with Geant4.
//      
//  To run aliroot with Geant4 using the g4menu.C:
//  aliroot
//  root [0] .x g4menu.C
//  --> First select "Geometry" to build geometry.root file
//  --> Then re-run aliroot and select "Run" button
//	    
// The Init button is kept for debugging purposes,
// it itializes MonteCarlo but it does not itialize
// completely ALICE framework. That's why to run simulation,
// you have to re-run aliroot and select Run button.
//
// The menu enables to start Geant4 interactive session:
//  --> Select "Geant4UI" button and use Geant4 interactive commands;
// To go back to Root UI, type exit.
//
// By I. Hrivnacova, IPN Orsay


#include <iostream>

void g4menu()
{

  // Load Geant4 libraries 
  if (!gInterpreter->IsLoaded("$ALICE/geant4_vmc/examples/macro/g4libs.C"))
    gROOT->LoadMacro("$ALICE/geant4_vmc/examples/macro/g4libs.C");
  gInterpreter->ProcessLine("g4libs()");

  // Menu
  TControlBar* menu = new TControlBar("vertical","Alice Geant4 menu");
  
  menu->AddButton("Geometry", "MakeGeometry()",  "Generate Root geometry file");
  menu->AddButton("Run",      "RunSimulation()",  "Process Alice run");
  menu->AddButton("Init",     "Init()",  "Initialize Alice for simulation");
  menu->AddButton("Geant4UI", "StartGeant4UI()","Go to Geant4 Interactive session");
  menu->AddButton("AGDD",     "GenerateAGDD()","Generate XML (AGDD) file with geometry description");
  //menu->AddButton("GDML",     "GenerateGDML()","Generate XML (GDML) file with geometry description");
  menu->AddButton("Quit",     "Quit()", "Quit aliroot");
  gROOT->SaveContext();
  
  cout << endl
       << "**************************************************************" << endl
       << "  To run simulation:"                                           << endl
       << "  First select <Geometry> to build geometry.root file."         << endl
       << "  Then re-run aliroot and select <Run> button"                  << endl
       << endl
       << "  The <Init> button is kept for debugging purposes,"            << endl
       << "  it itializes MonteCarlo but it does not itialize"             << endl
       << "  completely ALICE framework. That's why to run simulation,"    << endl
       << "  you have to re-run aliroot and select Run button."            << endl
       << endl
       << "  The menu enables to start Geant4 interactive session:"        << endl
       << "  Select <Geant4UI> button and use Geant4 interactive commands" << endl
       << "  To go back to Root UI, type exit."                            << endl
       << "**************************************************************" << endl
       << endl;
  
  menu->Show();
}

void MakeGeometry()
{  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(1);
  gAlice->Init("$ALICE_ROOT/macros/g4ConfigGeometry.C");
  
  // Generate geometry file
  //
  gGeoManager->Export("geometry.root");
  
  cout << endl
       << "Geometry file geometry.root has been generated." << endl
       << "You have to re-run aliroot and choose Run in g4menu." << endl;
       
  exit(0);     
}    


void Init()
{  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  gAlice->Init("$ALICE_ROOT/macros/g4Config.C");

  cout << endl
       << "Only MonteCarlo initialization has been performed. " << endl
       << "To run simulation you have to re-run aliroot and choose Run in g4menu." << endl;
}    


void RunSimulation()
{  
  AliSimulation sim("$ALICE_ROOT/macros/g4Config.C");
  sim.SetMakeDigits("");
  sim.SetMakeSDigits("");
  sim.SetRunHLT("");
  sim.SetNumberOfEvents(1);
  TStopwatch timer;
  timer.Start();
  sim.Run(1);
  timer.Stop();
  timer.Print();
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

void GenerateAGDD()
{
  if (gMC) {
    // release Root terminal control

    // go into non-raw term mode
    //Getlinem(kCleanUp, 0);
    
    // add test if gMC is TGeant4
    TGeant4* g4 = (TGeant4*)gMC;
    
    g4->ProcessGeantCommand("/vgm/generateAGDD");

    // new Root prompt
    //Getlinem(kInit, ((TRint*)gROOT->GetApplication())->GetPrompt());  
  }
  else {  
    cout << "Monte Carlo has not been yet created." << endl;
  } 
}        
/*
void GenerateGDML()
{
  if (gMC) {
    // release Root terminal control

    // go into non-raw term mode
    //Getlinem(kCleanUp, 0);
    
    // add test if gMC is TGeant4
    TGeant4* g4 = (TGeant4*)gMC;
    
    g4->ProcessGeantCommand("/vgm/generateGDML");

    // new Root prompt
    //Getlinem(kInit, ((TRint*)gROOT->GetApplication())->GetPrompt());  
  }
  else {  
    cout << "Monte Carlo has not been yet created." << endl;
  } 
}        
*/
void Quit()
{
  delete AliRunLoader::Instance();
  delete gAlice;
  
  exit(0);
}  
