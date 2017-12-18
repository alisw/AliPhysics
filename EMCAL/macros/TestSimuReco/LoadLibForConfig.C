///
/// \file LoadLibForConfig.C
/// \ingroup EMCAL_TestSimRec
/// \brief Macro loading needed libraries for simulation.
///
/// Example macro to be executed right before TestEMCALSimulation.C in case of Root6
/// Cppy paste from LoadLibs() method in the Config.C
/// It loads the necessayr libraries for Geant3 or Geant4 simulation.
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS). 
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>

#include <TSystem.h>
#include <TString.h>
#include <TApplication.h>
#include <TROOT.h>
#include <string>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TThread.h>

#endif

#ifndef TRANSPORTMODEL
#define TRANSPORTMODEL

TString kTransportModel = "Geant3";

#endif //TRANSPORTMODEL

void loadg4libs()
{
  
  // Load Pythia related libraries                                                                
  gSystem->Load("liblhapdf");      // Parton density functions                                 
  gSystem->Load("libEGPythia6");   // TGenerator interface                                     
  gSystem->Load("libpythia6");     // Pythia                                                   
  gSystem->Load("libAliPythia6");  // ALICE specific
  
  /// Macro function for loading Geant4 libraries
  /// from list of libraries build by geant4-config --libs 
  
  FILE* pipe = gSystem->OpenPipe("geant4-config --libs", "r");
  std::string all_lines;
  char line[1000];
  while ( fgets(line, sizeof(line), pipe ) != NULL ) {
    all_lines += line;
  }
  
  TString all_lines_t(all_lines.data());
  all_lines_t.Remove(all_lines_t.First('\n'));
  //cout << all_lines_t.Data() << endl;
  TObjArray* libs = all_lines_t.Tokenize(" ");
  
  //TString dynamicPath = gSystem->GetDynamicPath();
  for (Int_t i=libs->GetEntriesFast()-1; i>=0; i-- ) {
    TString addPath = ((TObjString*)libs->At(i))->GetString();
    if (addPath.BeginsWith("-L")) {
      addPath.Remove(0,2);
      addPath.ReplaceAll("\"", "");
      //cout << "Adding dynamic path " << addPath.Data() << endl;
      gSystem->AddDynamicPath(addPath.Data());
    }
  }
  
  cout << libs->GetEntriesFast() << endl;
  for (Int_t i=libs->GetEntriesFast()-1; i>=0; i-- ) {
    TString lib = ((TObjString*)libs->At(i))->GetString();
    lib.ReplaceAll("-l", "lib");
    //cout << "Loading |" << lib.Data() << "|" << endl; 
    if(lib.BeginsWith("lib"))
      gSystem->Load(lib.Data());
  } 
  
  gSystem->SetFPEMask(0); 
}   

Bool_t isLibrary(const char* libName)
{
  /// Helper function which testes the existence of the given library
  /// \param libName  The library name
  
  if (TString(gSystem->DynamicPathName(libName, kTRUE)) != TString(""))
    return kTRUE;
  else  
    return kFALSE;
}    

Bool_t isBatch()
{
  /// Helper function which testes if Root was started in batch mode
  
  for ( Int_t i=0; i<gApplication->Argc(); ++i ) 
    if ( TString(gROOT->GetApplication()->Argv(i)) == "-b" ) return true;
  
  return false;
}    

Bool_t isSet(const char* variable)
{
  /// Helper function which checks if the specified environment variable 
  /// is set.
  /// \param variable  The environment variable name
  
  TString value = gSystem->Getenv(variable);
  if ( value != "") return true;
  
  return false;
}  

Bool_t isMT()
{
  /// Macro function for detecting if Geant4 libraries
  /// are built in multi-threading mode via 
  /// geant4-config --has-feature multithreading
  
  FILE* pipe = gSystem->OpenPipe("geant4-config  --has-feature multithreading", "r");
  char line[10];
  fgets(line, sizeof(line), pipe);
  TString answer = line;
  answer.Remove(answer.First('\n'));
  
  return ( answer == "yes");
}  

void vgmlibs()
{ 
  /// Function for loading VGM libraries.
  
  if ( isSet("USE_VGM") ) { 
    cout << "Loading VGM libraries ... " << endl;
    gSystem->Load("libClhepVGM");
    gSystem->Load("libBaseVGM");
    gSystem->Load("libGeant4GM");
    gSystem->Load("libRootGM");
    gSystem->Load("libXmlVGM");
  }  
}

///
/// Main execution method
/// \param transport define here the transport model Geant3 or Geant4
///
void LoadLibForConfig(TString transport = "Geant4")
{
  
  kTransportModel = transport;
  
  // Load Pythia related libraries                                                                
  gSystem->Load("liblhapdf");      // Parton density functions                                 
  gSystem->Load("libEGPythia6");   // TGenerator interface                                     
  gSystem->Load("libpythia6");     // Pythia                                                   
  gSystem->Load("libAliPythia6");  // ALICE specific
  
  // Load Geant3 related libraries   
  if(kTransportModel == "Geant3")
  {
    gSystem->Load("libgeant321");
  }
  else if(kTransportModel == "Geant4")
  {
    // Load Geant4 libraries
    cout << "Loading Geant4 libraries (using geant4-config) ..." << endl;
    loadg4libs();
    
    // VGM librares
    vgmlibs();
    
    // G4Root library (if available)
    cout << "Loading g4root library ..." << endl;
    gSystem->Load("libg4root");
    
    // Geant4 VMC library
    cout << "Loading geant4vmc library ..." << endl;
    gSystem->Load("libgeant4vmc");
    
    // Geant4 VMC GUI library 
    // (if available and Root is not running in batch mode)
    if ( isLibrary("libgeant4vmc_gui") && ! isBatch() ) {
      cout << "Loading geant4vmc_gui library ... " << endl;
      gSystem->Load("libgeant4vmc_gui");
    }  
    
    // mtroot library (make optional)
    cout << "Loading mtroot library ..." << endl;
    gSystem->Load("libmtroot");
    
    // initialize Root threading  
    if ( isMT() ) {
      TThread::Initialize();
      gInterpreter->SetProcessLineLock(false);
    }
  }
  else printf("Transport model <%s> has no extra known libs\n",kTransportModel.Data());
  
}
