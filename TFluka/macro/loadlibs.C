// Macro for loading Geant4, Geant4 VMC and Flugg libraries
// Based on the very similar file g4libs.C in the geant4_vmc package

#include <iostream>

void loadlibs()
{
// Loads granular Geant4 libraries.
// Change the comment if global libraries are used.
// ---

  g4libs_granular();
  // g4libs_global();
}   

Bool_t isSet(const char* variable)
{
// Checks if the specified environment variable is set
// ---

  TString value = gSystem->Getenv(variable);
  if ( value != "") return true;
  
  return false;
}  

void g4libs_graphics() 
{
// Loads G4 graphics libraries, 
// external packages: graphics drivers, .. used by G4
// ---
  
  // Graphics configuration
  Bool_t isXm = !isSet("G4UI_NONE") && 
                (isSet("G4VIS_BUILD_OPENGLXM_DRIVER") ||
		 isSet("G4UI_BUILD_XM_SESSION"));
  Bool_t isGAG = !isSet("G4UI_NONE") && isSet("G4UI_USE_GAG");
  Bool_t isDAWN = !isSet("G4VIS_NONE");
  Bool_t isOpenGL = !isSet("G4VIS_NONE") &&
                    (isSet("G4VIS_BUILD_OPENGLX_DRIVER") ||
                     isSet("G4VIS_BUILD_OPENGLXM_DRIVER"));
  Bool_t isVRML = !isSet("G4VIS_NONE");
  Bool_t isRayTracer = !isSet("G4VIS_NONE");

  // Geant4 interfaces
  //
  if (isXm) {
    gSystem->Load("libXt");
    gSystem->Load("libXm");
  }
  gSystem->Load("libG4UIcommon");
  gSystem->Load("libG4UIbasic");
  if (isGAG) 
    gSystem->Load("libG4UIGAG");

  // Geant4 visualization
  //
  if (isOpenGL) {
    gSystem->Load("libGLU");
    gSystem->Load("libGL");
  }  
  gSystem->Load("libG4modeling");
  gSystem->Load("libG4vis_management");
  if (isDAWN)
    gSystem->Load("libG4FR");
  if (isOpenGL)
    gSystem->Load("libG4OpenGL");
  if (isVRML)
    gSystem->Load("libG4VRML");
  if (isRayTracer)
    gSystem->Load("libG4RayTracer");
}


void g4libs_granular()
{
// Loads G4 granular libraries and G4 VMC library. 
// external packages: CLHEP, graphics drivers, .. used by G4
// ---

  cout << "Loading Geant4 granular libraries ..." << endl;

  // CLHEP
  gSystem->Load("libCLHEP");

  // G4 categories

  // global
  gSystem->Load("libG4globman");  
  gSystem->Load("libG4hepnumerics");

//    // graphics_reps
//    gSystem->Load("libG4graphics_reps");   

  // intercoms
  gSystem->Load("libG4intercoms");

  // materials
  gSystem->Load("libG4materials");

  // geometry
  gSystem->Load("libG4geomver");
  gSystem->Load("libG4volumes");
  gSystem->Load("libG4magneticfield");
  // I.G.C.
  gSystem->Load("libFlugg");
  gSystem->Load("libG4geometrymng");  
  gSystem->Load("libG4geomBoolean");  
  gSystem->Load("libG4csg");  
  gSystem->Load("libG4step");
  gSystem->Load("libG4brep"); 
  gSystem->Load("libG4specsolids"); 
  gSystem->Load("libG4stepinterface");


  
  // particles  
  gSystem->Load("libG4partman");
  gSystem->Load("libG4bosons");   
  gSystem->Load("libG4baryons");  
  gSystem->Load("libG4ions");
  gSystem->Load("libG4mesons");
  gSystem->Load("libG4leptons");
  gSystem->Load("libG4shortlived");

  // track
  gSystem->Load("libG4track");

  // processes
  gSystem->Load("libG4procman");
  gSystem->Load("libG4parameterisation");
  gSystem->Load("libG4decay");  
  gSystem->Load("libG4emutils");  
  gSystem->Load("libG4emstandard");   
  gSystem->Load("libG4emlowenergy");  
  gSystem->Load("libG4muons");
  gSystem->Load("libG4xrays");
  gSystem->Load("libG4hadronic_xsect");
  gSystem->Load("libG4hadronic_mgt");   
  gSystem->Load("libG4hadronic_proc");
  gSystem->Load("libG4hadronic_util");
  gSystem->Load("libG4hadronic_man_gen");   
  gSystem->Load("libG4hadronic_util_gen");
  gSystem->Load("libG4hadronic_string_common");
  gSystem->Load("libG4hadronic_diffstring");  
  gSystem->Load("libG4hadronic_stringfrag");
  gSystem->Load("libG4hadronic_HE_gen");  
  gSystem->Load("libG4hadronic_kinetic");   
  gSystem->Load("libG4hadronic_qgstring");
  gSystem->Load("libG4hadronic_HE");  
  gSystem->Load("libG4hadronic_LE");  
  gSystem->Load("libG4hadronic_deex");
  gSystem->Load("libG4hadronic_preequ");  
  gSystem->Load("libG4hadronic_stop");
  gSystem->Load("libG4hadronic_neu");   
  gSystem->Load("libG4hadronic_iso");   
  gSystem->Load("libG4optical");
  gSystem->Load("libG4photolepton_hadron");
  gSystem->Load("libG4transportation");

  // tracking
  gSystem->Load("libG4tracking");

  // digits+hits  
  gSystem->Load("libG4hits");
  gSystem->Load("libG4digits");   
  gSystem->Load("libG4detector");   

  // event
  gSystem->Load("libG4event");  

  // readout
  gSystem->Load("libG4readout");
  
  // run
  gSystem->Load("libG4run");
  
  // g3tog4
  gSystem->Load("libG3toG4");   

  // interfaces and graphics
  g4libs_graphics();
  
  // geant4 mc
  gSystem->Load("libgeant4vmc");

  cout << "Loading Geant4 granular libraries ... finished" << endl;
}

void g4libs_global()
{
// Loads G4 global libraries, 
// external packages: CLHEP, graphics drivers, .. used by G4
// and Alice G4 libraries: AliGeant4, TGeant4
// ---

  cout << "Loading Geant4 global libraries ..." << endl;
 
   // CLHEP
  gSystem->Load("$(CLHEP_BASE_DIR)/lib/libCLHEP");

  // Geant4
  gSystem->Load("libG4global");
  gSystem->Load("libG4graphics_reps");
  gSystem->Load("libG4intercoms");
  gSystem->Load("libG4materials");
  gSystem->Load("libG4geometry");
  gSystem->Load("libG4particles");
  gSystem->Load("libG4track");
  gSystem->Load("libG4processes");
  gSystem->Load("libG4tracking");
  gSystem->Load("libG4digits+hits");
  gSystem->Load("libG4event");
  gSystem->Load("libG4readout");
  gSystem->Load("libG4run");
  gSystem->Load("libG3toG4");

  // interfaces and graphics
  g4libs_graphics();
 
  // geant4 mc
  gSystem->Load("libgeant4_vmc");

  cout << "Loading Geant4 global libraries ... finished" << endl;
}

