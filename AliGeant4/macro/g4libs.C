#include <iostream.h>

static Bool_t isGeant4 = false;
static Bool_t isSteer = false;
static Bool_t isDetector = false;

void g4libs()
{
  // g4libs_global();
  g4libs_granular();
}   

void g4libs_global()
{
// Loads G4 global libraries, 
// external packages: CLHEP, graphics drivers, .. used by G4
// and Alice G4 libraries: AliGeant4, TGeant4
// ---
  if (!isGeant4) {

    // CLHEP
    gSystem->Load("$(CLHEP_BASE_DIR)/lib/libCLHEP");

    // Geant4
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4global");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4graphics_reps");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4intercoms");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4materials");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4geometry");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4particles");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4track");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4processes");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4tracking");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4digits+hits");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4event");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4readout");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4run");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG3toG4");

    // Geant4 interfaces
    //gSystem->Load("/usr/X11R6/lib/libXt");
    //gSystem->Load("/usr/local/lib/libXm");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4UIcommon");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4UIbasic");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4UIGAG");

    // Geant4 visualization
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4vis_management");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4modeling");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4FR");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4OpenGL");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4VRML");
    
    // TGeant4, AliGeant4
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTGeant4");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libAliGeant4");
    
    isGeant4 = true;
    cout << "Geant4 global libraries have been loaded." << endl;
  }    
}

void g4libs_granular()
{
// Loads G4 granular libraries, 
// external packages: CLHEP, graphics drivers, .. used by G4
// and Alice G4 libraries: AliGeant4, TGeant4
// ---
  if (!isGeant4) {

    // CLHEP
    gSystem->Load("$(CLHEP_BASE_DIR)/lib/libCLHEP");

    // G4 categories

    // global
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4globman");                      
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hepnumerics");

    // graphics_reps
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4graphics_reps");       

    // intercoms
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4intercoms");

    // materials
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4materials");

    // geometry
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4magneticfield");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4volumes");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4geometrymng");    
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4geomBoolean");    
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4csg");                  
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4step");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4brep"); 
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4specsolids"); 
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4stepinterface");
  
    // particles          
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4partman");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4bosons");             
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4baryons");              
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4ions");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4mesons");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4leptons");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4shortlived");

    // track
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4track");

    // processes
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4procman");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4parameterisation");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4decay");                  
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4emutils");              
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4emstandard");           
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4emlowenergy");          
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4muons");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4xrays");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_xsect");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_mgt");         
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_proc");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_util");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_man_gen");     
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_util_gen");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_string_common");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_diffstring");  
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_stringfrag");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_HE_gen");      
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_kinetic");     
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_preequ");      
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_qgstring");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_HE");          
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_LE");          
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_deex");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_stop");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_neu");         
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hadronic_iso");         
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4optical");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4photolepton_hadron");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4transportation");

    // tracking
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4tracking");

    // digits+hits  
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4hits");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4digits");               
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4detector");             

    // event
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4event");                

    // readout
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4readout");
  
    // run
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4run");
  
    // g3tog4
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG3toG4");                     

    // interfaces
    //gSystem->Load("/usr/X11R6/lib/libXt");
    //gSystem->Load("/usr/local/lib/libXm");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4UIcommon");             
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4UIbasic");              
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4UIGAG");                

    // visualisation
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4modeling");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4vis_management");
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4FR");                   
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4OpenGL");               
    gSystem->Load("$(G4INSTALL)/lib/$(G4SYSTEM)/libG4VRML");                 
  
    // TGeant4, AliGeant4
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTGeant4");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libAliGeant4");

    isGeant4 = true;
    cout << "Geant4 granular libraries have been loaded." << endl;
  }
}

void steerlibs() {
// Loads AliRoot steer libraries
// ---
  if (!isSteer) {
  
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTEER");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libdummypythia6");
    gSystem->Load("$(ROOTSYS)/lib/libEGPythia6");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libEVGEN");

    isSteer = true;
    cout << "AliRoot steer libraries have been loaded." << endl;
  }  
}  

void detlibs() {
// Load AliRoot modules libraries
// ---
  if (!isDetector) {

    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libminicern");
    // minicern required by MUON
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTRUCT");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libFMD");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libMUON");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libPHOS");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libPMD");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libRICH");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTOF");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPC");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTRD");
    //gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libZDC");
       // requires symbols from geant3 lib
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libITS");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libCASTOR");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTART");
    
    isDetector = true;
    cout << "AliRoot detectors libraries have been loaded." << endl;
  }  
}  
