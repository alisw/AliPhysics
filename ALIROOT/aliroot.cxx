//////////////////////////////////////////////////////////////////////////
//                                                                      //
// aliroot                                                              //
//                                                                      //
// Main program used to create aliroot application.                     //
//                                                                      //
//                                                                      //
// To be able to communicate between the FORTRAN code of GEANT and the  //
// ROOT data structure, a number of interface routines have been        //
// developed.                                                           //
//Begin_Html
/*
<img src="gif/newg.gif">
*/
//End_Html
//////////////////////////////////////////////////////////////////////////

//Standard Root includes
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <AliRun.h>

#if defined __linux
//On linux Fortran wants this, so we give to it!
int xargv=0;
int xargc=0;
#endif

#if defined WIN32 
  extern "C" int __fastflag=0; 
  extern "C" int _pctype=0; 
  extern "C" int __mb_cur_max=0; 
#endif 

int gcbank_[3000000];

//Initialise the Root environment
 extern void InitGui();
 VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
 TROOT root("galice","The Alice/ROOT Interface", initfuncs);

//_____________________________________________________________________________
int main(int argc, char **argv)
{
  //
  // gAlice main program.
  // It creates the main Geant3 and AliRun objects.
  //
  // The Hits are written out after each track in a ROOT Tree structure TreeH,
  // of the file galice.root. There is one such tree per event. The kinematics
  // of all the particles that produce hits, together with their genealogy up
  // to the primary tracks is stared in the galice.root file in an other tree
  // TreeK of which exists one per event. Finally the information of the events
  // in the run is stored in the same file in the tree TreeE, containing the
  // run and event number, the number of vertices, tracks and primary tracks
  // in the event.
  //
  new AliRun("gAlice"," The Alice Geant3-based MonteCarlo");
  
  // Start interactive geant
  
  TRint *theApp = new TRint("aliroot", &argc, argv, 0, 0);
  
  // --- Initialisation of the GALICE package ---
  theApp->Run();
  
  // --- Simulation of all events ---
  return(0);
}


