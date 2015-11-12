/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

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
<img src="picts/newg.gif">
*/
//End_Html
//////////////////////////////////////////////////////////////////////////

//Standard Root includes
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <AliRun.h>
#include "Riostream.h"
#include "ARVersion.h"
// STD
#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <functional>
#include <AliLog.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <string>

#if defined __linux
//On linux Fortran wants this, so we give to it!
int xargv=0;
int xargc=0;
#endif

#ifdef FORTRAN_G95
extern "C" void g95_runtime_start();
#endif

#if defined WIN32 
  extern "C" int __fastflag=0; 
  extern "C" int _pctype=0; 
  extern "C" int __mb_cur_max=0; 
#endif 

using std::cout;
using std::endl;
using std::exception;
using std::cerr;

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

  //RS: use maximum stack size
  const long kMB = 1024L * 1024L;
  rlim_t newStackSize = 8L;   // new default stack size
  // check if it was overriden by env.var
  const char* envStack = getenv("ALIROOT_STACK_SIZE"); // expect size in MB
  if (envStack) {
    rlim_t envSiz = atoll(envStack);
    if (envSiz<8) {
      AliFatalGeneralF("AliRoot","ALIROOT_STACK_SIZE=%s must request stack size in MB, 8MB at least",envStack);
    }
    newStackSize = envSiz;
  }
  newStackSize *= kMB;
  struct rlimit rl;
  int result;
  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0) {
    rlim_t oldss = rl.rlim_cur;
    if (newStackSize > rl.rlim_max) {
      AliWarningGeneralF("AliRoot","Requestested new stack size %lld > hard limit %lld MB",newStackSize/kMB,rl.rlim_max/kMB);
      newStackSize = rl.rlim_max;
    }
    if (rl.rlim_cur < newStackSize) {
      rl.rlim_cur = newStackSize;
      result = setrlimit(RLIMIT_STACK, &rl);
      if (result != 0)	fprintf(stderr, "setrlimit returned result = %d\n", result);
      else AliInfoGeneralF("AliRoot","Set stack size from %lld to %lld MB",oldss/kMB,rl.rlim_cur/kMB);
    }
  }
  
  for ( int i = 1; i < argc; ++i ) 
  {
    TString argument(argv[i]);
    
    if (argument=="--version")
    {      
      cout << "aliroot " << ALIROOT_REVISION << " " << ALIROOT_VERSION << endl;
      return 0;
    }    
  }
  
  // Create new configuration 
  
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  AliLog::GetRootLogger();  // force AliLog to initialize at start time
  // Start interactive geant
  
  TRint *theApp = new TRint("aliroot", &argc, argv);
#ifdef FORTRAN_G95
  g95_runtime_start();
#endif
  
  // --- Start the event loop ---
  //  run.Info("Start AliRoot");
  try{
    theApp->Run();  
  } catch(const exception &e){
    cerr << "Excception catched" << e.what() << endl;
    AliLog::Message(0, "Exception catched", "ROOT", NULL, "Exception catched", NULL, 0);
  }
  
  return(0);
}

