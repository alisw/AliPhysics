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

/*
$Log$
Revision 1.13  2002/12/03 09:03:06  hristov
Changes needed on Itanium (F.Carminati)

Revision 1.12  2002/10/14 14:55:33  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.10.6.1  2002/08/28 15:06:49  alibrary
Updating to v3-09-01

Revision 1.11  2002/08/18 20:20:48  hristov
InitGui removed (J.Chudoba)

Revision 1.10  2002/02/18 14:26:14  hristov
#include <AliConfig.h> removed

Revision 1.9  2002/01/09 17:05:03  morsch
Increase memory allocated for ZEBRA.

Revision 1.8  2001/10/04 15:32:36  hristov
Instantiation of AliConfig removed

Revision 1.7  2001/05/22 11:39:48  buncic
Restored proper name for top level AliRoot folder.

Revision 1.6  2001/05/21 17:22:50  buncic
Fixed problem with missing AliConfig while reading galice.root

Revision 1.5  2001/05/16 14:57:07  alibrary
New files for folders and Stack

Revision 1.4  2000/12/20 08:39:37  fca
Support for Cerenkov and process list in Virtual MC

Revision 1.3  1999/09/29 09:24:07  fca
Introduction of the Copyright and cvs Log

*/

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
  
  // Create new configuration 
  
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");
    
  // Start interactive geant
  
  TRint *theApp = new TRint("aliroot", &argc, argv);
  
  // --- Start the event loop ---
  theApp->Run();
  
  return(0);
}


