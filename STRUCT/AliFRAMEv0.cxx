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
Revision 1.10  1999/10/06 20:56:55  fca
Introducing new frame 1099

Revision 1.9  1999/09/29 09:24:30  fca
Introduction of the Copyright and cvs Log

*/

////////////////////////////////////////////////
//  space frame class                            /
////////////////////////////////////////////////

#include "AliFRAMEv0.h"
#include "AliRun.h"
#include "TSystem.h"
 
ClassImp(AliFRAMEv0)
 
//_____________________________________________________________________________
AliFRAMEv0::AliFRAMEv0()
{
}

//_____________________________________________________________________________
AliFRAMEv0::AliFRAMEv0(const char *name, const char *title)
  : AliFRAME(name,title)
{
  printf("Create FRAMEv0 object\n");  
  fEuclidGeometry="$(ALICE_ROOT)/Euclid/frame1099h.euc";
  fEuclidMaterial="$(ALICE_ROOT)/Euclid/frame.tme";
}

 
//___________________________________________
void AliFRAMEv0::CreateGeometry()
{
//Begin_Html
/*
<img src="picts/frame.gif">
*/
//End_Html


//Begin_Html
/*
<img src="picts/tree_frame.gif">
*/
//End_Html

  char *filetmp;
  char topvol[5];
  
//
// The Space frame
  filetmp = gSystem->ExpandPathName(fEuclidGeometry.Data());
  FILE *file = fopen(filetmp,"r");
  delete [] filetmp;
  if(file) {
    fclose(file);
    printf(" Reading FRAME geometry\n");
    gAlice->ReadEuclid(fEuclidGeometry.Data(),this,topvol);
  } else {
    Warning("CreateGeometry","The Euclid file %s does not exist!\n",
	    fEuclidGeometry.Data());
    exit(1);
  }
//
// --- Place the FRAME ghost volume (B010) in its mother volume (ALIC)
//    and make it invisible
// 
//  AliMatrix(idrotm[2001],90.,0.,90.,90.,180.,0.);

  gMC->Gspos(topvol,1,"ALIC",0,0,0,0,"ONLY");

  gMC->Gsatt(topvol, "SEEN", 0);
}

 
//___________________________________________
void AliFRAMEv0::CreateMaterials()
{
  char *filetmp;
  printf("Create FRAMEv0 materials\n");
  filetmp = gSystem->ExpandPathName(fEuclidMaterial.Data());
  FILE *file = fopen(filetmp,"r");
  delete [] filetmp;
  if(file) {
    fclose(file);
    gAlice->ReadEuclidMedia(fEuclidMaterial.Data(),this);
  } else {
    Warning("CreateMaterials","The material file %s does not exist!\n",
	    fEuclidMaterial.Data());
    exit(1);
  }
}

//_____________________________________________________________________________
void AliFRAMEv0::Init()
{
  //
  // Initialise the module after the geometry has been defined
  //

  printf("**************************************"
	 " FRAME "
	 "**************************************\n");
  printf("\n     Version 0 of FRAME initialised, "
	 "with openings for PHOS and RICH\n\n");
  printf("**************************************"
	 " FRAME "
	 "**************************************\n");

}













