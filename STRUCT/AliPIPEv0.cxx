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
Revision 1.7  2000/02/23 16:25:24  fca
AliVMC and AliGeant3 classes introduced
ReadEuclid moved from AliRun to AliModule

Revision 1.6  1999/09/29 09:24:30  fca
Introduction of the Copyright and cvs Log

*/

////////////////////////////////////////////////
//  Beam pipe class                            /
////////////////////////////////////////////////

#include "AliPIPEv0.h"
#include "AliRun.h"
#include "TSystem.h"
 
ClassImp(AliPIPEv0)
 
//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0()
{
// Constructor
}

//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0(const char *name, const char *title)
  : AliPIPE(name,title)
{
// Constructor
}

 
//___________________________________________
void AliPIPEv0::CreateGeometry()
{
  printf("Create PIPEv0 geometry\n ");
//Begin_Html
/*
<img src="picts/pipe.gif">
*/
//End_Html


//Begin_Html
/*
<img src="picts/tree_pipe.gif">
*/
//End_Html

  const char *kPipeName = "$(ALICE_ROOT)/Euclid/bpipeb.euc";
  const char *kPumpName = "$(ALICE_ROOT)/Euclid/bpumpa.euc";
  char *filtmp;
  char topvol[5];
  printf("Create PIPEv0 geometry ");
  
  Int_t idrotm[2099];

//
// The peam pipe up to the Front Absorber
  filtmp=gSystem->ExpandPathName(kPipeName);
  FILE *file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    printf(" Reading PIPE \n");
    ReadEuclid(kPipeName,topvol);
  } else {
    printf(" THE GEOM FILE %s DOES NOT EXIST !\n",kPipeName);
    exit(1);
  }
//
// The Ion Pump
  filtmp=gSystem->ExpandPathName(kPumpName);
  file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    printf(" Reading PUMP \n");
    ReadEuclid(kPumpName,topvol);
  } else {
    printf(" THE GEOM FILE %s DOES NOT EXIST !\n",kPumpName);
    exit(1);
  }
//
// --- Place the PIPE ghost volume (QBPM) in its mother volume (ALIC)
//    and make it invisible
// 
  AliMatrix(idrotm[2001],90.,0.,90.,90.,180.,0.);
  
  gMC->Gspos("QBPM",1,"ALIC",0,0,0,idrotm[2001],"ONLY");
//
//    PLACE ION PUMP (QIPM) AT Z=-385.
//
  gMC->Gspos("QIPM",1,"ALIC",0,0,-385,idrotm[2001],"ONLY");

  gMC->Gsatt("QIPM", "SEEN", 0);
  gMC->Gsatt("QBPM", "SEEN", 0);
  gMC->Gsatt("QB20", "SEEN", 0);
}

 
//___________________________________________
void AliPIPEv0::DrawModule()
{
// Set drawing options
    ;
}

//___________________________________________
void AliPIPEv0::CreateMaterials()
{
// Create materials and media from Euclid file
  printf("Create PIPEv0 materials\n");
  const char *kName = "$(ALICE_ROOT)/Euclid/pipe.tme";
  char *filtmp;
  filtmp=gSystem->ExpandPathName(kName);
  FILE *file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    ReadEuclidMedia(kName);
  } else {
    printf(" THE MEDIA FILE %s DOES NOT EXIST !\n",kName);
    exit(1);
  }
}














