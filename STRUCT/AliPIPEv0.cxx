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
}

//_____________________________________________________________________________
AliPIPEv0::AliPIPEv0(const char *name, const char *title)
  : AliPIPE(name,title)
{
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

  const char *pipename = "$(ALICE_ROOT)/Euclid/bpipeb.euc";
  const char *pumpname = "$(ALICE_ROOT)/Euclid/bpumpa.euc";
  char *filtmp;
  char topvol[5];
  printf("Create PIPEv0 geometry ");
  
  Int_t idrotm[2099];

//
// The peam pipe up to the Front Absorber
  filtmp=gSystem->ExpandPathName(pipename);
  FILE *file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    printf(" Reading PIPE \n");
    gAlice->ReadEuclid(pipename,this,topvol);
  } else {
    printf(" THE GEOM FILE %s DOES NOT EXIST !\n",pipename);
    exit(1);
  }
//
// The Ion Pump
  filtmp=gSystem->ExpandPathName(pumpname);
  file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    printf(" Reading PUMP \n");
    gAlice->ReadEuclid(pumpname,this,topvol);
  } else {
    printf(" THE GEOM FILE %s DOES NOT EXIST !\n",pumpname);
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
}

//___________________________________________
void AliPIPEv0::CreateMaterials()
{
  printf("Create PIPEv0 materials\n");
  const char *name = "$(ALICE_ROOT)/Euclid/pipe.tme";
  char *filtmp;
  filtmp=gSystem->ExpandPathName(name);
  FILE *file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    gAlice->ReadEuclidMedia(name,this);
  } else {
    printf(" THE MEDIA FILE %s DOES NOT EXIST !\n",name);
    exit(1);
  }
}














