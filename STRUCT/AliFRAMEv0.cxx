////////////////////////////////////////////////
//  space frame class                            /
////////////////////////////////////////////////

#include <stdio.h> 
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include "AliFRAMEv0.h"
#include "AliRun.h"
#include "stdlib.h"
#include "AliMC.h"
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
  printf("Create FRAMEv0 object");  
}

 
//___________________________________________
void AliFRAMEv0::CreateGeometry()
{
  printf("Create FRAMEv0 geometry ");
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

  AliMC* pMC=AliMC::GetMC();
  char *filetmp;
  const char *framename = "$(ALICE_ROOT)/Euclid/frame0399.euc";
  char topvol[5];
  printf("Create FRAMEv0 geometry ");
  
//
// The Space frame
  filetmp = gSystem->ExpandPathName(framename);
  FILE *file = fopen(filetmp,"r");
  delete [] filetmp;
  if(file) {
    fclose(file);
    printf(" Reading FRAME \n");
    gAlice->ReadEuclid(framename,12,topvol);
  } else {
    printf(" THE GEOM FILE %s DOES NOT EXIST !\n",framename);
    exit(1);
  }
//
// --- Place the FRAME ghost volume (B010) in its mother volume (ALIC)
//    and make it invisible
// 
//  AliMatrix(idrotm[2001],90.,0.,90.,90.,180.,0.);
  
  pMC->Gspos("B010",1,"ALIC",0,0,0,0,"ONLY");

  pMC->Gsatt("B010", "SEEN", 0);
}

 
//___________________________________________
void AliFRAMEv0::CreateMaterials()
{
  char *filetmp;
  printf("Create FRAMEv0 materials");
  const char *name = "$(ALICE_ROOT)/Euclid/frame.tme";
  filetmp = gSystem->ExpandPathName(name);
  FILE *file = fopen(filetmp,"r");
  delete [] filetmp;
  if(file) {
    fclose(file);
    gAlice->ReadEuclidMedia(name,12);
  } else {
    printf(" THE MEDIA FILE %s DOES NOT EXIST !\n",name);
    exit(1);
  }
}














