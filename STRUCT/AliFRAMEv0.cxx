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
  printf("Create FRAMEv0 object\n");  
  fEuclidGeometry="$(ALICE_ROOT)/Euclid/frame0799.euc";
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














