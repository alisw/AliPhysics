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
<img src="gif/pipe.gif">
*/
//End_Html


//Begin_Html
/*
<img src="gif/tree_pipe.gif">
*/
//End_Html

  const char *pipename = "$(ALICE_ROOT)/Euclid/bpipeb.euc";
  const char *pumpname = "$(ALICE_ROOT)/Euclid/bpumpa.euc";
  char *filtmp;
  char topvol[5];
  printf("Create PIPEv0 geometry ");
  
  Int_t idrotm[2099];

  AliMC* pMC = AliMC::GetMC();
  
//
// The peam pipe up to the Front Absorber
  filtmp=gSystem->ExpandPathName(pipename);
  FILE *file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    printf(" Reading PIPE \n");
    gAlice->ReadEuclid(pipename,20,topvol);
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
    gAlice->ReadEuclid(pumpname,20,topvol);
  } else {
    printf(" THE GEOM FILE %s DOES NOT EXIST !\n",pumpname);
    exit(1);
  }
//
// --- Place the PIPE ghost volume (QBPM) in its mother volume (ALIC)
//    and make it invisible
// 
  AliMatrix(idrotm[2001],90.,0.,90.,90.,180.,0.);
  
  pMC->Gspos("QBPM",1,"ALIC",0,0,0,idrotm[2001],"ONLY");
//
//    PLACE ION PUMP (QIPM) AT Z=-385.
//
  pMC->Gspos("QIPM",1,"ALIC",0,0,-385,idrotm[2001],"ONLY");

  pMC->Gsatt("QIPM", "SEEN", 0);
  pMC->Gsatt("QBPM", "SEEN", 0);
  pMC->Gsatt("QB20", "SEEN", 0);
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
    gAlice->ReadEuclidMedia(name,20);
  } else {
    printf(" THE MEDIA FILE %s DOES NOT EXIST !\n",name);
    exit(1);
  }
}














