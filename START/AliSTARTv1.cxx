/////////////////////////////////////////////////////////////////////
//                                                                 //
// START ( T-zero) detector  version 0                        //
//
//Begin Html       
/*
<img src="gif/AliSTARTv0Class.gif">
*/
//End Html
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TGeometry.h>
#include <TTUBE.h>
#include <TNode.h>

#include "AliSTARTv0.h"
#include "AliRun.h"
#include "AliMC.h"

#include <iostream.h>
#include <fstream.h>

#include "TGeant3.h"
#include <stdlib.h>

ClassImp(AliSTARTv0)

//--------------------------------------------------------------------
AliSTARTv1::AliSTARTv1(const char *name, const char *title):
 AliSTART(name,title)
{
  //
  // Standart constructor for START Detector version 0
  //
  fIdSens1=0;
//  setBufferSize(128000);
}
//-------------------------------------------------------------------------
void AliSTARTv1::CreateGeometry()
{
  //
  // Create the geometry of START Detector version 0
  //
  // begin Html
  /*
   <img src="gif/AliSTARTv0.gif">
  */
  //

  AliMC* pMC = AliMC::GetMC();

  Int_t *idtmed = gAlice->Idtmed();
  
  Int_t is;
  Int_t idrotm[999];
  Float_t x,y,z;

  Float_t pstart[3]={5.,10.2,5.3};
  //  Float_t pscin[3]={0.,2.54/2.,1.5};
  Float_t ppmt[3]={0.,1.3,3.5};
  Float_t pdivider[3]={0.,1.2,1.75};
  Float_t pdiv2[3]={0.,1.2,1.25};
  Float_t pdiv1[3]={0.6,1.2,0.5};
  Float_t ptop[3]={0.,1.3,1.5};
  Float_t pbot[3]={0.6,1.2,0.1};
  Float_t pglass[3]={1.2,1.3,2.};
  Float_t pcer[3]={0.9,1.1,0.09};
  Float_t psteel[3]={0.9,1.1,0.01};
  Float_t ppins[3]={0.6,1.2,0.014};
  Float_t phole[3]={0.6,1.2,0.015};
  Float_t pknob[3]={0.5,0.6,0.4};
  Float_t pknob_vac[3]={0.,0.5,0.4};
  Float_t pknob_bot[3]={0.,0.6,0.05};
  Float_t pribber[3] = {0.,1.2,2.413/2.};
  Float_t presist[3] = {0.,1.2,0.087/2.};

  Float_t zdet=75.;
 //-------------------------------------------------------------------
 //  START volume 
 //-------------------------------------------------------------------

    AliMatrix(idrotm[901], 90., 0., 90., 90., 180., 0.);
    pMC->Gsvolu("STRT","TUBE",idtmed[2101-1],pstart,3);
    pMC->Gspos("STRT",1,"ALIC",0.,0.,zdet,0,"ONLY");
    pMC->Gspos("STRT",2,"ALIC",0.,0.,-zdet,idrotm[901],"ONLY");

//START interior
    //    pMC->Gsvolu("SCIN","TUBE",idtmed[2102-1],pscin,3);
    pMC->Gsvolu("PMT ","TUBE",idtmed[2103-1],ppmt,3);     
    pMC->Gsvolu("DIVI","TUBE",idtmed[2103-1],pdivider,3);     

// first ring: 13 units of Scintillator+PMT+divider
    for (is=1; is<=13; is++)
    {  
      x=6*TMath::Sin(is*2*3.1415/13);
      y=6*TMath::Cos(is*2*3.1415/13);
      z=-pstart[2]+ppmt[2];
      printf(" is %d, z PMT %f\n",is,z);
      //      pMC->Gspos("SCIN",is,"STRT",x,y,z,0,"ONLY");
      //      z=z+pscin[2]+ppmt[2];
      pMC->Gspos("PMT ",is,"STRT",x,y,z,0,"ONLY");
      z=ppmt[2]+pdivider[2];
      printf(" is %d, z Divider %f\n",is,z);
      pMC->Gspos("DIVI",is,"STRT",x,y,z,0,"ONLY");
    }
//second ring: 20 units of Scintillator+PMT+divider
    for (is=14; is<=33;is++)  
    {  
      x=8.8*TMath::Sin(2.*3.1415/26+(is-13)*2*3.1415/20);
      y=8.8*TMath::Cos(2.*3.1315/26+(is-13)*2*3.1415/20);
      z=-pstart[2]+ppmt[2];
      //      pMC->Gspos("SCIN",is,"STRT",x,y,z,0,"ONLY");
      //z=z+pscin[2]+ppmt[2];
      pMC->Gspos("PMT ",is,"STRT",x,y,z,0,"ONLY");
      z=ppmt[2]+pdiv2[2];
      pMC->Gspos("DIVI",is,"STRT",x,y,z,0,"ONLY");
    }
// PMT
      
// Entry window (glass)
      pMC->Gsvolu("PTOP","TUBE",idtmed[2106-1],ptop,3);
      z=-ppmt[2]+ptop[2];
      pMC->Gspos("PTOP",1,"PMT ",0,0,z,0,"ONLY");
      printf("Z PTOP %f -ppmt[2] %f ptop[2] %f\n",z,-ppmt[2],ptop[2]);
// Bottom glass
      pMC->Gsvolu("PBOT","TUBE",idtmed[2106-1],pbot,3);
      z=ppmt[2]-pbot[2];
      printf("Z bottom %f\n",z);
      pMC->Gspos("PBOT",1,"PMT ",0,0,z,0,"ONLY");
// Side cylinder glass
      pMC->Gsvolu("POUT","TUBE",idtmed[2106-1],pglass,3);
      z=ppmt[2]-pglass[2];
      printf("Z glass %f\n",z);
      pMC->Gspos("POUT",1,"PMT ",0,0,z,0,"ONLY");
//PMT electrodes support structure
      pMC->Gsvolu("PCER","TUBE",idtmed[2104-1],pcer,3);
      pMC->Gsvolu("PSTE","TUBE",idtmed[2108-1],psteel,3);
      z=-ppmt[2]+2*ptop[2]+0.3;;
      printf("Z Cer 1 %f\n",z);
      for (is=1; is<=15; is++)
      {
         z=z+psteel[2]+pcer[2];
         pMC->Gspos("PCER",is,"PMT",0,0,z,0,"ONLY");
         z=z+psteel[2]+pcer[2];
         pMC->Gspos("PSTE",is,"PMT",0,0,z,0,"ONLY");
       }

// Divider
// Knob at the bottom of PMT baloon
      pMC->Gsvolu("KNOB","TUBE",idtmed[2106-1],pknob,3);
      z=-pdivider[2]+pknob[2];
//      printf("zknob %f\n",z);
      pMC->Gspos("KNOB",1,"DIVI",0,0,z,0,"ONLY");
      pMC->Gsvolu("KNBO","TUBE",idtmed[2106-1],pknob_bot,3);
      z=-pdivider[2]+2*pknob[2]+pknob_bot[2];
//      printf("knobbot %f\n",z);
      pMC->Gspos("KNBO",1,"DIVI ",0,0,z,0,"ONLY");
      pMC->Gsvolu("KNVA","TUBE",idtmed[2106-1],pknob_vac,3);
      z=-pdivider[2]+pknob_vac[2];
//      printf("knobvac %f\n",z);
      pMC->Gspos("KNVA",1,"DIVI",0,0,z,0,"ONLY");
 //Steel pins + pin holes
      pMC->Gsvolu("PINS","TUBE",idtmed[2108-1],ppins,3);
      z=-pdivider[2]+ppins[2];
      pMC->Gspos("PINS",1,"DIVI",0,0,z,0,"ONLY");
      pMC->Gsvolu("HOLE","TUBE",idtmed[2111-1],phole,3);
      z=-pdivider[2]+2*ppins[2]+phole[2];
      pMC->Gspos("HOLE",1,"DIVI",0,0,z,0,"ONLY");

//Socket
      pMC->Gsvolu("DIV1","TUBE",idtmed[2104-1],pdiv1,3);
      z=-pdivider[2]+pdiv1[2];
      pMC->Gspos("DIV1",1,"DIVI",0,0,z,0,"ONLY");
//Resistors
      pMC->Gsvolu("DIV2","TUBE",idtmed[2101-1],pdiv2,3);
      z=pdivider[2]-pdiv2[2];
      pMC->Gspos("DIV2",1,"DIVI",0,0,z,0,"ONLY");
      pMC->Gsvolu("DRES","TUBE",idtmed[2104-1],presist,3);
      z=-pdiv2[2]+presist[2];
      pMC->Gspos("DRES",1,"DIV2",0,0,z,0,"ONLY");
      pMC->Gsvolu("DRIB","TUBE",idtmed[2109-1],pribber,3);
      z=pdiv2[2]-pribber[2];
      pMC->Gspos("DRIB",1,"DIV2",0,0,z,0,"ONLY");
      printf("z DRIB %f\n",z);

}    
//------------------------------------------------------------------------
void AliSTARTv1::CreateMaterials()
{
   Int_t ISXFLD   = gAlice->Field()->Integ();
   Float_t SXMGMX = gAlice->Field()->Max();
   Float_t a,z,d,radl,absl,buf[1];
   Int_t nbuf;

// Scintillator CH
   //    Float_t ascin[2]={1.01,12.01};
   // Float_t zscin[2]={1,6};
   // Float_t wscin[2]={1,1};
   // Float_t denscin=1.03;
// PMT glass SiO2
    Float_t aglass[2]={28.0855,15.9994};
    Float_t zglass[2]={14.,8.};
    Float_t wglass[2]={1.,2.};
    Float_t dglass=2.65;
// Ceramic   97.2% Al2O3 , 2.8% SiO2
    Float_t acer[2],zcer[2],wcer[2]={0.972,0.028};
    Float_t aal2o3[2]  = { 26.981539,15.9994 };
    Float_t zal2o3[2]  = { 13.,8. };
    Float_t wal2o3[2]  = { 2.,3. };
    Float_t denscer  = 3.6;

// Brass 80% Cu, 20% Zn
   Float_t abrass[2] = {63.546,65.39};
   Float_t zbrass[2] = {29,30};
   Float_t wbrass[2] = {0.8,0.2};
   Float_t denbrass=8.96;

//Ribber C6H12S
   Float_t aribber[3] = {12.,1.,32.};
   Float_t zribber[3] = {6.,1.,16.};
   Float_t wribber[3] = {6.,12.,1.};
   Float_t denribber=0.8;
 
   
   AliMC* pMC = AliMC::GetMC();

   Int_t *idtmed = gAlice->Idtmed();
   Int_t imat;

//*** Definition Of avaible START materials ***
 AliMaterial(0, "START Steel$", 55.850,26.,7.87,1.76,999);
 AliMaterial(1, "START Vacuum$", 1.e-16,1.e-16,1.e-16,1.e16,999);
 AliMaterial(2, "START Air$", 14.61, 7.3, .001205, 30423.,999); 

 AliMixture( 3, "Al2O3   $", aal2o3, zal2o3, denscer, -2, wcer);
 AliMixture( 4, "PMT glass   $",aglass,zglass,dglass,-2,wglass);
  char namate[21];
  pMC->Gfmate((*fIdmate)[3], namate, a, z, d, radl, absl, buf, nbuf);
  acer[0]=a;
  zcer[0]=z;
  pMC->Gfmate((*fIdmate)[4], namate, a, z, d, radl, absl, buf, nbuf);
  acer[1]=a;
  zcer[1]=z;
  
 AliMixture( 9, "Ceramic    $", acer, zcer, denscer, 2, wcer);
 // AliMixture( 5, "Scintillator$",ascin,zscin,denscin,-2,wscin);
 AliMixture( 6, "Brass    $", abrass, zbrass, denbrass, 2, wbrass);

 AliMixture( 7, "Ribber $",aribber,zribber,denribber,-3,wribber);
 
 

//**
 AliMedium(2101, "START Air$", 2, 0, ISXFLD, SXMGMX, 10., .1, 1., .003, .003);
 // AliMedium(2102, "Scintillator$", 5, 1, ISXFLD, SXMGMX, 10., .01, 1., .003, .003);
 AliMedium(2103, "Vacuum$", 1, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
 AliMedium(2104, "Ceramic$", 9, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
 AliMedium(2106, "Glass$", 4, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
 AliMedium(2108, "Steel$", 0, 0, ISXFLD, SXMGMX, 1., .001, 1., .001, .001);
 AliMedium(2111, "Brass  $", 6, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
 AliMedium(2109, "Ribber  $", 7, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);


}
//---------------------------------------------------------------------
void AliSTARTv1::DrawDetector()
{
//
// Draw a shaded view of the Forward multiplicity detector version 0
//

AliMC* pMC = AliMC::GetMC();

//Set ALIC mother transparent
pMC->Gsatt("ALIC","SEEN",0);
//
//Set volumes visible
pMC->Gsatt("STRT","SEEN",0);
//pMC->Gsatt("SCIN","SEEN",1);
pMC->Gsatt("PMT ","SEEN",1);
pMC->Gsatt("DIVI","SEEN",1);
//
pMC->Gdopt("hide","on");
pMC->Gdopt("shad","on");
pMC->SetClipBox(".");
pMC->SetClipBox("*",0,1000,-1000,1000,-1000,1000);
pMC->DefaultRange();
pMC->Gdraw("alic",40,30,0,12,9.5,.2,0.2);
pMC->Gdhead(1111,"T-Zero detector");
pMC->Gdopt("hide","off");
}
//-------------------------------------------------------------------
void AliSTARTv1::Init()
{
// Initialises version 0 of the Forward Multiplicity Detector
//
AliMC* pMC=AliMC::GetMC();
//Int_t *idtmed  = gAlice->Idtmed();
AliSTART::Init();
fIdSens1=pMC->VolId("PTOP");
printf("*** START version 0 initialized ***\n");
 
}

//-------------------------------------------------------------------

void AliSTARTv1::StepManager()
{
  //
  // Called for every step in the START Detector
  //
  Int_t id,copy,copy1,i;
  static Float_t hits[7];
  static Float_t edep;
  static Int_t vol[2];

  TClonesArray &lhits = *fHits;
  AliMC* pMC=AliMC::GetMC();
  TGeant3* geant3= (TGeant3*)pMC; 
  
  if(!pMC->TrackAlive()) return; // particle has disappeared
  Float_t charge = pMC->TrackCharge();
  if(TMath::Abs(charge)<=0.) return; //take only charged particles
  //    geant3->Gpcxyz();
    id=pMC->CurrentVol(0,copy);

 
//  printf("igauto %d\n",pMC->ctrak->igauto);
//  printf("pMC->ckine->ipart %d",pMC->ckine->ipart);
// Check the sensetive volume
   if(id==fIdSens1 )
   {
     if(pMC->TrackEntering())
     {
       pMC->CurrentVolOff(2,0,copy);
       vol[0]=copy;
       pMC->CurrentVolOff(1,0,copy1);
       vol[1]=copy1;
       pMC->TrackPosition(hits);
       Float_t etot=pMC->Etot();
//    printf("pMC->ctrak->ekin %f\n",pMC->ctrak->Ekin);
       hits[4]=etot;
       Int_t part= pMC->TrackPid();
       hits[5]=part;
       Float_t ttime=geant3->Gctrak()->tofg;
       hits[6]=ttime*1e9;
       edep=0;
       //       Float_t xV=geant3->Gckine()->vert[0];
       //Float_t yV=geant3->Gckine()->vert[1];
       //Float_t zV=geant3->Gckine()->vert[2];
       //Float_t tl=pMC -> TrackLength();
       //      if(hits[6]<2.4){     
       //	    for (i=0; i<=6; i++){
       //	    printf(" HITS on START entr %f\n",hits[i]);} 
       //}
     }
     if(pMC->TrackInside())
     {
       Float_t de=pMC->Edep(); 
       edep=edep+de;
//       printf ("E deposition %f\n",edep);
//    for (i=0; i<=6; i++){
//    printf(" HITS on START inside %f\n",hits[i]); } 
     }
     if(pMC->TrackExiting())
    {
       Float_t de=pMC->Edep(); 
       edep=edep+de;
       hits[3]=edep*1e3;
       //       for (i=0; i<=6; i++){
       //	 printf(" HITS on START Exit %f\n",hits[i]); } 
       //for (i=0; i<=1; i++) { printf("START vol %d\n",vol[i]);}
    
    new(lhits[fNhits++]) AliSTARThit(fIshunt,gAlice->CurrentTrack(),vol,hits);      }
    }
//---------------------------------------------------------------------
    }
 //}












