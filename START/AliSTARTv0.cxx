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

//#include "TGeant3.h"
#include <stdlib.h>

ClassImp(AliSTARTv0)

//--------------------------------------------------------------------
AliSTARTv0::AliSTARTv0(const char *name, const char *title):
 AliSTART(name,title)
{
  //
  // Standart constructor for START Detector version 0
  //
  fIdSens1=0;
//  setBufferSize(128000);
}
//-------------------------------------------------------------------------
void AliSTARTv0::CreateGeometry()
{
  //
  // Create the geometry of START Detector version 0
  //
  // begin Html
  /*
   <img src="gif/AliSTARTv0.gif">
  */
  //


  Int_t *idtmed = fIdtmed->GetArray()-2099;
  
  Int_t is;
  Int_t idrotm[999];
  Float_t x,y,z;

  Float_t pstart[3]={4.5,10.7,5.3};
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
    gMC->Gsvolu("STRT","TUBE",idtmed[2101-1],pstart,3);
    gMC->Gspos("STRT",1,"ALIC",0.,0.,zdet,0,"ONLY");
    gMC->Gspos("STRT",2,"ALIC",0.,0.,-zdet,idrotm[901],"ONLY");

//START interior
    //    gMC->Gsvolu("SCIN","TUBE",idtmed[2102-1],pscin,3);
    gMC->Gsvolu("PMT ","TUBE",idtmed[2103-1],ppmt,3);     
    gMC->Gsvolu("DIVI","TUBE",idtmed[2103-1],pdivider,3);     

// first ring: 13 units of Scintillator+PMT+divider
    for (is=1; is<=13; is++)
    {  
      x=6.5*TMath::Sin(is*2*3.1415/13);
      y=6.5*TMath::Cos(is*2*3.1415/13);
      z=-pstart[2]+ppmt[2];
      printf(" is %d x %f y %f r%f\n",is,x,y,sqrt(x*x+y*y));
      //      gMC->Gspos("SCIN",is,"STRT",x,y,z,0,"ONLY");
      //      z=z+pscin[2]+ppmt[2];
      gMC->Gspos("PMT ",is,"STRT",x,y,z,0,"ONLY");
      z=ppmt[2]+pdivider[2];
      printf(" is %d, z Divider %f\n",is,z);
      gMC->Gspos("DIVI",is,"STRT",x,y,z,0,"ONLY");
    }
//second ring: 20 units of Scintillator+PMT+divider
    for (is=14; is<=33;is++)  
    {  
      x=9.3*TMath::Sin(2.*3.1415/26+(is-13)*2*3.1415/20);
      y=9.3*TMath::Cos(2.*3.1315/26+(is-13)*2*3.1415/20);
      z=-pstart[2]+ppmt[2];
      //      gMC->Gspos("SCIN",is,"STRT",x,y,z,0,"ONLY");
      //z=z+pscin[2]+ppmt[2];
      gMC->Gspos("PMT ",is,"STRT",x,y,z,0,"ONLY");
      z=ppmt[2]+pdiv2[2];
      gMC->Gspos("DIVI",is,"STRT",x,y,z,0,"ONLY");
    }
// PMT
      
// Entry window (glass)
      gMC->Gsvolu("PTOP","TUBE",idtmed[2106-1],ptop,3);
      z=-ppmt[2]+ptop[2];
      gMC->Gspos("PTOP",1,"PMT ",0,0,z,0,"ONLY");
 //     printf("Z PTOP %f -ppmt[2] %f ptop[2] %f\n",z,-ppmt[2],ptop[2]);
// Bottom glass
      gMC->Gsvolu("PBOT","TUBE",idtmed[2106-1],pbot,3);
      z=ppmt[2]-pbot[2];
      printf("Z bottom %f\n",z);
      gMC->Gspos("PBOT",1,"PMT ",0,0,z,0,"ONLY");
// Side cylinder glass
      gMC->Gsvolu("POUT","TUBE",idtmed[2106-1],pglass,3);
      z=ppmt[2]-pglass[2];
//      printf("Z glass %f\n",z);
      gMC->Gspos("POUT",1,"PMT ",0,0,z,0,"ONLY");
//PMT electrodes support structure
      gMC->Gsvolu("PCER","TUBE",idtmed[2104-1],pcer,3);
      gMC->Gsvolu("PSTE","TUBE",idtmed[2108-1],psteel,3);
      z=-ppmt[2]+2*ptop[2]+0.3;;
//      printf("Z Cer 1 %f\n",z);
      for (is=1; is<=15; is++)
      {
         z=z+psteel[2]+pcer[2];
         gMC->Gspos("PCER",is,"PMT",0,0,z,0,"ONLY");
         z=z+psteel[2]+pcer[2];
         gMC->Gspos("PSTE",is,"PMT",0,0,z,0,"ONLY");
       }

// Divider
// Knob at the bottom of PMT baloon
      gMC->Gsvolu("KNOB","TUBE",idtmed[2106-1],pknob,3);
      z=-pdivider[2]+pknob[2];
//      printf("zknob %f\n",z);
      gMC->Gspos("KNOB",1,"DIVI",0,0,z,0,"ONLY");
      gMC->Gsvolu("KNBO","TUBE",idtmed[2106-1],pknob_bot,3);
      z=-pdivider[2]+2*pknob[2]+pknob_bot[2];
//      printf("knobbot %f\n",z);
      gMC->Gspos("KNBO",1,"DIVI ",0,0,z,0,"ONLY");
      gMC->Gsvolu("KNVA","TUBE",idtmed[2106-1],pknob_vac,3);
      z=-pdivider[2]+pknob_vac[2];
//      printf("knobvac %f\n",z);
      gMC->Gspos("KNVA",1,"DIVI",0,0,z,0,"ONLY");
 //Steel pins + pin holes
      gMC->Gsvolu("PINS","TUBE",idtmed[2108-1],ppins,3);
      z=-pdivider[2]+ppins[2];
      gMC->Gspos("PINS",1,"DIVI",0,0,z,0,"ONLY");
      gMC->Gsvolu("HOLE","TUBE",idtmed[2111-1],phole,3);
      z=-pdivider[2]+2*ppins[2]+phole[2];
      gMC->Gspos("HOLE",1,"DIVI",0,0,z,0,"ONLY");

//Socket
      gMC->Gsvolu("DIV1","TUBE",idtmed[2104-1],pdiv1,3);
      z=-pdivider[2]+pdiv1[2];
      gMC->Gspos("DIV1",1,"DIVI",0,0,z,0,"ONLY");
//Resistors
      gMC->Gsvolu("DIV2","TUBE",idtmed[2101-1],pdiv2,3);
      z=pdivider[2]-pdiv2[2];
      gMC->Gspos("DIV2",1,"DIVI",0,0,z,0,"ONLY");
      gMC->Gsvolu("DRES","TUBE",idtmed[2104-1],presist,3);
      z=-pdiv2[2]+presist[2];
      gMC->Gspos("DRES",1,"DIV2",0,0,z,0,"ONLY");
      gMC->Gsvolu("DRIB","TUBE",idtmed[2109-1],pribber,3);
      z=pdiv2[2]-pribber[2];
      gMC->Gspos("DRIB",1,"DIV2",0,0,z,0,"ONLY");
//      printf("z DRIB %f\n",z);

}    
//------------------------------------------------------------------------
void AliSTARTv0::CreateMaterials()
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
 
   
//   AliMC* gMC = AliMC::GetMC();

 //  Int_t *idtmed = gAlice->Idtmed();
 //  Int_t imat;

//*** Definition Of avaible START materials ***
 AliMaterial(0, "START Steel$", 55.850,26.,7.87,1.76,999);
 AliMaterial(1, "START Vacuum$", 1.e-16,1.e-16,1.e-16,1.e16,999);
 AliMaterial(2, "START Air$", 14.61, 7.3, .001205, 30423.,999); 

 AliMixture( 3, "Al2O3   $", aal2o3, zal2o3, denscer, -2, wal2o3);
 AliMixture( 4, "PMT glass   $",aglass,zglass,dglass,-2,wglass);
  char namate[21];
  gMC->Gfmate((*fIdmate)[3], namate, a, z, d, radl, absl, buf, nbuf);
  acer[0]=a;
  zcer[0]=z;
  gMC->Gfmate((*fIdmate)[4], namate, a, z, d, radl, absl, buf, nbuf);
  acer[1]=a;
  zcer[1]=z;
  
 AliMixture( 9, "Ceramic    $", acer, zcer, denscer, 2, wcer);
 // AliMixture( 5, "Scintillator$",ascin,zscin,denscin,-2,wscin);
 AliMixture( 6, "Brass    $", abrass, zbrass, denbrass, 2, wbrass);

 AliMixture( 7, "Ribber $",aribber,zribber,denribber,-3,wribber);
 
 

//**
 AliMedium(1, "START Air$", 2, 0, ISXFLD, SXMGMX, 10., .1, 1., .003, .003);
 // AliMedium(2102, "Scintillator$", 5, 1, ISXFLD, SXMGMX, 10., .01, 1., .003, .003);
 AliMedium(3, "Vacuum$", 1, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
 AliMedium(4, "Ceramic$", 9, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
 AliMedium(6, "Glass$", 4, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
 AliMedium(8, "Steel$", 0, 0, ISXFLD, SXMGMX, 1., .001, 1., .001, .001);
 AliMedium(11, "Brass  $", 6, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
 AliMedium(9, "Ribber  $", 7, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);


}
//---------------------------------------------------------------------
void AliSTARTv0::DrawModule()
{
//
// Draw a shaded view of the Forward multiplicity detector version 0
//

//Set ALIC mother transparent
gMC->Gsatt("ALIC","SEEN",0);
//
//Set volumes visible
gMC->Gsatt("STRT","SEEN",0);
//gMC->Gsatt("SCIN","SEEN",1);
gMC->Gsatt("PMT ","SEEN",1);
gMC->Gsatt("DIVI","SEEN",1);
//
gMC->Gdopt("hide","on");
gMC->Gdopt("shad","on");
gMC->SetClipBox(".");
gMC->SetClipBox("*",0,1000,-1000,1000,-1000,1000);
gMC->DefaultRange();
gMC->Gdraw("alic",40,30,0,12,9.5,.2,0.2);
gMC->Gdhead(1111,"T-Zero detector");
gMC->Gdopt("hide","off");
}
//-------------------------------------------------------------------
void AliSTARTv0::Init()
{
  // Initialises version 0 of the Forward Multiplicity Detector
  //
  
  //Int_t *idtmed  = gAlice->Idtmed();
  AliSTART::Init();
  fIdSens1=gMC->VolId("PTOP");
  printf("*** START version 0 initialized ***\n");
 
}

//-------------------------------------------------------------------

void AliSTARTv0::StepManager()
{
  //
  // Called for every step in the START Detector
  //
  Int_t id,copy,copy1;
  static Float_t hits[7];
  static Float_t edep;
  static Int_t vol[2];
  TLorentzVector pos;

  TClonesArray &lhits = *fHits;
  //AliMC* gMC=AliMC::GetMC();
  //  TGeant3* geant3= (TGeant3*)gMC; 
  
  if(!gMC->IsTrackAlive()) return; // particle has disappeared
  Float_t charge = gMC->TrackCharge();
  if(TMath::Abs(charge)<=0.) return; //take only charged particles
 //    gMC->Gpcxyz();
    id=gMC->CurrentVolID(copy);

 
//  printf("igauto %d\n",gMC->ctrak->igauto);
//  printf("gMC->ckine->ipart %d",gMC->ckine->ipart);
// Check the sensetive volume
   if(id==fIdSens1 )
   {
     if(gMC->IsTrackEntering())
     {
       gMC->CurrentVolOffID(2,copy);
       vol[0]=copy;
       gMC->CurrentVolOffID(1,copy1);
       vol[1]=copy1;
       gMC->TrackPosition(pos);
       hits[0] = pos[0];
       hits[1] = pos[1];
       hits[2] = pos[2];
       Float_t etot=gMC->Etot();
       hits[4]=etot;
       Int_t part= gMC->TrackPid();
       hits[5]=part;
       Float_t ttime=gMC->TrackTime();
       hits[6]=ttime*1e9;
       edep=0;
       //       Float_t xV=geant3->Gckine()->vert[0];
       //Float_t yV=geant3->Gckine()->vert[1];
       //Float_t zV=geant3->Gckine()->vert[2];
       //Float_t tl=gMC -> TrackLength();
       //      if(hits[6]<2.4){     
       //	    for (i=0; i<=6; i++){
       //	    printf(" HITS on START entr %f\n",hits[i]);} 
       //}
     }
     if(gMC->IsTrackInside())
     {
       Float_t de=gMC->Edep(); 
       edep=edep+de;
//       printf ("E deposition %f\n",edep);
//    for (i=0; i<=6; i++){
//    printf(" HITS on START inside %f\n",hits[i]); } 
     }
     if(gMC->IsTrackExiting())
    {
       Float_t de=gMC->Edep(); 
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






