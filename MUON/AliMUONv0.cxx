/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 

#include "AliMUONv0.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h"
#include "AliConst.h" 

#define trig trig_

extern "C" void type_of_call trig(float (*)[4], float (*)[4], int& iflag);

ClassImp(AliMUONv0)
 
//___________________________________________
AliMUONv0::AliMUONv0() : AliMUON()
{
    fChambers = 0;
}
 
//___________________________________________
AliMUONv0::AliMUONv0(const char *name, const char *title)
       : AliMUON(name,title)
{
//
//  z-Positions of Chambers
    const Float_t zch[7]={511., 686., 971., 1245., 1445., 1610, 1710.};
//
//  inner diameter
    const Float_t dmi[7]={ 35.,  47.,  67.,   86.,  100., 96., 101.96};
//
//  outer diameter
    const Float_t dma[7]={183., 245., 346.,  442.,  513., 824., 874.};
//
    Int_t k;

    fChambers = new TObjArray(14);
    
    for (Int_t i=0; i<7; i++) {
	for (Int_t j=0; j< 2; j++) {
//
//    
//    Default Parameters for Muon Tracking Stations
	    k=2*i+j;
//
	    (*fChambers)[k] = new AliMUONchamber();
	    AliMUONchamber* chamber = (AliMUONchamber*) (*fChambers)[k];
	    chamber->SetGid(0);
	    chamber->SetZPOS(zch[i]);
//
	    chamber->InitGeo(zch[i]);
	    chamber->frMin=dmi[i]/2.;
	    chamber->frMax=dma[i]/2.;
//
	} // Chamber j in 
    }     // Station i
    fMaxStepGas=0.01; 
    fMaxStepAlu=0.1; 
    fMaxDestepGas=-1;
    fMaxDestepAlu=-1;
}

//___________________________________________
void AliMUONv0::Trigger(Float_t (*x)[4], Float_t (*y)[4], Int_t& iflag)
{
  trig(x,y,iflag);
}

//___________________________________________
void AliMUONv0::CreateGeometry()
{
    Int_t *idtmed = fIdtmed->GetArray()-1099;
//
//   Note: all chambers have the same structure, which could be 
//   easily parameterised. This was intentionally not done in order
//   to give a starting point for the implementation of the actual 
//   design of each station. 
//
//   Distance between Stations
     const Float_t dstation = 8.;
//
     Float_t bpar[3];
     Float_t tpar[3];
     Float_t tspar[5], pgpar[10];
     Float_t zpos1, zpos2, zfpos;
     Float_t dframep=3.;
     Float_t dframez=0.9;
     Float_t dr, rMin;
//
//   Rotation matrices in the x-y plane  
     Int_t idrotm[1199];
//   phi=   0 deg
     AliMatrix(idrotm[1100],  90.,   0., 90.,  90., 0., 0.);
//   phi=  90 deg
     AliMatrix(idrotm[1101],  90.,  90., 90., 180., 0., 0.);
//   phi= 180 deg
     AliMatrix(idrotm[1102],  90., 180., 90., 270., 0., 0.);
//   phi= 270 deg
     AliMatrix(idrotm[1103],  90., 270., 90.,   0., 0., 0.);
//
     Float_t phi=2*TMath::Pi()/12/2;

//
//   pointer to the current chamber
     AliMUONchamber *iChamber;
//********************************************************************
//                            Station 1                             **
//********************************************************************
//  CONCENTRIC
     iChamber=(AliMUONchamber*) (*fChambers)[0];
     zpos1=iChamber->ZPosition()-dstation/2; 
     zpos2=zpos1+dstation; 
     zfpos=-(iChamber->fdGas+dframez)/2;
     
//
//   Mother volume
     tpar[0] = iChamber->frMin-dframep; 
     tpar[1] = (iChamber->frMax+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/2;

     gMC->Gsvolu("C01M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gsvolu("C02M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gspos("C01M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C02M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
// Aluminium frames
// Outer frames
     pgpar[0] = 360/12/2;
     pgpar[1] = 360.;
     pgpar[2] = 12.;
     pgpar[3] =   2;
     pgpar[4] = -dframez/2;
     pgpar[5] = iChamber->frMax;
     pgpar[6] = pgpar[5]+dframep;
     pgpar[7] = +dframez/2;
     pgpar[8] = pgpar[5];
     pgpar[9] = pgpar[6];
     gMC->Gsvolu("C01O", "PGON", idtmed[1103], pgpar, 10);
     gMC->Gsvolu("C02O", "PGON", idtmed[1103], pgpar, 10);
     gMC->Gspos("C01O",1,"C01M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C01O",2,"C01M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C02O",1,"C02M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C02O",2,"C02M", 0.,0.,+zfpos,  0,"ONLY");
//
// Inner frame
     tpar[0]= iChamber->frMin-dframep;
     tpar[1]= iChamber->frMin;
     tpar[2]= dframez/2;
     gMC->Gsvolu("C01I", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gsvolu("C02I", "TUBE", idtmed[1103], tpar, 3);

     gMC->Gspos("C01I",1,"C01M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C01I",2,"C01M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C02I",1,"C02M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C02I",2,"C02M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     bpar[0] = (iChamber->frMax - iChamber->frMin)/2;
     bpar[1] = dframep/2;
     bpar[2] = dframez/2;
     gMC->Gsvolu("C01B", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C02B", "BOX", idtmed[1103], bpar, 3);

     gMC->Gspos("C01B",1,"C01M", +iChamber->frMin+bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C01B",2,"C01M", -iChamber->frMin-bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C01B",3,"C01M", 0, +iChamber->frMin+bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C01B",4,"C01M", 0, -iChamber->frMin-bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C01B",5,"C01M", +iChamber->frMin+bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C01B",6,"C01M", -iChamber->frMin-bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C01B",7,"C01M", 0, +iChamber->frMin+bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C01B",8,"C01M", 0, -iChamber->frMin-bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");

     gMC->Gspos("C02B",1,"C02M", +iChamber->frMin+bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C02B",2,"C02M", -iChamber->frMin-bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C02B",3,"C02M", 0, +iChamber->frMin+bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C02B",4,"C02M", 0, -iChamber->frMin-bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C02B",5,"C02M", +iChamber->frMin+bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C02B",6,"C02M", -iChamber->frMin-bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C02B",7,"C02M", 0, +iChamber->frMin+bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C02B",8,"C02M", 0, -iChamber->frMin-bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");

//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->frMin+dframep*2;
     tpar[1]= iChamber->frMax-dframep*2;
     tpar[2] = (iChamber->fdGas+iChamber->fdAlu)/2;
     gMC->Gsvolu("C01A", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gsvolu("C02A", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gspos("C01A", 1, "C01M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C02A", 1, "C02M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->fdGas;
     tpar[2] = iChamber->fdGas/2;
     gMC->Gsvolu("C01G", "TUBE", idtmed[1105], tpar, 3);
     gMC->Gsvolu("C02G", "TUBE", idtmed[1105], tpar, 3);
     gMC->Gspos("C01G", 1, "C01A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C02G", 1, "C02A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     dr = (iChamber->frMax - iChamber->frMin);
     bpar[0] = TMath::Sqrt(dr*dr-dframep*dframep/4)/2;
     bpar[1] = dframep/2;
     bpar[2] = iChamber->fdGas/2;
     gMC->Gsvolu("C01F", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C02F", "BOX", idtmed[1103], bpar, 3);

     gMC->Gspos("C01F",1,"C01G", +iChamber->frMin+bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C01F",2,"C01G", -iChamber->frMin-bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C01F",3,"C01G", 0, +iChamber->frMin+bpar[0] , 0, 
idrotm[1101],"ONLY");
     gMC->Gspos("C01F",4,"C01G", 0, -iChamber->frMin-bpar[0] , 0, 
idrotm[1101],"ONLY");

     gMC->Gspos("C02F",1,"C02G", +iChamber->frMin+bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C02F",2,"C02G", -iChamber->frMin-bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C02F",3,"C02G", 0, +iChamber->frMin+bpar[0] , 0, 
idrotm[1101],"ONLY");
     gMC->Gspos("C02F",4,"C02G", 0, -iChamber->frMin-bpar[0] , 0, 
idrotm[1101],"ONLY");
//
//
//********************************************************************
//                            Station 2                             **
//********************************************************************
     iChamber=(AliMUONchamber*) (*fChambers)[2];
     zpos1=iChamber->ZPosition()-dstation/2; 
     zpos2=zpos1+dstation; 
     zfpos=-(iChamber->fdGas+dframez)/2;
     
//
//   Mother volume
     tpar[0] = iChamber->frMin-dframep; 
     tpar[1] = (iChamber->frMax+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/2;

     gMC->Gsvolu("C03M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gsvolu("C04M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gspos("C03M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C04M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
// Aluminium frames
// Outer frames
     pgpar[0] = 360/12/2;
     pgpar[1] = 360.;
     pgpar[2] = 12.;
     pgpar[3] =   2;
     pgpar[4] = -dframez/2;
     pgpar[5] = iChamber->frMax;
     pgpar[6] = pgpar[5]+dframep;
     pgpar[7] = +dframez/2;
     pgpar[8] = pgpar[5];
     pgpar[9] = pgpar[6];
     gMC->Gsvolu("C03O", "PGON", idtmed[1103], pgpar, 10);
     gMC->Gsvolu("C04O", "PGON", idtmed[1103], pgpar, 10);
     gMC->Gspos("C03O",1,"C03M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C03O",2,"C03M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C04O",1,"C04M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C04O",2,"C04M", 0.,0.,+zfpos,  0,"ONLY");
//
// Inner frame
     tpar[0]= iChamber->frMin-dframep;
     tpar[1]= iChamber->frMin;
     tpar[2]= dframez/2;
     gMC->Gsvolu("C03I", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gsvolu("C04I", "TUBE", idtmed[1103], tpar, 3);

     gMC->Gspos("C03I",1,"C03M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C03I",2,"C03M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C04I",1,"C04M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C04I",2,"C04M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     bpar[0] = (iChamber->frMax - iChamber->frMin)/2;
     bpar[1] = dframep/2;
     bpar[2] = dframez/2;
     gMC->Gsvolu("C03B", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C04B", "BOX", idtmed[1103], bpar, 3);

     gMC->Gspos("C03B",1,"C03M", +iChamber->frMin+bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C03B",2,"C03M", -iChamber->frMin-bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C03B",3,"C03M", 0, +iChamber->frMin+bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C03B",4,"C03M", 0, -iChamber->frMin-bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C03B",5,"C03M", +iChamber->frMin+bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C03B",6,"C03M", -iChamber->frMin-bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C03B",7,"C03M", 0, +iChamber->frMin+bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C03B",8,"C03M", 0, -iChamber->frMin-bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");

     gMC->Gspos("C04B",1,"C04M", +iChamber->frMin+bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C04B",2,"C04M", -iChamber->frMin-bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C04B",3,"C04M", 0, +iChamber->frMin+bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C04B",4,"C04M", 0, -iChamber->frMin-bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C04B",5,"C04M", +iChamber->frMin+bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C04B",6,"C04M", -iChamber->frMin-bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C04B",7,"C04M", 0, +iChamber->frMin+bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C04B",8,"C04M", 0, -iChamber->frMin-bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");

//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->frMin+dframep*2;
     tpar[1]= iChamber->frMax-dframep*2;
     tpar[2] = (iChamber->fdGas+iChamber->fdAlu)/2;
     gMC->Gsvolu("C03A", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gsvolu("C04A", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gspos("C03A", 1, "C03M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C04A", 1, "C04M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->fdGas;
     tpar[2] = iChamber->fdGas/2;
     gMC->Gsvolu("C03G", "TUBE", idtmed[1105], tpar, 3);
     gMC->Gsvolu("C04G", "TUBE", idtmed[1105], tpar, 3);
     gMC->Gspos("C03G", 1, "C03A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C04G", 1, "C04A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     dr = (iChamber->frMax - iChamber->frMin);
     bpar[0] = TMath::Sqrt(dr*dr-dframep*dframep/4)/2;
     bpar[1] = dframep/2;
     bpar[2] = iChamber->fdGas/2;
     gMC->Gsvolu("C03F", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C04F", "BOX", idtmed[1103], bpar, 3);

     gMC->Gspos("C03F",1,"C03G", +iChamber->frMin+bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C03F",2,"C03G", -iChamber->frMin-bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C03F",3,"C03G", 0, +iChamber->frMin+bpar[0] , 0, 
idrotm[1101],"ONLY");
     gMC->Gspos("C03F",4,"C03G", 0, -iChamber->frMin-bpar[0] , 0, 
idrotm[1101],"ONLY");

     gMC->Gspos("C04F",1,"C04G", +iChamber->frMin+bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C04F",2,"C04G", -iChamber->frMin-bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C04F",3,"C04G", 0, +iChamber->frMin+bpar[0] , 0, 
idrotm[1101],"ONLY");
     gMC->Gspos("C04F",4,"C04G", 0, -iChamber->frMin-bpar[0] , 0, 
idrotm[1101],"ONLY");

//********************************************************************
//                            Station 3                             **
//********************************************************************
//  CONCENTRIC
     iChamber=(AliMUONchamber*) (*fChambers)[4];
     zpos1=iChamber->ZPosition(); // 975-13.75
     zpos2=zpos1 // +dstation;
                    +27.5;
//
//   Mother volume
     tpar[0] = iChamber->frMin; 
     tpar[1]= TMath::Sqrt(iChamber->frMax*iChamber->frMax + dframep*dframep) ;
     tpar[2] = // 3.; 
            5.325*2;
     gMC->Gsvolu("C05M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gsvolu("C06M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gspos("C05M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C06M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
//
//   Mother volume for one quadrant
     tspar[0]= iChamber->frMin;
     tspar[1]= TMath::Sqrt(iChamber->frMax*iChamber->frMax + dframep*dframep) ;
     tspar[2]= // dframez; 
                 5.325;
     tspar[3] = 0.-TMath::ATan2(dframep,iChamber->frMin)*180/kPI;
     tspar[4] = 90.+TMath::ATan2(dframep,iChamber->frMin)*180/kPI;
     gMC->Gsvolu("C05Q", "TUBS", idtmed[1100], tspar, 5);
     gMC->Gsvolu("C06Q", "TUBS", idtmed[1100], tspar, 5);
//   Position the four quadrants
     gMC->Gspos("C05Q",1,"C05M", 0., 0., 5.325, idrotm[1100], "ONLY");
     gMC->Gspos("C05Q",2,"C05M", 0., 0.,-5.325, idrotm[1101], "ONLY");
     gMC->Gspos("C05Q",3,"C05M", 0., 0., 5.325, idrotm[1102], "ONLY");
     gMC->Gspos("C05Q",4,"C05M", 0., 0.,-5.325, idrotm[1103], "ONLY");

     gMC->Gspos("C06Q",1,"C06M", 0., 0., 5.325, idrotm[1100], "ONLY");
     gMC->Gspos("C06Q",2,"C06M", 0., 0.,-5.325, idrotm[1101], "ONLY");
     gMC->Gspos("C06Q",3,"C06M", 0., 0., 5.325, idrotm[1102], "ONLY");
     gMC->Gspos("C06Q",4,"C06M", 0., 0.,-5.325, idrotm[1103], "ONLY");
// Aluminium frames
// Outer frame
     tspar[0]= iChamber->frMax-dframep*2;
     tspar[1]= iChamber->frMax;
     tspar[3] = 0.;
     tspar[4] = 90.;
     gMC->Gsvolu("C05O", "TUBS", idtmed[1100], tspar, 5);
     gMC->Gsvolu("C06O", "TUBS", idtmed[1100], tspar, 5);
     gMC->Gspos("C05O",1,"C05Q", 0.,0.,0.,  0,"ONLY");
     gMC->Gspos("C06O",1,"C06Q", 0.,0.,0.,  0,"ONLY");
//
// Inner frame
     tspar[0]= iChamber->frMin;
     tspar[1]= iChamber->frMin+dframep*2;
     gMC->Gsvolu("C05I", "TUBS", idtmed[1100], tspar, 5);
     gMC->Gsvolu("C06I", "TUBS", idtmed[1100], tspar, 5);
     gMC->Gspos("C05I",1,"C05Q", 0.,0.,0.,  0,"ONLY");
     gMC->Gspos("C06I",1,"C06Q", 0.,0.,0.,  0,"ONLY");
//
// Boundary half frame
     bpar[0] = (iChamber->frMax - iChamber->frMin)/2;
     bpar[1] = dframep/2;
     bpar[2] = 5.325;
     gMC->Gsvolu("C05B", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C06B", "BOX", idtmed[1103], bpar, 3);
//place 2 boudaries
     gMC->Gspos("C05B",1,"C05Q", iChamber->frMin+bpar[0] ,-bpar[1],0.,  idrotm[1100],"ONLY");
     gMC->Gspos("C05B",2,"C05Q", -bpar[1],iChamber->frMin+bpar[0] ,0.,  idrotm[1101],"ONLY");
     gMC->Gspos("C06B",1,"C06Q", iChamber->frMin+bpar[0] ,-bpar[1],0.,  idrotm[1100],"ONLY");
     gMC->Gspos("C06B",2,"C06Q", -bpar[1],iChamber->frMin+bpar[0] ,0.,  idrotm[1101],"ONLY");
//
// Boundary second half frame (should not overlapp with sensitive surface, nor frames)
//          Effective inner radius due to circle effect
     rMin =  TMath::Sqrt(
        (iChamber->frMin+2*dframep)*(iChamber->frMin+2*dframep) - dframep*dframep );
     bpar[0] = (iChamber->frMax - 2*dframep - rMin ) /2;
     bpar[2] = (5.325- (0.055 + 0.325)) / 2;
     gMC->Gsvolu("C05H", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C06H", "BOX", idtmed[1103], bpar, 3);
//place 2 boudaries
     gMC->Gspos("C05H",1,"C05Q", rMin+bpar[0],bpar[1],  0.055+0.325+bpar[2] , idrotm[1100],"ONLY");
     gMC->Gspos("C05H",2,"C05Q", rMin+bpar[0],bpar[1],-(0.055+0.325+bpar[2]), idrotm[1100],"ONLY");
     gMC->Gspos("C05H",3,"C05Q", bpar[1],rMin+bpar[0],  0.055+0.325+bpar[2] , idrotm[1101],"ONLY");
     gMC->Gspos("C05H",4,"C05Q", bpar[1],rMin+bpar[0],-(0.055+0.325+bpar[2]), idrotm[1101],"ONLY");
     gMC->Gspos("C06H",1,"C06Q", rMin+bpar[0],bpar[1],  0.055+0.325+bpar[2] , idrotm[1100],"ONLY");
     gMC->Gspos("C06H",2,"C06Q", rMin+bpar[0],bpar[1],-(0.055+0.325+bpar[2]), idrotm[1100],"ONLY");
     gMC->Gspos("C06H",3,"C06Q", bpar[1],rMin+bpar[0],  0.055+0.325+bpar[2] , idrotm[1101],"ONLY");
     gMC->Gspos("C06H",4,"C06Q", bpar[1],rMin+bpar[0],-(0.055+0.325+bpar[2]), idrotm[1101],"ONLY");
//
//   Chamber Material represented by Alu sheet
     //     tspar[2] = (iChamber->fdAlu)+(iChamber->fdGas);
     tspar[0]= iChamber->frMin+dframep*2;
     tspar[1]= iChamber->frMax-dframep*2;
     tspar[2] = 0.055 + 0.325;
     gMC->Gsvolu("C05A", "TUBS", idtmed[1103], tspar, 5);
     gMC->Gsvolu("C06A", "TUBS", idtmed[1103], tspar, 5);
     gMC->Gspos("C05A", 1, "C05Q", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C06A", 1, "C06Q", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->fdGas;
     tspar[2] = 0.325;
     gMC->Gsvolu("C05G", "TUBS", idtmed[1105], tspar, 5);
     gMC->Gsvolu("C06G", "TUBS", idtmed[1105], tspar, 5);
     gMC->Gspos("C05G", 1, "C05A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C06G", 1, "C06A", 0., 0., 0.,  0, "ONLY");
//
//   Overwrite sensitive volume with ALU
// Overwrite Gaz volume
     bpar[2] = 0.325;
     gMC->Gsvolu("C05Z", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C06Z", "BOX", idtmed[1103], bpar, 3);
     gMC->Gspos("C05Z",1,"C05G", rMin+bpar[0] ,bpar[1],0.,  idrotm[1100],"ONLY");
     gMC->Gspos("C05Z",2,"C05G", bpar[1], rMin+bpar[0] ,0., idrotm[1101],"ONLY");
     gMC->Gspos("C06Z",1,"C06G", rMin+bpar[0] ,bpar[1],0.,  idrotm[1100],"ONLY");
     gMC->Gspos("C06Z",2,"C06G", bpar[1], rMin+bpar[0] ,0., idrotm[1101],"ONLY");

//********************************************************************
//                            Station 4                             **
//********************************************************************
     iChamber=(AliMUONchamber*) (*fChambers)[6];
     zpos1=iChamber->ZPosition()-dstation/2; 
     zpos2=zpos1+dstation; 
     zfpos=-(iChamber->fdGas+dframez)/2;
     
//
//   Mother volume
     tpar[0] = iChamber->frMin-dframep; 
     tpar[1] = (iChamber->frMax+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/2;

     gMC->Gsvolu("C07M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gsvolu("C08M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gspos("C07M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C08M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
// Aluminium frames
// Outer frames
     pgpar[0] = 360/12/2;
     pgpar[1] = 360.;
     pgpar[2] = 12.;
     pgpar[3] =   2;
     pgpar[4] = -dframez/2;
     pgpar[5] = iChamber->frMax;
     pgpar[6] = pgpar[5]+dframep;
     pgpar[7] = +dframez/2;
     pgpar[8] = pgpar[5];
     pgpar[9] = pgpar[6];
     gMC->Gsvolu("C07O", "PGON", idtmed[1103], pgpar, 10);
     gMC->Gsvolu("C08O", "PGON", idtmed[1103], pgpar, 10);
     gMC->Gspos("C07O",1,"C07M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C07O",2,"C07M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C08O",1,"C08M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C08O",2,"C08M", 0.,0.,+zfpos,  0,"ONLY");
//
// Inner frame
     tpar[0]= iChamber->frMin-dframep;
     tpar[1]= iChamber->frMin;
     tpar[2]= dframez/2;
     gMC->Gsvolu("C07I", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gsvolu("C08I", "TUBE", idtmed[1103], tpar, 3);

     gMC->Gspos("C07I",1,"C07M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C07I",2,"C07M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C08I",1,"C08M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C08I",2,"C08M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     bpar[0] = (iChamber->frMax - iChamber->frMin)/2;
     bpar[1] = dframep/2;
     bpar[2] = dframez/2;
     gMC->Gsvolu("C07B", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C08B", "BOX", idtmed[1103], bpar, 3);

     gMC->Gspos("C07B",1,"C07M", +iChamber->frMin+bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C07B",2,"C07M", -iChamber->frMin-bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C07B",3,"C07M", 0, +iChamber->frMin+bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C07B",4,"C07M", 0, -iChamber->frMin-bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C07B",5,"C07M", +iChamber->frMin+bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C07B",6,"C07M", -iChamber->frMin-bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C07B",7,"C07M", 0, +iChamber->frMin+bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C07B",8,"C07M", 0, -iChamber->frMin-bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");

     gMC->Gspos("C08B",1,"C08M", +iChamber->frMin+bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C08B",2,"C08M", -iChamber->frMin-bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C08B",3,"C08M", 0, +iChamber->frMin+bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C08B",4,"C08M", 0, -iChamber->frMin-bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C08B",5,"C08M", +iChamber->frMin+bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C08B",6,"C08M", -iChamber->frMin-bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C08B",7,"C08M", 0, +iChamber->frMin+bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C08B",8,"C08M", 0, -iChamber->frMin-bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");

//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->frMin+dframep*2;
     tpar[1]= iChamber->frMax-dframep*2;
     tpar[2] = (iChamber->fdGas+iChamber->fdAlu)/2;
     gMC->Gsvolu("C07A", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gsvolu("C08A", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gspos("C07A", 1, "C07M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C08A", 1, "C08M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->fdGas;
     tpar[2] = iChamber->fdGas/2;
     gMC->Gsvolu("C07G", "TUBE", idtmed[1105], tpar, 3);
     gMC->Gsvolu("C08G", "TUBE", idtmed[1105], tpar, 3);
     gMC->Gspos("C07G", 1, "C07A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C08G", 1, "C08A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     dr = (iChamber->frMax - iChamber->frMin);
     bpar[0] = TMath::Sqrt(dr*dr-dframep*dframep/4)/2;
     bpar[1] = dframep/2;
     bpar[2] = iChamber->fdGas/2;
     gMC->Gsvolu("C07F", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C08F", "BOX", idtmed[1103], bpar, 3);

     gMC->Gspos("C07F",1,"C07G", +iChamber->frMin+bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C07F",2,"C07G", -iChamber->frMin-bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C07F",3,"C07G", 0, +iChamber->frMin+bpar[0] , 0, 
idrotm[1101],"ONLY");
     gMC->Gspos("C07F",4,"C07G", 0, -iChamber->frMin-bpar[0] , 0, 
idrotm[1101],"ONLY");

     gMC->Gspos("C08F",1,"C08G", +iChamber->frMin+bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C08F",2,"C08G", -iChamber->frMin-bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C08F",3,"C08G", 0, +iChamber->frMin+bpar[0] , 0, 
idrotm[1101],"ONLY");
     gMC->Gspos("C08F",4,"C08G", 0, -iChamber->frMin-bpar[0] , 0, 
idrotm[1101],"ONLY");

//********************************************************************
//                            Station 5                             **
//********************************************************************
     iChamber=(AliMUONchamber*) (*fChambers)[8];
     zpos1=iChamber->ZPosition()-dstation/2; 
     zpos2=zpos1+dstation; 
     zfpos=-(iChamber->fdGas+dframez)/2;
     
//
//   Mother volume
     tpar[0] = iChamber->frMin-dframep; 
     tpar[1] = (iChamber->frMax+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/2;

     gMC->Gsvolu("C09M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gsvolu("C10M", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gspos("C09M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C10M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
// Aluminium frames
// Outer frames
     pgpar[0] = 360/12/2;
     pgpar[1] = 360.;
     pgpar[2] = 12.;
     pgpar[3] =   2;
     pgpar[4] = -dframez/2;
     pgpar[5] = iChamber->frMax;
     pgpar[6] = pgpar[5]+dframep;
     pgpar[7] = +dframez/2;
     pgpar[8] = pgpar[5];
     pgpar[9] = pgpar[6];
     gMC->Gsvolu("C09O", "PGON", idtmed[1103], pgpar, 10);
     gMC->Gsvolu("C10O", "PGON", idtmed[1103], pgpar, 10);
     gMC->Gspos("C09O",1,"C09M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C09O",2,"C09M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C10O",1,"C10M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C10O",2,"C10M", 0.,0.,+zfpos,  0,"ONLY");
//
// Inner frame
     tpar[0]= iChamber->frMin-dframep;
     tpar[1]= iChamber->frMin;
     tpar[2]= dframez/2;
     gMC->Gsvolu("C09I", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gsvolu("C10I", "TUBE", idtmed[1103], tpar, 3);

     gMC->Gspos("C09I",1,"C09M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C09I",2,"C09M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C10I",1,"C10M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C10I",2,"C10M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     bpar[0] = (iChamber->frMax - iChamber->frMin)/2;
     bpar[1] = dframep/2;
     bpar[2] = dframez/2;
     gMC->Gsvolu("C09B", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C10B", "BOX", idtmed[1103], bpar, 3);

     gMC->Gspos("C09B",1,"C09M", +iChamber->frMin+bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C09B",2,"C09M", -iChamber->frMin-bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C09B",3,"C09M", 0, +iChamber->frMin+bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C09B",4,"C09M", 0, -iChamber->frMin-bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C09B",5,"C09M", +iChamber->frMin+bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C09B",6,"C09M", -iChamber->frMin-bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C09B",7,"C09M", 0, +iChamber->frMin+bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C09B",8,"C09M", 0, -iChamber->frMin-bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");

     gMC->Gspos("C10B",1,"C10M", +iChamber->frMin+bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C10B",2,"C10M", -iChamber->frMin-bpar[0] , 0,-zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C10B",3,"C10M", 0, +iChamber->frMin+bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C10B",4,"C10M", 0, -iChamber->frMin-bpar[0] ,-zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C10B",5,"C10M", +iChamber->frMin+bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C10B",6,"C10M", -iChamber->frMin-bpar[0] , 0,+zfpos, 
idrotm[1100],"ONLY");
     gMC->Gspos("C10B",7,"C10M", 0, +iChamber->frMin+bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");
     gMC->Gspos("C10B",8,"C10M", 0, -iChamber->frMin-bpar[0] ,+zfpos, 
idrotm[1101],"ONLY");

//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->frMin+dframep*2;
     tpar[1]= iChamber->frMax-dframep*2;
     tpar[2] = (iChamber->fdGas+iChamber->fdAlu)/2;
     gMC->Gsvolu("C09A", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gsvolu("C10A", "TUBE", idtmed[1103], tpar, 3);
     gMC->Gspos("C09A", 1, "C09M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C10A", 1, "C10M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->fdGas;
     tpar[2] = iChamber->fdGas/2;
     gMC->Gsvolu("C09G", "TUBE", idtmed[1105], tpar, 3);
     gMC->Gsvolu("C10G", "TUBE", idtmed[1105], tpar, 3);
     gMC->Gspos("C09G", 1, "C09A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C10G", 1, "C10A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     dr = (iChamber->frMax - iChamber->frMin);
     bpar[0] = TMath::Sqrt(dr*dr-dframep*dframep/4)/2;
     bpar[1] = dframep/2;
     bpar[2] = iChamber->fdGas/2;
     gMC->Gsvolu("C09F", "BOX", idtmed[1103], bpar, 3);
     gMC->Gsvolu("C10F", "BOX", idtmed[1103], bpar, 3);

     gMC->Gspos("C09F",1,"C09G", +iChamber->frMin+bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C09F",2,"C09G", -iChamber->frMin-bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C09F",3,"C09G", 0, +iChamber->frMin+bpar[0] , 0, 
idrotm[1101],"ONLY");
     gMC->Gspos("C09F",4,"C09G", 0, -iChamber->frMin-bpar[0] , 0, 
idrotm[1101],"ONLY");

     gMC->Gspos("C10F",1,"C10G", +iChamber->frMin+bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C10F",2,"C10G", -iChamber->frMin-bpar[0] , 0, 0, 
idrotm[1100],"ONLY");
     gMC->Gspos("C10F",3,"C10G", 0, +iChamber->frMin+bpar[0] , 0, 
idrotm[1101],"ONLY");
     gMC->Gspos("C10F",4,"C10G", 0, -iChamber->frMin-bpar[0] , 0, 
idrotm[1101],"ONLY");

///////////////////////////////////////
// GEOMETRY FOR THE TRIGGER CHAMBERS //
///////////////////////////////////////
     
//  Distance between planes inside each trigger station
    const Float_t DTPLANES = 15.;

//  Parameters of the Trigger Chambers
    //Station 1
    		
    const Float_t X_MC1_MIN=38.;       
    const Float_t X_MC1_MED=51.;                                
    const Float_t X_MC1_MAX=272.;                               
    const Float_t Y_MC1_MIN=34.;                              
    const Float_t Y_MC1_MAX=51.;                              
    const Float_t R_MIN1=48.;
    const Float_t R_MAX1=64.;

// Station 1
     iChamber=(AliMUONchamber*) (*fChambers)[10];
     zpos1=iChamber->ZPosition();
     zpos2=zpos1+DTPLANES;

// Mother volume definition     
     tpar[0] = iChamber->frMin; 
     tpar[1] = iChamber->frMax;
     tpar[2] = 0.4;    
     gMC->Gsvolu("CM11", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gsvolu("CM12", "TUBE", idtmed[1100], tpar, 3);
     
// Definition of the flange between the beam shielding and the RPC 
     tpar[0]= R_MIN1;
     tpar[1]= R_MAX1;
     tpar[2]= 0.4;
   
     gMC->Gsvolu("CF1A", "TUBE", idtmed[1103], tpar, 3);     //Al
     gMC->Gspos("CF1A", 1, "CM11", 0., 0., 0., 0, "MANY");
     gMC->Gspos("CF1A", 2, "CM12", 0., 0., 0., 0, "MANY");

// Definition of prototype for chambers in the first plane     
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC1A", "BOX ", idtmed[1103], tpar, 0);     //Al      
     gMC->Gsvolu("CB1A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG1A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t X_MC1A=X_MC1_MED+(X_MC1_MAX-X_MC1_MED)/2.;
     const Float_t Y_MC1A=0.;
     const Float_t Z_MC1A=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG1A", 1, "CB1A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB1A", 1, "CC1A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.4;
     tpar[0] = (X_MC1_MAX-X_MC1_MED)/2.;
     tpar[1] = Y_MC1_MIN;
     gMC->Gsposp("CC1A", 1, "CM11",X_MC1A,Y_MC1A,Z_MC1A, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 2, "CM11",-X_MC1A,Y_MC1A,Z_MC1A, 0, "ONLY", tpar, 3);
     
//  chamber type B    
     tpar[0] = (X_MC1_MAX-X_MC1_MIN)/2.;
     tpar[1] = (Y_MC1_MAX-Y_MC1_MIN)/2.;
     
     const Float_t X_MC1B=X_MC1_MIN+tpar[0];
     const Float_t Y_MC1B=Y_MC1_MIN+tpar[1];
     const Float_t Z_MC1B=0.;

     gMC->Gsposp("CC1A", 3, "CM11",X_MC1B,Y_MC1B,Z_MC1B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 4, "CM11",-X_MC1B,Y_MC1B,Z_MC1B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 5, "CM11",X_MC1B,-Y_MC1B,Z_MC1B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 6, "CM11",-X_MC1B,-Y_MC1B,Z_MC1B, 0, "ONLY", tpar, 3);
     
//  chamber type C        
     tpar[0] = X_MC1_MAX/2;
     tpar[1] = Y_MC1_MAX/2;
     
     const Float_t X_MC1C=tpar[0];
     const Float_t Y_MC1C=Y_MC1_MAX+tpar[1];
     const Float_t Z_MC1C=0.;
     
     gMC->Gsposp("CC1A", 7, "CM11",X_MC1C,Y_MC1C,Z_MC1C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 8, "CM11",-X_MC1C,Y_MC1C,Z_MC1C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 9, "CM11",X_MC1C,-Y_MC1C,Z_MC1C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 10, "CM11",-X_MC1C,-Y_MC1C,Z_MC1C, 0, "ONLY", tpar, 3);
     
//  chamber type D        
     tpar[0] = X_MC1_MAX/2.;
     tpar[1] = Y_MC1_MIN;
     
     const Float_t X_MC1D=tpar[0];
     const Float_t Z_MC1D=0.;
     
     Float_t Y_MC1D=4.*Y_MC1_MIN;
     gMC->Gsposp("CC1A", 11, "CM11",X_MC1D,Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 12, "CM11",X_MC1D,-Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 13, "CM11",-X_MC1D,Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 14, "CM11",-X_MC1D,-Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);

     Y_MC1D=6.*Y_MC1_MIN;
     gMC->Gsposp("CC1A", 15, "CM11",X_MC1D,Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 16, "CM11",X_MC1D,-Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 17, "CM11",-X_MC1D,Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 18, "CM11",-X_MC1D,-Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);

     Y_MC1D=8.*Y_MC1_MIN;
     gMC->Gsposp("CC1A", 19, "CM11",X_MC1D,Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 20, "CM11",X_MC1D,-Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 21, "CM11",-X_MC1D,Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 22, "CM11",-X_MC1D,-Y_MC1D,Z_MC1D, 0, "ONLY", tpar, 3);

// Positioning first plane in ALICE     
     gMC->Gspos("CM11", 1, "ALIC", 0., 0., zpos1, 0, "ONLY");

// End of geometry definition for the first plane

// Station 1 - plan 2  -  same RPCs as plan 1 ==> small non covered area
// Y position moved (ratio zpos2/zpos1)
     const Float_t Z_1S2=zpos2/zpos1;
     
// Definition of prototype for chambers in the second plane     
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC2A", "BOX ", idtmed[1103], tpar, 0);     //Al      
     gMC->Gsvolu("CB2A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG2A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t X_MC2A=X_MC1A;
     const Float_t Y_MC2A=0.;
     const Float_t Z_MC2A=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG2A", 1, "CB2A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB2A", 1, "CC2A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.4;
     tpar[0] = (X_MC1_MAX-X_MC1_MED)/2.;
     tpar[1] = Y_MC1_MIN;
     gMC->Gsposp("CC2A", 1, "CM12",X_MC2A,Y_MC2A,Z_MC2A, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 2, "CM12",-X_MC2A,Y_MC2A,Z_MC2A, 0, "ONLY", tpar, 3);
     
//  chamber type B    
     tpar[0] = (X_MC1_MAX-X_MC1_MIN)/2.;
     tpar[1] = (Y_MC1_MAX-Y_MC1_MIN)/2.;
     
     const Float_t X_MC2B=X_MC1B;
     const Float_t Y_MC2B=2.*Y_MC1_MIN*Z_1S2-Y_MC1_MIN*1.5+Y_MC1_MAX*0.5;
     const Float_t Z_MC2B=0.;

     gMC->Gsposp("CC2A", 3, "CM12",X_MC2B,Y_MC2B,Z_MC2B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 4, "CM12",-X_MC2B,Y_MC2B,Z_MC2B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 5, "CM12",X_MC2B,-Y_MC2B,Z_MC2B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 6, "CM12",-X_MC2B,-Y_MC2B,Z_MC2B, 0, "ONLY", tpar, 3);
     
//  chamber type C        
     tpar[0] = X_MC1_MAX/2;
     tpar[1] = Y_MC1_MAX/2;
     
     const Float_t X_MC2C=X_MC1C;
     const Float_t Y_MC2C=2.*Y_MC1_MIN*Z_1S2-Y_MC1_MIN*2.+Y_MC1_MAX*1.5;
     const Float_t Z_MC2C=0.;
     
     gMC->Gsposp("CC2A", 7, "CM12",X_MC2C,Y_MC2C,Z_MC2C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 8, "CM12",-X_MC2C,Y_MC2C,Z_MC2C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 9, "CM12",X_MC2C,-Y_MC2C,Z_MC2C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 10, "CM12",-X_MC2C,-Y_MC2C,Z_MC2C, 0, "ONLY", tpar, 3);
     
//  chamber type D        
     tpar[0] = X_MC1_MAX/2.;
     tpar[1] = Y_MC1_MIN;
     
     const Float_t X_MC2D=X_MC1D;
     const Float_t Z_MC2D=0.;
     
     Float_t Y_MC2D=4.*Y_MC1_MIN*Z_1S2;
     gMC->Gsposp("CC2A", 11, "CM12",X_MC2D,Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 12, "CM12",X_MC2D,-Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 13, "CM12",-X_MC2D,Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 14, "CM12",-X_MC2D,-Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);

     Y_MC2D=6.*Y_MC1_MIN*Z_1S2;
     gMC->Gsposp("CC2A", 15, "CM12",X_MC2D,Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 16, "CM12",X_MC2D,-Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 17, "CM12",-X_MC2D,Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 18, "CM12",-X_MC2D,-Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);

     Y_MC2D=8.*Y_MC1_MIN*Z_1S2;
     gMC->Gsposp("CC2A", 19, "CM12",X_MC2D,Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 20, "CM12",X_MC2D,-Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 21, "CM12",-X_MC2D,Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 22, "CM12",-X_MC2D,-Y_MC2D,Z_MC2D, 0, "ONLY", tpar, 3);
     
     gMC->Gspos("CM12", 1, "ALIC", 0., 0., zpos2, 0, "ONLY");

// Station 2
     iChamber=(AliMUONchamber*) (*fChambers)[12];
     Float_t zpos3=iChamber->ZPosition();
     Float_t zpos4=zpos3+DTPLANES;

//  Parameters of the Trigger Chambers
                                              //Station 2
     const Float_t X_MC3_MIN=X_MC1_MIN*zpos3/zpos1;          
     const Float_t X_MC3_MED=X_MC1_MED*zpos3/zpos1;          
     const Float_t X_MC3_MAX=X_MC1_MAX*zpos3/zpos1;           
     const Float_t Y_MC3_MIN=Y_MC1_MIN*zpos3/zpos1;            
     const Float_t Y_MC3_MAX=Y_MC1_MAX*zpos3/zpos1;             
     const Float_t R_MIN3=R_MIN1*zpos3/zpos1;
     const Float_t R_MAX3=R_MAX1*zpos3/zpos1;

// Mother volume definition     
     tpar[0] = iChamber->frMin; 
     tpar[1] = iChamber->frMax;
     tpar[2] = 0.4;    
     gMC->Gsvolu("CM21", "TUBE", idtmed[1100], tpar, 3);
     gMC->Gsvolu("CM22", "TUBE", idtmed[1100], tpar, 3);
     
// Definition of the flange between the beam shielding and the RPC 
     tpar[0]= R_MIN3;
     tpar[1]= R_MAX3;
     tpar[2]= 0.4;
   
     gMC->Gsvolu("CF2A", "TUBE", idtmed[1103], tpar, 3);     //Al
     gMC->Gspos("CF2A", 1, "CM21", 0., 0., 0., 0, "MANY");
     gMC->Gspos("CF2A", 2, "CM22", 0., 0., 0., 0, "MANY");

// Definition of prototype for chambers in the third plane     
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC3A", "BOX ", idtmed[1103], tpar, 0);     //Al      
     gMC->Gsvolu("CB3A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG3A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t X_MC3A=X_MC3_MED+(X_MC3_MAX-X_MC3_MED)/2.;
     const Float_t Y_MC3A=0.;
     const Float_t Z_MC3A=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG3A", 1, "CB3A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB3A", 1, "CC3A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[0] = (X_MC3_MAX-X_MC3_MED)/2.;
     tpar[1] = Y_MC3_MIN;
     tpar[2] = 0.4;
     gMC->Gsposp("CC3A", 1, "CM21",X_MC3A,Y_MC3A,Z_MC3A, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 2, "CM21",-X_MC3A,Y_MC3A,Z_MC3A, 0, "ONLY", tpar, 3);
     
//  chamber type B    
     tpar[0] = (X_MC3_MAX-X_MC3_MIN)/2.;
     tpar[1] = (Y_MC3_MAX-Y_MC3_MIN)/2.;
     
     const Float_t X_MC3B=X_MC3_MIN+tpar[0];
     const Float_t Y_MC3B=Y_MC3_MIN+tpar[1];
     const Float_t Z_MC3B=0.;

     gMC->Gsposp("CC3A", 3, "CM21",X_MC3B,Y_MC3B,Z_MC3B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 4, "CM21",-X_MC3B,Y_MC3B,Z_MC3B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 5, "CM21",X_MC3B,-Y_MC3B,Z_MC3B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 6, "CM21",-X_MC3B,-Y_MC3B,Z_MC3B, 0, "ONLY", tpar, 3);
     
//  chamber type C        
     tpar[0] = X_MC3_MAX/2.;
     tpar[1] = Y_MC3_MAX/2.;
     
     const Float_t X_MC3C=tpar[0];
     const Float_t Y_MC3C=Y_MC3_MAX+tpar[1];
     const Float_t Z_MC3C=0.;
     
     gMC->Gsposp("CC3A", 7, "CM21",X_MC3C,Y_MC3C,Z_MC3C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 8, "CM21",-X_MC3C,Y_MC3C,Z_MC3C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 9, "CM21",X_MC3C,-Y_MC3C,Z_MC3C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 10, "CM21",-X_MC3C,-Y_MC3C,Z_MC3C, 0, "ONLY", tpar, 3);
     
//  chamber type D        
     tpar[0] = X_MC3_MAX/2.;
     tpar[1] = Y_MC3_MIN;
     
     const Float_t X_MC3D=tpar[0];
     const Float_t Z_MC3D=0.;
     
     Float_t Y_MC3D=4.*Y_MC3_MIN;
     gMC->Gsposp("CC3A", 11, "CM21",X_MC3D,Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 12, "CM21",X_MC3D,-Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 13, "CM21",-X_MC3D,Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 14, "CM21",-X_MC3D,-Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);

     Y_MC3D=6.*Y_MC3_MIN;
     gMC->Gsposp("CC3A", 15, "CM21",X_MC3D,Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 16, "CM21",X_MC3D,-Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 17, "CM21",-X_MC3D,Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 18, "CM21",-X_MC3D,-Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);

     Y_MC3D=8.*Y_MC3_MIN;
     gMC->Gsposp("CC3A", 19, "CM21",X_MC3D,Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 20, "CM21",X_MC3D,-Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 21, "CM21",-X_MC3D,Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 22, "CM21",-X_MC3D,-Y_MC3D,Z_MC3D, 0, "ONLY", tpar, 3);

// Positioning third plane in ALICE     
     gMC->Gspos("CM21", 1, "ALIC", 0., 0., zpos3, 0, "ONLY");

// End of geometry definition for the third plane

// Station 2 - plan 4  -  same RPCs as plan 3 ==> small non covered area
// Y position moved (ratio zpos4/zpos3)
     const Float_t Z_3S4=zpos4/zpos3;
     
// Definition of prototype for chambers in the fourth plane     
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC4A", "BOX ", idtmed[1103], tpar, 0);     //Al      
     gMC->Gsvolu("CB4A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG4A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t X_MC4A=X_MC3A;
     const Float_t Y_MC4A=0.;
     const Float_t Z_MC4A=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG4A", 1, "CB4A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB4A", 1, "CC4A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.4;
     tpar[0] = (X_MC3_MAX-X_MC3_MED)/2.;
     tpar[1] = Y_MC3_MIN;
     gMC->Gsposp("CC4A", 1, "CM22",X_MC4A,Y_MC4A,Z_MC4A, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 2, "CM22",-X_MC4A,Y_MC4A,Z_MC4A, 0, "ONLY", tpar, 3);
     
//  chamber type B    
     tpar[0] = (X_MC3_MAX-X_MC3_MIN)/2.;
     tpar[1] = (Y_MC3_MAX-Y_MC3_MIN)/2.;
     
     const Float_t X_MC4B=X_MC3B;
     const Float_t Y_MC4B=2.*Y_MC3_MIN*Z_3S4-Y_MC3_MIN*1.5+Y_MC3_MAX*0.5;
     const Float_t Z_MC4B=0.;

     gMC->Gsposp("CC4A", 3, "CM22",X_MC4B,Y_MC4B,Z_MC4B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 4, "CM22",-X_MC4B,Y_MC4B,Z_MC4B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 5, "CM22",X_MC4B,-Y_MC4B,Z_MC4B, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 6, "CM22",-X_MC4B,-Y_MC4B,Z_MC4B, 0, "ONLY", tpar, 3);
     
//  chamber type C        
     tpar[0] = X_MC3_MAX/2;
     tpar[1] = Y_MC3_MAX/2;
     
     const Float_t X_MC4C=X_MC3C;
     const Float_t Y_MC4C=2.*Y_MC3_MIN*Z_3S4-Y_MC3_MIN*2.+Y_MC3_MAX*1.5;
     const Float_t Z_MC4C=0.;
     
     gMC->Gsposp("CC4A", 7, "CM22",X_MC4C,Y_MC4C,Z_MC4C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 8, "CM22",-X_MC4C,Y_MC4C,Z_MC4C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 9, "CM22",X_MC4C,-Y_MC4C,Z_MC4C, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 10, "CM22",-X_MC4C,-Y_MC4C,Z_MC4C, 0, "ONLY", tpar, 3);
     
//  chamber type D        
     tpar[0] = X_MC3_MAX/2.;
     tpar[1] = Y_MC3_MIN;
     
     const Float_t X_MC4D=X_MC3D;
     const Float_t Z_MC4D=0.;
     
     Float_t Y_MC4D=4.*Y_MC3_MIN*Z_3S4;
     gMC->Gsposp("CC4A", 11, "CM22",X_MC4D,Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 12, "CM22",X_MC4D,-Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 13, "CM22",-X_MC4D,Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 14, "CM22",-X_MC4D,-Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);

     Y_MC4D=6.*Y_MC3_MIN*Z_3S4;
     gMC->Gsposp("CC4A", 15, "CM22",X_MC4D,Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 16, "CM22",X_MC4D,-Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 17, "CM22",-X_MC4D,Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 18, "CM22",-X_MC4D,-Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);

     Y_MC4D=8.*Y_MC3_MIN*Z_3S4;
     gMC->Gsposp("CC4A", 19, "CM22",X_MC4D,Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 20, "CM22",X_MC4D,-Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 21, "CM22",-X_MC4D,Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 22, "CM22",-X_MC4D,-Y_MC4D,Z_MC4D, 0, "ONLY", tpar, 3);
     
     gMC->Gspos("CM22", 1, "ALIC", 0., 0., zpos4, 0, "ONLY");
}

 
//___________________________________________
void AliMUONv0::CreateMaterials()
{
  // *** DEFINITION OF AVAILABLE MUON MATERIALS *** 
  //
  //     Ar-CO2 gas 
    Float_t ag1[3]   = { 39.95,12.01,16. };
    Float_t zg1[3]   = { 18.,6.,8. };
    Float_t wg1[3]   = { .8,.0667,.13333 };
    Float_t dg1      = .001821;
    //
    //     Ar-buthane-freon gas -- trigger chambers 
    Float_t atr1[4]  = { 39.95,12.01,1.01,19. };
    Float_t ztr1[4]  = { 18.,6.,1.,9. };
    Float_t wtr1[4]  = { .56,.1262857,.2857143,.028 };
    Float_t dtr1     = .002599;
    //
    //     Ar-CO2 gas 
    Float_t agas[3]  = { 39.95,12.01,16. };
    Float_t zgas[3]  = { 18.,6.,8. };
    Float_t wgas[3]  = { .74,.086684,.173316 };
    Float_t dgas     = .0018327;
    //
    //     Ar-Isobutane gas (80%+20%) -- tracking 
    Float_t ag[3]    = { 39.95,12.01,1.01 };
    Float_t zg[3]    = { 18.,6.,1. };
    Float_t wg[3]    = { .8,.057,.143 };
    Float_t dg       = .0019596;
    //
    //     Ar-Isobutane-Forane-SF6 gas (49%+7%+40%+4%) -- trigger 
    Float_t atrig[5] = { 39.95,12.01,1.01,19.,32.066 };
    Float_t ztrig[5] = { 18.,6.,1.,9.,16. };
    Float_t wtrig[5] = { .49,1.08,1.5,1.84,0.04 };
    Float_t dtrig    = .0031463;
    //
    //     bakelite 

    Float_t abak[3] = {12.01 , 1.01 , 16.};
    Float_t zbak[3] = {6.     , 1.   , 8.};
    Float_t wbak[3] = {6.     , 6.   , 1.}; 
    Float_t dbak = 1.4;

    Float_t epsil, stmin, deemax, tmaxfd, stemax;

    Int_t ISXFLD   = gAlice->Field()->Integ();
    Float_t SXMGMX = gAlice->Field()->Max();
    //
    // --- Define the various materials for GEANT --- 
    AliMaterial(9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
    AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
    AliMixture(19, "Bakelite$", abak, zbak, dbak, -3, wbak);
    AliMixture(20, "ArC4H10 GAS$", ag, zg, dg, 3, wg);
    AliMixture(21, "TRIG GAS$", atrig, ztrig, dtrig, -5, wtrig);
    AliMixture(22, "ArCO2 80%$", ag1, zg1, dg1, 3, wg1);
    AliMixture(23, "Ar-freon $", atr1, ztr1, dtr1, 4, wtr1);
    AliMixture(24, "ArCO2 GAS$", agas, zgas, dgas, 3, wgas);

    epsil  = .001; // Tracking precision, 
    stemax = -1.;  // Maximum displacement for multiple scat 
    tmaxfd = -20.; // Maximum angle due to field deflection 
    deemax = -.3;  // Maximum fractional energy loss, DLS 
    stmin  = -.8;
    //
    //    Air 
    AliMedium(1, "AIR_CH_US         ", 15, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    //
    //    Aluminum 

    AliMedium(4, "ALU_CH_US          ", 9, 0, ISXFLD, SXMGMX, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);
    //
    //    Ar-isoC4H10 gas 

    AliMedium(6, "AR_CH_US          ", 20, 1, ISXFLD, SXMGMX, tmaxfd, fMaxStepGas, 
	    fMaxDestepGas, epsil, stmin);
//
    //    Ar-Isobuthane-Forane-SF6 gas 

    AliMedium(7, "GAS_CH_TRIGGER    ", 21, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);

    AliMedium(8, "BAKE_CH_TRIGGER   ", 19, 0, ISXFLD, SXMGMX, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);

}

//___________________________________________

void AliMUONv0::Init()
{
   printf("\n\n\n Start Init for version 0 - CPC chamber type\n\n\n");

   // 
   // Initialize Tracking Chambers
   //
   for (Int_t i=0; i<NCH; i++) {
       ( (AliMUONchamber*) (*fChambers)[i])->Init();
   }
   //
   // Set the chamber (sensitive region) GEANT identifier
   ((AliMUONchamber*)(*fChambers)[0])->SetGid(gMC->VolId("C01G"));
   ((AliMUONchamber*)(*fChambers)[1])->SetGid(gMC->VolId("C02G"));
   ((AliMUONchamber*)(*fChambers)[2])->SetGid(gMC->VolId("C03G"));
   ((AliMUONchamber*)(*fChambers)[3])->SetGid(gMC->VolId("C04G"));
   ((AliMUONchamber*)(*fChambers)[4])->SetGid(gMC->VolId("C05G"));
   ((AliMUONchamber*)(*fChambers)[5])->SetGid(gMC->VolId("C06G"));
   ((AliMUONchamber*)(*fChambers)[6])->SetGid(gMC->VolId("C07G"));
   ((AliMUONchamber*)(*fChambers)[7])->SetGid(gMC->VolId("C08G"));
   ((AliMUONchamber*)(*fChambers)[8])->SetGid(gMC->VolId("C09G"));
   ((AliMUONchamber*)(*fChambers)[9])->SetGid(gMC->VolId("C10G"));
   ((AliMUONchamber*)(*fChambers)[10])->SetGid(gMC->VolId("CG1A"));
   ((AliMUONchamber*)(*fChambers)[11])->SetGid(gMC->VolId("CG2A"));
   ((AliMUONchamber*)(*fChambers)[12])->SetGid(gMC->VolId("CG3A"));
   ((AliMUONchamber*)(*fChambers)[13])->SetGid(gMC->VolId("CG4A"));

   printf("\n\n\n Finished Init for version 0 - CPC chamber type\n\n\n");
}

//___________________________________________
void AliMUONv0::StepManager()
{
  Int_t          copy, id;
  static Int_t   idvol;
  static Int_t   vol[2];
  Int_t          ipart;
  static Float_t hits[10];
  Float_t        pos[3];
  Float_t        mom[4];
  Float_t        theta,phi;
  Float_t        destep, step;
  static Float_t eloss, xhit, yhit, tlength;
  const  Float_t big=1.e10;
  
  TClonesArray &lhits = *fHits;

  //
  // Set maximum step size for gas
  // numed=gMC->GetMedium();
  //
  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  idvol=-1;
  id=gMC->CurrentVol(0,copy);
  
    for (Int_t i=1; i<=NCH; i++) {
      if(id==((AliMUONchamber*)(*fChambers)[i-1])->GetGid()){ 
	  vol[0]=i; 
	  idvol=i-1;
      }
    }
    if (idvol == -1) return;
  //
  // Get current particle id (ipart), track position (pos)  and momentum (mom) 
  gMC->TrackPosition(pos);
  gMC->TrackMomentum(mom);

  ipart  = gMC->TrackPid();
  //
  // momentum loss and steplength in last step
  destep = gMC->Edep();
  step   = gMC->TrackStep();
  
  //
  // record hits when track enters ...
  if( gMC->TrackEntering()) {
      gMC->SetMaxStep(fMaxStepGas);
      Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
      Double_t rt = TMath::Sqrt(tc);
      theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
      phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
      hits[0] = Float_t(ipart);         // Geant3 particle type
      hits[1] = pos[0];                 // X-position for hit
      hits[2] = pos[1];                 // Y-position for hit
      hits[3] = pos[2];                 // Z-position for hit
      hits[4] = theta;                  // theta angle of incidence
      hits[5] = phi;                    // phi angle of incidence 
      hits[8] = (Float_t) fNclusters;   // first padhit
      hits[9] = -1;                     // last pad hit
      // phi angle of incidence
      tlength = 0;
      eloss   = 0;
      xhit    = pos[0];
      yhit    = pos[1];      
      // Only if not trigger chamber
      if(idvol<10) {
	  //
	  //  Initialize hit position (cursor) in the segmentation model 
	  ((AliMUONchamber*) (*fChambers)[idvol])
	      ->SigGenInit(pos[0], pos[1], pos[2]);
      } else {
	  //geant3->Gpcxyz();
	  //printf("In the Trigger Chamber #%d\n",idvol-9);
      }
  }
  
  // 
  // Calculate the charge induced on a pad (disintegration) in case 
  //
  // Mip left chamber ...
  if( gMC->TrackExiting() || gMC->TrackStop() || gMC->TrackDisappear()){
      gMC->SetMaxStep(big);
      eloss   += destep;
      tlength += step;
      
      // Only if not trigger chamber
      if(idvol<10) {
	  if (eloss > 0) MakePadHits(xhit,yhit,eloss,idvol);
      }
	  
      hits[6]=tlength;
      hits[7]=eloss;
      if (fNclusters > (Int_t)hits[8]) {
	  hits[8]= hits[8]+1;
	  hits[9]= (Float_t) fNclusters;
      }
    
      new(lhits[fNhits++]) 
	  AliMUONhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
      eloss = 0; 
      //
      // Check additional signal generation conditions 
      // defined by the segmentation
      // model (boundary crossing conditions) 
  } else if 
      (((AliMUONchamber*) (*fChambers)[idvol])
       ->SigGenCond(pos[0], pos[1], pos[2]))
  {
      ((AliMUONchamber*) (*fChambers)[idvol])
	  ->SigGenInit(pos[0], pos[1], pos[2]);
//      printf("\n-> MakePadHits, reason special %d",ipart);
      if (eloss > 0) MakePadHits(xhit,yhit,eloss,idvol);
      xhit     = pos[0];
      yhit     = pos[1]; 
      eloss    = destep;
      tlength += step ;
      //
      // nothing special  happened, add up energy loss
  } else {        
      eloss   += destep;
      tlength += step ;
  }
}

//___________________________________________
void AliMUON::MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol)
{
//
//  Calls the charge disintegration method of the current chamber and adds
//  the simulated cluster to the root treee 
//
    Int_t clhits[7];
    Float_t newclust[6][500];
    Int_t nnew;
    
    
//
//  Integrated pulse height on chamber

    
    clhits[0]=fNhits+1;
//
//
    ((AliMUONchamber*) (*fChambers)[idvol])->DisIntegration(eloss, xhit, yhit, nnew, newclust);
//    printf("\n Add new clusters %d %f", nnew, eloss*1.e9);
    Int_t ic=0;
    
//
//  Add new clusters
    for (Int_t i=0; i<nnew; i++) {
	if (Int_t(newclust[3][i]) > 0) {
	    ic++;
// Cathode plane
	    clhits[1] = Int_t(newclust[5][i]);
//  Cluster Charge
	    clhits[2] = Int_t(newclust[0][i]);
//  Pad: ix
	    clhits[3] = Int_t(newclust[1][i]);
//  Pad: iy 
	    clhits[4] = Int_t(newclust[2][i]);
//  Pad: charge
	    clhits[5] = Int_t(newclust[3][i]);
//  Pad: chamber sector
	    clhits[6] = Int_t(newclust[4][i]);
	    
	    AddCluster(clhits);
	}
    }
//    printf("\n %d new clusters added", ic);
}

ClassImp(AliMUONchamber)	
    AliMUONchamber::AliMUONchamber() 
{
    fSegmentation = new TObjArray(2);
    fResponse=0;
    fnsec=1;
}

void AliMUONchamber::Init()
{
    
    ((AliMUONsegmentation *) (*fSegmentation)[0])->Init(this);
    if (fnsec==2) {
	((AliMUONsegmentation *) (*fSegmentation)[1])->Init(this);
    }
    
}

void AliMUONchamber::DisIntegration(Float_t eloss, Float_t xhit, Float_t yhit,
				    Int_t& nnew,Float_t newclust[6][500]) 
{
//    
//  Generates pad hits (simulated cluster) 
//  using the segmentation and the response model 
    Float_t dx, dy;
    //
    // Width of the integration area
    //
    dx=fResponse->Nsigma()*fResponse->ChwX();
    dy=fResponse->Nsigma()*fResponse->ChwY();
    //
    // Get pulse height from energy loss
    Float_t qtot = fResponse->IntPH(eloss);
    //
    // Loop Over Pads
    
    Float_t qcheck=0, qp;
    nnew=0;
    for (Int_t i=1; i<=fnsec; i++) {
	qcheck=0;
	AliMUONsegmentation * segmentation=(AliMUONsegmentation *) (*fSegmentation)[i-1];
	for (segmentation->FirstPad(xhit, yhit, dx, dy); 
	     segmentation->MorePads(); 
	     segmentation->NextPad()) 
	{
	    qp=fResponse->IntXY(segmentation);
	    qp=TMath::Abs(qp);
	    
//
//
	    if (qp > 1.e-4) {
		qcheck+=qp;
	    //
	    // --- store signal information
		newclust[0][nnew]=qtot;
		newclust[1][nnew]=segmentation->Ix();
		newclust[2][nnew]=segmentation->Iy();
		newclust[3][nnew]=qp * qtot;
		newclust[4][nnew]=segmentation->ISector();
		newclust[5][nnew]=(Float_t) i;
//		printf("\n pad hit %d %d %f %f ",nnew,i,newclust[1][nnew],newclust[2][nnew]);
		nnew++;

		
	    }
	} // Pad loop
//	printf("\n check sum is %f %f %f %f %d",qcheck,qtot,xhit,yhit,fGid);
    } // Cathode plane loop
}


ClassImp(AliMUONsegmentation)
ClassImp(AliMUONresponse)	
//___________________________________________
ClassImp(AliMUONsegmentationV0)
    void AliMUONsegmentationV0::Init(AliMUONchamber* Chamber)
{
    fNpx=(Int_t) (Chamber->frMax/fDpx+1);
    fNpy=(Int_t) (Chamber->frMax/fDpy+1);
}


Float_t AliMUONsegmentationV0::GetAnod(Float_t xhit)
{
    Float_t wire= (xhit<0)? Int_t(xhit/fWireD)+0.5:Int_t(xhit/fWireD)-0.5;
    return fWireD*wire;
}

void AliMUONsegmentationV0::SetPADSIZ(Float_t p1, Float_t p2)
{
    fDpx=p1;
    fDpy=p2;
}
void AliMUONsegmentationV0::
    GetPadIxy(Float_t x, Float_t y, Int_t &ix, Int_t &iy)
{
//  returns pad coordinates (ix,iy) for given real coordinates (x,y)
//
    ix = (x>0)? Int_t(x/fDpx)+1 : Int_t(x/fDpx)-1;
    iy = (y>0)? Int_t(y/fDpy)+1 : Int_t(y/fDpy)-1;
    if (iy >  fNpy) iy= fNpy;
    if (iy < -fNpy) iy=-fNpy;
    if (ix >  fNpx) ix= fNpx;
    if (ix < -fNpx) ix=-fNpx;
}
void AliMUONsegmentationV0::
GetPadCxy(Int_t ix, Int_t iy, Float_t &x, Float_t &y)
{
//  returns real coordinates (x,y) for given pad coordinates (ix,iy)
//
    x = (ix>0) ? Float_t(ix*fDpx)-fDpx/2. : Float_t(ix*fDpx)+fDpx/2.;
    y = (iy>0) ? Float_t(iy*fDpy)-fDpy/2. : Float_t(iy*fDpy)+fDpy/2.;
}

void AliMUONsegmentationV0::
FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy)
{
    //
    // Find the wire position (center of charge distribution)
    Float_t x0a=GetAnod(xhit);
    //
    // and take fNsigma*sigma around this center
    Float_t x01=x0a  - dx;
    Float_t x02=x0a  + dx;
    Float_t y01=yhit - dy;
    Float_t y02=yhit + dy;
    //
    // find the pads over which the charge distributes
    GetPadIxy(x01,y01,fixmin,fiymin);
    GetPadIxy(x02,y02,fixmax,fiymax);    
//    printf("\n %f %f %d %d",x02,y02,fixmax,fiymax);
//    printf("\n FirstPad called %f %f ", fDpx, fDpy);
//    printf("\n Hit Position %f %f",xhit,yhit);
//    printf("\n Integration limits: %i %i %i %i",fixmin,fixmax,fiymin,fiymax);
//    printf("\n Integration limits: %f %f %f %f \n",x01,x02,y01,y02);
    // 
    // Set current pad to lower left corner
    fix=fixmin;
    fiy=fiymin;
    GetPadCxy(fix,fiy,fx,fy);
}

void AliMUONsegmentationV0::NextPad()
{
  // 
  // Step to next pad in integration region
    if (fix != fixmax) {
	fix++;
    } else if (fiy != fiymax) {
	fix=fixmin;
	fiy++;
    } else {
	printf("\n Error: Stepping outside integration region\n ");
    }
    GetPadCxy(fix,fiy,fx,fy);
}

Int_t AliMUONsegmentationV0::MorePads()
//
// Are there more pads in the integration region
{
    if (fix == fixmax && fiy == fiymax) {
	return 0;
    } else {
	return 1;
	
    }
}

void AliMUONsegmentationV0::SigGenInit(Float_t x,Float_t y,Float_t)
{
//
//  Initialises pad and wire position during stepping
    fxt =x;
    fyt =y;
    GetPadIxy(x,y,fixt,fiyt);
    fiwt=Int_t(x/fWireD)+1;
}

Int_t AliMUONsegmentationV0::SigGenCond(Float_t x,Float_t y,Float_t)
{
//
//  Signal will be generated if particle crosses pad boundary or
//  boundary between two wires. 
    Int_t ixt, iyt;
    GetPadIxy(x,y,ixt,iyt);
    Int_t iwt=Int_t(x/fWireD)+1;
    
    if ((ixt != fixt) || (iyt !=fiyt) || (iwt != fiwt)) {
	return 1;
    } else {
	return 0;
    }
}
void AliMUONsegmentationV0::
IntegrationLimits(Float_t& x1,Float_t& x2,Float_t& y1, Float_t& y2)
{
    x1=fxt-fx-fDpx/2.;
    x2=x1+fDpx;
    y1=fyt-fy-fDpy/2.;
    y2=y1+fDpy;    
}

void AliMUONsegmentationV0::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[7], Int_t Ylist[7])
{
    *Nlist=4;Xlist[0]=Xlist[1]=iX;Xlist[2]=iX-1;Xlist[3]=iX+1;
    Ylist[0]=iY-1;Ylist[1]=iY+1;Ylist[2]=Ylist[3]=iY;
}

void AliMUONsegmentationV0::
FitXY(AliMUONRecCluster* ,TClonesArray* )
    // Default : Centre of gravity method
{
    ;
}


//___________________________________________
ClassImp(AliMUONresponseV0)	
Float_t AliMUONresponseV0::IntPH(Float_t eloss)
{
  // Get number of electrons and return charge
     
  Int_t nel;
  nel= Int_t(eloss*1.e9/26.);
  Float_t charge=0;
  if (nel == 0) nel=1;
  for (Int_t i=1;i<=nel;i++) {
    charge -= fChslope*TMath::Log(gRandom->Rndm());    
  }
  return charge;
}
// -------------------------------------------
Float_t AliMUONresponseV0::IntXY(AliMUONsegmentation * segmentation)
{

    const Float_t invpitch = 1/fPitch;
//
//  Integration limits defined by segmentation model
//  
    Float_t xi1, xi2, yi1, yi2;
    segmentation->IntegrationLimits(xi1,xi2,yi1,yi2);
    xi1=xi1*invpitch;
    xi2=xi2*invpitch;
    yi1=yi1*invpitch;
    yi2=yi2*invpitch;
//
// The Mathieson function 
    Double_t ux1=fSqrtKx3*TMath::TanH(fKx2*xi1);
    Double_t ux2=fSqrtKx3*TMath::TanH(fKx2*xi2);

    Double_t uy1=fSqrtKy3*TMath::TanH(fKy2*yi1);
    Double_t uy2=fSqrtKy3*TMath::TanH(fKy2*yi2);

    
    return Float_t(4.*fKx4*(TMath::ATan(ux2)-TMath::ATan(ux1))*
		      fKy4*(TMath::ATan(uy2)-TMath::ATan(uy1)));
}

// -------------------------------------------
ClassImp(AliMUONgeometry)	
    void AliMUONgeometry::InitGeo(Float_t)
{
      fdGas= 0.5;
      fdAlu= 2.5/100*8.9;
}
















