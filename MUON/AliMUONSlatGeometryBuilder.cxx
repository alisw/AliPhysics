// $Id$
//
// Class AliMUONSlatGeometryBuilder
// -------------------------------
// Abstract base class for geometry construction per chamber.
//
// Author: Eric Dumonteil

#include <TVirtualMC.h>
#include <TArrayI.h>
#include <TGeoMatrix.h>
#include "AliRun.h"

#include "AliMUONSlatGeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONChamberGeometry.h"

ClassImp(AliMUONSlatGeometryBuilder)

    Int_t   ConvertSlatNum(Int_t numslat, Int_t quadnum, Int_t fspq);



//______________________________________________________________________________
AliMUONSlatGeometryBuilder::AliMUONSlatGeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder(&muon->Chamber(4), &muon->Chamber(5), &muon->Chamber(6), &muon->Chamber(7), &muon->Chamber(8), &muon->Chamber(9)),
//  : AliMUONVGeometryBuilder(&muon->Chamber(4), &muon->Chamber(5)),
   fMUON(muon)
{
// Standard constructor

}

//______________________________________________________________________________
AliMUONSlatGeometryBuilder::AliMUONSlatGeometryBuilder() 
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONSlatGeometryBuilder::AliMUONSlatGeometryBuilder(const AliMUONSlatGeometryBuilder& rhs)
  : AliMUONVGeometryBuilder(rhs)
{
  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONSlatGeometryBuilder::~AliMUONSlatGeometryBuilder() {
//
}

//______________________________________________________________________________
AliMUONSlatGeometryBuilder& 
AliMUONSlatGeometryBuilder::operator = (const AliMUONSlatGeometryBuilder& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  Fatal("operator=", 
        "Assignment operator is not implemented.");
    
  return *this;  
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::CreateGeometry()
{
 
     Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;

     Float_t angle;
     Float_t *dum=0;

      // define the id of tracking media:
     Int_t idCopper = idtmed[1110];
     Int_t idGlass  = idtmed[1111];
     Int_t idCarbon = idtmed[1112];
     Int_t idRoha   = idtmed[1113];
     Int_t idGas=idtmed[1108]; // medium 9 = Ar-CO2 gas (80%+20%)
     Int_t idAir= idtmed[1100]; // medium 1

      // sensitive area: 40*40 cm**2
     const Float_t sensLength = 40.; 
     const Float_t sensHeight = 40.; 
     const Float_t sensWidth  = 0.5; // according to TDR fig 2.120 
     const Int_t sensMaterial = idGas;
     const Float_t yOverlap   = 1.5; 

     // PCB dimensions in cm; width: 30 mum copper   
     const Float_t pcbLength  = sensLength; 
     const Float_t pcbHeight  = 60.; 
     const Float_t pcbWidth   = 0.003;   
     const Int_t pcbMaterial  = idCopper;

     // Insulating material: 200 mum glass fiber glued to pcb  
     const Float_t insuLength = pcbLength; 
     const Float_t insuHeight = pcbHeight; 
     const Float_t insuWidth  = 0.020;   
     const Int_t insuMaterial = idGlass;

     // Carbon fiber panels: 200mum carbon/epoxy skin   
     const Float_t panelLength = sensLength; 
     const Float_t panelHeight = sensHeight; 
     const Float_t panelWidth  = 0.020;      
     const Int_t panelMaterial = idCarbon;

     // rohacell between the two carbon panels   
     const Float_t rohaLength = sensLength; 
     const Float_t rohaHeight = sensHeight; 
     const Float_t rohaWidth  = 0.5;
     const Int_t rohaMaterial = idRoha;

     // Frame around the slat: 2 sticks along length,2 along height  
     // H: the horizontal ones 
     const Float_t hFrameLength = pcbLength; 
     const Float_t hFrameHeight = 1.5; 
     const Float_t hFrameWidth  = sensWidth; 
     const Int_t hFrameMaterial = idGlass;

     // V: the vertical ones 
     const Float_t vFrameLength = 4.0; 
     const Float_t vFrameHeight = sensHeight + hFrameHeight; 
     const Float_t vFrameWidth  = sensWidth;
     const Int_t vFrameMaterial = idGlass;

     // B: the horizontal border filled with rohacell 
     const Float_t bFrameLength = hFrameLength; 
     const Float_t bFrameHeight = (pcbHeight - sensHeight)/2. - hFrameHeight; 
     const Float_t bFrameWidth  = hFrameWidth;
     const Int_t bFrameMaterial = idRoha;

     // NULOC: 30 mum copper + 200 mum vetronite (same radiation length as 14mum copper)
     const Float_t nulocLength = 2.5; 
     const Float_t nulocHeight = 7.5; 
     const Float_t nulocWidth  = 0.0030 + 0.0014; // equivalent copper width of vetronite; 
     const Int_t   nulocMaterial = idCopper;

     const Float_t slatHeight = pcbHeight; 
     const Float_t slatWidth = sensWidth + 2.*(pcbWidth + insuWidth + 
					       2.* panelWidth + rohaWidth);
     const Int_t slatMaterial = idAir;
     const Float_t dSlatLength = vFrameLength; // border on left and right 

     Float_t spar[3];  
     Int_t i, j;

     // the panel volume contains the rohacell

     Float_t twidth = 2 * panelWidth + rohaWidth; 
     Float_t panelpar[3] = { panelLength/2., panelHeight/2., twidth/2. }; 
     Float_t rohapar[3] = { rohaLength/2., rohaHeight/2., rohaWidth/2. }; 

     // insulating material contains PCB-> gas-> 2 borders filled with rohacell

     twidth = 2*(insuWidth + pcbWidth) + sensWidth;  
     Float_t insupar[3] = { insuLength/2., insuHeight/2., twidth/2. }; 
     twidth -= 2 * insuWidth; 
     Float_t pcbpar[3] = { pcbLength/2., pcbHeight/2., twidth/2. }; 
     Float_t senspar[3] = { sensLength/2., sensHeight/2., sensWidth/2. }; 
     Float_t theight = 2*hFrameHeight + sensHeight;
     Float_t hFramepar[3]={hFrameLength/2., theight/2., hFrameWidth/2.}; 
     Float_t bFramepar[3]={bFrameLength/2., bFrameHeight/2., bFrameWidth/2.}; 
     Float_t vFramepar[3]={vFrameLength/2., vFrameHeight/2., vFrameWidth/2.};
     Float_t nulocpar[3]={nulocLength/2., nulocHeight/2., nulocWidth/2.}; 
     Float_t xx;
     Float_t xxmax = (bFrameLength - nulocLength)/2.; 
     Int_t index=0;
      
    AliMUONChamber *iChamber, *iChamber1, *iChamber2;

    Int_t* fStations = new Int_t[5];
    for (Int_t i=0; i<5; i++) fStations[i] = 1;

    if (fStations[2])
    {
	
//********************************************************************
//                            Station 3                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers

     iChamber = GetChamber(4);
     iChamber1 = iChamber;
     iChamber2 = GetChamber(5);
     
     //iChamber1->GetGeometry()->SetDebug(kTRUE);
     //iChamber2->GetGeometry()->SetDebug(kTRUE);

     if (gAlice->GetModule("DIPO")) {
       // if DIPO is preset, the whole station will be placed in DDIP volume
       iChamber1->GetGeometry()->SetMotherVolume("DDIP");
       iChamber2->GetGeometry()->SetMotherVolume("DDIP");
     }

//      if (gAlice->GetModule("DIPO")) {
//        slats5Mother="DDIP";
//        slats6Mother="DDIP";

//        zoffs5 = zpos1;
//        zoffs6 = zpos2;
//      }
//      else {
//        gMC->Gsvolu("S05M", "TUBE", idAir, tpar, 3);
//        gMC->Gsvolu("S06M", "TUBE", idAir, tpar, 3);
//        gMC->Gspos("S05M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");

//        gMC->Gspos("S06M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
//      }

     // volumes for slat geometry (xx=5,..,10 chamber id): 
     // Sxx0 Sxx1 Sxx2 Sxx3  -->   Slat Mother volumes 
     // SxxG                          -->   Sensitive volume (gas)
     // SxxP                          -->   PCB (copper) 
     // SxxI                          -->   Insulator (vetronite) 
     // SxxC                          -->   Carbon panel 
     // SxxR                          -->   Rohacell
     // SxxH, SxxV                    -->   Horizontal and Vertical frames (vetronite)
     // SB5x                          -->   Volumes for the 35 cm long PCB
     // slat dimensions: slat is a MOTHER volume!!! made of air

     // only for chamber 5: slat 1 has a PCB shorter by 5cm!

     Float_t tlength = 35.;
     Float_t panelpar2[3]  = { tlength/2., panelpar[1],  panelpar[2]}; 
     Float_t rohapar2[3]   = { tlength/2., rohapar[1],   rohapar[2]}; 
     Float_t insupar2[3]   = { tlength/2., insupar[1],   insupar[2]}; 
     Float_t pcbpar2[3]    = { tlength/2., pcbpar[1],    pcbpar[2]}; 
     Float_t senspar2[3]   = { tlength/2., senspar[1],   senspar[2]}; 
     Float_t hFramepar2[3] = { tlength/2., hFramepar[1], hFramepar[2]}; 
     Float_t bFramepar2[3] = { tlength/2., bFramepar[1], bFramepar[2]}; 
     Float_t *dum=0;

     const Int_t nSlats3 = 5;  // number of slats per quadrant
     const Int_t nPCB3[nSlats3] = {3,4,4,3,2}; // n PCB per slat
     const Float_t xpos3[nSlats3] = {31., 0., 0., 0., 0.};
     Float_t slatLength3[nSlats3]; 

     // create and position the slat (mother) volumes 

//      char volNam5[5];
//      char volNam6[5];
     char idSlatCh5[5];
     char idSlatCh6[5];
     Float_t xSlat3;
     Float_t angle = 0.;
     
     Float_t spar2[3];
     for (i = 0; i<nSlats3; i++){
       slatLength3[i] = pcbLength * nPCB3[i] + 2. * dSlatLength; 
       xSlat3 = slatLength3[i]/2. - vFrameLength/2. + xpos3[i]; 
       if (i==1 || i==0) slatLength3[i] -=  2. *dSlatLength; // frame out in PCB with circular border 
       Float_t ySlat31 =  sensHeight * i - yOverlap * i; 
       Float_t ySlat32 = -sensHeight * i + yOverlap * i; 
       spar[0] = slatLength3[i]/2.; 
       spar[1] = slatHeight/2.;
       spar[2] = slatWidth/2. * 1.01; 
       // take away 5 cm from the first slat in chamber 5
       Float_t xSlat32 = 0;
       if (i==1 || i==2) { // 1 pcb is shortened by 5cm
	 spar2[0] = spar[0]-5./2.;
	 xSlat32 = xSlat3 - 5/2.;
       }
       else {
	 spar2[0] = spar[0];
	 xSlat32 = xSlat3;
       }
       spar2[1] = spar[1];
       spar2[2] = spar[2]; 
       Float_t dzCh3=spar[2] * 1.01;
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? -spar[2] : spar[2]; 
//        sprintf(volNam5,"S05%d",i);
//        gMC->Gsvolu(volNam5,"BOX",slatMaterial,spar2,3);
//        gMC->Gspos(volNam5, i*4+1,slats5Mother, xSlat32, ySlat31, zoffs5+zSlat+2.*dzCh3, 0, "ONLY");
//        gMC->Gspos(volNam5, i*4+2,slats5Mother,-xSlat32, ySlat31, zoffs5+zSlat-2.*dzCh3, 0, "ONLY");

       sprintf(idSlatCh5,"LA%d",nSlats3-1+i);
       gMC->Gsvolu(idSlatCh5,"BOX",slatMaterial,spar2,3);
       GetChamber(4)->GetGeometry()->AddEnvelope(idSlatCh5, true, TGeoTranslation(xSlat32, ySlat31, zSlat+2.*dzCh3) ,TGeoRotation("rot1",90,angle,90,90+angle,0,0)
	   );

       sprintf(idSlatCh5,"LA%d",3*nSlats3-2+i);
       gMC->Gsvolu(idSlatCh5,"BOX",slatMaterial,spar2,3);
       GetChamber(4)->GetGeometry()->AddEnvelope(idSlatCh5, true, TGeoTranslation(-xSlat32, ySlat31, zSlat-2.*dzCh3) ,TGeoRotation("rot2",90,180+angle,90,90+angle,180,0)
	   );

       if (i>0) { 

       sprintf(idSlatCh5,"LA%d",nSlats3-1-i);
       gMC->Gsvolu(idSlatCh5,"BOX",slatMaterial,spar2,3);
       GetChamber(4)->GetGeometry()->AddEnvelope(idSlatCh5, true, TGeoTranslation(xSlat32, ySlat32, zSlat+2.*dzCh3) ,TGeoRotation("rot3",90,angle,90,270+angle,180,0)
	   );

       sprintf(idSlatCh5,"LA%d",3*nSlats3-2-i);
       gMC->Gsvolu(idSlatCh5,"BOX",slatMaterial,spar2,3);
       GetChamber(4)->GetGeometry()->AddEnvelope(idSlatCh5, true, TGeoTranslation(-xSlat32, ySlat32, zSlat-2.*dzCh3) ,TGeoRotation("rot4",90,180+angle,90,270+angle,0,0)
	   );
       }

       sprintf(idSlatCh6,"LB%d",nSlats3-1+i);       
       gMC->Gsvolu(idSlatCh6,"BOX",slatMaterial,spar2,3);
       GetChamber(5)->GetGeometry()->AddEnvelope(idSlatCh6, true, TGeoTranslation(xSlat3, ySlat31, zSlat+2.*dzCh3) ,TGeoRotation("rot5",90,angle,90,90+angle,0,0)
	   );
       sprintf(idSlatCh6,"LB%d",3*nSlats3-2+i);
       gMC->Gsvolu(idSlatCh6,"BOX",slatMaterial,spar2,3);
       GetChamber(5)->GetGeometry()->AddEnvelope(idSlatCh6, true, TGeoTranslation(-xSlat3, ySlat31, zSlat-2.*dzCh3) ,TGeoRotation("rot6",90,180+angle,90,90+angle,180,0)
	   );

     if (i>0) { 
       sprintf(idSlatCh6,"LB%d",nSlats3-1-i);
       gMC->Gsvolu(idSlatCh6,"BOX",slatMaterial,spar2,3);
       GetChamber(5)->GetGeometry()->AddEnvelope(idSlatCh6, true, TGeoTranslation(xSlat3, ySlat32, zSlat+2.*dzCh3) ,TGeoRotation("rot7",90,angle,90,270+angle,180,0)
	   );

       sprintf(idSlatCh6,"LB%d",3*nSlats3-2-i);
       gMC->Gsvolu(idSlatCh6,"BOX",slatMaterial,spar2,3);
       GetChamber(5)->GetGeometry()->AddEnvelope(idSlatCh6, true, TGeoTranslation(-xSlat3, ySlat32, zSlat-2.*dzCh3) ,TGeoRotation("rot8",90,180+angle,90,270+angle,0,0)
	   );
      }
	 }
     
     // create the panel volume 
 
     gMC->Gsvolu("S05C","BOX",panelMaterial,panelpar,3);
     gMC->Gsvolu("SB5C","BOX",panelMaterial,panelpar2,3);
     gMC->Gsvolu("S06C","BOX",panelMaterial,panelpar,3);

     // create the rohacell volume 

     gMC->Gsvolu("S05R","BOX",rohaMaterial,rohapar,3);
     gMC->Gsvolu("SB5R","BOX",rohaMaterial,rohapar2,3);
     gMC->Gsvolu("S06R","BOX",rohaMaterial,rohapar,3);

     // create the insulating material volume 

     gMC->Gsvolu("S05I","BOX",insuMaterial,insupar,3);
     gMC->Gsvolu("SB5I","BOX",insuMaterial,insupar2,3);
     gMC->Gsvolu("S06I","BOX",insuMaterial,insupar,3);

     // create the PCB volume 

     gMC->Gsvolu("S05P","BOX",pcbMaterial,pcbpar,3);
     gMC->Gsvolu("SB5P","BOX",pcbMaterial,pcbpar2,3);
     gMC->Gsvolu("S06P","BOX",pcbMaterial,pcbpar,3);
 
     // create the sensitive volumes,
     gMC->Gsvolu("S05G","BOX",sensMaterial,dum,0);
     gMC->Gsvolu("S06G","BOX",sensMaterial,dum,0);


     // create the vertical frame volume 

     gMC->Gsvolu("S05V","BOX",vFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S06V","BOX",vFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 


     gMC->Gsvolu("S05H","BOX",hFrameMaterial,hFramepar,3);
     gMC->Gsvolu("SB5H","BOX",hFrameMaterial,hFramepar2,3);
     gMC->Gsvolu("S06H","BOX",hFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S05B","BOX",bFrameMaterial,bFramepar,3);
     gMC->Gsvolu("SB5B","BOX",bFrameMaterial,bFramepar2,3);
     gMC->Gsvolu("S06B","BOX",bFrameMaterial,bFramepar,3);

     index=0; 
     for (i = 0; i<nSlats3; i++){
	 for (Int_t quadrant=1; quadrant<=4; quadrant++) {

	     if (i==0&&quadrant==2) continue;
	     if (i==0&&quadrant==4) continue;

       sprintf(idSlatCh5,"LA%d",ConvertSlatNum(i,quadrant,4));
       sprintf(idSlatCh6,"LB%d",ConvertSlatNum(i,quadrant,4));
       Float_t xvFrame  = (slatLength3[i] - vFrameLength)/2.;
       Float_t xvFrame2  = xvFrame;

       if ( i==1 || i ==2 ) xvFrame2 -= 5./2.;

       // position the vertical frames 
       if (i!=1 && i!=0) { 
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("S05V", idSlatCh5, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame2,0.,0.));
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("S05V", idSlatCh5, (2*i)*10+quadrant,TGeoTranslation(-xvFrame2,0.,0.));
       GetChamber(5)->GetGeometry()->AddEnvelopeConstituent("S06V", idSlatCh6, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
       GetChamber(5)->GetGeometry()->AddEnvelopeConstituent("S06V", idSlatCh6, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));

       }       
       // position the panels and the insulating material 
       for (j=0; j<nPCB3[i]; j++){
	 if (i==1&&j==0) continue;
	 index++;
	 Float_t xx = sensLength * (-nPCB3[i]/2.+j+.5); 
	 Float_t xx2 = xx + 5/2.; 
	 
	 Float_t zPanel = spar[2] - panelpar[2]; 
	 if ( (i==1 || i==2) && j == nPCB3[i]-1) { // 1 pcb is shortened by 5cm 
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("SB5C", idSlatCh5, 2*index-1,TGeoTranslation(xx,0.,zPanel));
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("SB5C", idSlatCh5, 2*index,TGeoTranslation(xx,0.,-zPanel));
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("SB5I", idSlatCh5, index,TGeoTranslation(xx,0.,0.));
	 }
	 else if ( (i==1 || i==2) && j < nPCB3[i]-1) {
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("S05C", idSlatCh5, 2*index-1,TGeoTranslation(xx2,0.,zPanel));
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("S05C", idSlatCh5, 2*index,TGeoTranslation(xx2,0.,-zPanel));
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("S05I", idSlatCh5, index,TGeoTranslation(xx2,0.,0.));
	 }
	 else {
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("S05C", idSlatCh5, 2*index-1,TGeoTranslation(xx,0.,zPanel));
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("S05C", idSlatCh5, 2*index,TGeoTranslation(xx,0.,-zPanel));
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituent("S05I", idSlatCh5, index,TGeoTranslation(xx,0.,0.));
	 }
       GetChamber(5)->GetGeometry()->AddEnvelopeConstituent("S06C", idSlatCh6, 2*index-1,TGeoTranslation(xx,0.,zPanel));
       GetChamber(5)->GetGeometry()->AddEnvelopeConstituent("S06C", idSlatCh6, 2*index,TGeoTranslation(xx,0.,-zPanel));
       GetChamber(5)->GetGeometry()->AddEnvelopeConstituent("S06I", idSlatCh6, index,TGeoTranslation(xx,0.,0.));
 
       } 
	 }
     }
     
     // position the rohacell volume inside the panel volume
     gMC->Gspos("S05R",1,"S05C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("SB5R",1,"SB5C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06R",1,"S06C",0.,0.,0.,0,"ONLY"); 

     // position the PCB volume inside the insulating material volume
     gMC->Gspos("S05P",1,"S05I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("SB5P",1,"SB5I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06P",1,"S06I",0.,0.,0.,0,"ONLY"); 
     // position the horizontal frame volume inside the PCB volume
     gMC->Gspos("S05H",1,"S05P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("SB5H",1,"SB5P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06H",1,"S06P",0.,0.,0.,0,"ONLY"); 
     // position the sensitive volume inside the horizontal frame volume
     gMC->Gsposp("S05G",1,"S05H",0.,0.,0.,0,"ONLY",senspar,3); 
     gMC->Gsposp("S05G",1,"SB5H",0.,0.,0.,0,"ONLY",senspar2,3); 
     gMC->Gsposp("S06G",1,"S06H",0.,0.,0.,0,"ONLY",senspar,3); 
     // position the border volumes inside the PCB volume
     Float_t yborder = ( pcbHeight - bFrameHeight ) / 2.; 
     gMC->Gspos("S05B",1,"S05P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S05B",2,"S05P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("SB5B",1,"SB5P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("SB5B",2,"SB5P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("S06B",1,"S06P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S06B",2,"S06P",0.,-yborder,0.,0,"ONLY"); 

     // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S05N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S06N","BOX",nulocMaterial,nulocpar,3);
     index = 0;
     Float_t xxmax2 = xxmax - 5./2.;
     for (xx = -xxmax; xx<=xxmax; xx+=2*nulocLength) { 
       index++; 
       gMC->Gspos("S05N",2*index-1,"S05B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S05N",2*index  ,"S05B", xx, 0., bFrameWidth/4., 0, "ONLY");
       if (xx > -xxmax2 && xx< xxmax2) {
	 gMC->Gspos("S05N",2*index-1,"SB5B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
	 gMC->Gspos("S05N",2*index  ,"SB5B", xx, 0., bFrameWidth/4., 0, "ONLY");
       }
       gMC->Gspos("S06N",2*index-1,"S06B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S06N",2*index  ,"S06B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }

     // position the volumes approximating the circular section of the pipe
     Float_t yoffs = sensHeight/2.-yOverlap; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Double_t divpar[3];
     Double_t dydiv= sensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv;
     Int_t imax=0; 
     imax = 1; 
     Float_t rmin = 33.; 
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (pcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = sensWidth/2.; 
       Float_t xvol=(pcbLength+xdiv)/2.;
       Float_t yvol=ydiv + dydiv/2.; 

       for (Int_t quadrant=1; quadrant<=4; quadrant++)
       {
	   sprintf(idSlatCh5,"LA%d",ConvertSlatNum(1,quadrant,4));
	   sprintf(idSlatCh6,"LB%d",ConvertSlatNum(1,quadrant,4));
	   
       GetChamber(4)->GetGeometry()->AddEnvelopeConstituentParam("S05G", idSlatCh5, quadrant*100+imax+4*idiv+1,TGeoTranslation(xvol-(pcbLength * (nPCB3[1]-1)/2. + 35./2.),yvol-pcbLength+yOverlap,0.),3,divpar);
       GetChamber(5)->GetGeometry()->AddEnvelopeConstituentParam("S06G", idSlatCh6,  quadrant*100+imax+4*idiv+1,TGeoTranslation(xvol-pcbLength * nPCB3[1]/2.,yvol-pcbLength+yOverlap,0.),3,divpar);
       }
       
     }
     cout << "Geometry for Station 3...... done" << endl;
    }
    
    if (fStations[3]) {


// //********************************************************************
// //                            Station 4                             **
// //********************************************************************
//      // indices 1 and 2 for first and second chambers in the station
//      // iChamber (first chamber) kept for other quanties than Z,
//      // assumed to be the same in both chambers
 
     iChamber = GetChamber(6);
     iChamber1 = iChamber;
     iChamber2 = GetChamber(7);

     const Int_t nSlats4 = 6;  // number of slats per quadrant
     const Int_t nPCB4[nSlats4] = {4,4,5,5,4,3}; // n PCB per slat
     const Float_t xpos4[nSlats4] = {38.5, 40., 0., 0., 0., 0.};
     Float_t slatLength4[nSlats4];     

//      // create and position the slat (mother) volumes 

     char idSlatCh7[5];
     char idSlatCh8[5];
     Float_t xSlat4;
     Float_t ySlat4;
     angle = 0.;

     for (i = 0; i<nSlats4; i++){
       slatLength4[i] = pcbLength * nPCB4[i] + 2. * dSlatLength; 
       xSlat4 = slatLength4[i]/2. - vFrameLength/2. + xpos4[i]; 
       if (i==1) slatLength4[i] -=  2. *dSlatLength; // frame out in PCB with circular border 
       ySlat4 =  sensHeight * i - yOverlap *i;
       
       spar[0] = slatLength4[i]/2.; 
       spar[1] = slatHeight/2.;
       spar[2] = slatWidth/2.*1.01; 
       Float_t dzCh4=spar[2]*1.01;
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? spar[2] : -spar[2]; 

       sprintf(idSlatCh7,"LC%d",nSlats4-1+i);
       gMC->Gsvolu(idSlatCh7,"BOX",slatMaterial,spar,3);
       GetChamber(6)->GetGeometry()->AddEnvelope(idSlatCh7, true, TGeoTranslation(xSlat4, ySlat4, zSlat+2.*dzCh4));

       sprintf(idSlatCh7,"LC%d",3*nSlats4-2+i);
       gMC->Gsvolu(idSlatCh7,"BOX",slatMaterial,spar,3);
       GetChamber(6)->GetGeometry()->AddEnvelope(idSlatCh7, true, TGeoTranslation(-xSlat4, ySlat4, zSlat-2.*dzCh4));
 
       if (i>0) { 

       sprintf(idSlatCh7,"LC%d",nSlats4-1-i);
       gMC->Gsvolu(idSlatCh7,"BOX",slatMaterial,spar,3);
       GetChamber(6)->GetGeometry()->AddEnvelope(idSlatCh7, true, TGeoTranslation(xSlat4, -ySlat4, zSlat+2.*dzCh4) ,TGeoRotation("rot3",90,angle,90,270+angle,180,0)
	   );

       sprintf(idSlatCh7,"LC%d",3*nSlats4-2-i);
       gMC->Gsvolu(idSlatCh7,"BOX",slatMaterial,spar,3);
       GetChamber(6)->GetGeometry()->AddEnvelope(idSlatCh7, true, TGeoTranslation(-xSlat4, -ySlat4, zSlat-2.*dzCh4) ,TGeoRotation("rot3",90,angle,90,270+angle,180,0)
	   );
       }

       sprintf(idSlatCh8,"LD%d",nSlats4-1+i);
       gMC->Gsvolu(idSlatCh8,"BOX",slatMaterial,spar,3);
       GetChamber(7)->GetGeometry()->AddEnvelope(idSlatCh8, true, TGeoTranslation(xSlat4, ySlat4, zSlat+2.*dzCh4) ,TGeoRotation("rot5",90,angle,90,90+angle,0,0)
	   );
       sprintf(idSlatCh8,"LD%d",3*nSlats4-2+i);
       gMC->Gsvolu(idSlatCh8,"BOX",slatMaterial,spar,3);
       GetChamber(7)->GetGeometry()->AddEnvelope(idSlatCh8, true, TGeoTranslation(-xSlat4, ySlat4, zSlat-2.*dzCh4) ,TGeoRotation("rot6",90,180+angle,90,90+angle,180,0)
	   );
       if (i>0) { 
       sprintf(idSlatCh8,"LD%d",nSlats4-1-i);
       gMC->Gsvolu(idSlatCh8,"BOX",slatMaterial,spar,3);
       GetChamber(7)->GetGeometry()->AddEnvelope(idSlatCh8, true, TGeoTranslation(xSlat4, -ySlat4, zSlat+2.*dzCh4) ,TGeoRotation("rot7",90,angle,90,270+angle,180,0)
	   );
       sprintf(idSlatCh8,"LD%d",3*nSlats4-2-i);
       gMC->Gsvolu(idSlatCh8,"BOX",slatMaterial,spar,3);
       GetChamber(7)->GetGeometry()->AddEnvelope(idSlatCh8, true, TGeoTranslation(-xSlat4, -ySlat4, zSlat-2.*dzCh4) ,TGeoRotation("rot8",90,180+angle,90,270+angle,0,0)
						 );
       }
     }
     

     // create the panel volume 
 
     gMC->Gsvolu("S07C","BOX",panelMaterial,panelpar,3);
     gMC->Gsvolu("S08C","BOX",panelMaterial,panelpar,3);

     // create the rohacell volume 

     gMC->Gsvolu("S07R","BOX",rohaMaterial,rohapar,3);
     gMC->Gsvolu("S08R","BOX",rohaMaterial,rohapar,3);

     // create the insulating material volume 

     gMC->Gsvolu("S07I","BOX",insuMaterial,insupar,3);
     gMC->Gsvolu("S08I","BOX",insuMaterial,insupar,3);

     // create the PCB volume 

     gMC->Gsvolu("S07P","BOX",pcbMaterial,pcbpar,3);
     gMC->Gsvolu("S08P","BOX",pcbMaterial,pcbpar,3);
 
     // create the sensitive volumes,

     gMC->Gsvolu("S07G","BOX",sensMaterial,dum,0);
     gMC->Gsvolu("S08G","BOX",sensMaterial,dum,0);

     // create the vertical frame volume 

     gMC->Gsvolu("S07V","BOX",vFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S08V","BOX",vFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 

     gMC->Gsvolu("S07H","BOX",hFrameMaterial,hFramepar,3);
     gMC->Gsvolu("S08H","BOX",hFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S07B","BOX",bFrameMaterial,bFramepar,3);
     gMC->Gsvolu("S08B","BOX",bFrameMaterial,bFramepar,3);

     index=0; 
     for (i = 0; i<nSlats4; i++){
	 for (Int_t quadrant=1; quadrant<=4; quadrant++) {

	     if (i==0&&quadrant==2) continue;
	     if (i==0&&quadrant==4) continue;

       sprintf(idSlatCh7,"LC%d",ConvertSlatNum(i,quadrant,5));
       sprintf(idSlatCh8,"LD%d",ConvertSlatNum(i,quadrant,5));
       Float_t xvFrame  = (slatLength4[i] - vFrameLength)/2.;

       // position the vertical frames 
       if (i!=1 && i!=0) { 
       GetChamber(6)->GetGeometry()->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
       GetChamber(6)->GetGeometry()->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
       GetChamber(7)->GetGeometry()->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
       GetChamber(7)->GetGeometry()->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
       }
       // position the panels and the insulating material 
       for (j=0; j<nPCB4[i]; j++){
	 index++;
	 Float_t xx = sensLength * (-nPCB4[i]/2.+j+.5); 

	 Float_t zPanel = spar[2] - panelpar[2]; 
       GetChamber(6)->GetGeometry()->AddEnvelopeConstituent("S07C", idSlatCh7, 2*index-1,TGeoTranslation(xx,0.,zPanel));
       GetChamber(6)->GetGeometry()->AddEnvelopeConstituent("S07C", idSlatCh7, 2*index,TGeoTranslation(xx,0.,-zPanel));
       GetChamber(6)->GetGeometry()->AddEnvelopeConstituent("S07I", idSlatCh7, index,TGeoTranslation(xx,0.,0.));
       GetChamber(7)->GetGeometry()->AddEnvelopeConstituent("S08C", idSlatCh8, 2*index-1,TGeoTranslation(xx,0.,zPanel));
       GetChamber(7)->GetGeometry()->AddEnvelopeConstituent("S08C", idSlatCh8, 2*index,TGeoTranslation(xx,0.,-zPanel));
       GetChamber(7)->GetGeometry()->AddEnvelopeConstituent("S08I", idSlatCh8, index,TGeoTranslation(xx,0.,0.));
       }
	 } 
     }

     // position the rohacell volume inside the panel volume
     gMC->Gspos("S07R",1,"S07C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S08R",1,"S08C",0.,0.,0.,0,"ONLY"); 

     // position the PCB volume inside the insulating material volume
     gMC->Gspos("S07P",1,"S07I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S08P",1,"S08I",0.,0.,0.,0,"ONLY"); 
     // position the horizontal frame volume inside the PCB volume
     gMC->Gspos("S07H",1,"S07P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S08H",1,"S08P",0.,0.,0.,0,"ONLY"); 
     // position the sensitive volume inside the horizontal frame volume
     gMC->Gsposp("S07G",1,"S07H",0.,0.,0.,0,"ONLY",senspar,3); 
     gMC->Gsposp("S08G",1,"S08H",0.,0.,0.,0,"ONLY",senspar,3); 
     // position the border volumes inside the PCB volume
     Float_t yborder = ( pcbHeight - bFrameHeight ) / 2.; 
     gMC->Gspos("S07B",1,"S07P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S07B",2,"S07P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("S08B",1,"S08P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S08B",2,"S08P",0.,-yborder,0.,0,"ONLY"); 

//      // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S07N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S08N","BOX",nulocMaterial,nulocpar,3);
     index = 0;
     for (xx = -xxmax; xx<=xxmax; xx+=2*nulocLength) { 
       index++; 
       gMC->Gspos("S07N",2*index-1,"S07B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S07N",2*index  ,"S07B", xx, 0., bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S08N",2*index-1,"S08B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S08N",2*index  ,"S08B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }

//      // position the volumes approximating the circular section of the pipe
     Float_t yoffs = sensHeight/2. - yOverlap; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Double_t divpar[3];
     Double_t dydiv= sensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv;
     Int_t imax=0; 
     imax = 1; 
     Float_t rmin = 40.;
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (pcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = sensWidth/2.; 
       Float_t xvol=(pcbLength+xdiv)/2.+1.999;
       Float_t yvol=ydiv + dydiv/2.;

       for (Int_t quadrant=1; quadrant<=4; quadrant++)
       {
	   sprintf(idSlatCh7,"LC%d",ConvertSlatNum(1,quadrant,5));
	   sprintf(idSlatCh8,"LD%d",ConvertSlatNum(1,quadrant,5));

       GetChamber(6)->GetGeometry()->AddEnvelopeConstituentParam("S07G", idSlatCh7, quadrant*100+imax+4*idiv+1,TGeoTranslation(xvol-pcbLength * nPCB4[1]/2.,yvol-pcbLength+yOverlap,0.),3,divpar);
       GetChamber(7)->GetGeometry()->AddEnvelopeConstituentParam("S08G", idSlatCh8,  quadrant*100+imax+4*idiv+1,TGeoTranslation(xvol-pcbLength * nPCB4[1]/2.,yvol-pcbLength+yOverlap,0.),3,divpar);
       }
     }
     cout << "Geometry for Station 4...... done" << endl;

    }
    
    if (fStations[4]) {
	

// //********************************************************************
// //                            Station 5                             **
// //********************************************************************
//      // indices 1 and 2 for first and second chambers in the station
//      // iChamber (first chamber) kept for other quanties than Z,
//      // assumed to be the same in both chambers

     iChamber = GetChamber(8);
     iChamber1 = iChamber;
     iChamber2 = GetChamber(9);
 
     const Int_t nSlats5 = 7;  // number of slats per quadrant
     const Int_t nPCB5[nSlats5] = {5,5,6,6,5,4,3}; // n PCB per slat
     const Float_t xpos5[nSlats5] = {38.5, 40., 0., 0., 0., 0., 0.};
     Float_t slatLength5[nSlats5]; 

//      // create and position the slat (mother) volumes 

     char idSlatCh9[5];
     char idSlatCh10[5];
     Float_t xSlat5;
     Float_t ySlat5;
     angle = 0.;

     for (i = 0; i<nSlats5; i++){
       slatLength5[i] = pcbLength * nPCB5[i] + 2. * dSlatLength; 
       xSlat5 = slatLength5[i]/2. - vFrameLength/2. +xpos5[i]; 
       if (i==1 || i==0) slatLength5[i] -=  2. *dSlatLength; // frame out in PCB with circular border 
       ySlat5 = sensHeight * i - yOverlap * i;
 
       spar[0] = slatLength5[i]/2.; 
       spar[1] = slatHeight/2.;
       spar[2] = slatWidth/2. * 1.01; 
       Float_t dzCh5=spar[2]*1.01;
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? -spar[2] : spar[2]; 

       sprintf(idSlatCh9,"LE%d",nSlats5-1+i);
       gMC->Gsvolu(idSlatCh9,"BOX",slatMaterial,spar,3);
       GetChamber(8)->GetGeometry()->AddEnvelope(idSlatCh9, true, TGeoTranslation(xSlat5, ySlat5, zSlat+2.*dzCh5));

       sprintf(idSlatCh9,"LE%d",3*nSlats5-2+i);
       gMC->Gsvolu(idSlatCh9,"BOX",slatMaterial,spar,3);
       GetChamber(8)->GetGeometry()->AddEnvelope(idSlatCh9, true, TGeoTranslation(-xSlat5, ySlat5, zSlat-2.*dzCh5));
 
       if (i>0) { 

       sprintf(idSlatCh9,"LE%d",nSlats5-1-i);
       gMC->Gsvolu(idSlatCh9,"BOX",slatMaterial,spar,3);
       GetChamber(8)->GetGeometry()->AddEnvelope(idSlatCh9, true, TGeoTranslation(xSlat5, -ySlat5, zSlat+2.*dzCh5) ,TGeoRotation("rot3",90,angle,90,270+angle,180,0)
	   );

       sprintf(idSlatCh9,"LE%d",3*nSlats5-2-i);
       gMC->Gsvolu(idSlatCh9,"BOX",slatMaterial,spar,3);
       GetChamber(8)->GetGeometry()->AddEnvelope(idSlatCh9, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat-2.*dzCh5) ,TGeoRotation("rot3",90,angle,90,270+angle,180,0)
	   );
       }

       sprintf(idSlatCh10,"LF%d",nSlats5-1+i);
       gMC->Gsvolu(idSlatCh10,"BOX",slatMaterial,spar,3);
       GetChamber(9)->GetGeometry()->AddEnvelope(idSlatCh10, true, TGeoTranslation(xSlat5, ySlat5, zSlat+2.*dzCh5) ,TGeoRotation("rot5",90,angle,90,90+angle,0,0)
	   );

       sprintf(idSlatCh10,"LF%d",3*nSlats5-2+i);
       gMC->Gsvolu(idSlatCh10,"BOX",slatMaterial,spar,3);
       GetChamber(9)->GetGeometry()->AddEnvelope(idSlatCh10, true, TGeoTranslation(-xSlat5, ySlat5, zSlat-2.*dzCh5) ,TGeoRotation("rot6",90,180+angle,90,90+angle,180,0)
       	   );

	     if (i>0) { 

       sprintf(idSlatCh10,"LF%d",nSlats5-1-i);
       gMC->Gsvolu(idSlatCh10,"BOX",slatMaterial,spar,3);
       GetChamber(9)->GetGeometry()->AddEnvelope(idSlatCh10, true, TGeoTranslation(xSlat5, -ySlat5, zSlat+2.*dzCh5) ,TGeoRotation("rot7",90,angle,90,270+angle,180,0)
	   );
       sprintf(idSlatCh10,"LF%d",3*nSlats5-2-i);
       gMC->Gsvolu(idSlatCh10,"BOX",slatMaterial,spar,3);
       GetChamber(9)->GetGeometry()->AddEnvelope(idSlatCh10, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat-2.*dzCh5) ,TGeoRotation("rot8",90,180+angle,90,270+angle,0,0)
	   );
	     }
     }
//      // create the panel volume 
 
     gMC->Gsvolu("S09C","BOX",panelMaterial,panelpar,3);
     gMC->Gsvolu("S10C","BOX",panelMaterial,panelpar,3);

     // create the rohacell volume 

     gMC->Gsvolu("S09R","BOX",rohaMaterial,rohapar,3);
     gMC->Gsvolu("S10R","BOX",rohaMaterial,rohapar,3);

     // create the insulating material volume 

     gMC->Gsvolu("S09I","BOX",insuMaterial,insupar,3);
     gMC->Gsvolu("S10I","BOX",insuMaterial,insupar,3);

     // create the PCB volume 

     gMC->Gsvolu("S09P","BOX",pcbMaterial,pcbpar,3);
     gMC->Gsvolu("S10P","BOX",pcbMaterial,pcbpar,3);
 
     // create the sensitive volumes,

     gMC->Gsvolu("S09G","BOX",sensMaterial,dum,0);
     gMC->Gsvolu("S10G","BOX",sensMaterial,dum,0);

     // create the vertical frame volume 

     gMC->Gsvolu("S09V","BOX",vFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S10V","BOX",vFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 

     gMC->Gsvolu("S09H","BOX",hFrameMaterial,hFramepar,3);
     gMC->Gsvolu("S10H","BOX",hFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S09B","BOX",bFrameMaterial,bFramepar,3);
     gMC->Gsvolu("S10B","BOX",bFrameMaterial,bFramepar,3);

     index=0; 
     for (i = 0; i<nSlats5; i++){
	 for (Int_t quadrant=1; quadrant<=4; quadrant++) {

	     if (i==0&&quadrant==2) continue;
	     if (i==0&&quadrant==4) continue;

       sprintf(idSlatCh9,"LE%d",ConvertSlatNum(i,quadrant,6));
       sprintf(idSlatCh10,"LF%d",ConvertSlatNum(i,quadrant,6));
       Float_t xvFrame  = (slatLength5[i] - vFrameLength)/2.;

       // position the vertical frames 
       if (i!=1 && i!=0) { 
       GetChamber(8)->GetGeometry()->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
       GetChamber(8)->GetGeometry()->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
       GetChamber(9)->GetGeometry()->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
       GetChamber(9)->GetGeometry()->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
       }
       
       // position the panels and the insulating material 
       for (j=0; j<nPCB5[i]; j++){
	 index++;
	 Float_t xx = sensLength * (-nPCB5[i]/2.+j+.5); 

	 Float_t zPanel = spar[2] - panelpar[2]; 
       GetChamber(8)->GetGeometry()->AddEnvelopeConstituent("S09C", idSlatCh9, 2*index-1,TGeoTranslation(xx,0.,zPanel));
       GetChamber(8)->GetGeometry()->AddEnvelopeConstituent("S09C", idSlatCh9, 2*index,TGeoTranslation(xx,0.,-zPanel));
       GetChamber(8)->GetGeometry()->AddEnvelopeConstituent("S09I", idSlatCh9, index,TGeoTranslation(xx,0.,0.));
       GetChamber(9)->GetGeometry()->AddEnvelopeConstituent("S10C", idSlatCh10, 2*index-1,TGeoTranslation(xx,0.,zPanel));
       GetChamber(9)->GetGeometry()->AddEnvelopeConstituent("S10C", idSlatCh10, 2*index,TGeoTranslation(xx,0.,-zPanel));
       GetChamber(9)->GetGeometry()->AddEnvelopeConstituent("S10I", idSlatCh10, index,TGeoTranslation(xx,0.,0.));
       }
	 } 
     }

     // position the rohacell volume inside the panel volume
     gMC->Gspos("S09R",1,"S09C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S10R",1,"S10C",0.,0.,0.,0,"ONLY"); 

     // position the PCB volume inside the insulating material volume
     gMC->Gspos("S09P",1,"S09I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S10P",1,"S10I",0.,0.,0.,0,"ONLY"); 
     // position the horizontal frame volume inside the PCB volume
     gMC->Gspos("S09H",1,"S09P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S10H",1,"S10P",0.,0.,0.,0,"ONLY"); 
     // position the sensitive volume inside the horizontal frame volume
     gMC->Gsposp("S09G",1,"S09H",0.,0.,0.,0,"ONLY",senspar,3); 
     gMC->Gsposp("S10G",1,"S10H",0.,0.,0.,0,"ONLY",senspar,3); 
     // position the border volumes inside the PCB volume
     Float_t yborder = ( pcbHeight - bFrameHeight ) / 2.; 
     gMC->Gspos("S09B",1,"S09P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S09B",2,"S09P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("S10B",1,"S10P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S10B",2,"S10P",0.,-yborder,0.,0,"ONLY"); 

//      // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S09N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S10N","BOX",nulocMaterial,nulocpar,3);
     index = 0;
     for (xx = -xxmax; xx<=xxmax; xx+=2*nulocLength) { 
       index++; 
       gMC->Gspos("S09N",2*index-1,"S09B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S09N",2*index  ,"S09B", xx, 0., bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S10N",2*index-1,"S10B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S10N",2*index  ,"S10B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }

//      // position the volumes approximating the circular section of the pipe
     Float_t yoffs = sensHeight/2. - yOverlap; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Double_t divpar[3];
     Double_t dydiv= sensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv;
     Int_t imax=0; 
     //     for (Int_t islat=0; islat<nSlats3; islat++) imax += nPCB3[islat]; 
     imax = 1; 
     Float_t rmin = 40.; 
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (pcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = sensWidth/2.; 
       Float_t xvol=(pcbLength+xdiv)/2. + 1.999;
       Float_t yvol=ydiv + dydiv/2.;

       for (Int_t quadrant=1; quadrant<=4; quadrant++)
       {
	   sprintf(idSlatCh9,"LE%d",ConvertSlatNum(1,quadrant,6));
	   sprintf(idSlatCh10,"LF%d",ConvertSlatNum(1,quadrant,6));

       GetChamber(8)->GetGeometry()->AddEnvelopeConstituentParam("S09G", idSlatCh9, quadrant*100+imax+4*idiv+1,TGeoTranslation(xvol-pcbLength * nPCB5[1]/2.,yvol-pcbLength+yOverlap,0.),3,divpar);
       GetChamber(9)->GetGeometry()->AddEnvelopeConstituentParam("S10G", idSlatCh10,  quadrant*100+imax+4*idiv+1,TGeoTranslation(xvol-pcbLength * nPCB5[1]/2.,yvol-pcbLength+yOverlap,0.),3,divpar);
       }
     }
     cout << "Geometry for Station 5...... done" << endl;

    }
}


//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::SetTransformations()
{
// Defines the transformations for the station2 chambers.
// ---

  AliMUONChamber* iChamber1 = GetChamber(4);
  Double_t zpos1 = - iChamber1->Z(); 
  iChamber1->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos1));

  AliMUONChamber* iChamber2 = GetChamber(5);
  Double_t zpos2 = - iChamber2->Z(); 
  iChamber2->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos2));

 iChamber1 = GetChamber(6);
  zpos1 = - iChamber1->Z(); 
  iChamber1->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos1));

  iChamber2 = GetChamber(7);
  zpos2 = - iChamber2->Z(); 
  iChamber2->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos2));

 iChamber1 = GetChamber(8);
  zpos1 = - iChamber1->Z(); 
  iChamber1->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos1));

  iChamber2 = GetChamber(9);
  zpos2 = - iChamber2->Z(); 
  iChamber2->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos2));

}

//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::SetSensitiveVolumes()
{
// Defines the sensitive volumes for slat stations chambers.
// ---

  GetChamber(4)->GetGeometry()->SetSensitiveVolume("S05G");
  GetChamber(5)->GetGeometry()->SetSensitiveVolume("S06G");
  GetChamber(6)->GetGeometry()->SetSensitiveVolume("S07G");
  GetChamber(7)->GetGeometry()->SetSensitiveVolume("S08G");
  GetChamber(8)->GetGeometry()->SetSensitiveVolume("S09G");
  GetChamber(9)->GetGeometry()->SetSensitiveVolume("S10G");
}

//______________________________________________________________________________
Int_t  AliMUONSlatGeometryBuilder::ConvertSlatNum(Int_t numslat, Int_t quadnum, Int_t fspq) const
{
	      numslat=numslat+1;
	      if (quadnum==2||quadnum==3) numslat=numslat+fspq;
	      else numslat=fspq+2-numslat;
	      numslat=numslat-1;
	      
	      if (quadnum==3||quadnum==4) numslat=numslat+2*fspq+1;
	      return numslat;
}
