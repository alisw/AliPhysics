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
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
*/

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////

#include <TLorentzVector.h> 
#include "AliMUONvTemp.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliMUONChamber.h"
#include "AliMUONResponseV0.h"
#include "AliMUONResponseTrigger.h"
#include "AliMUONSegmentationV0.h"
#include "AliMUONSegmentationV01.h"
#include "AliMUONSegmentationV02.h"
#include "AliMUONSegmentationV04.h"
#include "AliMUONSegmentationV05.h"
#include "AliMUONSegmentationSlat.h"
#include "AliMUONSegmentationSlatN.h"
#include "AliMUONSegmentationTrigger.h"
#include "AliMUONSegmentationTriggerX.h"
#include "AliMUONSegmentationTriggerY.h"
#include "AliMUONConstants.h"

ClassImp(AliMUONvTemp)
AliMUONvTemp::AliMUONvTemp(const char *name, const char *title)
       : AliMUONv1(name,title)
{

    SetIshunt(1);
    SetMaxStepGas(0.1);
    SetMaxStepAlu(0.1);
//
// Version 0
//
// First define the number of planes that are segmented (1 or 2) by a call
// to SetNsec. 
// Then chose for each chamber (chamber plane) the segmentation 
// and response model.
// They should be equal for the two chambers of each station. In a future
// version this will be enforced.
//
//  
    Int_t chamber;
    // Default response: 5 mm of gas
    AliMUONResponseV0* response0 = new AliMUONResponseV0;
    response0->SetSqrtKx3AndDeriveKx2Kx4(0.7131); // sqrt(0.5085)
    response0->SetSqrtKy3AndDeriveKy2Ky4(0.7642); // sqrt(0.5840)
    response0->SetPitch(0.25); // anode-cathode distance
    response0->SetSigmaIntegration(10.);
    response0->SetChargeSlope(50);
    response0->SetChargeSpread(0.18, 0.18);
    response0->SetMaxAdc(4096);
    response0->SetZeroSuppression(6);
    
    // Response for 4 mm of gas (station 1)
    // automatic consistency with width of sensitive medium in CreateGeometry ????
    AliMUONResponseV0* responseSt1 = new AliMUONResponseV0;
    // Mathieson parameters from L.Kharmandarian's thesis, page 190
    responseSt1->SetSqrtKx3AndDeriveKx2Kx4(0.7000); // sqrt(0.4900)
    responseSt1->SetSqrtKy3AndDeriveKy2Ky4(0.7550); // sqrt(0.5700)
    responseSt1->SetPitch(0.20); // anode-cathode distance
    responseSt1->SetSigmaIntegration(10.);
    // ChargeSlope larger to compensate for the smaller anode-cathode distance
    // and keep the same most probable ADC channel for mip's
    responseSt1->SetChargeSlope(62.5); 
    // assumed proportionality to anode-cathode distance for ChargeSpread
    responseSt1->SetChargeSpread(0.144, 0.144);
    responseSt1->SetMaxAdc(4096);
    responseSt1->SetZeroSuppression(6);
    
//--------------------------------------------------------
// Configuration for Chamber TC1/2  (Station 1) ----------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Float_t rseg1[4]={17.5, 55.2, 71.3, 95.5};
    Int_t   nseg1[4]={4, 4, 2, 1};
//
    chamber=1;
//^^^^^^^^^
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationV01 *seg11=new AliMUONSegmentationV01(4);
    
    seg11->SetSegRadii(rseg1);
    seg11->SetPadSize(2.4, 0.4); // smaller pad size
    seg11->SetDAnod(0.20); // smaller distance between anode wires
    seg11->SetPadDivision(nseg1);
    
    SetSegmentationModel(chamber-1, 1, seg11);
//
    AliMUONSegmentationV02 *seg12=new AliMUONSegmentationV02(4);
    seg12->SetSegRadii(rseg1); 
    seg12->SetPadSize(0.6, 1.6); // smaller pad size
    seg12->SetDAnod(0.20); // smaller distance between anode wires
    seg12->SetPadDivision(nseg1);

    SetSegmentationModel(chamber-1, 2, seg12);

    SetResponseModel(chamber-1, responseSt1); // special response	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
    chamber=2;
//^^^^^^^^^
//
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationV01 *seg21=new AliMUONSegmentationV01(4);
    seg21->SetSegRadii(rseg1);
    seg21->SetPadSize(2.4, 0.4); // smaller pad size
    seg21->SetDAnod(0.20); // smaller distance between anode wires
    seg21->SetPadDivision(nseg1);
    SetSegmentationModel(chamber-1, 1, seg21);
//
    AliMUONSegmentationV02 *seg22=new AliMUONSegmentationV02(4);
    seg22->SetSegRadii(rseg1); 
    seg22->SetPadSize(0.6, 1.6); // smaller pad size
    seg22->SetDAnod(0.20); // smaller distance between anode wires
    seg22->SetPadDivision(nseg1);
    SetSegmentationModel(chamber-1, 2, seg22);
    
    SetResponseModel(chamber-1, responseSt1); // special response
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
	    
//
//--------------------------------------------------------
// Configuration for Chamber TC3/4 (Station 2) -----------
///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Float_t rseg2[4]={23.5, 47.1, 87.7, 122.5};
    Int_t   nseg2[4]={4, 4, 2, 1};
//
    chamber=3;
//^^^^^^^^^
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationV01 *seg31=new AliMUONSegmentationV01(4);
    seg31->SetSegRadii(rseg2);
    seg31->SetPadSize(3.0, 0.5);
    seg31->SetDAnod(3.0/3./4);
    seg31->SetPadDivision(nseg2);
    SetSegmentationModel(chamber-1, 1, seg31);
//
    AliMUONSegmentationV02 *seg32=new AliMUONSegmentationV02(4);
    seg32->SetSegRadii(rseg2); 
    seg32->SetPadSize(0.75, 2.0);
    seg32->SetPadDivision(nseg2);
    seg32->SetDAnod(3.0/3./4);
    
    SetSegmentationModel(chamber-1, 2, seg32);
    
    SetResponseModel(chamber-1, response0);	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
    chamber=4;
//^^^^^^^^^
//
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationV01 *seg41=new AliMUONSegmentationV01(4);
    seg41->SetSegRadii(rseg2);
    seg41->SetPadSize(3.0, 0.5);
    seg41->SetDAnod(3.0/3./4);
    seg41->SetPadDivision(nseg2);
    SetSegmentationModel(chamber-1, 1, seg41);
//
    AliMUONSegmentationV02 *seg42=new AliMUONSegmentationV02(4);
    seg42->SetSegRadii(rseg2); 
    seg42->SetPadSize(0.75, 2.0);
    seg42->SetPadDivision(nseg2);
    seg42->SetDAnod(3.0/3./4);
    
    SetSegmentationModel(chamber-1, 2, seg42);
    
    SetResponseModel(chamber-1, response0);	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
    
//--------------------------------------------------------
// Configuration for Chamber TC5/6  (Station 3) ----------          
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Int_t   nseg3[3]={4, 4, 2};
    Int_t   npcb5[28] = {0,0,1,0,
			 0,1,1,0,
			 0,2,1,0,
			 0,1,1,0, 
			 0,2,1,0, 
			 0,1,1,0, 
			 0,0,1,0};
    
    Float_t shift = 1.5/2.;
    Float_t xpos5[7]    = {2., 2., 2., 34., 2., 2., 2.};
    Float_t ypos5       = -(20.+3.*(40.-2.*shift));

    chamber=5;
    SetNsec(chamber-1,2);
    AliMUONSegmentationSlat *seg51=new AliMUONSegmentationSlat(4);
    seg51->SetNSlats(7); 
    seg51->SetShift(shift);  
    seg51->SetNPCBperSector(npcb5); 
    seg51->SetSlatXPositions(xpos5);
    seg51->SetSlatYPosition(ypos5);
    seg51->SetPadSize(10.,0.5);
    seg51->SetDAnod(0.25);
    seg51->SetPadDivision(nseg3);
    SetSegmentationModel(chamber-1, 1, seg51);
    
    AliMUONSegmentationSlatN *seg52=new AliMUONSegmentationSlatN(4);
    seg52->SetNSlats(7); 
    seg52->SetShift(shift);  
    seg52->SetNPCBperSector(npcb5); 
    seg52->SetSlatXPositions(xpos5);
    seg52->SetSlatYPosition(ypos5);
    seg52->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
    seg52->SetDAnod(0.25);
    seg52->SetPadDivision(nseg3);
    SetSegmentationModel(chamber-1, 2, seg52);
    SetResponseModel(chamber-1, response0);	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
    chamber=6;
    SetNsec(chamber-1,2);
    AliMUONSegmentationSlat *seg61=new AliMUONSegmentationSlat(4);
    seg61->SetNSlats(7); 
    seg61->SetShift(shift);  
    seg61->SetNPCBperSector(npcb5); 
    seg61->SetSlatXPositions(xpos5);
    seg61->SetSlatYPosition(ypos5);
    seg61->SetPadSize(10.,0.5);
    seg61->SetDAnod(0.25);
    seg61->SetPadDivision(nseg3);
    SetSegmentationModel(chamber-1, 1, seg61);
    
    AliMUONSegmentationSlatN *seg62=new AliMUONSegmentationSlatN(4);
    seg62->SetNSlats(7); 
    seg62->SetShift(shift);  
    seg62->SetNPCBperSector(npcb5); 
    seg62->SetSlatXPositions(xpos5);
    seg62->SetSlatYPosition(ypos5);
    seg62->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
    seg62->SetDAnod(0.25);
    seg62->SetPadDivision(nseg3);
    SetSegmentationModel(chamber-1, 2, seg62);
    SetResponseModel(chamber-1, response0);	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
//--------------------------------------------------------
// Configuration for Chamber TC7/8  (Station 4) ----------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    Int_t   nseg4[4]={4, 4, 2, 1};
    
    chamber=7;
//^^^^^^^^^
    
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationSlat *seg71=new AliMUONSegmentationSlat(4);
    Int_t npcb7[44] = {0,0,0,3,
		       0,0,2,2,
		       0,0,3,2,
		       0,2,2,1,
		       0,2,2,1,
		       0,1,2,1, 
		       0,2,2,1, 
		       0,2,2,1, 
		       0,0,3,2, 
		       0,0,2,2, 
		       0,0,0,3};
    Float_t xpos7[11]   = {2., 2., 2., 2., 2., 39.5, 2., 2., 2., 2., 2.};
    Float_t ypos7       = -(20.+5.*(40.-2.*shift));
    
    seg71->SetNSlats(11);  
    seg71->SetShift(shift);  
    seg71->SetNPCBperSector(npcb7); 
    seg71->SetSlatXPositions(xpos7);
    seg71->SetSlatYPosition(ypos7);
    
    seg71->SetPadSize(10.,0.5);
    seg71->SetDAnod(0.25);
    seg71->SetPadDivision(nseg4);
    SetSegmentationModel(chamber-1, 1, seg71);
    
    AliMUONSegmentationSlatN *seg72=new AliMUONSegmentationSlatN(4);
    
    SetSegmentationModel(chamber-1, 2, seg72);
    seg72->SetNSlats(11);  
    seg72->SetShift(shift);   
    seg72->SetNPCBperSector(npcb7); 
    seg72->SetSlatXPositions(xpos7);
    seg72->SetSlatYPosition(ypos7);
    seg72->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
    seg72->SetDAnod(0.25);
    seg72->SetPadDivision(nseg4);
    
    SetResponseModel(chamber-1, response0);	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
    chamber=8;
//^^^^^^^^^
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationSlat *seg81=new AliMUONSegmentationSlat(4);
    
    seg81->SetNSlats(11);  
    seg81->SetShift(shift);  
    seg81->SetNPCBperSector(npcb7); 
    seg81->SetSlatXPositions(xpos7);
    seg81->SetSlatYPosition(ypos7);
    seg81->SetPadSize(10.,0.5);
    seg81->SetDAnod(0.25);
    seg81->SetPadDivision(nseg4);
    SetSegmentationModel(chamber-1, 1, seg81);

    AliMUONSegmentationSlat *seg82=new AliMUONSegmentationSlatN(4);

    SetSegmentationModel(chamber-1, 2, seg82);
    seg82->SetNSlats(11);  
    seg82->SetShift(shift);  
    seg82->SetNPCBperSector(npcb7); 
    seg82->SetSlatXPositions(xpos7);
    seg82->SetSlatYPosition(ypos7);
    seg82->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
    seg82->SetDAnod(0.25);
    seg82->SetPadDivision(nseg4);
    
    SetResponseModel(chamber-1, response0);	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    

//--------------------------------------------------------
// Configuration for Chamber TC9/10  (Station 5) ---------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    chamber=9;
//^^^^^^^^^

    SetNsec(chamber-1,2);
//
    AliMUONSegmentationSlat *seg91=new AliMUONSegmentationSlat(4);
    Int_t   npcb9[52] = {0,0,0,3,
			 0,0,0,4,
			 0,0,2,3,
			 0,0,3,3,
			 0,2,2,2,
			 0,2,2,2,
			 0,1,2,2, 
			 0,2,2,2, 
			 0,2,2,2, 
			 0,0,3,3, 
			 0,0,2,3, 
			 0,0,0,4, 
			 0,0,0,3};   
    Float_t xpos9[13]   = {2., 2., 2., 2., 2., 2., 39.5, 2., 2., 2., 2., 2., 2.};
    Float_t ypos9       = -(20.+6.*(40.-2.*shift));
    
    seg91->SetNSlats(13);  
    seg91->SetShift(shift);  
    seg91->SetNPCBperSector(npcb9); 
    seg91->SetSlatXPositions(xpos9);
    seg91->SetSlatYPosition(ypos9);
    seg91->SetPadSize(10.,0.5);
    seg91->SetDAnod(0.25);
    seg91->SetPadDivision(nseg4);
    SetSegmentationModel(chamber-1, 1, seg91);
    
    AliMUONSegmentationSlatN *seg92=new AliMUONSegmentationSlatN(4);
    
    SetSegmentationModel(chamber-1, 2, seg92);
    seg92->SetNSlats(13);  
    seg92->SetShift(shift);   
    seg92->SetNPCBperSector(npcb9); 
    seg92->SetSlatXPositions(xpos9);
    seg92->SetSlatYPosition(ypos9);
    seg92->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
    seg92->SetDAnod(0.25);
    seg92->SetPadDivision(nseg4);
    
    SetResponseModel(chamber-1, response0);	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
    chamber=10;
//^^^^^^^^^
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationSlat *seg101=new AliMUONSegmentationSlat(4);
    
    seg101->SetNSlats(13);  
    seg101->SetShift(shift);  
    seg101->SetNPCBperSector(npcb9); 
    seg101->SetSlatXPositions(xpos9);
    seg101->SetSlatYPosition(ypos9);
    seg101->SetPadSize(10.,0.5);
    seg101->SetDAnod(0.25);
    seg101->SetPadDivision(nseg4);
    SetSegmentationModel(chamber-1, 1, seg101);
    
    AliMUONSegmentationSlatN *seg102=new AliMUONSegmentationSlatN(4);
    
    SetSegmentationModel(chamber-1, 2, seg102);
    seg102->SetNSlats(13);  
    seg102->SetShift(shift);   
    seg102->SetNPCBperSector(npcb9); 
    seg102->SetSlatXPositions(xpos9);
    seg102->SetSlatYPosition(ypos9);
    seg102->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
    seg102->SetDAnod(0.25);
    seg102->SetPadDivision(nseg4);
    
    SetResponseModel(chamber-1, response0);	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
//--------------------------------------------------------
// Configuration for Trigger Stations -------------------- 
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// Cluster-size off
    AliMUONResponseTrigger* responseTrigger0 =  new AliMUONResponseTrigger;
// Cluster-size on  
// AliMUONResponseTriggerV1* responseTrigger0 =  new AliMUONResponseTriggerV1;
    
    chamber=11;
    SetNsec(chamber-1,2);
    AliMUONSegmentationTriggerX *seg111=new AliMUONSegmentationTriggerX;
    SetSegmentationModel(chamber-1, 1, seg111);
    AliMUONSegmentationTriggerY *seg112=new AliMUONSegmentationTriggerY;
    SetSegmentationModel(chamber-1, 2, seg112);
    
    SetResponseModel(chamber-1, responseTrigger0);      
    Chamber(chamber-1).SetChargeCorrel(0); // same charge on cathodes
    
 
    chamber=12;
    SetNsec(chamber-1,2);
    AliMUONSegmentationTriggerX *seg121=new AliMUONSegmentationTriggerX;
    SetSegmentationModel(chamber-1, 1, seg121);
    AliMUONSegmentationTriggerY *seg122=new AliMUONSegmentationTriggerY;
    SetSegmentationModel(chamber-1, 2, seg122);
    
    SetResponseModel(chamber-1, responseTrigger0);
    Chamber(chamber-1).SetChargeCorrel(0); // same charge on cathodes
    
    chamber=13;
    SetNsec(chamber-1,2);
    AliMUONSegmentationTriggerX *seg131=new AliMUONSegmentationTriggerX;
    SetSegmentationModel(chamber-1, 1, seg131);
    AliMUONSegmentationTriggerY *seg132=new AliMUONSegmentationTriggerY;
    SetSegmentationModel(chamber-1, 2, seg132);
    SetResponseModel(chamber-1, responseTrigger0);      
    Chamber(chamber-1).SetChargeCorrel(0); // same charge on cathodes
    
    chamber=14;
    SetNsec(chamber-1,2);
    AliMUONSegmentationTriggerX *seg141=new AliMUONSegmentationTriggerX;
    SetSegmentationModel(chamber-1, 1, seg141);
    AliMUONSegmentationTriggerY *seg142=new AliMUONSegmentationTriggerY;
    SetSegmentationModel(chamber-1, 2, seg142);
    
    SetResponseModel(chamber-1, responseTrigger0); 
    Chamber(chamber-1).SetChargeCorrel(0); // same charge on cathodes
    
}
//___________________________________________
void AliMUONvTemp::CreateGeometry()
{
//
//   Note: all chambers have the same structure, which could be 
//   easily parameterised. This was intentionally not done in order
//   to give a starting point for the implementation of the actual 
//   design of each station. 
  Int_t *idtmed = fIdtmed->GetArray()-1099;

//   Distance between Stations
//
     Float_t bpar[3];
     Float_t tpar[3];
//      Float_t pgpar[10];
     Float_t zpos1, zpos2, zfpos;
     // Outer excess and inner recess for mother volume radius
     // with respect to ROuter and RInner
     Float_t dframep=.001; // Value for station 3 should be 6 ...
     // Width (RdPhi) of the frame crosses for stations 1 and 2 (cm)
//      Float_t dframep1=.001;
     Float_t dframep1 = 11.0;
//      Bool_t frameCrosses=kFALSE;     
     Bool_t frameCrosses=kTRUE;     
     
//      Float_t dframez=0.9;
     // Half of the total thickness of frame crosses (including DAlu)
     // for each chamber in stations 1 and 2:
     // 3% of X0 of composite material,
     // but taken as Aluminium here, with same thickness in number of X0
     Float_t dframez = 3. * 8.9 / 100;
//      Float_t dr;
     Float_t dstation;

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
//   pointer to the current chamber
     Int_t idAlu1=idtmed[1103]; // medium 4
     Int_t idAlu2=idtmed[1104]; // medium 5
//     Int_t idAlu1=idtmed[1100];
//     Int_t idAlu2=idtmed[1100];
     Int_t idAir=idtmed[1100]; // medium 1
//      Int_t idGas=idtmed[1105]; // medium 6 = Ar-isoC4H10 gas
     Int_t idGas=idtmed[1108]; // medium 9 = Ar-CO2 gas (80%+20%)
     

     AliMUONChamber *iChamber, *iChamber1, *iChamber2;
     Int_t stations[5] = {1, 1, 1, 1, 1};
     
     if (stations[0]) {
	 
//********************************************************************
//                            Station 1                             **
//********************************************************************
//  CONCENTRIC
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[0];
     iChamber2 =(AliMUONChamber*) (*fChambers)[1];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = zpos2 - zpos1;
     // DGas decreased from standard one (0.5)
     iChamber->SetDGas(0.4); iChamber2->SetDGas(0.4);
     // DAlu increased from standard one (3% of X0),
     // because more electronics with smaller pads
     iChamber->SetDAlu(3.5 * 8.9 / 100.); iChamber2->SetDAlu(3.5 * 8.9 / 100.);
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/8;

     gMC->Gsvolu("C01M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C02M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C01M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C02M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");     
// // Aluminium frames
// // Outer frames
//      pgpar[0] = 360/12/2;
//      pgpar[1] = 360.;
//      pgpar[2] = 12.;
//      pgpar[3] =   2;
//      pgpar[4] = -dframez/2;
//      pgpar[5] = iChamber->ROuter();
//      pgpar[6] = pgpar[5]+dframep1;
//      pgpar[7] = +dframez/2;
//      pgpar[8] = pgpar[5];
//      pgpar[9] = pgpar[6];
//      gMC->Gsvolu("C01O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gsvolu("C02O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gspos("C01O",1,"C01M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("C01O",2,"C01M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("C02O",1,"C02M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("C02O",2,"C02M", 0.,0.,+zfpos,  0,"ONLY");
// //
// // Inner frame
//      tpar[0]= iChamber->RInner()-dframep1;
//      tpar[1]= iChamber->RInner();
//      tpar[2]= dframez/2;
//      gMC->Gsvolu("C01I", "TUBE", idAlu1, tpar, 3);
//      gMC->Gsvolu("C02I", "TUBE", idAlu1, tpar, 3);

//      gMC->Gspos("C01I",1,"C01M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("C01I",2,"C01M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("C02I",1,"C02M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("C02I",2,"C02M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     if (frameCrosses) {
         // outside gas
         // security for inside mother volume
	 bpar[0] = (iChamber->ROuter() - iChamber->RInner())
	   * TMath::Cos(TMath::ASin(dframep1 /
				   (iChamber->ROuter() - iChamber->RInner())))
	   / 2.0;
	 bpar[1] = dframep1/2;
	 // total thickness will be (4 * bpar[2]) for each chamber,
	 // which has to be equal to (2 * dframez) - DAlu
	 bpar[2] = (2.0 * dframez - iChamber->DAlu()) / 4.0;
	 gMC->Gsvolu("C01B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("C02B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("C01B",1,"C01M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01B",2,"C01M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01B",3,"C01M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C01B",4,"C01M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C01B",5,"C01M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01B",6,"C01M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01B",7,"C01M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C01B",8,"C01M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 
	 gMC->Gspos("C02B",1,"C02M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02B",2,"C02M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02B",3,"C02M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C02B",4,"C02M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C02B",5,"C02M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02B",6,"C02M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02B",7,"C02M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C02B",8,"C02M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->RInner();
     tpar[1]= iChamber->ROuter();
     tpar[2] = (iChamber->DGas()+iChamber->DAlu())/2;
     gMC->Gsvolu("C01A", "TUBE",  idAlu2, tpar, 3);
     gMC->Gsvolu("C02A", "TUBE",idAlu2, tpar, 3);
     gMC->Gspos("C01A", 1, "C01M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C02A", 1, "C02M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->DGas();
     tpar[2] = iChamber->DGas()/2;
     gMC->Gsvolu("C01G", "TUBE", idGas, tpar, 3);
     gMC->Gsvolu("C02G", "TUBE", idGas, tpar, 3);
     gMC->Gspos("C01G", 1, "C01A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C02G", 1, "C02A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas
     // NONE: chambers are sensitive everywhere
//      if (frameCrosses) {

// 	 dr = (iChamber->ROuter() - iChamber->RInner());
// 	 bpar[0] = TMath::Sqrt(dr*dr-dframep1*dframep1/4)/2;
// 	 bpar[1] = dframep1/2;
// 	 bpar[2] = iChamber->DGas()/2;
// 	 gMC->Gsvolu("C01F", "BOX", idAlu1, bpar, 3);
// 	 gMC->Gsvolu("C02F", "BOX", idAlu1, bpar, 3);
	 
// 	 gMC->Gspos("C01F",1,"C01G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("C01F",2,"C01G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("C01F",3,"C01G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
// 	 gMC->Gspos("C01F",4,"C01G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
	 
// 	 gMC->Gspos("C02F",1,"C02G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("C02F",2,"C02G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("C02F",3,"C02G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
// 	 gMC->Gspos("C02F",4,"C02G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
//      }
     }
     if (stations[1]) {
	 
//********************************************************************
//                            Station 2                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[2];
     iChamber2 =(AliMUONChamber*) (*fChambers)[3];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = zpos2 - zpos1;
     // DGas and DAlu not changed from standard values
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/10;

     gMC->Gsvolu("C03M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C04M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C03M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C04M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");

// // Aluminium frames
// // Outer frames
//      pgpar[0] = 360/12/2;
//      pgpar[1] = 360.;
//      pgpar[2] = 12.;
//      pgpar[3] =   2;
//      pgpar[4] = -dframez/2;
//      pgpar[5] = iChamber->ROuter();
//      pgpar[6] = pgpar[5]+dframep;
//      pgpar[7] = +dframez/2;
//      pgpar[8] = pgpar[5];
//      pgpar[9] = pgpar[6];
//      gMC->Gsvolu("C03O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gsvolu("C04O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gspos("C03O",1,"C03M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("C03O",2,"C03M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("C04O",1,"C04M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("C04O",2,"C04M", 0.,0.,+zfpos,  0,"ONLY");
// //
// // Inner frame
//      tpar[0]= iChamber->RInner()-dframep;
//      tpar[1]= iChamber->RInner();
//      tpar[2]= dframez/2;
//      gMC->Gsvolu("C03I", "TUBE", idAlu1, tpar, 3);
//      gMC->Gsvolu("C04I", "TUBE", idAlu1, tpar, 3);

//      gMC->Gspos("C03I",1,"C03M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("C03I",2,"C03M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("C04I",1,"C04M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("C04I",2,"C04M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     if (frameCrosses) {
         // outside gas
         // security for inside mother volume
	 bpar[0] = (iChamber->ROuter() - iChamber->RInner())
	   * TMath::Cos(TMath::ASin(dframep1 /
				   (iChamber->ROuter() - iChamber->RInner())))
	   / 2.0;
	 bpar[1] = dframep1/2;
	 // total thickness will be (4 * bpar[2]) for each chamber,
	 // which has to be equal to (2 * dframez) - DAlu
	 bpar[2] = (2.0 * dframez - iChamber->DAlu()) / 4.0;
	 gMC->Gsvolu("C03B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("C04B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("C03B",1,"C03M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03B",2,"C03M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03B",3,"C03M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C03B",4,"C03M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C03B",5,"C03M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03B",6,"C03M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03B",7,"C03M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C03B",8,"C03M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 
	 gMC->Gspos("C04B",1,"C04M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04B",2,"C04M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04B",3,"C04M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C04B",4,"C04M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C04B",5,"C04M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04B",6,"C04M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04B",7,"C04M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C04B",8,"C04M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->RInner();
     tpar[1]= iChamber->ROuter();
     tpar[2] = (iChamber->DGas()+iChamber->DAlu())/2;
     gMC->Gsvolu("C03A", "TUBE", idAlu2, tpar, 3);
     gMC->Gsvolu("C04A", "TUBE", idAlu2, tpar, 3);
     gMC->Gspos("C03A", 1, "C03M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C04A", 1, "C04M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->DGas();
     tpar[2] = iChamber->DGas()/2;
     gMC->Gsvolu("C03G", "TUBE", idGas, tpar, 3);
     gMC->Gsvolu("C04G", "TUBE", idGas, tpar, 3);
     gMC->Gspos("C03G", 1, "C03A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C04G", 1, "C04A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     // NONE: chambers are sensitive everywhere
//      if (frameCrosses) {

// 	 dr = (iChamber->ROuter() - iChamber->RInner());
// 	 bpar[0] = TMath::Sqrt(dr*dr-dframep1*dframep1/4)/2;
// 	 bpar[1] = dframep1/2;
// 	 bpar[2] = iChamber->DGas()/2;
// 	 gMC->Gsvolu("C03F", "BOX", idAlu1, bpar, 3);
// 	 gMC->Gsvolu("C04F", "BOX", idAlu1, bpar, 3);
	 
// 	 gMC->Gspos("C03F",1,"C03G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("C03F",2,"C03G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("C03F",3,"C03G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
// 	 gMC->Gspos("C03F",4,"C03G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
	 
// 	 gMC->Gspos("C04F",1,"C04G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("C04F",2,"C04G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("C04F",3,"C04G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
// 	 gMC->Gspos("C04F",4,"C04G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
//      }
     }
     // define the id of tracking media:
     Int_t idCopper = idtmed[1110];
     Int_t idGlass  = idtmed[1111];
     Int_t idCarbon = idtmed[1112];
     Int_t idRoha   = idtmed[1113];

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
     
     if (stations[2]) {
	 
//********************************************************************
//                            Station 3                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[4];
     iChamber2 =(AliMUONChamber*) (*fChambers)[5];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = zpos2 - zpos1;

//      zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2; // not used any more
//
//   Mother volume
//     tpar[0] = iChamber->RInner()-vFrameLength; 
     tpar[0] = 30.; 
//     tpar[1] = (iChamber->ROuter()+dframep)*TMath::Sqrt(2.);
     tpar[1] = 160.;
     tpar[2] = dstation/4;
     gMC->Gsvolu("C05M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C06M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C05M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C06M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
 
     // volumes for slat geometry (xx=5,..,10 chamber id): 
     // Sxx0 Sxx1 Sxx2 Sxx3  -->   Slat Mother volumes 
     // SxxG                          -->   Sensitive volume (gas)
     // SxxP                          -->   PCB (copper) 
     // SxxI                          -->   Insulator (vetronite) 
     // SxxC                          -->   Carbon panel 
     // SxxR                          -->   Rohacell
     // SxxH, SxxV                    -->   Horizontal and Vertical frames (vetronite)

     // slat dimensions: slat is a MOTHER volume!!! made of air

     const Int_t nSlats3 = 4;  // number of slats per quadrant
     const Int_t nPCB3[nSlats3] = {2, 2, 2, 1}; // n PCB per slat
     const Float_t xpos3[nSlats3] = {32., 40., 0., 0.};
     Float_t slatLength3[nSlats3]; 

     // create and position the slat (mother) volumes 

     char volNam5[5];
     char volNam6[5];
     Float_t xSlat3;

     for (i = 0; i<nSlats3; i++){
       slatLength3[i] = pcbLength * nPCB3[i] + 2. * dSlatLength; 
       xSlat3 = slatLength3[i]/2. - vFrameLength/2. + xpos3[i]; 
       if (i==1) slatLength3[i] -=  2. *dSlatLength; // frame out in PCB with circular border 
       Float_t ySlat31 =  sensHeight * i - yOverlap * i; 
       Float_t ySlat32 = -sensHeight * i + yOverlap * i; 
       spar[0] = slatLength3[i]/2.; 
       spar[1] = slatHeight/2.;
       spar[2] = slatWidth/2. * 1.01; 
       Float_t dzCh3=spar[2] * 1.01;
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? spar[2] : -spar[2]; 
       sprintf(volNam5,"S05%d",i);
       gMC->Gsvolu(volNam5,"BOX",slatMaterial,spar,3);
       gMC->Gspos(volNam5, i*4+1,"C05M", xSlat3, ySlat31, zSlat+2.*dzCh3, 0, "ONLY");
       gMC->Gspos(volNam5, i*4+2,"C05M",-xSlat3, ySlat31, zSlat-2.*dzCh3, 0, "ONLY");
       
       if (i>0) { 
	 gMC->Gspos(volNam5, i*4+3,"C05M", xSlat3, ySlat32, zSlat+2.*dzCh3, 0, "ONLY");
	 gMC->Gspos(volNam5, i*4+4,"C05M",-xSlat3, ySlat32, zSlat-2.*dzCh3, 0, "ONLY");
       }

       sprintf(volNam6,"S06%d",i);
       gMC->Gsvolu(volNam6,"BOX",slatMaterial,spar,3);
       gMC->Gspos(volNam6, i*4+1,"C06M", xSlat3, ySlat31, zSlat+2.*dzCh3, 0, "ONLY");
       gMC->Gspos(volNam6, i*4+2,"C06M",-xSlat3, ySlat31, zSlat-2.*dzCh3, 0, "ONLY");
       if (i>0) { 
	 gMC->Gspos(volNam6, i*4+3,"C06M", xSlat3, ySlat32, zSlat+2.*dzCh3, 0, "ONLY");
	 gMC->Gspos(volNam6, i*4+4,"C06M",-xSlat3, ySlat32, zSlat-2.*dzCh3, 0, "ONLY");
       }
     }

     // create the panel volume 
 
     gMC->Gsvolu("S05C","BOX",panelMaterial,panelpar,3);
     gMC->Gsvolu("S06C","BOX",panelMaterial,panelpar,3);

     // create the rohacell volume 

     gMC->Gsvolu("S05R","BOX",rohaMaterial,rohapar,3);
     gMC->Gsvolu("S06R","BOX",rohaMaterial,rohapar,3);

     // create the insulating material volume 

     gMC->Gsvolu("S05I","BOX",insuMaterial,insupar,3);
     gMC->Gsvolu("S06I","BOX",insuMaterial,insupar,3);

     // create the PCB volume 

     gMC->Gsvolu("S05P","BOX",pcbMaterial,pcbpar,3);
     gMC->Gsvolu("S06P","BOX",pcbMaterial,pcbpar,3);
 
     // create the sensitive volumes,
     gMC->Gsvolu("S05G","BOX",sensMaterial,0,0);
     gMC->Gsvolu("S06G","BOX",sensMaterial,0,0);


     // create the vertical frame volume 

     gMC->Gsvolu("S05V","BOX",vFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S06V","BOX",vFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 

     gMC->Gsvolu("S05H","BOX",hFrameMaterial,hFramepar,3);
     gMC->Gsvolu("S06H","BOX",hFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S05B","BOX",bFrameMaterial,bFramepar,3);
     gMC->Gsvolu("S06B","BOX",bFrameMaterial,bFramepar,3);

     index=0; 
     for (i = 0; i<nSlats3; i++){
       sprintf(volNam5,"S05%d",i);
       sprintf(volNam6,"S06%d",i);
       Float_t xvFrame  = (slatLength3[i] - vFrameLength)/2.;
       // position the vertical frames 
       if (i!=1) { 
	 gMC->Gspos("S05V",2*i-1,volNam5, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S05V",2*i  ,volNam5,-xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S06V",2*i-1,volNam6, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S06V",2*i  ,volNam6,-xvFrame, 0., 0. , 0, "ONLY");
       }       
       // position the panels and the insulating material 
       for (j=0; j<nPCB3[i]; j++){
	 index++;
	 Float_t xx = sensLength * (-nPCB3[i]/2.+j+.5); 
	 
	 Float_t zPanel = spar[2] - panelpar[2]; 
	 gMC->Gspos("S05C",2*index-1,volNam5, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S05C",2*index  ,volNam5, xx, 0.,-zPanel , 0, "ONLY");
	 gMC->Gspos("S06C",2*index-1,volNam6, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S06C",2*index  ,volNam6, xx, 0.,-zPanel , 0, "ONLY");

	 gMC->Gspos("S05I",index,volNam5, xx, 0., 0 , 0, "ONLY");
	 gMC->Gspos("S06I",index,volNam6, xx, 0., 0 , 0, "ONLY");
       } 
     }

     // position the rohacell volume inside the panel volume
     gMC->Gspos("S05R",1,"S05C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06R",1,"S06C",0.,0.,0.,0,"ONLY"); 

     // position the PCB volume inside the insulating material volume
     gMC->Gspos("S05P",1,"S05I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06P",1,"S06I",0.,0.,0.,0,"ONLY"); 
     // position the horizontal frame volume inside the PCB volume
     gMC->Gspos("S05H",1,"S05P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06H",1,"S06P",0.,0.,0.,0,"ONLY"); 
     // position the sensitive volume inside the horizontal frame volume
     gMC->Gsposp("S05G",1,"S05H",0.,0.,0.,0,"ONLY",senspar,3); 
     gMC->Gsposp("S06G",1,"S06H",0.,0.,0.,0,"ONLY",senspar,3); 
     // position the border volumes inside the PCB volume
     Float_t yborder = ( pcbHeight - bFrameHeight ) / 2.; 
     gMC->Gspos("S05B",1,"S05P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S05B",2,"S05P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("S06B",1,"S06P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S06B",2,"S06P",0.,-yborder,0.,0,"ONLY"); 

     // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S05N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S06N","BOX",nulocMaterial,nulocpar,3);
     index = 0;
     for (xx = -xxmax; xx<=xxmax; xx+=3*nulocLength) { 
       index++; 
       gMC->Gspos("S05N",2*index-1,"S05B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S05N",2*index  ,"S05B", xx, 0., bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S06N",2*index-1,"S06B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S06N",2*index  ,"S06B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }
     
     // position the volumes approximating the circular section of the pipe
     Float_t yoffs = sensHeight/2. - yOverlap; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Float_t divpar[3];
     Double_t dydiv= sensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv - yOverlap/2.;
     Int_t imax=0; 
     //     for (Int_t islat=0; islat<nSlats3; islat++) imax += nPCB3[islat]; 
     imax = 1; 
     Float_t rmin = 35.; 
     Float_t z1 = -spar[2], z2=2*spar[2]*1.01; 
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (pcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = sensWidth/2.; 
       Float_t xvol=(pcbLength+xdiv)/2.+1.999;
       Float_t yvol=ydiv + dydiv/2.; 
       gMC->Gsposp("S05G",imax+4*idiv+1,"C05M", xvol, yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S06G",imax+4*idiv+1,"C06M", xvol, yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S05G",imax+4*idiv+2,"C05M", xvol,-yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S06G",imax+4*idiv+2,"C06M", xvol,-yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S05G",imax+4*idiv+3,"C05M",-xvol, yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S06G",imax+4*idiv+3,"C06M",-xvol, yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S05G",imax+4*idiv+4,"C05M",-xvol,-yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S06G",imax+4*idiv+4,"C06M",-xvol,-yvol, z1-z2, 0, "ONLY",divpar,3);
     }
     }
     

 if (stations[3]) {

//********************************************************************
//                            Station 4                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[6];
     iChamber2 =(AliMUONChamber*) (*fChambers)[7];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = zpos2 - zpos1;
//      zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2; // not used any more
     
//
//   Mother volume
     tpar[0] = 37.5-vFrameLength/2.-0.1; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = 3.252;

     gMC->Gsvolu("C07M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C08M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C07M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C08M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
     

     const Int_t nSlats4 = 6;  // number of slats per quadrant
     const Int_t nPCB4[nSlats4] = {4,4,5,5,4,3}; // n PCB per slat
     const Float_t xpos4[nSlats4] = {37.5, 40., 0., 0., 0., 0.};
     Float_t slatLength4[nSlats4];     

     // create and position the slat (mother) volumes 

     char volNam7[5];
     char volNam8[5];
     Float_t xSlat4;
     Float_t ySlat4;

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
       sprintf(volNam7,"S07%d",i);
       gMC->Gsvolu(volNam7,"BOX",slatMaterial,spar,3);
       gMC->Gspos(volNam7, i*4+1,"C07M", xSlat4, ySlat4, zSlat+2.*dzCh4, 0, "ONLY");
       gMC->Gspos(volNam7, i*4+2,"C07M",-xSlat4, ySlat4, zSlat-2.*dzCh4, 0, "ONLY");
       if (i>0) { 
	 gMC->Gspos(volNam7, i*4+3,"C07M", xSlat4,-ySlat4, zSlat+2.*dzCh4, 0, "ONLY");
	 gMC->Gspos(volNam7, i*4+4,"C07M",-xSlat4,-ySlat4, zSlat-2.*dzCh4, 0, "ONLY");
       }
       sprintf(volNam8,"S08%d",i);
       gMC->Gsvolu(volNam8,"BOX",slatMaterial,spar,3);
       gMC->Gspos(volNam8, i*4+1,"C08M", xSlat4, ySlat4, zSlat+2.*dzCh4, 0, "ONLY");
       gMC->Gspos(volNam8, i*4+2,"C08M",-xSlat4, ySlat4, zSlat-2.*dzCh4, 0, "ONLY");
       if (i>0) { 
	 gMC->Gspos(volNam8, i*4+3,"C08M", xSlat4,-ySlat4, zSlat+2.*dzCh4, 0, "ONLY");
	 gMC->Gspos(volNam8, i*4+4,"C08M",-xSlat4,-ySlat4, zSlat-2.*dzCh4, 0, "ONLY");
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

     gMC->Gsvolu("S07G","BOX",sensMaterial,0,0);
     gMC->Gsvolu("S08G","BOX",sensMaterial,0,0);

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
       sprintf(volNam7,"S07%d",i);
       sprintf(volNam8,"S08%d",i);
       Float_t xvFrame  = (slatLength4[i] - vFrameLength)/2.;
       // position the vertical frames 
       if (i!=1) { 
	 gMC->Gspos("S07V",2*i-1,volNam7, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S07V",2*i  ,volNam7,-xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S08V",2*i-1,volNam8, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S08V",2*i  ,volNam8,-xvFrame, 0., 0. , 0, "ONLY");
       }
       // position the panels and the insulating material 
       for (j=0; j<nPCB4[i]; j++){
	 index++;
	 Float_t xx = sensLength * (-nPCB4[i]/2.+j+.5); 

	 Float_t zPanel = spar[2] - panelpar[2]; 
	 gMC->Gspos("S07C",2*index-1,volNam7, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S07C",2*index  ,volNam7, xx, 0.,-zPanel , 0, "ONLY");
	 gMC->Gspos("S08C",2*index-1,volNam8, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S08C",2*index  ,volNam8, xx, 0.,-zPanel , 0, "ONLY");

	 gMC->Gspos("S07I",index,volNam7, xx, 0., 0 , 0, "ONLY");
	 gMC->Gspos("S08I",index,volNam8, xx, 0., 0 , 0, "ONLY");
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

     // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S07N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S08N","BOX",nulocMaterial,nulocpar,3);
     index = 0;
     for (xx = -xxmax; xx<=xxmax; xx+=3*nulocLength) { 
       index++; 
       gMC->Gspos("S07N",2*index-1,"S07B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S07N",2*index  ,"S07B", xx, 0., bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S08N",2*index-1,"S08B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S08N",2*index  ,"S08B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }

     // position the volumes approximating the circular section of the pipe
     Float_t yoffs = sensHeight/2. - yOverlap/2.; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Float_t divpar[3];
     Double_t dydiv= sensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv - yOverlap/2.;
     Int_t imax=0; 
     //     for (Int_t islat=0; islat<nSlats3; islat++) imax += nPCB3[islat]; 
     imax = 1; 
     Float_t rmin = 40.; 
     Float_t z1 = -spar[2], z2=2*spar[2]*1.01; 
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (pcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = sensWidth/2.; 
       Float_t xvol=(pcbLength+xdiv)/2.+1.999;
       Float_t yvol=ydiv + dydiv/2.;
       gMC->Gsposp("S07G",imax+4*idiv+1,"C07M", xvol, yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S08G",imax+4*idiv+1,"C08M", xvol, yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S07G",imax+4*idiv+2,"C07M", xvol,-yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S08G",imax+4*idiv+2,"C08M", xvol,-yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S07G",imax+4*idiv+3,"C07M",-xvol, yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S08G",imax+4*idiv+3,"C08M",-xvol, yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S07G",imax+4*idiv+4,"C07M",-xvol,-yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S08G",imax+4*idiv+4,"C08M",-xvol,-yvol, z1-z2, 0, "ONLY",divpar,3);
     }

 }

 if (stations[4]) {
     

//********************************************************************
//                            Station 5                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[8];
     iChamber2 =(AliMUONChamber*) (*fChambers)[9];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = zpos2 - zpos1;
//      zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2; // not used any more
     
//
//   Mother volume
     tpar[0] = 37.5-vFrameLength/2.-0.1; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5.;

     gMC->Gsvolu("C09M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C10M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C09M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C10M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");


     const Int_t nSlats5 = 7;  // number of slats per quadrant
     const Int_t nPCB5[nSlats5] = {5,5,6,6,5,4,3}; // n PCB per slat
     const Float_t xpos5[nSlats5] = {37.5, 40., 0., 0., 0., 0., 0.};
     Float_t slatLength5[nSlats5]; 
     char volNam9[5];
     char volNam10[5];
     Float_t xSlat5;
     Float_t ySlat5;

     for (i = 0; i<nSlats5; i++){
       slatLength5[i] = pcbLength * nPCB5[i] + 2. * dSlatLength; 
       xSlat5 = slatLength5[i]/2. - vFrameLength/2. +xpos5[i]; 
       if (i==1) slatLength5[i] -=  2. *dSlatLength; // frame out in PCB with circular border 
       ySlat5 = sensHeight * i - yOverlap * i; 
       spar[0] = slatLength5[i]/2.; 
       spar[1] = slatHeight/2.;
       spar[2] = slatWidth/2. * 1.01; 
       Float_t dzCh5=spar[2]*1.01;
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? -spar[2] : spar[2]; 
       sprintf(volNam9,"S09%d",i);
       gMC->Gsvolu(volNam9,"BOX",slatMaterial,spar,3);
       gMC->Gspos(volNam9, i*4+1,"C09M", xSlat5, ySlat5, zSlat+2.*dzCh5, 0, "ONLY");
       gMC->Gspos(volNam9, i*4+2,"C09M",-xSlat5, ySlat5, zSlat-2.*dzCh5, 0, "ONLY");
       if (i>0) { 
	   gMC->Gspos(volNam9, i*4+3,"C09M", xSlat5,-ySlat5, zSlat+2.*dzCh5, 0, "ONLY");
	   gMC->Gspos(volNam9, i*4+4,"C09M",-xSlat5,-ySlat5, zSlat-2.*dzCh5, 0, "ONLY");
       }
       sprintf(volNam10,"S10%d",i);
       gMC->Gsvolu(volNam10,"BOX",slatMaterial,spar,3);
       gMC->Gspos(volNam10, i*4+1,"C10M", xSlat5, ySlat5, zSlat+2.*dzCh5, 0, "ONLY");
       gMC->Gspos(volNam10, i*4+2,"C10M",-xSlat5, ySlat5, zSlat-2.*dzCh5, 0, "ONLY");
       if (i>0) { 
	   gMC->Gspos(volNam10, i*4+3,"C10M", xSlat5,-ySlat5, zSlat+2.*dzCh5, 0, "ONLY");
	   gMC->Gspos(volNam10, i*4+4,"C10M",-xSlat5,-ySlat5, zSlat-2.*dzCh5, 0, "ONLY");
       }
     }

     // create the panel volume 
 
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

     gMC->Gsvolu("S09G","BOX",sensMaterial,0,0);
     gMC->Gsvolu("S10G","BOX",sensMaterial,0,0);

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
       sprintf(volNam9,"S09%d",i);
       sprintf(volNam10,"S10%d",i);
       Float_t xvFrame  = (slatLength5[i] - vFrameLength)/2.;
       // position the vertical frames 
       if (i!=1) { 
	 gMC->Gspos("S09V",2*i-1,volNam9, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S09V",2*i  ,volNam9,-xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S10V",2*i-1,volNam10, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S10V",2*i  ,volNam10,-xvFrame, 0., 0. , 0, "ONLY");
       }
       
       // position the panels and the insulating material 
       for (j=0; j<nPCB5[i]; j++){
	 index++;
	 Float_t xx = sensLength * (-nPCB5[i]/2.+j+.5); 

	 Float_t zPanel = spar[2] - panelpar[2]; 
	 gMC->Gspos("S09C",2*index-1,volNam9, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S09C",2*index  ,volNam9, xx, 0.,-zPanel , 0, "ONLY");
	 gMC->Gspos("S10C",2*index-1,volNam10, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S10C",2*index  ,volNam10, xx, 0.,-zPanel , 0, "ONLY");

	 gMC->Gspos("S09I",index,volNam9, xx, 0., 0 , 0, "ONLY");
	 gMC->Gspos("S10I",index,volNam10, xx, 0., 0 , 0, "ONLY");
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

     // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S09N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S10N","BOX",nulocMaterial,nulocpar,3);
     index = 0;
     for (xx = -xxmax; xx<=xxmax; xx+=3*nulocLength) { 
       index++; 
       gMC->Gspos("S09N",2*index-1,"S09B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S09N",2*index  ,"S09B", xx, 0., bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S10N",2*index-1,"S10B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S10N",2*index  ,"S10B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }
     // position the volumes approximating the circular section of the pipe
     Float_t yoffs = sensHeight/2. - yOverlap/2.; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Float_t divpar[3];
     Double_t dydiv= sensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv - yOverlap/2.;
     Int_t imax=0; 
     //     for (Int_t islat=0; islat<nSlats3; islat++) imax += nPCB3[islat]; 
     imax = 1; 
     Float_t rmin = 40.; 
     Float_t z1 = spar[2], z2=2*spar[2]*1.01; 
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (pcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = sensWidth/2.; 
       Float_t xvol=(pcbLength+xdiv)/2. + 1.999;
       Float_t yvol=ydiv + dydiv/2.;
       gMC->Gsposp("S09G",imax+4*idiv+1,"C09M", xvol, yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S10G",imax+4*idiv+1,"C10M", xvol, yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S09G",imax+4*idiv+2,"C09M", xvol,-yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S10G",imax+4*idiv+2,"C10M", xvol,-yvol, z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S09G",imax+4*idiv+3,"C09M",-xvol, yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S10G",imax+4*idiv+3,"C10M",-xvol, yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S09G",imax+4*idiv+4,"C09M",-xvol,-yvol, z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S10G",imax+4*idiv+4,"C10M",-xvol,-yvol, z1-z2, 0, "ONLY",divpar,3);
     }

 }
 

///////////////////////////////////////
// GEOMETRY FOR THE TRIGGER CHAMBERS //
///////////////////////////////////////

// 03/00 P. Dupieux : introduce a slighly more realistic  
//                    geom. of the trigger readout planes with
//                    2 Zpos per trigger plane (alternate
//                    between left and right of the trigger)  

//  Parameters of the Trigger Chambers

    		
     const Float_t kXMC1MIN=34.;       
     const Float_t kXMC1MED=51.;                                
     const Float_t kXMC1MAX=272.;                               
     const Float_t kYMC1MIN=34.;                              
     const Float_t kYMC1MAX=51.;                              
     const Float_t kRMIN1=50.;
     const Float_t kRMAX1=62.;
     const Float_t kRMIN2=50.;
     const Float_t kRMAX2=66.;

//   zposition of the middle of the gas gap in mother vol 
     const Float_t kZMCm=-3.6;
     const Float_t kZMCp=+3.6;


// TRIGGER STATION 1 - TRIGGER STATION 1 - TRIGGER STATION 1

     // iChamber 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[10];
     iChamber2 =(AliMUONChamber*) (*fChambers)[11]; 

     // 03/00 
     // zpos1 and zpos2 are now the middle of the first and second
     // plane of station 1 : 
     // zpos1=(16075+15995)/2=16035 mm, thick/2=40 mm
     // zpos2=(16225+16145)/2=16185 mm, thick/2=40 mm
     //
     // zpos1m=15999 mm , zpos1p=16071 mm (middles of gas gaps)
     // zpos2m=16149 mm , zpos2p=16221 mm (middles of gas gaps)
     // rem : the total thickness accounts for 1 mm of al on both 
     // side of the RPCs (see zpos1 and zpos2), as previously

     zpos1=iChamber1->Z();
     zpos2=iChamber2->Z();


// Mother volume definition     
     tpar[0] = iChamber->RInner(); 
     tpar[1] = iChamber->ROuter();
     tpar[2] = 4.0;    
     gMC->Gsvolu("CM11", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("CM12", "TUBE", idAir, tpar, 3);
     
// Definition of the flange between the beam shielding and the RPC 
     tpar[0]= kRMIN1;
     tpar[1]= kRMAX1;
     tpar[2]= 4.0;
   
     gMC->Gsvolu("CF1A", "TUBE", idAlu1, tpar, 3);     //Al
     gMC->Gspos("CF1A", 1, "CM11", 0., 0., 0., 0, "MANY");
     gMC->Gspos("CF1A", 2, "CM12", 0., 0., 0., 0, "MANY");


// FIRST PLANE OF STATION 1

//   ratios of zpos1m/zpos1p and inverse for first plane
     Float_t zmp=(zpos1-3.6)/(zpos1+3.6);
     Float_t zpm=1./zmp;
   

// Definition of prototype for chambers in the first plane     
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC1A", "BOX ", idAlu1, tpar, 0);           //Al    
     gMC->Gsvolu("CB1A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG1A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t kXMC1A=kXMC1MED+(kXMC1MAX-kXMC1MED)/2.;
     const Float_t kYMC1Am=0.;
     const Float_t kYMC1Ap=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG1A", 1, "CB1A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB1A", 1, "CC1A", 0., 0., 0., 0, "ONLY",tpar,3);

     tpar[2] = 0.4;
     tpar[0] = (kXMC1MAX-kXMC1MED)/2.;
     tpar[1] = kYMC1MIN;

     gMC->Gsposp("CC1A", 1, "CM11",kXMC1A,kYMC1Am,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 2, "CM11",-kXMC1A,kYMC1Ap,kZMCp, 0, "ONLY", tpar, 3);
     
//  chamber type B    
     Float_t tpar1save=tpar[1];
     Float_t y1msave=kYMC1Am;
     Float_t y1psave=kYMC1Ap;
 
     tpar[0] = (kXMC1MAX-kXMC1MIN)/2.;
     tpar[1] = (kYMC1MAX-kYMC1MIN)/2.;
     
     const Float_t kXMC1B=kXMC1MIN+tpar[0];
     const Float_t kYMC1Bp=(y1msave+tpar1save)*zpm+tpar[1];
     const Float_t kYMC1Bm=(y1psave+tpar1save)*zmp+tpar[1];

     gMC->Gsposp("CC1A", 3, "CM11",kXMC1B,kYMC1Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 4, "CM11",-kXMC1B,kYMC1Bm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 5, "CM11",kXMC1B,-kYMC1Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 6, "CM11",-kXMC1B,-kYMC1Bm,kZMCm, 0, "ONLY", tpar, 3);
     
//  chamber type C  (end of type B !!)      
     tpar1save=tpar[1];
     y1msave=kYMC1Bm;
     y1psave=kYMC1Bp;

     tpar[0] = kXMC1MAX/2;
     tpar[1] = kYMC1MAX/2;
     
     const Float_t kXMC1C=tpar[0];
// warning : same Z than type B
     const Float_t kYMC1Cp=(y1psave+tpar1save)*1.+tpar[1];
     const Float_t kYMC1Cm=(y1msave+tpar1save)*1.+tpar[1];
     
     gMC->Gsposp("CC1A", 7, "CM11",kXMC1C,kYMC1Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 8, "CM11",-kXMC1C,kYMC1Cm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 9, "CM11",kXMC1C,-kYMC1Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 10, "CM11",-kXMC1C,-kYMC1Cm,kZMCm, 0, "ONLY", tpar, 3);
     
//  chamber type D, E and F (same size)        
     tpar1save=tpar[1];
     y1msave=kYMC1Cm;
     y1psave=kYMC1Cp;

     tpar[0] = kXMC1MAX/2.;
     tpar[1] = kYMC1MIN;
     
     const Float_t kXMC1D=tpar[0];
     const Float_t kYMC1Dp=(y1msave+tpar1save)*zpm+tpar[1];
     const Float_t kYMC1Dm=(y1psave+tpar1save)*zmp+tpar[1];
     
     gMC->Gsposp("CC1A", 11, "CM11",kXMC1D,kYMC1Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 12, "CM11",-kXMC1D,kYMC1Dp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 13, "CM11",kXMC1D,-kYMC1Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 14, "CM11",-kXMC1D,-kYMC1Dp,kZMCp, 0, "ONLY", tpar, 3);


     tpar1save=tpar[1];
     y1msave=kYMC1Dm;
     y1psave=kYMC1Dp;
     const Float_t kYMC1Ep=(y1msave+tpar1save)*zpm+tpar[1];
     const Float_t kYMC1Em=(y1psave+tpar1save)*zmp+tpar[1];
     
     gMC->Gsposp("CC1A", 15, "CM11",kXMC1D,kYMC1Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 16, "CM11",-kXMC1D,kYMC1Em,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 17, "CM11",kXMC1D,-kYMC1Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 18, "CM11",-kXMC1D,-kYMC1Em,kZMCm, 0, "ONLY", tpar, 3);

     tpar1save=tpar[1];
     y1msave=kYMC1Em;
     y1psave=kYMC1Ep;
     const Float_t kYMC1Fp=(y1msave+tpar1save)*zpm+tpar[1];
     const Float_t kYMC1Fm=(y1psave+tpar1save)*zmp+tpar[1];
    
     gMC->Gsposp("CC1A", 19, "CM11",kXMC1D,kYMC1Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 20, "CM11",-kXMC1D,kYMC1Fp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 21, "CM11",kXMC1D,-kYMC1Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 22, "CM11",-kXMC1D,-kYMC1Fp,kZMCp, 0, "ONLY", tpar, 3);

// Positioning first plane in ALICE     
     gMC->Gspos("CM11", 1, "ALIC", 0., 0., zpos1, 0, "ONLY");

// End of geometry definition for the first plane of station 1



// SECOND PLANE OF STATION 1 : proj ratio = zpos2/zpos1

     const Float_t kZ12=zpos2/zpos1;
      
// Definition of prototype for chambers in the second plane of station 1    
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC2A", "BOX ", idAlu1, tpar, 0);           //Al    
     gMC->Gsvolu("CB2A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG2A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t kXMC2A=kXMC1A*kZ12;
     const Float_t kYMC2Am=0.;
     const Float_t kYMC2Ap=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG2A", 1, "CB2A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB2A", 1, "CC2A", 0., 0., 0., 0, "ONLY",tpar,3);

     tpar[2] = 0.4;
     tpar[0] = ((kXMC1MAX-kXMC1MED)/2.)*kZ12;
     tpar[1] = kYMC1MIN*kZ12;

     gMC->Gsposp("CC2A", 1, "CM12",kXMC2A,kYMC2Am,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 2, "CM12",-kXMC2A,kYMC2Ap,kZMCp, 0, "ONLY", tpar, 3);
     

//  chamber type B    

     tpar[0] = ((kXMC1MAX-kXMC1MIN)/2.)*kZ12;
     tpar[1] = ((kYMC1MAX-kYMC1MIN)/2.)*kZ12;
     
     const Float_t kXMC2B=kXMC1B*kZ12;
     const Float_t kYMC2Bp=kYMC1Bp*kZ12;
     const Float_t kYMC2Bm=kYMC1Bm*kZ12;
     gMC->Gsposp("CC2A", 3, "CM12",kXMC2B,kYMC2Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 4, "CM12",-kXMC2B,kYMC2Bm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 5, "CM12",kXMC2B,-kYMC2Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 6, "CM12",-kXMC2B,-kYMC2Bm,kZMCm, 0, "ONLY", tpar, 3);

     
//  chamber type C   (end of type B !!)     

     tpar[0] = (kXMC1MAX/2)*kZ12;
     tpar[1] = (kYMC1MAX/2)*kZ12;
     
     const Float_t kXMC2C=kXMC1C*kZ12;
     const Float_t kYMC2Cp=kYMC1Cp*kZ12;
     const Float_t kYMC2Cm=kYMC1Cm*kZ12;     
     gMC->Gsposp("CC2A", 7, "CM12",kXMC2C,kYMC2Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 8, "CM12",-kXMC2C,kYMC2Cm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 9, "CM12",kXMC2C,-kYMC2Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 10, "CM12",-kXMC2C,-kYMC2Cm,kZMCm, 0, "ONLY", tpar, 3);
     
//  chamber type D, E and F (same size)        

     tpar[0] = (kXMC1MAX/2.)*kZ12;
     tpar[1] = kYMC1MIN*kZ12;
     
     const Float_t kXMC2D=kXMC1D*kZ12;
     const Float_t kYMC2Dp=kYMC1Dp*kZ12;
     const Float_t kYMC2Dm=kYMC1Dm*kZ12;     
     gMC->Gsposp("CC2A", 11, "CM12",kXMC2D,kYMC2Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 12, "CM12",-kXMC2D,kYMC2Dp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 13, "CM12",kXMC2D,-kYMC2Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 14, "CM12",-kXMC2D,-kYMC2Dp,kZMCp, 0, "ONLY", tpar, 3);

     const Float_t kYMC2Ep=kYMC1Ep*kZ12;
     const Float_t kYMC2Em=kYMC1Em*kZ12;
     gMC->Gsposp("CC2A", 15, "CM12",kXMC2D,kYMC2Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 16, "CM12",-kXMC2D,kYMC2Em,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 17, "CM12",kXMC2D,-kYMC2Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 18, "CM12",-kXMC2D,-kYMC2Em,kZMCm, 0, "ONLY", tpar, 3);


     const Float_t kYMC2Fp=kYMC1Fp*kZ12;
     const Float_t kYMC2Fm=kYMC1Fm*kZ12;
     gMC->Gsposp("CC2A", 19, "CM12",kXMC2D,kYMC2Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 20, "CM12",-kXMC2D,kYMC2Fp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 21, "CM12",kXMC2D,-kYMC2Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 22, "CM12",-kXMC2D,-kYMC2Fp,kZMCp, 0, "ONLY", tpar, 3);

// Positioning second plane of station 1 in ALICE     
     
     gMC->Gspos("CM12", 1, "ALIC", 0., 0., zpos2, 0, "ONLY");

// End of geometry definition for the second plane of station 1



// TRIGGER STATION 2 - TRIGGER STATION 2 - TRIGGER STATION 2    

     // 03/00 
     // zpos3 and zpos4 are now the middle of the first and second
     // plane of station 2 : 
     // zpos3=(17075+16995)/2=17035 mm, thick/2=40 mm
     // zpos4=(17225+17145)/2=17185 mm, thick/2=40 mm
     //
     // zpos3m=16999 mm , zpos3p=17071 mm (middles of gas gaps)
     // zpos4m=17149 mm , zpos4p=17221 mm (middles of gas gaps)
     // rem : the total thickness accounts for 1 mm of al on both 
     // side of the RPCs (see zpos3 and zpos4), as previously
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[12];
     iChamber2 =(AliMUONChamber*) (*fChambers)[13];
     Float_t zpos3=iChamber1->Z();
     Float_t zpos4=iChamber2->Z();


// Mother volume definition     
     tpar[0] = iChamber->RInner(); 
     tpar[1] = iChamber->ROuter();
     tpar[2] = 4.0;    
 
     gMC->Gsvolu("CM21", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("CM22", "TUBE", idAir, tpar, 3);
     
// Definition of the flange between the beam shielding and the RPC 
//  ???? interface shielding

     tpar[0]= kRMIN2;
     tpar[1]= kRMAX2;
     tpar[2]= 4.0;
   
     gMC->Gsvolu("CF2A", "TUBE", idAlu1, tpar, 3);            //Al
     gMC->Gspos("CF2A", 1, "CM21", 0., 0., 0., 0, "MANY");
     gMC->Gspos("CF2A", 2, "CM22", 0., 0., 0., 0, "MANY");
    


// FIRST PLANE OF STATION 2 : proj ratio = zpos3/zpos1

     const Float_t kZ13=zpos3/zpos1; 

// Definition of prototype for chambers in the first plane of station 2       
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC3A", "BOX ", idAlu1, tpar, 0);           //Al  
     gMC->Gsvolu("CB3A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG3A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer


// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t kXMC3A=kXMC1A*kZ13;
     const Float_t kYMC3Am=0.;
     const Float_t kYMC3Ap=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG3A", 1, "CB3A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB3A", 1, "CC3A", 0., 0., 0., 0, "ONLY",tpar,3);

     tpar[2] = 0.4;
     tpar[0] = ((kXMC1MAX-kXMC1MED)/2.)*kZ13;
     tpar[1] = kYMC1MIN*kZ13;
     gMC->Gsposp("CC3A", 1, "CM21",kXMC3A,kYMC3Am,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 2, "CM21",-kXMC3A,kYMC3Ap,kZMCp, 0, "ONLY", tpar, 3);

     
//  chamber type B    
     tpar[0] = ((kXMC1MAX-kXMC1MIN)/2.)*kZ13;
     tpar[1] = ((kYMC1MAX-kYMC1MIN)/2.)*kZ13;
     
     const Float_t kXMC3B=kXMC1B*kZ13;
     const Float_t kYMC3Bp=kYMC1Bp*kZ13;
     const Float_t kYMC3Bm=kYMC1Bm*kZ13;
     gMC->Gsposp("CC3A", 3, "CM21",kXMC3B,kYMC3Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 4, "CM21",-kXMC3B,kYMC3Bm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 5, "CM21",kXMC3B,-kYMC3Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 6, "CM21",-kXMC3B,-kYMC3Bm,kZMCm, 0, "ONLY", tpar, 3);

     
//  chamber type C  (end of type B !!)      
     tpar[0] = (kXMC1MAX/2)*kZ13;
     tpar[1] = (kYMC1MAX/2)*kZ13;
     
     const Float_t kXMC3C=kXMC1C*kZ13;
     const Float_t kYMC3Cp=kYMC1Cp*kZ13;
     const Float_t kYMC3Cm=kYMC1Cm*kZ13;     
     gMC->Gsposp("CC3A", 7, "CM21",kXMC3C,kYMC3Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 8, "CM21",-kXMC3C,kYMC3Cm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 9, "CM21",kXMC3C,-kYMC3Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 10, "CM21",-kXMC3C,-kYMC3Cm,kZMCm, 0, "ONLY", tpar, 3);
     

//  chamber type D, E and F (same size)         

     tpar[0] = (kXMC1MAX/2.)*kZ13;
     tpar[1] = kYMC1MIN*kZ13;
     
     const Float_t kXMC3D=kXMC1D*kZ13;
     const Float_t kYMC3Dp=kYMC1Dp*kZ13;
     const Float_t kYMC3Dm=kYMC1Dm*kZ13;          
     gMC->Gsposp("CC3A", 11, "CM21",kXMC3D,kYMC3Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 12, "CM21",-kXMC3D,kYMC3Dp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 13, "CM21",kXMC3D,-kYMC3Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 14, "CM21",-kXMC3D,-kYMC3Dp,kZMCp, 0, "ONLY", tpar, 3);

     const Float_t kYMC3Ep=kYMC1Ep*kZ13;
     const Float_t kYMC3Em=kYMC1Em*kZ13;
     gMC->Gsposp("CC3A", 15, "CM21",kXMC3D,kYMC3Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 16, "CM21",-kXMC3D,kYMC3Em,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 17, "CM21",kXMC3D,-kYMC3Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 18, "CM21",-kXMC3D,-kYMC3Em,kZMCm, 0, "ONLY", tpar, 3);

     const Float_t kYMC3Fp=kYMC1Fp*kZ13;
     const Float_t kYMC3Fm=kYMC1Fm*kZ13;
     gMC->Gsposp("CC3A", 19, "CM21",kXMC3D,kYMC3Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 20, "CM21",-kXMC3D,kYMC3Fp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 21, "CM21",kXMC3D,-kYMC3Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 22, "CM21",-kXMC3D,-kYMC3Fp,kZMCp, 0, "ONLY", tpar, 3);
       

// Positioning first plane of station 2 in ALICE
     
     gMC->Gspos("CM21", 1, "ALIC", 0., 0., zpos3, 0, "ONLY");

// End of geometry definition for the first plane of station 2




// SECOND PLANE OF STATION 2 : proj ratio = zpos4/zpos1

     const Float_t kZ14=zpos4/zpos1;
     
// Definition of prototype for chambers in the second plane of station 2    
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC4A", "BOX ", idAlu1, tpar, 0);           //Al      
     gMC->Gsvolu("CB4A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG4A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t kXMC4A=kXMC1A*kZ14;
     const Float_t kYMC4Am=0.;
     const Float_t kYMC4Ap=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG4A", 1, "CB4A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB4A", 1, "CC4A", 0., 0., 0., 0, "ONLY",tpar,3);

     tpar[2] = 0.4;
     tpar[0] = ((kXMC1MAX-kXMC1MED)/2.)*kZ14;
     tpar[1] = kYMC1MIN*kZ14;
     gMC->Gsposp("CC4A", 1, "CM22",kXMC4A,kYMC4Am,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 2, "CM22",-kXMC4A,kYMC4Ap,kZMCp, 0, "ONLY", tpar, 3);
     

//  chamber type B    
     tpar[0] = ((kXMC1MAX-kXMC1MIN)/2.)*kZ14;
     tpar[1] = ((kYMC1MAX-kYMC1MIN)/2.)*kZ14;
     
     const Float_t kXMC4B=kXMC1B*kZ14;
     const Float_t kYMC4Bp=kYMC1Bp*kZ14;
     const Float_t kYMC4Bm=kYMC1Bm*kZ14;
     gMC->Gsposp("CC4A", 3, "CM22",kXMC4B,kYMC4Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 4, "CM22",-kXMC4B,kYMC4Bm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 5, "CM22",kXMC4B,-kYMC4Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 6, "CM22",-kXMC4B,-kYMC4Bm,kZMCm, 0, "ONLY", tpar, 3);

     
//  chamber type C   (end of type B !!)      
     tpar[0] =(kXMC1MAX/2)*kZ14;
     tpar[1] =  (kYMC1MAX/2)*kZ14;
     
     const Float_t kXMC4C=kXMC1C*kZ14;
     const Float_t kYMC4Cp=kYMC1Cp*kZ14;
     const Float_t kYMC4Cm=kYMC1Cm*kZ14;     
     gMC->Gsposp("CC4A", 7, "CM22",kXMC4C,kYMC4Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 8, "CM22",-kXMC4C,kYMC4Cm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 9, "CM22",kXMC4C,-kYMC4Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 10, "CM22",-kXMC4C,-kYMC4Cm,kZMCm, 0, "ONLY", tpar, 3);

     
//  chamber type D, E and F (same size)      
     tpar[0] = (kXMC1MAX/2.)*kZ14;
     tpar[1] =  kYMC1MIN*kZ14;
     
     const Float_t kXMC4D=kXMC1D*kZ14;
     const Float_t kYMC4Dp=kYMC1Dp*kZ14;
     const Float_t kYMC4Dm=kYMC1Dm*kZ14;          
     gMC->Gsposp("CC4A", 11, "CM22",kXMC4D,kYMC4Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 12, "CM22",-kXMC4D,kYMC4Dp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 13, "CM22",kXMC4D,-kYMC4Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 14, "CM22",-kXMC4D,-kYMC4Dp,kZMCp, 0, "ONLY", tpar, 3);

     const Float_t kYMC4Ep=kYMC1Ep*kZ14;
     const Float_t kYMC4Em=kYMC1Em*kZ14;          
     gMC->Gsposp("CC4A", 15, "CM22",kXMC4D,kYMC4Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 16, "CM22",-kXMC4D,kYMC4Em,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 17, "CM22",kXMC4D,-kYMC4Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 18, "CM22",-kXMC4D,-kYMC4Em,kZMCm, 0, "ONLY", tpar, 3);

     const Float_t kYMC4Fp=kYMC1Fp*kZ14;
     const Float_t kYMC4Fm=kYMC1Fm*kZ14;          
     gMC->Gsposp("CC4A", 19, "CM22",kXMC4D,kYMC4Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 20, "CM22",-kXMC4D,kYMC4Fp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 21, "CM22",kXMC4D,-kYMC4Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 22, "CM22",-kXMC4D,-kYMC4Fp,kZMCp, 0, "ONLY", tpar, 3);
     

// Positioning second plane of station 2 in ALICE
     
     gMC->Gspos("CM22", 1, "ALIC", 0., 0., zpos4, 0, "ONLY");

// End of geometry definition for the second plane of station 2

// End of trigger geometry definition

}
