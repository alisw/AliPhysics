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
Revision 1.13  2000/10/06 10:03:38  morsch
Call to gMC->VolId() moved to Init()

Revision 1.12  2000/10/02 21:28:09  fca
Removal of useless dependecies via forward declarations

Revision 1.11  2000/06/27 07:31:07  morsch
fChambers = 0; deleted from constructor.

Revision 1.10  2000/06/26 14:02:38  morsch
Add class AliMUONConstants with MUON specific constants using static memeber data and access methods.

Revision 1.9  2000/06/15 07:58:49  morsch
Code from MUON-dev joined

Revision 1.8.4.9  2000/06/12 19:20:49  morsch
Constructor sets default geometry, segmentation and response parameters.

Revision 1.8.4.8  2000/06/09 21:55:28  morsch
Most coding rule violations corrected.

Revision 1.8.4.7  2000/05/02 13:15:18  morsch
Coding rule violations RS3, RN13 corected

Revision 1.8.4.6  2000/05/02 10:24:26  morsch
Public access to fdGas and fdAlu of AliMUONChamber replaced by getters.

Revision 1.8.4.5  2000/04/26 19:58:47  morsch
Obsolete reference to trig_ removed.

Revision 1.8.4.4  2000/04/19 19:42:47  morsch
change NCH to kNCH

Revision 1.8.4.3  2000/02/17 08:17:43  morsch
Gammas and neutrons are also scored in the stepmanager
*/

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 
#include <TLorentzVector.h> 
#include <iostream.h>

#include "AliMUONv0.h"
#include "AliMUONChamber.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliCallf77.h"
#include "AliConst.h" 
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

ClassImp(AliMUONv0)
 
//___________________________________________
AliMUONv0::AliMUONv0() : AliMUON()
{
// Constructor
    fChambers = 0;
}
 
//___________________________________________
AliMUONv0::AliMUONv0(const char *name, const char *title)
       : AliMUON(name,title)
{
// Constructor
    SetIshunt(0);
    SetMaxStepGas(0.1);
    SetMaxStepAlu(0.1);
    
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
    Int_t station;
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
    AliMUONSegmentationV01 *seg11=new AliMUONSegmentationV01;
    
    seg11->SetSegRadii(rseg1);
//  seg11->SetPadSize(3, 0.5);
    seg11->SetPadSize(2.4, 0.4); // smaller pad size
//  seg11->SetDAnod(3.0/3./4);
    seg11->SetDAnod(0.20); // smaller distance between anode wires
    seg11->SetPadDivision(nseg1);
    
    SetSegmentationModel(chamber-1, 1, seg11);
//
    AliMUONSegmentationV02 *seg12=new AliMUONSegmentationV02;
    seg12->SetSegRadii(rseg1); 
//  seg12->SetPadSize(0.75, 2.0);
    seg12->SetPadSize(0.6, 1.6); // smaller pad size
//  seg12->SetDAnod(3.0/3./4);
    seg12->SetDAnod(0.20); // smaller distance between anode wires
    seg12->SetPadDivision(nseg1);
    
    SetSegmentationModel(chamber-1, 2, seg12);
    
//  SetResponseModel(chamber-1, response0);	    
    SetResponseModel(chamber-1, responseSt1); // special response	    
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
    chamber=2;
//^^^^^^^^^
//
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationV01 *seg21=new AliMUONSegmentationV01;
    seg21->SetSegRadii(rseg1);
//  seg21->SetPadSize(3, 0.5);
    seg21->SetPadSize(2.4, 0.4); // smaller pad size
//  seg21->SetDAnod(3.0/3./4);
    seg21->SetDAnod(0.20); // smaller distance between anode wires
    seg21->SetPadDivision(nseg1);
    SetSegmentationModel(chamber-1, 1, seg21);
//
    AliMUONSegmentationV02 *seg22=new AliMUONSegmentationV02;
    seg22->SetSegRadii(rseg1); 
//  seg22->SetPadSize(0.75, 2.);
    seg22->SetPadSize(0.6, 1.6); // smaller pad size
//  seg22->SetDAnod(3.0/3./4);
    seg22->SetDAnod(0.20); // smaller distance between anode wires
    seg22->SetPadDivision(nseg1);
    SetSegmentationModel(chamber-1, 2, seg22);
    
//  SetResponseModel(chamber-1, response0);	    
    SetResponseModel(chamber-1, responseSt1); // special response
    Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
    
//
//--------------------------------------------------------
// Configuration for Chamber TC3/4 (Station 2) -----------
///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// Float_t rseg2[4]={23.5, 87.7, 122.4, 122.5};
    Float_t rseg2[4]={23.5, 47.1, 87.7, 122.5};
    Int_t   nseg2[4]={4, 4, 2, 1};
//
    chamber=3;
//^^^^^^^^^
    SetNsec(chamber-1,2);
//
    AliMUONSegmentationV01 *seg31=new AliMUONSegmentationV01;
    seg31->SetSegRadii(rseg2);
    seg31->SetPadSize(3.0, 0.5);
    seg31->SetDAnod(3.0/3./4);
    seg31->SetPadDivision(nseg2);
    SetSegmentationModel(chamber-1, 1, seg31);
//
    AliMUONSegmentationV02 *seg32=new AliMUONSegmentationV02;
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
    AliMUONSegmentationV01 *seg41=new AliMUONSegmentationV01;
    seg41->SetSegRadii(rseg2);
    seg41->SetPadSize(3.0, 0.5);
    seg41->SetDAnod(3.0/3./4);
    seg41->SetPadDivision(nseg2);
    SetSegmentationModel(chamber-1, 1, seg41);
//
    AliMUONSegmentationV02 *seg42=new AliMUONSegmentationV02;
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
    Int_t   nseg3[4]={4, 4, 2, 1};
    Int_t   npcb5[36] = {0,0,2,0,
			 0,0,3,0,
			 0,1,3,0,
			 0,2,2,0,
			 0,1,2,0, 
			 0,2,2,0, 
			 0,1,3,0, 
			 0,0,3,0,
			 0,0,2,0};

    Float_t shift = 1.5/2.;
    // Float_t xpos5[8]    = {2., 2., 2., 42., 42., 2., 2., 2.};
    Float_t xpos5[9]    = {2., 2., 2., 2.,32., 2., 2., 2., 2.};
    Float_t ypos5       = -(20.+4.*(40.-2.*shift));
    
    chamber=5;
    SetNsec(chamber-1,2);
    AliMUONSegmentationSlat *seg51=new AliMUONSegmentationSlat;
    seg51->SetNSlats(9); 
    seg51->SetShift(shift);  
    seg51->SetNPCBperSector(npcb5); 
    seg51->SetSlatXPositions(xpos5);
    seg51->SetSlatYPosition(ypos5);
    seg51->SetPadSize(10.,0.5);
    seg51->SetDAnod(0.25);
    seg51->SetPadDivision(nseg3);
    SetSegmentationModel(chamber-1, 1, seg51);
    
    AliMUONSegmentationSlatN *seg52=new AliMUONSegmentationSlatN;
    seg52->SetNSlats(9); 
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
    AliMUONSegmentationSlat *seg61=new AliMUONSegmentationSlat;
    seg61->SetNSlats(9); 
    seg61->SetShift(shift);  
    seg61->SetNPCBperSector(npcb5); 
    seg61->SetSlatXPositions(xpos5);
    seg61->SetSlatYPosition(ypos5);
    seg61->SetPadSize(10.,0.5);
    seg61->SetDAnod(0.25);
    seg61->SetPadDivision(nseg3);
    SetSegmentationModel(chamber-1, 1, seg61);
    
    AliMUONSegmentationSlatN *seg62=new AliMUONSegmentationSlatN;
    seg62->SetNSlats(9); 
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
    AliMUONSegmentationSlat *seg71=new AliMUONSegmentationSlat;
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
    
    AliMUONSegmentationSlatN *seg72=new AliMUONSegmentationSlatN;
    
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
    AliMUONSegmentationSlat *seg81=new AliMUONSegmentationSlat;
    
    seg81->SetNSlats(11);  
    seg81->SetShift(shift);  
    seg81->SetNPCBperSector(npcb7); 
    seg81->SetSlatXPositions(xpos7);
    seg81->SetSlatYPosition(ypos7);
    seg81->SetPadSize(10.,0.5);
    seg81->SetDAnod(0.25);
    seg81->SetPadDivision(nseg4);
    SetSegmentationModel(chamber-1, 1, seg81);
    
    AliMUONSegmentationSlat *seg82=new AliMUONSegmentationSlatN;

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
    AliMUONSegmentationSlat *seg91=new AliMUONSegmentationSlat;
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
    
    // Float_t xpos9[13]   = {2., 2., 2., 2., 2., 2., 39.5 , 2., 2., 2., 2., 2., 2.};
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
    
    AliMUONSegmentationSlatN *seg92=new AliMUONSegmentationSlatN;
    
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
    AliMUONSegmentationSlat *seg101=new AliMUONSegmentationSlat;
    
    seg101->SetNSlats(13);  
    seg101->SetShift(shift);  
    seg101->SetNPCBperSector(npcb9); 
    seg101->SetSlatXPositions(xpos9);
    seg101->SetSlatYPosition(ypos9);
    seg101->SetPadSize(10.,0.5);
    seg101->SetDAnod(0.25);
    seg101->SetPadDivision(nseg4);
    SetSegmentationModel(chamber-1, 1, seg101);
    
    AliMUONSegmentationSlatN *seg102=new AliMUONSegmentationSlatN;
    
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

void AliMUONv0::CreateGeometry()
{
// Creates coarse geometry for hit density simulations
    Int_t *idtmed = fIdtmed->GetArray()-1099;
//
     Float_t zpos, dAlu, tpar[3];
     Int_t idAir=idtmed[1100];
     Int_t idAlu=idtmed[1103];     

     AliMUONChamber *iChamber;
     // Loop over all chambers (tracking and trigger)
     for (Int_t ch = 0; ch < AliMUONConstants::NCh(); ch++) {
	 char alu[8];
	 char gas[8];
     
	 iChamber=(AliMUONChamber*) (*fChambers)[ch];
	 // Z of the chamber
	 zpos=iChamber->Z(); 
	 dAlu=iChamber->DAlu();
	 if (ch < AliMUONConstants::NTrackingCh()) {
	   // tracking chambers
	     sprintf(alu,"CA0%1d",ch);
	     sprintf(gas,"CG0%1d",ch);	 
	 } else {
	   // trigger chambers
	     sprintf(alu,"CA%2d",ch);
	     sprintf(gas,"CG%2d",ch);	 
	 }
//
	 tpar[0] = iChamber->RInner(); 
	 tpar[1] = iChamber->ROuter();
	 tpar[2] = (dAlu+0.2)/2.;
	 if (ch !=4 && ch !=5) {
	     gMC->Gsvolu(alu, "TUBE", idAlu, tpar, 3);
	     tpar[2] = 0.1;
	     gMC->Gsvolu(gas, "TUBE", idAir, tpar, 3);
	 } else {
	     gMC->Gsvolu(alu, "TUBE", idAlu, tpar, 3);
	     tpar[2] = 0.1;
	     gMC->Gsvolu(gas, "TUBE", idAir, tpar, 3);
	 }
	 gMC->Gspos(gas, 1, alu,  0., 0., 0., 0, "ONLY");
	 gMC->Gspos(alu, 1, "ALIC", 0., 0., zpos, 0, "ONLY");
     }
}

//___________________________________________
void AliMUONv0::CreateMaterials()
{
// Creates materials for coarse geometry
    AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
    AliMaterial(9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);

    Float_t epsil  = .001; // Tracking precision, 
    Float_t stemax = -1.;  // Maximum displacement for multiple scat 
    Float_t tmaxfd = -20.; // Maximum angle due to field deflection 
    Float_t deemax = -.3;  // Maximum fractional energy loss, DLS 
    Float_t stmin  = -.8;
    Int_t isxfld   = gAlice->Field()->Integ();
    Float_t sxmgmx = gAlice->Field()->Max();

    //
    //    Air 
    AliMedium(1, "AIR_CH_US         ", 15, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(4, "ALU_CH_US          ", 9, 0, isxfld, sxmgmx, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);

}

void AliMUONv0::Init()
{
   // 
   // Initialize Tracking Chambers
   //
    char vName[8];
    printf("\n\n\n Start Init for version 0 - CPC chamber type\n\n\n");
    for (Int_t i=0; i<AliMUONConstants::NCh(); i++) {
// Initialise chamber
	((AliMUONChamber*) (*fChambers)[i])->Init();
// Set sensitive volume Id
	if (i < AliMUONConstants::NTrackingCh()) {
	    // tracking chambers
	    sprintf(vName,"CG0%1d",i);	 
	} else {
	    // trigger chambers
	    sprintf(vName,"CG%2d",i);	 
	}
	((AliMUONChamber*) (*fChambers)[i])->SetGid(gMC->VolId(vName));
    }
}

void AliMUONv0::StepManager()
{
//
// Step manager for hit density simulations
  Int_t          copy, id;
  static Int_t   idvol;
  static Int_t   vol[2];
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;
  Float_t        theta,phi;
  
  //  modifs perso
  static Float_t hits[15];

  TClonesArray &lhits = *fHits;
  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  idvol=-1;
  id=gMC->CurrentVolID(copy);
  
    for (Int_t i=1; i<=AliMUONConstants::NCh(); i++) {
      if(id==((AliMUONChamber*)(*fChambers)[i-1])->GetGid()){ 
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
  // record hits when track enters ...
  if( !(gMC->TrackCharge()) ) return; 
  if( gMC->IsTrackEntering()) {
      Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
      Double_t rt = TMath::Sqrt(tc);
      theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
      phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
      hits[0] = Float_t(ipart);             // Geant3 particle type
      hits[1] = pos[0];                     // X-position for hit
      hits[2] = pos[1];                     // Y-position for hit
      hits[3] = pos[2];                     // Z-position for hit
      hits[4] = theta;                      // theta angle of incidence
      hits[5] = phi;                        // phi angle of incidence 
      hits[8] = -1;                         // first padhit
      hits[9] = -1;                         // last pad hit

      // modifs personel
      hits[10] = mom[3]; // hit Energy
      hits[11] = mom[0]; // Px
      hits[12] = mom[1]; // Py
      hits[13] = mom[2]; // Pz
      hits[14] = gMC->TrackTime();
      
      // fin modifs perso
      new(lhits[fNhits++]) 
	  AliMUONHit(fIshunt,gAlice->CurrentTrack(),vol,hits);

  }
}









