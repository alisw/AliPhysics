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

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response //
////////////////////////////////////////////////////////////

/* $Id$ */

#include "AliMUONFactory.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONResponseV0.h"
#include "AliMUONResponseTrigger.h"
#include "AliMUONSegmentationV01.h"
#include "AliMUONSegmentationV02.h"
#include "AliMUONSegmentationSlat.h"
#include "AliMUONSegmentationSlatN.h"
#include "AliMUONSegmentationTriggerX.h"
#include "AliMUONSegmentationTriggerY.h"
#include "AliLog.h"

ClassImp(AliMUONFactory)

//__________________________________________________________________________
AliMUONFactory::AliMUONFactory()
  : TObject(),
    fMUON(0),
    fResponse0(0)
{
//
}

//__________________________________________________________________________
AliMUONFactory::AliMUONFactory(const AliMUONFactory& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//__________________________________________________________________________

AliMUONFactory::~AliMUONFactory()
{
//
}

//__________________________________________________________________________
AliMUONFactory&  AliMUONFactory::operator=(const AliMUONFactory& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    
          
//__________________________________________________________________________
void AliMUONFactory::BuildCommon() 
{
//
// Construct the default response.
//

	// Default response: 5 mm of gas
	fResponse0 = new AliMUONResponseV0;
	fResponse0->SetSqrtKx3AndDeriveKx2Kx4(0.7131); // sqrt(0.5085)
	fResponse0->SetSqrtKy3AndDeriveKy2Ky4(0.7642); // sqrt(0.5840)
	fResponse0->SetPitch(0.25); // anode-cathode distance
	fResponse0->SetSigmaIntegration(10.);
	fResponse0->SetChargeSlope(10);
	fResponse0->SetChargeSpread(0.18, 0.18);
	fResponse0->SetMaxAdc(4096);
	fResponse0->SetSaturation(3000);
	fResponse0->SetZeroSuppression(6);
}	
	
//__________________________________________________________________________
void AliMUONFactory::BuildStation1() 
{
//--------------------------------------------------------
// Configuration for Chamber TC1/2  (Station 1) ----------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
	responseSt1->SetSaturation(3000);
	responseSt1->SetZeroSuppression(6);
	
        //--------------------------------------------------------
        // Configuration for Chamber TC1/2  (Station 1) ----------           

	Float_t rseg1[4]={17.5, 55.2, 71.3, 95.5};
	Int_t   nseg1[4]={4, 4, 2, 1};
//
	Int_t chamber=1;
//      ^^^^^^^^^^^^^^^^
	fMUON->SetNsec(chamber-1,2);
//
	AliMUONSegmentationV01 *seg11=new AliMUONSegmentationV01(4);
	
	seg11->SetSegRadii(rseg1);
	seg11->SetPadSize(2.4, 0.4); // smaller pad size
	seg11->SetDAnod(0.20); // smaller distance between anode wires
	seg11->SetPadDivision(nseg1);
	
	fMUON->SetSegmentationModel(chamber-1, 1, seg11);
	
	AliMUONSegmentationV02 *seg12=new AliMUONSegmentationV02(4);
	seg12->SetSegRadii(rseg1); 
	seg12->SetPadSize(0.6, 1.6); // smaller pad size
	seg12->SetDAnod(0.20); // smaller distance between anode wires
	seg12->SetPadDivision(nseg1);
	
	fMUON->SetSegmentationModel(chamber-1, 2, seg12);
	
	fMUON->SetResponseModel(chamber-1, responseSt1); // special response	    
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
	
	chamber=2;
//      ^^^^^^^^^
//
	fMUON->SetNsec(chamber-1,2);
//
	AliMUONSegmentationV01 *seg21=new AliMUONSegmentationV01(4);
	seg21->SetSegRadii(rseg1);
	seg21->SetPadSize(2.4, 0.4); // smaller pad size
	seg21->SetDAnod(0.20); // smaller distance between anode wires
	seg21->SetPadDivision(nseg1);
	fMUON->SetSegmentationModel(chamber-1, 1, seg21);
//
	AliMUONSegmentationV02 *seg22=new AliMUONSegmentationV02(4);
	seg22->SetSegRadii(rseg1); 
	seg22->SetPadSize(0.6, 1.6); // smaller pad size
	seg22->SetDAnod(0.20); // smaller distance between anode wires
	seg22->SetPadDivision(nseg1);
	fMUON->SetSegmentationModel(chamber-1, 2, seg22);
	
	fMUON->SetResponseModel(chamber-1, responseSt1); // special response
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
	
}

//__________________________________________________________________________
void AliMUONFactory::BuildStation2() 
{
//
//--------------------------------------------------------
// Configuration for Chamber TC3/4 (Station 2) -----------
///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// Float_t rseg2[4]={23.5, 87.7, 122.4, 122.5};

	Float_t rseg2[4]={23.5, 53.5, 90.5, 122.5};       
	Int_t   nseg2[4]={4, 4, 2, 1};
//
	Int_t chamber=3;
//      ^^^^^^^^^^^^^^^^
	fMUON->SetNsec(chamber-1,2);
//
	AliMUONSegmentationV01 *seg31=new AliMUONSegmentationV01(4);
	seg31->SetSegRadii(rseg2);
	seg31->SetPadSize(3.0, 0.5);
	seg31->SetDAnod(3.0/3./4);
	seg31->SetPadDivision(nseg2);
	fMUON->SetSegmentationModel(chamber-1, 1, seg31);
//
	AliMUONSegmentationV02 *seg32=new AliMUONSegmentationV02(4);
	seg32->SetSegRadii(rseg2); 
	seg32->SetPadSize(0.75, 2.0);
	seg32->SetPadDivision(nseg2);
	seg32->SetDAnod(3.0/3./4);
	
	fMUON->SetSegmentationModel(chamber-1, 2, seg32);
	
	fMUON->SetResponseModel(chamber-1, fResponse0);	    
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
	
	chamber=4;
//      ^^^^^^^^^
//
	fMUON->SetNsec(chamber-1,2);
//
	AliMUONSegmentationV01 *seg41=new AliMUONSegmentationV01(4);
	seg41->SetSegRadii(rseg2);
	seg41->SetPadSize(3.0, 0.5);
	seg41->SetDAnod(3.0/3./4);
	seg41->SetPadDivision(nseg2);
	fMUON->SetSegmentationModel(chamber-1, 1, seg41);
//
	AliMUONSegmentationV02 *seg42=new AliMUONSegmentationV02(4);
	seg42->SetSegRadii(rseg2); 
	seg42->SetPadSize(0.75, 2.0);
	seg42->SetPadDivision(nseg2);
	seg42->SetDAnod(3.0/3./4);
	
	fMUON->SetSegmentationModel(chamber-1, 2, seg42);
	
	fMUON->SetResponseModel(chamber-1, fResponse0);	    
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
}	
	
	
//__________________________________________________________________________
void AliMUONFactory::BuildStation3() 
{
//--------------------------------------------------------
// Configuration for Chamber TC5/6  (Station 3) ----------          
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	Int_t   nseg3[4]={4, 4, 2, 1};
	Int_t   npcb5[36] = {0,0,2,0,
			     0,0,3,0,
			     0,1,3,0,
			     0,2,2,0,
			     0,2,2,0, 
			     0,2,2,0, 
			     0,1,3,0, 
			     0,0,3,0,
			     0,0,2,0};
	
	Float_t shift = 1.5/2.;
	Float_t xpos5[9]    = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5};
	Float_t ypos5       = -(20.+4.*(40.-2.*shift));
	
	Int_t chamber=5;
	fMUON->SetNsec(chamber-1,2);
	AliMUONSegmentationSlat *seg51=new AliMUONSegmentationSlat(4);
	seg51->SetNSlats(9); 
	seg51->SetShift(shift);  
	seg51->SetNPCBperSector(npcb5); 
	seg51->SetSlatXPositions(xpos5);
	seg51->SetSlatYPosition(ypos5);
	seg51->SetPadSize(10.,0.5);
	seg51->SetDAnod(0.25);
	seg51->SetPadDivision(nseg3);
	fMUON->SetSegmentationModel(chamber-1, 1, seg51);
	
	AliMUONSegmentationSlatN *seg52=new AliMUONSegmentationSlatN(4);
	seg52->SetNSlats(9); 
	seg52->SetShift(shift);  
	seg52->SetNPCBperSector(npcb5); 
	seg52->SetSlatXPositions(xpos5);
	seg52->SetSlatYPosition(ypos5);
	seg52->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
	seg52->SetDAnod(0.25);
	seg52->SetPadDivision(nseg3);
	fMUON->SetSegmentationModel(chamber-1, 2, seg52);
	fMUON->SetResponseModel(chamber-1, fResponse0);      
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
	
	chamber=6;
	fMUON->SetNsec(chamber-1,2);
	AliMUONSegmentationSlat *seg61=new AliMUONSegmentationSlat(4);
	seg61->SetNSlats(9); 
	seg61->SetShift(shift);  
	seg61->SetNPCBperSector(npcb5); 
	seg61->SetSlatXPositions(xpos5);
	seg61->SetSlatYPosition(ypos5);
	seg61->SetPadSize(10.,0.5);
	seg61->SetDAnod(0.25);
	seg61->SetPadDivision(nseg3);
	fMUON->SetSegmentationModel(chamber-1, 1, seg61);
	
	AliMUONSegmentationSlatN *seg62=new AliMUONSegmentationSlatN(4);
	seg62->SetNSlats(9); 
	seg62->SetShift(shift);  
	seg62->SetNPCBperSector(npcb5); 
	seg62->SetSlatXPositions(xpos5);
	seg62->SetSlatYPosition(ypos5);
	seg62->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
	seg62->SetDAnod(0.25);
	seg62->SetPadDivision(nseg3);
	fMUON->SetSegmentationModel(chamber-1, 2, seg62);
	fMUON->SetResponseModel(chamber-1, fResponse0);      
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
}

	
//__________________________________________________________________________
void AliMUONFactory::BuildStation4() 
{
//--------------------------------------------------------
// Configuration for Chamber TC7/8  (Station 4) ----------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	Int_t   nseg4[4]={4, 4, 2, 1};
	
	Int_t chamber=7;
//      ^^^^^^^^^^^^^^^^
	
	fMUON->SetNsec(chamber-1,2);
//
	AliMUONSegmentationSlat *seg71=new AliMUONSegmentationSlat(4);
	Float_t shift = 1.5/2.;
	Int_t npcb7[52] = {0,0,0,2,
			   0,0,0,3,
			   0,0,2,2,
			   0,0,3,2,
			   0,2,2,1,
			   0,2,2,2,
			   0,1,2,2, 
			   0,2,2,2, 
			   0,2,2,1, 
			   0,0,3,2, 
			   0,0,2,2, 
			   0,0,0,3,
			   0,0,0,2};
	Float_t xpos7[13]   = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 44., 4.5, 4.5, 4.5, 4.5, 4.5, 4.5};
	Float_t ypos7       = -(20.+6.*(40.-2.*shift));  
	
	seg71->SetNSlats(13);  
	seg71->SetShift(shift);  
	seg71->SetNPCBperSector(npcb7); 
	seg71->SetSlatXPositions(xpos7);
	seg71->SetSlatYPosition(ypos7);
	
	seg71->SetPadSize(10.,0.5);
	seg71->SetDAnod(0.25);
	seg71->SetPadDivision(nseg4);
	fMUON->SetSegmentationModel(chamber-1, 1, seg71);
	
	AliMUONSegmentationSlatN *seg72=new AliMUONSegmentationSlatN(4);
	
	fMUON->SetSegmentationModel(chamber-1, 2, seg72);
	seg72->SetNSlats(13);  
	seg72->SetShift(shift);   
	seg72->SetNPCBperSector(npcb7); 
	seg72->SetSlatXPositions(xpos7);
	seg72->SetSlatYPosition(ypos7);
	seg72->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
	seg72->SetDAnod(0.25);
	seg72->SetPadDivision(nseg4);
	
	fMUON->SetResponseModel(chamber-1, fResponse0);	    
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
	
	chamber=8;
//      ^^^^^^^^^
	fMUON->SetNsec(chamber-1,2);
//
	AliMUONSegmentationSlat *seg81=new AliMUONSegmentationSlat(4);
	
	seg81->SetNSlats(13);  
	seg81->SetShift(shift);  
	seg81->SetNPCBperSector(npcb7); 
	seg81->SetSlatXPositions(xpos7);
	seg81->SetSlatYPosition(ypos7);
	seg81->SetPadSize(10.,0.5);
	seg81->SetDAnod(0.25);
	seg81->SetPadDivision(nseg4);
	fMUON->SetSegmentationModel(chamber-1, 1, seg81);
	
	AliMUONSegmentationSlat *seg82=new AliMUONSegmentationSlatN(4);
	
	fMUON->SetSegmentationModel(chamber-1, 2, seg82);
	seg82->SetNSlats(13);  
	seg82->SetShift(shift);  
	seg82->SetNPCBperSector(npcb7); 
	seg82->SetSlatXPositions(xpos7);
	seg82->SetSlatYPosition(ypos7);
	seg82->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
	seg82->SetDAnod(0.25);
	seg82->SetPadDivision(nseg4);
	
	fMUON->SetResponseModel(chamber-1, fResponse0);	    
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
}

//__________________________________________________________________________
void AliMUONFactory::BuildStation5() 
{	
//--------------------------------------------------------
// Configuration for Chamber TC9/10  (Station 5) ---------           
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	Int_t chamber=9;
//      ^^^^^^^^^^^^^^^^

	fMUON->SetNsec(chamber-1,2);
//
	AliMUONSegmentationSlat *seg91=new AliMUONSegmentationSlat(4);

	Int_t   nseg4[4]={4, 4, 2, 1};
	Float_t shift = 1.5/2.;
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
	
	Float_t xpos9[13]   = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 44., 4.5, 4.5, 4.5, 4.5, 4.5, 4.5};
	Float_t ypos9       = -(20.+6.*(40.-2.*shift));
	
	seg91->SetNSlats(13);  
	seg91->SetShift(shift);  
	seg91->SetNPCBperSector(npcb9); 
	seg91->SetSlatXPositions(xpos9);
	seg91->SetSlatYPosition(ypos9);
	seg91->SetPadSize(10.,0.5);
	seg91->SetDAnod(0.25);
	seg91->SetPadDivision(nseg4);
	fMUON->SetSegmentationModel(chamber-1, 1, seg91);
	
	AliMUONSegmentationSlatN *seg92=new AliMUONSegmentationSlatN(4);
	
	fMUON->SetSegmentationModel(chamber-1, 2, seg92);
	seg92->SetNSlats(13);  
	seg92->SetShift(shift);   
	seg92->SetNPCBperSector(npcb9); 
	seg92->SetSlatXPositions(xpos9);
	seg92->SetSlatYPosition(ypos9);
	seg92->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
	seg92->SetDAnod(0.25);
	seg92->SetPadDivision(nseg4);
	
	fMUON->SetResponseModel(chamber-1, fResponse0);	    
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
	
	chamber=10;
//      ^^^^^^^^^
	fMUON->SetNsec(chamber-1,2);
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
	fMUON->SetSegmentationModel(chamber-1, 1, seg101);
	
	AliMUONSegmentationSlatN *seg102=new AliMUONSegmentationSlatN(4);
	
	fMUON->SetSegmentationModel(chamber-1, 2, seg102);
	seg102->SetNSlats(13);  
	seg102->SetShift(shift);   
	seg102->SetNPCBperSector(npcb9); 
	seg102->SetSlatXPositions(xpos9);
	seg102->SetSlatYPosition(ypos9);
	seg102->SetPadSize(1., 10.); // DeltaX(non bending) = 2 * DeltaY(bending)
	seg102->SetDAnod(0.25);
	seg102->SetPadDivision(nseg4);
	
	fMUON->SetResponseModel(chamber-1, fResponse0);	    
	fMUON->Chamber(chamber-1).SetChargeCorrel(0.11); // 11% charge spread
}

//__________________________________________________________________________
void AliMUONFactory::BuildStation6() 
{	
//--------------------------------------------------------
// Configuration for Trigger Stations -------------------- 
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

// Cluster-size off
	AliMUONResponseTrigger* responseTrigger0 =  new AliMUONResponseTrigger;
// Cluster-size on  
// AliMUONResponseTriggerV1* responseTrigger0 =  new AliMUONResponseTriggerV1;
 
	Int_t chamber=11;
	fMUON->SetNsec(chamber-1,2);
	AliMUONSegmentationTriggerX *seg111=new AliMUONSegmentationTriggerX;
	fMUON->SetSegmentationModel(chamber-1, 1, seg111);
	AliMUONSegmentationTriggerY *seg112=new AliMUONSegmentationTriggerY;
	fMUON->SetSegmentationModel(chamber-1, 2, seg112);
	
	fMUON->SetResponseModel(chamber-1, responseTrigger0);      
	fMUON->Chamber(chamber-1).SetChargeCorrel(0); // same charge on cathodes
	
 
	chamber=12;
	fMUON->SetNsec(chamber-1,2);
	AliMUONSegmentationTriggerX *seg121=new AliMUONSegmentationTriggerX;
	fMUON->SetSegmentationModel(chamber-1, 1, seg121);
	AliMUONSegmentationTriggerY *seg122=new AliMUONSegmentationTriggerY;
	fMUON->SetSegmentationModel(chamber-1, 2, seg122);
	
	fMUON->SetResponseModel(chamber-1, responseTrigger0);
	fMUON->Chamber(chamber-1).SetChargeCorrel(0); // same charge on cathodes
	
	chamber=13;
	fMUON->SetNsec(chamber-1,2);
	AliMUONSegmentationTriggerX *seg131=new AliMUONSegmentationTriggerX;
	fMUON->SetSegmentationModel(chamber-1, 1, seg131);
	AliMUONSegmentationTriggerY *seg132=new AliMUONSegmentationTriggerY;
	fMUON->SetSegmentationModel(chamber-1, 2, seg132);
	fMUON->SetResponseModel(chamber-1, responseTrigger0);      
	fMUON->Chamber(chamber-1).SetChargeCorrel(0); // same charge on cathodes
	
	chamber=14;
	fMUON->SetNsec(chamber-1,2);
	AliMUONSegmentationTriggerX *seg141=new AliMUONSegmentationTriggerX;
	fMUON->SetSegmentationModel(chamber-1, 1, seg141);
	AliMUONSegmentationTriggerY *seg142=new AliMUONSegmentationTriggerY;
	fMUON->SetSegmentationModel(chamber-1, 2, seg142);
	
	fMUON->SetResponseModel(chamber-1, responseTrigger0); 
	fMUON->Chamber(chamber-1).SetChargeCorrel(0); // same charge on cathodes
}	

//__________________________________________________________________________
void AliMUONFactory::Build(AliMUON* where, const char* what) 
{
//
// Construct MUON from chambers, segmentation and responses
//

    fMUON = where;
    char tmp[20];
    strcpy(tmp, what);

    if (strcmp(tmp, "default")==0) {
      // Set default parameters
      fMUON->SetIshunt(0);
      fMUON->SetMaxStepGas(0.1);
      fMUON->SetMaxStepAlu(0.1);

      // Build all stations
      BuildCommon();
      BuildStation1();
      BuildStation2();
      BuildStation3();
      BuildStation4();
      BuildStation5();
      BuildStation6();
    } 
    else
	  AliDebug(0,"Non default version of MUON selected. You have to construct yourself the MUON elements !!");
}

//__________________________________________________________________________
void AliMUONFactory::BuildStation(AliMUON* where, Int_t stationNumber) 
{
//
// Construct MUON from chambers, segmentation and responses
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

    fMUON = where;
    if (!fResponse0) BuildCommon(); 
    
    switch (stationNumber) {    
      case 1:  BuildStation1(); break;
      case 2:  BuildStation2(); break;
      case 3:  BuildStation3(); break;
      case 4:  BuildStation4(); break;
      case 5:  BuildStation5(); break;
      case 6:  BuildStation6(); break;
    
      default: AliFatal("Wrong station number");
    }  
}         
