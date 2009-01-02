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

/* $Id: AliACORDEv1.cxx,v 1.2 2007/12/03 08:40:00 hristov Exp $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// ALICE Cosmic Ray Trigger                                                  //
//                                                                           //
//  This class contains the functions for final version of the ALICE Cosmic  //
//  Ray Trigger. This version will be used to simulation comic rays in alice //
//  with all the detectors. 						     //
//  It include geometry and hits (position and momentum)                     //
//                                                                           //
//                  Send comments to:                                        //
//									     //
//      Arturo Fernandez Tellez   	<afernand@fcfm.buap.mx>              //
//      Eleazar Cuautle Flores    	<ecuautle@nucleares.unam.mx>         //
//	Mario Rodriguez	Cahuantzi 	<mrodrigu@mail.cern.ch>		     //	
//									     //
//			Puebla, Pue. Mexico December 2007                    //
//									     //
//	Last Update: Aug. 4th 2008				             //
///////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <TGeoMatrix.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>

#include "AliConst.h"
#include "AliRun.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoPcon.h"
#include "TGeoPgon.h"
#include "TGeoTrd1.h"
#include "TGeoCompositeShape.h"
#include "TGeoPara.h"

#include "AliACORDEv1.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVirtualMC.h>
#include <TPDGCode.h>

#include "AliRun.h"
#include "AliConst.h"
#include "AliACORDEhit.h"
#include "AliACORDEConstants.h"
#include "AliMC.h"
#include "AliLog.h"

ClassImp(AliACORDEv1)
 
//_____________________________________________________________________________
AliACORDEv1::AliACORDEv1()
  : AliACORDE()
{
  //
  // Default constructor
  fIshunt = 0;
  fHits = 0;
  //
} 
//_____________________________________________________________________________
AliACORDEv1::AliACORDEv1(const char *name, const char *title)
  : AliACORDE(name, title)
{
  //
  // Standard constructor
  //
  fIshunt = 1; // All hits are associated with primary particles 
  fHits =  new TClonesArray("AliACORDEhit",400);
  gAlice->GetMCApp()->AddHitList(fHits);
}
//_____________________________________________________________________________
AliACORDEv1::~AliACORDEv1()
{
  //
  // Default destructor
  //
}
//_____________________________________________________________________________
void AliACORDEv1::CreateGeometry()
{
  CreateAcorde();
}
void AliACORDEv1::CreateAcorde()
{  


  //  _______________________________________________________________________________
  // |										     |	
  // |										     |	
  // |	**** Acorde's Geometry using the TGeo Class....January 2008 ****             |
  // |										     |	
  // |	 	ACORDE--> Volume for ACORDE in Alice's Magnet                        |
  // |		ACORDE1_a--> Volume for frame of Acorde's Module		     |
  // |		ACORDE10--> Volume for frame of Acorde's Module			     |
  // |		ACORDE2--> Volume for scintillators in Acorde's Module  	     |
  // |		ACORDE7--> Volume for bars					     |
  // |		ACORDE7_1--> Volume for side's bars				     |	
  // |	        ACORDE8--> Volume for supports					     |	
  // |          ACORDE9--> Volume for supports					     |		
  // |		ACORDE_SUPPORT--> Volume that contains a full Acorde's support	     |
  // |		ALL_ACORDE_MODULES--> Volume that contains ALL Acorde's module	     |	
  // |          ACORDE_MODULE--> Volume that represents ONE Acorde-Module            |
  // |		ACORDE_1--> Volume that contains the bars&supports in-face	     |
  // |		ACORDE_2--> Volume that contains the bars&supports up-face	     | 	
  // |		ACORDE_3--> Volume that contains the bars&supports out-face	     |
  // |										     |	
  // |_______________________________________________________________________________|


	// Call the global constants for the Modules
	
	AliACORDEConstants* constants = AliACORDEConstants::Instance();

	// Get the Alice Volume

	TGeoVolume *alice = gGeoManager->GetVolume("ALIC");

	// Define some materials & medium

	//*** Support & Bars***

	TGeoMedium* al    = gGeoManager->GetMedium("ACORDE_ALU_C0");
	TGeoMedium* med6  = gGeoManager->GetMedium("ACORDE_CPV scint.1");
	
	//Define a FULL-ACORDE-VOLUME

	TGeoVolume *aCORDE = new TGeoVolumeAssembly("ACORDE");


	// Define 6 master volumes for ACORDE
	
	TGeoVolume *inFace = new TGeoVolumeAssembly("ACORDE_1");
	TGeoVolume *upFace = new TGeoVolumeAssembly("ACORDE_2");
	TGeoVolume *outFace = new TGeoVolumeAssembly("ACORDE_3");
	TGeoVolume *modules = new TGeoVolumeAssembly("ALL_ACORDE_MODULES");

	// Define global variables

	Float_t box[3];
	Int_t count;
	Float_t dy=10;//-->displacement of the support and bars of ACORDE
	Float_t dy2=66.5;//-->displacement of the support and bars of ACORDE
	Float_t placedAt;
	Float_t small=0.05;

	// Define the position of support and bars for the sides faces

	Float_t des = 22*0.7071;
  
	// Define rotation Matrix for Side's faces in Alice
	
	TGeoRotation *idrotm231 = new TGeoRotation("idrotm231",90, 45, 90, 135, 0, 0);
	TGeoRotation *idrotm232 = new TGeoRotation("idrotm232",90, 315, 90, 45, 0, 0);

	// Begin the Geometry of the structure for ACORDE

	// *** Definition of ACORDE's Modules ***

	// Define Measures of ACORDE's Modules

	box[0] = constants->ModuleLength()/2;
	box[1] = constants->ModuleHeight()/2;
	box[2] = constants->ModuleWidth()/2;

	// Define Measures of Scintillators
	
	Float_t pbox[3];	
	pbox[0] = constants->PlasticLength()/2;
	pbox[1] = constants->PlasticHeight()/2;
	pbox[2] = constants->PlasticWidth()/2;


	// Create the Modules, Scintillators & Metallic Frame

	//*** Aluminium frame ***

	TGeoBBox *acorde1 = new TGeoBBox("acorde1",box[0],box[1],26/20+2);
	TGeoBBox *acorde10 = new TGeoBBox("acorde10",26/20,box[1],box[2]+3);
	TGeoVolume *aCORDE1qa = new TGeoVolume("ACORDE1_a",acorde1,al);
	TGeoVolume *aCORDE10 = new TGeoVolume("ACORDE10",acorde10,al);

	//*** Scintillators ***

	TGeoBBox *acorde2 = new TGeoBBox("acorde2",pbox[0],pbox[1],pbox[2]);
	TGeoVolume *aCORDE2 = new TGeoVolume("ACORDE2",acorde2,med6);


	// Here I define & construct a Master Volume ("ACORDE_MODULE") for one Module in ACORDE

	TGeoVolume *acomodule = new TGeoVolumeAssembly("ACORDE_MODULE");
	acomodule->AddNode(aCORDE1qa,1,new TGeoTranslation("aco1",0,0,13));
	acomodule->AddNode(aCORDE1qa,2,new TGeoTranslation("aco10",0,0,-13));
	acomodule->AddNode(aCORDE10,3,new TGeoTranslation("aco10",293/2+5,0,0));
	acomodule->AddNode(aCORDE10,4,new TGeoTranslation("aco10",-293/2-5,0,0));
        placedAt = pbox[1]+constants->ProfileThickness()-constants->ModuleHeight()/2+small;
	acomodule->AddNode(aCORDE2,5,new TGeoTranslation("aco2",placedAt,0,0));
        placedAt = placedAt + 2.0*pbox[1]+small;
	acomodule->AddNode(aCORDE2,6,new TGeoTranslation("aco2",placedAt,-1,0));

	// Put the Modules of In-Face
	
	count=1;
	for(Int_t i=1;i<9;i++){

		Float_t posx = constants->ModulePositionX(i);
		Float_t posy = constants->ModulePositionY(i);
		Float_t posz = constants->ModulePositionZ(i);	
                Int_t moduleElectronicID = constants->ModuleElectronicChannel(i); 
		
		modules->AddNode(acomodule,moduleElectronicID,
			new TGeoCombiTrans("aco01",posx,posy,posz,idrotm232));
		count++;

	}

	count=9;
	for(Int_t i=10;i<20;i++){
		Float_t posx = constants->ModulePositionX(i);
		Float_t posy = constants->ModulePositionY(i);
		Float_t posz = constants->ModulePositionZ(i);
		Int_t moduleElectronicID = constants->ModuleElectronicChannel(i); 	

		modules->AddNode(acomodule,moduleElectronicID,
			new TGeoCombiTrans("aco01",posx,posy,posz,idrotm232));
	}

	// Put he Modules of Up-Face

	count=1;
	for(Int_t i=20;i<40;i++){
		Float_t posx = constants->ModulePositionX(i);
		Float_t posy = constants->ModulePositionY(i);
		Float_t posz = constants->ModulePositionZ(i);	
		Int_t moduleElectronicID = constants->ModuleElectronicChannel(i); 

		modules->AddNode(acomodule,moduleElectronicID,new TGeoTranslation("aco01",posx,posy,posz));
		count++;
	}

	// Put the Modules of Out-Face

	count=1;
	for(Int_t i=40;i<50;i++){
		Float_t posx = constants->ModulePositionX(i);
		Float_t posy = constants->ModulePositionY(i);
		Float_t posz = constants->ModulePositionZ(i);	
		Int_t moduleElectronicID = constants->ModuleElectronicChannel(i); 

		modules->AddNode(acomodule,moduleElectronicID,
			new TGeoCombiTrans("aco01",posx,posy,posz,idrotm231));
		count++;
	}

	// Put the Modules of Out-Face

	count=11;
	for(Int_t i=51;i<59;i++){
		Float_t posx = constants->ModulePositionX(i);
		Float_t posy = constants->ModulePositionY(i);
		Float_t posz = constants->ModulePositionZ(i);	
		Int_t moduleElectronicID = constants->ModuleElectronicChannel(i); 

	if ((i==57) || (i==56)){
		 modules->AddNode(acomodule,moduleElectronicID,
					new TGeoCombiTrans("aco01",posx,posy,posz,idrotm231));
	}else{
		modules->AddNode(acomodule,moduleElectronicID,
			new TGeoCombiTrans("aco01",posx,posy,posz,idrotm231));
		}count++;
	}


	// Put th Modules ITS-ACORDE

	if (GetITSGeometry()) {

		modules->AddNode(acomodule,constants->ModuleElectronicChannel(50),new TGeoTranslation("ITS-3",
				constants->ExtraModulePositionX(),
				constants->ExtraModulePositionY(),
				constants->ExtraModulePositionZ(0)));

		modules->AddNode(acomodule,constants->ModuleElectronicChannel(59),new TGeoTranslation("ITS-4",
				constants->ExtraModulePositionX(),
				constants->ExtraModulePositionY(),
				constants->ExtraModulePositionZ(1)));

		modules->AddNode(acomodule,constants->ModuleElectronicChannel(0),new TGeoTranslation("ITS-1",
				constants->ExtraModulePositionX(),
				constants->ExtraModulePositionY(),
				constants->ExtraModulePositionZ(2)));

		modules->AddNode(acomodule,constants->ModuleElectronicChannel(9),new TGeoTranslation("ITS-2",
				constants->ExtraModulePositionX(),
				constants->ExtraModulePositionY(),
				constants->ExtraModulePositionZ(3)));


		} 
	else {


		modules->AddNode(acomodule,61,new TGeoTranslation("its1",
				constants->ModulePositionX(0),
				constants->ModulePositionY(0),
				constants->ModulePositionZ(0)));

		modules->AddNode(acomodule,62,new TGeoTranslation("its2",
				constants->ModulePositionX(9),
				constants->ModulePositionY(9),
				constants->ModulePositionZ(9)));

		modules->AddNode(acomodule,63,new TGeoTranslation("its3",
				constants->ModulePositionX(50),
				constants->ModulePositionY(50),
				constants->ModulePositionZ(50)));

		modules->AddNode(acomodule,64,new TGeoTranslation("its4",
				constants->ModulePositionX(59),
				constants->ModulePositionY(59),
				constants->ModulePositionZ(59)));

		} // end if (fITSGeometry)



	//*** Begin the structure of support & bars for ACORDE ***

	// Define a volume for the bars (up-face)

	box[0]=5;
//	box[1]=40;
	box[1]=33;
	box[2]=5;
	Float_t z1 = 21 ;
	TGeoBBox *acorde00 = new TGeoBBox("acorde00",box[0],box[1],box[2]);

	TGeoVolume *aCORDE00 = new TGeoVolume("ACORDE00",acorde00,al);

	count=25;
	for (Int_t ma=20;ma<=24;ma++)
	{
		TGeoTranslation *aco00=new TGeoTranslation("aco00",
					constants->SupportModulePositionX(ma)-0.5*293+dy2,
					constants->SupportModulePositionY(ma)-box[1]-z1,
					constants->SupportModulePositionZ(ma));

		upFace->AddNode(aCORDE00,count,aco00);

		TGeoTranslation *aco00q1=new TGeoTranslation("aco00q1",
					-(constants->SupportModulePositionX(ma)-0.5*293+dy2),
					constants->SupportModulePositionY(ma)-box[1]-z1,
					constants->SupportModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+1,aco00q1);

		TGeoTranslation *aco00q2=new TGeoTranslation("aco00q2",
					constants->SupportModulePositionX(ma)+0.5*293-dy2,
					constants->SupportModulePositionY(ma)-box[1]-z1,
					constants->SupportModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+2,aco00q2);

		TGeoTranslation *aco00q3=new TGeoTranslation("aco00q3",
					-(constants->SupportModulePositionX(ma)+0.5*293-dy2),
					constants->SupportModulePositionY(ma)-box[1]-z1,
					constants->SupportModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+3,aco00q3);
		count=count+4;

		ma++;
	}

	count=41;
	for(Int_t ma=25;ma<=29;ma++)
	{
		TGeoTranslation *aco00=new TGeoTranslation("aco00",
					constants->SupportModulePositionX(ma)-0.5*293+dy2,
					constants->SupportModulePositionY(ma)-box[1]-z1,
					constants->SupportModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count,aco00);

		TGeoTranslation *aco00q1=new TGeoTranslation("aco00q1",
					-(constants->SupportModulePositionX(ma)-0.5*293+dy2),
					constants->SupportModulePositionY(ma)-box[1]-z1,
					constants->SupportModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+1,aco00q1);

		TGeoTranslation *aco00q2=new TGeoTranslation("aco00q2",
					constants->SupportModulePositionX(ma)+0.5*293-dy2,
					constants->SupportModulePositionY(ma)-box[1]-z1,
					constants->SupportModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+2,aco00q2);

		TGeoTranslation *aco00q3=new TGeoTranslation("aco00q3",
					-(constants->SupportModulePositionX(ma)+0.5*293-dy2),
					constants->SupportModulePositionY(ma)-box[1]-z1,
					constants->SupportModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+3,aco00q3);
		count=count+4;
		ma++;
	}

	TGeoTranslation *c1 = new TGeoTranslation ("c1",
				constants->SupportModulePositionX(20)-0.5*293,
				constants->SupportModulePositionY(20)-box[1]-z1,
				constants->SupportModulePositionZ(20)-40);
	TGeoTranslation *c2 = new TGeoTranslation ("c2",
				constants->SupportModulePositionX(23)-0.5*293,
				constants->SupportModulePositionY(23)-box[1]-z1,
				constants->SupportModulePositionZ(23)-40);
	TGeoTranslation *c3 = new TGeoTranslation ("c3",
				constants->SupportModulePositionX(24)-0.5*293,
				constants->SupportModulePositionY(24)-box[1]-z1,
				constants->SupportModulePositionZ(25)-40);
	TGeoTranslation *c4 = new TGeoTranslation ("c4",
				constants->SupportModulePositionX(27)-0.5*293,
				constants->SupportModulePositionY(27)-box[1]-z1,
				constants->SupportModulePositionZ(28)-40);
	upFace->AddNode(aCORDE00,57,c1);
	upFace->AddNode(aCORDE00,58,c2);
	upFace->AddNode(aCORDE00,59,c3);
	upFace->AddNode(aCORDE00,60,c4);


	// Construct Bars for lateral supports (up-face)

	TGeoTranslation *aco00=new TGeoTranslation("aco00",
				constants->SupportModulePositionX(20)+0.5*293-dy,
				constants->SupportModulePositionY(20)-box[1]-z1,
				constants->SupportModulePositionZ(20)-40);
	upFace->AddNode(aCORDE00,61,aco00);

	TGeoTranslation *aco00q1=new TGeoTranslation("aco00q1",
				constants->SupportModulePositionX(23)+0.5*293-dy,
				constants->SupportModulePositionY(23)-box[1]-z1,
				constants->SupportModulePositionZ(23)-40);
	upFace->AddNode(aCORDE00,62,aco00q1);

	TGeoTranslation *aco00q2=new TGeoTranslation("aco00q2",
				constants->SupportModulePositionX(24)+0.5*293-dy,
				constants->SupportModulePositionY(24)-box[1]-z1,
				constants->SupportModulePositionZ(25)-40);
	upFace->AddNode(aCORDE00,63,aco00q2);

	TGeoTranslation *aco00q3=new TGeoTranslation("aco00q3",
				constants->SupportModulePositionX(27)+0.5*293-dy,
				constants->SupportModulePositionY(27)-box[1]-z1,
				constants->SupportModulePositionZ(28)-40);
	upFace->AddNode(aCORDE00,64,aco00q3);


	TGeoTranslation *aco01=new TGeoTranslation("aco01",
				constants->SupportModulePositionX(30)-0.5*293+dy,
				constants->SupportModulePositionY(30)-box[1]-z1,
				constants->SupportModulePositionZ(30)-40);
	upFace->AddNode(aCORDE00,65,aco01);

	TGeoTranslation *aco01q1=new TGeoTranslation("aco01q1",
				constants->SupportModulePositionX(33)-0.5*293+dy,
				constants->SupportModulePositionY(33)-box[1]-z1,
				constants->SupportModulePositionZ(33)-40);
	upFace->AddNode(aCORDE00,66,aco01q1);

	TGeoTranslation *aco01q2=new TGeoTranslation("aco01q2",
				constants->SupportModulePositionX(34)-0.5*293+dy,
				constants->SupportModulePositionY(34)-box[1]-z1,
				constants->SupportModulePositionZ(35)-40);
	upFace->AddNode(aCORDE00,67,aco01q2);

	TGeoTranslation *aco01q3=new TGeoTranslation("aco01q3",
				constants->SupportModulePositionX(37)-0.5*293+dy,
				constants->SupportModulePositionY(37)-box[1]-z1,
				constants->SupportModulePositionZ(38)-40);
	upFace->AddNode(aCORDE00,68,aco01q3);



	// Acorde's support bars (side's faces)

	//*** In Face ***

//	box[0]=39;
	box[0]=27;
	box[1]=5;
	box[2]=5;
	Float_t kro=3;
	Float_t q1=0;
	Float_t posx=constants->SupportModulePositionX(0)+0.5*293*0.7071-56*0.7071-18;
	Float_t posy=constants->SupportModulePositionY(0)-0.5*293*0.7071-56*0.7071+3-q1+kro;
	Float_t posz=constants->SupportModulePositionZ(0);

	TGeoBBox *acorde7 = new TGeoBBox("acorde7",box[0],box[1],box[2]);

	TGeoVolume *aCORDE7 = new TGeoVolume("ACORDE7",acorde7,al);

	TGeoCombiTrans *aco7 = new TGeoCombiTrans("aco7",posx,posy,posz-4*dy,idrotm231);
	TGeoCombiTrans *aco7q1 = new TGeoCombiTrans("aco7q1",posx,posy,
					constants->ModulePositionZ(3)-4*dy,idrotm231);
	TGeoCombiTrans *aco7q2 = new TGeoCombiTrans("aco7q2",posx,posy,
					constants->ModulePositionZ(5)-4*dy,idrotm231);
	TGeoCombiTrans *aco7q3 = new TGeoCombiTrans("aco7q3",posx,posy,
					constants->ModulePositionZ(8)-4*dy,idrotm231);

	inFace->AddNode(aCORDE7,20,aco7);
	inFace->AddNode(aCORDE7,21,aco7q1);
	inFace->AddNode(aCORDE7,22,aco7q2);
	inFace->AddNode(aCORDE7,23,aco7q3);


	count=24;
	for(Int_t dyA=0;dyA<=4;dyA++)
	{

		Float_t posx1=constants->SupportModulePositionX(dyA)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->SupportModulePositionY(dyA)-0.1*293*0.7071-56*0.7071+3-des-q1+kro;
		Float_t posza=constants->SupportModulePositionZ(dyA);
		Float_t posx2=constants->SupportModulePositionX(dyA)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->SupportModulePositionY(dyA)+0.27*293*0.7071-56*0.7071+3-des-q1+kro;
		TGeoCombiTrans *aco7q4 = new TGeoCombiTrans("aco7q4",posx1,posy1,posza,idrotm231);
		TGeoCombiTrans *aco7q5 = new TGeoCombiTrans("aco7q5",posx2,posy2,posza,idrotm231);
		inFace->AddNode(aCORDE7,count,aco7q4);
		inFace->AddNode(aCORDE7,count+1,aco7q5);
		count=count+2;
		dyA++;
	}	


	count=34;
	for(Int_t dyb=5;dyb<=9;dyb++)
	{

		Float_t posx1=constants->SupportModulePositionX(dyb)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->SupportModulePositionY(dyb)-0.1*293*0.7071-56*0.7071+3-des-q1+kro;
		Float_t poszb=constants->SupportModulePositionZ(dyb+10);
		Float_t posx2=constants->SupportModulePositionX(dyb)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->SupportModulePositionY(dyb)+0.27*293*0.7071-56*0.7071+3-des-q1+kro;
		TGeoCombiTrans *aco7q6 = new TGeoCombiTrans("aco7q6",posx1,posy1,poszb,idrotm231);
		TGeoCombiTrans *aco7q7 = new TGeoCombiTrans("aco7q7",posx2,posy2,poszb,idrotm231);
		inFace->AddNode(aCORDE7,count,aco7q6);
		inFace->AddNode(aCORDE7,count+1,aco7q7);
		count=count+2;
		dyb++;
	}	



	Float_t posxq1=constants->SupportModulePositionX(10)+0.5*293*0.7071-56*0.7071-18;
	Float_t posyq1=constants->SupportModulePositionY(10)-0.5*293*0.7071-56*0.7071+3-q1+kro;
	Float_t poszq1=constants->SupportModulePositionZ(10);
	TGeoCombiTrans *aco7q8 = new TGeoCombiTrans("aco7q8",posxq1,posyq1,poszq1-4*dy,idrotm231);
	TGeoCombiTrans *aco7q9 = new TGeoCombiTrans("aco7q9",posxq1,posyq1,
					constants->ModulePositionZ(13)-4*dy,idrotm231);
	TGeoCombiTrans *aco7q10 = new TGeoCombiTrans("aco7q10",posxq1,posyq1,
					constants->ModulePositionZ(15)-4*dy,idrotm231);
	TGeoCombiTrans *aco7q11 = new TGeoCombiTrans("aco7q11",posxq1,posyq1,
					constants->ModulePositionZ(18)-4*dy,idrotm231);
	inFace->AddNode(aCORDE7,44,aco7q8);
	inFace->AddNode(aCORDE7,45,aco7q9);
	inFace->AddNode(aCORDE7,46,aco7q10);
	inFace->AddNode(aCORDE7,47,aco7q11);


	count=48;
	for(Int_t dyc=10;dyc<=14;dyc++)

	{

		Float_t posx1=constants->SupportModulePositionX(dyc)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->SupportModulePositionY(dyc)-0.1*293*0.7071-56*0.7071+3-des-0.8+kro;
		Float_t poszc=constants->SupportModulePositionZ(dyc);
		Float_t posx2=constants->SupportModulePositionX(dyc)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->SupportModulePositionY(dyc)+0.27*293*0.7071-56*0.7071+3-des-1.5-0.8+kro;
		TGeoRotation *rot1 = new TGeoRotation();
		rot1->RotateZ(70);
		TGeoCombiTrans *aco7q12 = new TGeoCombiTrans("aco7q12",posx1,posy1,poszc,idrotm231);
		TGeoCombiTrans *aco7q13 = new TGeoCombiTrans("aco7q13",posx2+15,posy2-10,poszc,rot1);
		inFace->AddNode(aCORDE7,count,aco7q12);
		inFace->AddNode(aCORDE7,count+1,aco7q13);// bars 25 grades
		count=count+2;
		dyc++;
	}


	count=57;
	for(Int_t dyd=15;dyd<=19;dyd++)

	{

		Float_t posx1=constants->SupportModulePositionX(dyd)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->SupportModulePositionY(dyd)-0.1*293*0.7071-56*0.7071+3-des-q1-0.8+kro;
		Float_t poszd=constants->SupportModulePositionZ(dyd);
		Float_t posx2=constants->SupportModulePositionX(dyd)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->SupportModulePositionY(dyd)+0.27*293*0.7071-56*0.7071+3-des-1.5-q1-0.8+kro;
		TGeoRotation *rot1 = new TGeoRotation();
		rot1->RotateZ(70);
		TGeoCombiTrans *aco7q14 = new TGeoCombiTrans("aco7q14",posx1,posy1,poszd,idrotm231);
		TGeoCombiTrans *aco7q15 = new TGeoCombiTrans("aco7q15",posx2+15,posy2-10,poszd,rot1);
		inFace->AddNode(aCORDE7,count,aco7q14);
		inFace->AddNode(aCORDE7,count+1,aco7q15);// bars 25 grades
		count=count+2;
		dyd++;
	}


	//*** Out Face ***

//	box[0]=39;
	box[0]=27;
	box[1]=5;
	box[2]=5;
	Float_t s1=2.5;
	Float_t posxqa=constants->SupportModulePositionX(50)-0.5*293*0.7071+56*0.7071+18;
	Float_t posyqa=constants->SupportModulePositionY(50)-0.5*293*0.7071-56*0.7071+3-s1+kro;
	Float_t poszqa=constants->SupportModulePositionZ(50);
	TGeoCombiTrans *aco7q16 = new TGeoCombiTrans("aco7q16",
					posxqa,posyqa,poszqa-4*dy,idrotm232);
	TGeoCombiTrans *aco7q17 = new TGeoCombiTrans("aco7q17",
					posxqa,posyqa,
					constants->ModulePositionZ(43)-4*dy,idrotm232);
	TGeoCombiTrans *aco7q18 = new TGeoCombiTrans("aco7q18",posxqa,posyqa,
					constants->ModulePositionZ(55)-4*dy,idrotm232);
	TGeoCombiTrans *aco7q19 = new TGeoCombiTrans("aco7q19",posxqa,posyqa,
					constants->ModulePositionZ(58)-4*dy,idrotm232);
	TGeoCombiTrans *aco7q20 = new TGeoCombiTrans("aco7q20",
					constants->SupportModulePositionX(50)-0.1*293*0.7071
					+56*0.7071+18-des,
					constants->SupportModulePositionY
					(50)-0.1*293*0.7071-56*0.7071+3-des-s1,
					constants->SupportModulePositionZ(45),idrotm232);
	TGeoCombiTrans *aco7q21 = new TGeoCombiTrans("aco7q21",
					constants->SupportModulePositionX(50)+0.27*293*0.7071
					+56*0.7071+18-des,
					constants->SupportModulePositionY(50)
					+0.27*293*0.7071-56*0.7071+3-des-s1,
					constants->SupportModulePositionZ(45),idrotm232);
	outFace->AddNode(aCORDE7,19,aco7q16);
	outFace->AddNode(aCORDE7,20,aco7q17);
	outFace->AddNode(aCORDE7,21,aco7q18);
	outFace->AddNode(aCORDE7,22,aco7q19);
	outFace->AddNode(aCORDE7,23,aco7q20);
	outFace->AddNode(aCORDE7,24,aco7q21);


	count=25;
	for(Int_t dye=50;dye<=54;dye++)
	{

		Float_t posx1=constants->SupportModulePositionX(dye)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->SupportModulePositionY(dye)-0.1*293*0.7071-56*0.7071+3-des-s1+kro;
		Float_t posze=constants->SupportModulePositionZ(dye);
		Float_t posx2=constants->SupportModulePositionX(dye)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->SupportModulePositionY(dye)+0.27*293*0.7071-56*0.7071+3-des-s1+kro;
		TGeoCombiTrans *aco7q22 = new TGeoCombiTrans("aco7q22",posx1,posy1,posze,idrotm232);
		TGeoCombiTrans *aco7q23 = new TGeoCombiTrans("aco7q23",posx2,posy2,posze,idrotm232);
		outFace->AddNode(aCORDE7,count,aco7q22);
		outFace->AddNode(aCORDE7,count+1,aco7q23);
		count=count+2;
		dye++;
	}


	count=35;
	for(Int_t dyf=57;dyf<=59;dyf++)
	{

		Float_t posx1=constants->SupportModulePositionX(dyf)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->SupportModulePositionY(dyf)-0.1*293*0.7071-56*0.7071+3-des-s1+kro;
		Float_t poszf=constants->SupportModulePositionZ(dyf-10);
		Float_t posx2=constants->SupportModulePositionX(dyf)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->SupportModulePositionY(dyf)+0.27*293*0.7071-56*0.7071+3-des-s1+kro;
		TGeoCombiTrans *aco7q24 = new TGeoCombiTrans("aco7q24",posx1,posy1,poszf,idrotm232);
		TGeoCombiTrans *aco7q25 = new TGeoCombiTrans("aco7q25",posx2,posy2,poszf,idrotm232);
		outFace->AddNode(aCORDE7,count,aco7q24);
		outFace->AddNode(aCORDE7,count+1,aco7q25);
		count=count+2;
		dyf++;
	}


	Float_t posxqb=constants->SupportModulePositionX(40)-0.5*293*0.7071+56*0.7071+18;
	Float_t posyqb=constants->SupportModulePositionY(40)-0.5*293*0.7071-56*0.7071+3-s1+kro;
	Float_t poszqb=constants->SupportModulePositionZ(40);
	TGeoCombiTrans *aco7q26 = new TGeoCombiTrans("aco7q26",
					posxqb,posyqb,poszqb-4*dy,idrotm232);
	TGeoCombiTrans *aco7q27 = new TGeoCombiTrans("aco7q27",
					posxqb,posyqb,
					constants->SupportModulePositionZ(43)-4*dy,idrotm232);
	TGeoCombiTrans *aco7q28 = new TGeoCombiTrans("aco7q28",
					posxqb,posyqb,
					constants->SupportModulePositionZ(45)-4*dy,idrotm232);
	TGeoCombiTrans *aco7q29 = new TGeoCombiTrans("aco7q29",posxqb,posyqb,
					constants->SupportModulePositionZ(48)-4*dy,idrotm232);
	outFace->AddNode(aCORDE7,41,aco7q26);
	outFace->AddNode(aCORDE7,42,aco7q27);
	outFace->AddNode(aCORDE7,43,aco7q28);
	outFace->AddNode(aCORDE7,44,aco7q29);

	count=45;
	for(Int_t dyg=40;dyg<=44;dyg++)
	{

		Float_t posx1=constants->SupportModulePositionX(dyg)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->SupportModulePositionY(dyg)-0.1*293*0.7071-56*0.7071+3-des-s1+kro;
		Float_t poszg=constants->SupportModulePositionZ(dyg);
		Float_t posx2=constants->SupportModulePositionX(dyg)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->SupportModulePositionY(dyg)+0.27*293*0.7071-56*0.7071+3-des-s1+kro;
		TGeoRotation *rot1 = new TGeoRotation();
		rot1->RotateZ(105);
		TGeoCombiTrans *aco7q30 = new TGeoCombiTrans("aco7q30",posx1,posy1,poszg,idrotm232);
		TGeoCombiTrans *aco7q31 = new TGeoCombiTrans("aco7q31",posx2-15,posy2-10,poszg,rot1);
		outFace->AddNode(aCORDE7,count,aco7q30);
		outFace->AddNode(aCORDE7,count+1,aco7q31);// bars 25 grades
		count=count+2;
		dyg++;
	}


	count=55;
	for(Int_t dyh=45;dyh<=49;dyh++)
	{

		Float_t posx1=constants->SupportModulePositionX(dyh)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->SupportModulePositionY(dyh)-0.1*293*0.7071-56*0.7071+3-des-s1+kro;
		Float_t poszh=constants->SupportModulePositionZ(dyh);
		Float_t posx2=constants->SupportModulePositionX(dyh)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->SupportModulePositionY(dyh)+0.27*293*0.7071-56*0.7071+3-des-s1+kro;
		TGeoRotation *rot1 = new TGeoRotation();
		rot1->RotateZ(105);
		TGeoCombiTrans *aco7q32 = new TGeoCombiTrans("aco7q32",posx1,posy1,poszh,idrotm232);
		TGeoCombiTrans *aco7q33 = new TGeoCombiTrans("aco7q33",posx2-15,posy2-10,poszh,rot1);
		outFace->AddNode(aCORDE7,count,aco7q32);
		outFace->AddNode(aCORDE7,count+1,aco7q33);// bars 25 grades
		count=count+2;
		dyh++;
	}



	// Set the bars non perpendicular at side faces

	//*** In-Face ***

	box[0]=5;
//	box[1]=55.15;
	box[1]=40;
	box[2]=5;
	Float_t sm=2;
	Float_t re=1;
	Float_t posx1=constants->SupportModulePositionX(0)+0.5*293*0.7071-4*box[0]-8+re;
	Float_t posy1=constants->SupportModulePositionY(0)-0.5*293*0.7071-box[1]-18-2+sm;
	Float_t posz1=constants->SupportModulePositionZ(0);

	TGeoBBox *acorde7q1 = new TGeoBBox("acorde7q1",box[0],box[1],box[2]);

	TGeoVolume *aCORDE7q1 = new TGeoVolume("ACORDE7_1",acorde7q1,al);
	TGeoTranslation *aco71 = new TGeoTranslation("aco71",posx1,posy1,posz1-4*dy);
	TGeoTranslation *aco72 = new TGeoTranslation("aco72",posx1,posy1,
					constants->SupportModulePositionZ(3)-4*dy);
	TGeoTranslation *aco73 = new TGeoTranslation("aco73",posx1,posy1,
					constants->SupportModulePositionZ(5)-4*dy);
	TGeoTranslation *aco74 = new TGeoTranslation("aco74",posx1,posy1,
					constants->SupportModulePositionZ(8)-4*dy);
	inFace->AddNode(aCORDE7q1,67,aco71);
	inFace->AddNode(aCORDE7q1,68,aco72);
	inFace->AddNode(aCORDE7q1,69,aco73);
	inFace->AddNode(aCORDE7q1,70,aco74);


	count=71;
	for(Int_t dyi=0;dyi<=4;dyi++)
	{

		Float_t posx1a=constants->SupportModulePositionX(dyi)+0.1*293*0.7071-4*box[0]-8+des+re;
		Float_t posy1a=constants->SupportModulePositionY(dyi)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1a=constants->SupportModulePositionZ(dyi);
		Float_t dyx2=constants->SupportModulePositionX(dyi)-0.27*293*0.7071-4*box[0]-8+des+re;
		Float_t dyy2=constants->SupportModulePositionY(dyi)+0.27*293*0.7071-box[1]-18-2-des+sm;
		TGeoTranslation *aco75=new TGeoTranslation("aco75",posx1a,posy1a,posz1a);
		TGeoTranslation *aco76=new TGeoTranslation("aco76",dyx2,dyy2,posz1a);
		inFace->AddNode(aCORDE7q1,count,aco75);
		inFace->AddNode(aCORDE7q1,count+1,aco76);
		count=count+2;
		dyi++;
	}


	count=81;
	for(Int_t dyj=5;dyj<=9;dyj++)
	{

		Float_t posx1b=constants->SupportModulePositionX(dyj)+0.1*293*0.7071-4*box[0]-8+des+re;
		Float_t posy1b=constants->SupportModulePositionY(dyj)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1b=constants->SupportModulePositionZ(dyj+10);
		Float_t dyx2=constants->SupportModulePositionX(dyj)-0.27*293*0.7071-4*box[0]-8+des+re;
		Float_t dyy2=constants->SupportModulePositionY(dyj)+0.27*293*0.7071-box[1]-18-2-des+sm;
		TGeoTranslation *aco75=new TGeoTranslation("aco75",posx1b,posy1b,posz1b);
		TGeoTranslation *aco76=new TGeoTranslation("aco76",dyx2,dyy2,posz1b);
		inFace->AddNode(aCORDE7q1,count,aco75);
		inFace->AddNode(aCORDE7q1,count+1,aco76);
		count=count+2;
		dyj++;
	}


	Float_t posx1q1=constants->SupportModulePositionX(10)+0.5*293*0.7071-4*box[0]-8+re;
	Float_t posy1q1=constants->SupportModulePositionY(10)-0.5*293*0.7071-box[1]-18-2+sm;
	Float_t posz1q1=constants->SupportModulePositionZ(10);
	TGeoTranslation *aco77=new TGeoTranslation("aco77",posx1q1,posy1q1,posz1q1-4*dy);
	TGeoTranslation *aco78=new TGeoTranslation("aco78",posx1q1,posy1q1,
					constants->SupportModulePositionZ(13)-4*dy);

	TGeoTranslation *aco79=new TGeoTranslation("aco79",posx1q1,posy1q1,
					constants->SupportModulePositionZ(15)-4*dy);
	TGeoTranslation *aco710=new TGeoTranslation("aco710",posx1q1,posy1q1,
					constants->SupportModulePositionZ(18)-4*dy);
	inFace->AddNode(aCORDE7q1,91,aco77);
	inFace->AddNode(aCORDE7q1,92,aco78);
	inFace->AddNode(aCORDE7q1,93,aco79);
	inFace->AddNode(aCORDE7q1,94,aco710);

	count=95;
	for(Int_t dyk=10;dyk<=14;dyk++)
	{

		Float_t posx1c=constants->SupportModulePositionX(dyk)+0.1*293*0.7071-4*box[0]-8+des+re+.83;
		Float_t posy1c=constants->SupportModulePositionY(dyk)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1c=constants->SupportModulePositionZ(dyk);
		Float_t dyx2=constants->SupportModulePositionX(dyk)-0.27*293*0.7071-4*box[0]-4+des+re+0.83;
		Float_t dyy2=constants->SupportModulePositionY(dyk)+0.27*293*0.7071-box[1]-18-5-des+sm;
		TGeoTranslation *aco711=new TGeoTranslation("aco711",posx1c,posy1c,posz1c);
		TGeoTranslation *aco712=new TGeoTranslation("aco712",dyx2,dyy2,posz1c);
		inFace->AddNode(aCORDE7q1,count,aco711);
		inFace->AddNode(aCORDE7q1,count+1,aco712);
		count=count+2;
		dyk++;
	}



	count=105;
	for(Int_t dyl=15;dyl<=19;dyl++)
	{

		Float_t posx1d=constants->SupportModulePositionX(dyl)+0.1*293*0.7071-4*box[0]-8+des+re+0.83;
		Float_t posy1d=constants->SupportModulePositionY(dyl)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1d=constants->SupportModulePositionZ(dyl);
		Float_t dyx2=constants->SupportModulePositionX(dyl)-0.27*293*0.7071-4*box[0]-4+des+re+0.83;
		Float_t dyy2=constants->SupportModulePositionY(dyl)+0.27*293*0.7071-box[1]-18-5-des;
		TGeoTranslation *aco713=new TGeoTranslation("aco713",posx1d,posy1d,posz1d);
		TGeoTranslation *aco714=new TGeoTranslation("aco714",dyx2,dyy2,posz1d);
		inFace->AddNode(aCORDE7q1,count,aco713);
		inFace->AddNode(aCORDE7q1,count+1,aco714);
		count=count+2;
		dyl++;
	}

		//*** Out-Face ***

	Float_t posx1qa=constants->SupportModulePositionX(50)-0.5*293*0.7071+4*box[0]+8-re-1;
	Float_t posy1qa=constants->SupportModulePositionY(50)-0.5*293*0.7071-box[1]-18-2+sm-2.5;
	Float_t posz1qa=constants->SupportModulePositionZ(50);
	TGeoTranslation *aco715=new TGeoTranslation("aco715",posx1qa,posy1qa,posz1qa-4*dy);
	TGeoTranslation *aco716=new TGeoTranslation("aco716",posx1qa,posy1qa,
				constants->SupportModulePositionZ(43)-4*dy);
	TGeoTranslation *aco717=new TGeoTranslation("aco717",posx1qa,posy1qa,
				constants->SupportModulePositionZ(55)-4*dy);
	TGeoTranslation *aco718=new TGeoTranslation("aco718",posx1qa,posy1qa,
				constants->SupportModulePositionZ(58)-4*dy);
	TGeoTranslation *aco719=new TGeoTranslation("aco719",
				constants->SupportModulePositionX(50)-0.1*293*0.7071+4*box[0]+8-des-re-1,		
				constants->SupportModulePositionY(50)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5,
				constants->SupportModulePositionZ(45));
	TGeoTranslation *aco720=new TGeoTranslation("aco720",
				constants->SupportModulePositionX(50)+0.27*293*0.7071+4*box[0]+8-des-re-1,
				constants->SupportModulePositionY(50)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5,
				constants->SupportModulePositionZ(45));


	outFace->AddNode(aCORDE7q1,115,aco715);
	outFace->AddNode(aCORDE7q1,116,aco716);
	outFace->AddNode(aCORDE7q1,117,aco717);
	outFace->AddNode(aCORDE7q1,118,aco718);
	outFace->AddNode(aCORDE7q1,119,aco719);
	outFace->AddNode(aCORDE7q1,120,aco720);




	count=65;
	for(Int_t dym=50;dym<=54;dym++)
	{

		Float_t posx1e=constants->SupportModulePositionX(dym)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1e=constants->SupportModulePositionY(dym)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1e=constants->SupportModulePositionZ(dym);
		Float_t dyx2=constants->SupportModulePositionX(dym)+0.27*293*0.7071+4*box[0]+8-des-re-1;
		Float_t dyy2=constants->SupportModulePositionY(dym)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5;
		TGeoTranslation *aco721=new TGeoTranslation("aco721",posx1e,posy1e,posz1e);
		TGeoTranslation *aco722=new TGeoTranslation("aco722",dyx2,dyy2,posz1e);
		outFace->AddNode(aCORDE7q1,count,aco721);
		outFace->AddNode(aCORDE7q1,count+1,aco722);
		count=count+2;
		dym++;
	}



	count=75;
	for(Int_t dyn=57;dyn<=59;dyn++)
	{

		Float_t posx1f=constants->SupportModulePositionX(dyn)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1f=constants->SupportModulePositionY(dyn)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1f=constants->SupportModulePositionZ(dyn-10);
		Float_t dyx2=constants->SupportModulePositionX(dyn)+0.27*293*0.7071+4*box[0]+8-des-re-1;
		Float_t dyy2=constants->SupportModulePositionY(dyn)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5;
		TGeoTranslation *aco723=new TGeoTranslation("aco723",posx1f,posy1f,posz1f);
		TGeoTranslation *aco724=new TGeoTranslation("aco724",dyx2,dyy2,posz1f);
		outFace->AddNode(aCORDE7q1,count,aco723);
		outFace->AddNode(aCORDE7q1,count+1,aco724);
		count=count+2;
		dyn++;
	}


	Float_t posx1qb=constants->SupportModulePositionX(40)-0.5*293*0.7071+4*box[0]+5;
	Float_t posy1qb=constants->SupportModulePositionY(40)-0.5*293*0.7071-box[1]-18-2;
	Float_t posz1qb=constants->SupportModulePositionZ(40);
	TGeoTranslation *aco725=new TGeoTranslation("aco725",posx1qb,posy1qb,posz1qb-4*dy);
	TGeoTranslation *aco726=new TGeoTranslation("aco726",posx1qb,posy1qb,
				constants->SupportModulePositionZ(43)-4*dy);
	TGeoTranslation *aco727=new TGeoTranslation("aco727",posx1qb,posy1qb,
				constants->SupportModulePositionZ(45)-4*dy);
	TGeoTranslation *aco728=new TGeoTranslation("aco728",posx1qb,posy1qb,
				constants->SupportModulePositionZ(48)-4*dy);
	outFace->AddNode(aCORDE7q1,85,aco725);
	outFace->AddNode(aCORDE7q1,86,aco726);
	outFace->AddNode(aCORDE7q1,87,aco727);
	outFace->AddNode(aCORDE7q1,88,aco728);



	count=89;
	for(Int_t dyo=40;dyo<=44;dyo++)
	{

		Float_t posx1g=constants->SupportModulePositionX(dyo)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1g=constants->SupportModulePositionY(dyo)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1g=constants->SupportModulePositionZ(dyo);
		Float_t dyx2=constants->SupportModulePositionX(dyo)+0.27*293*0.7071+4*box[0]+4-des-re-1+2.8;
		Float_t dyy2=constants->SupportModulePositionY(dyo)+0.27*293*0.7071-box[1]-18-5-des+sm-2.5+3;
		TGeoTranslation *aco729=new TGeoTranslation("aco729",posx1g,posy1g,posz1g);
		TGeoTranslation *aco730=new TGeoTranslation("aco730",dyx2,dyy2,posz1g);
		outFace->AddNode(aCORDE7q1,count,aco729);
		outFace->AddNode(aCORDE7q1,count+1,aco730);
		count=count+2;
		dyo++;
	}



	count=99;
	for(Int_t dyp=45;dyp<=49;dyp++)
	{

		Float_t posx1h=constants->SupportModulePositionX(dyp)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1h=constants->SupportModulePositionY(dyp)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1h=constants->SupportModulePositionZ(dyp);
		Float_t dyx2=constants->SupportModulePositionX(dyp)+0.27*293*0.7071+4*box[0]+4-des-re-1+2.8;
		Float_t dyy2=constants->SupportModulePositionY(dyp)+0.27*293*0.7071-box[1]-18-5-des+sm-2.5+3;
		TGeoTranslation *aco729=new TGeoTranslation("aco729",posx1h,posy1h,posz1h);
		TGeoTranslation *aco730=new TGeoTranslation("aco730",dyx2,dyy2,posz1h);
		outFace->AddNode(aCORDE7q1,count,aco729);
		outFace->AddNode(aCORDE7q1,count+1,aco730);
		count=count+2;
		dyp++;
	}


	// Here I define a master volume "ACORDE_SUPPORT" for Acorde's support

	//---> Set the support of ACORDE alice MODULES


	Float_t dy1=20;
	box[0]=10;
	box[1]=0.5;
	box[2]=500;

	Float_t sx=constants->SupportModulePositionX(24)-0.5*293;
	Float_t sy=constants->SupportModulePositionY(24)-box[1];
	Float_t sz=0;
	Float_t sx2=constants->SupportModulePositionX(24)+0.5*293-dy;
	Float_t sy2=constants->SupportModulePositionY(24)-box[1];
	Float_t sx4=constants->SupportModulePositionX(24)-0.5*293+dy2;
	Float_t sy4=constants->SupportModulePositionY(24)-box[1];
	Float_t sx5=constants->SupportModulePositionX(24)+0.5*293-dy2;
	Float_t sy5=constants->SupportModulePositionY(24)-box[1];

	Float_t dyx=constants->SupportModulePositionX(4)+0.5*293*0.7071-box[0];
	Float_t dyy=constants->SupportModulePositionY(4)-0.5*293*0.7071-box[1];
	Float_t dyz=0;
	Float_t dyx1=constants->SupportModulePositionX(4)+0.1*293*0.7071-box[0];
	Float_t dyy1=constants->SupportModulePositionY(4)-0.1*293*0.7071-box[1];
	Float_t dyx2=constants->SupportModulePositionX(4)-0.27*293*0.7071-box[0];
	Float_t dyy2=constants->SupportModulePositionY(4)+0.27*293*0.7071-box[1];


	Float_t dx1=constants->SupportModulePositionX(14)+0.5*293*0.7071-box[0];
	Float_t dy11=constants->SupportModulePositionY(14)-0.5*293*0.7071-box[1];
	Float_t dyx11=constants->SupportModulePositionX(14)+0.1*293*0.7071-box[0];
	Float_t dyy11=constants->SupportModulePositionY(14)-0.1*293*0.7071-box[1];
	Float_t dyx21=constants->SupportModulePositionX(14)-0.27*293*0.7071-box[0];
	Float_t dyy21=constants->SupportModulePositionY(14)+0.27*293*0.7071-box[1];


	Float_t tbox[3];
	tbox[0]=1;
	tbox[1]=7;
	tbox[2]=500;

	TGeoVolume *support = new TGeoVolumeAssembly("ACORDE_SUPPORT");

	TGeoBBox *acorde8 = new TGeoBBox("acorde8",box[0],box[1],box[2]);
	TGeoVolume *aCORDE8 = new TGeoVolume("ACORDE8",acorde8,al);

	TGeoBBox *acorde9 = new TGeoBBox("acorde9",tbox[0],tbox[1],tbox[2]);
	TGeoVolume *aCORDE9 = new TGeoVolume("ACORDE9",acorde9,al);

	support->AddNode(aCORDE8,1,new TGeoTranslation(0,-5,0));
	support->AddNode(aCORDE8,2,new TGeoTranslation(0,-dy1,0));
	support->AddNode(aCORDE9,3,new TGeoTranslation(0,-tbox[1]-5.5,0));


	// Put "support" on Up-Face

	upFace->AddNode(support,69,new TGeoTranslation("aco8",sx,sy,sz));
	upFace->AddNode(support,70,new TGeoTranslation("aco8_2",sx2,sy2,sz));
	upFace->AddNode(support,71,new TGeoTranslation("aco8_4",sx4,sy4,sz));
	upFace->AddNode(support,72,new TGeoTranslation("aco8_6",sx5,sy5,sz));
	upFace->AddNode(support,73,new TGeoTranslation("aco8_2",-sx2,sy2,sz));
	upFace->AddNode(support,74,new TGeoTranslation("aco8_4",-sx4,sy4,sz));
	upFace->AddNode(support,75,new TGeoTranslation("aco8_6",-sx5,sy5,sz));

	// Put "support" on In-Face
	Float_t ms = 1.3;
	inFace->AddNode(support,121,new TGeoCombiTrans("aco8_81",dyx,dyy+ms,dyz,idrotm232));
	inFace->AddNode(support,122,new TGeoCombiTrans("aco8_121",dyx1+des,ms+dyy1-des,dyz,idrotm232));
	inFace->AddNode(support,123,new TGeoCombiTrans("aco8_161",dyx2+des,ms+dyy2-des,dyz,idrotm232));
	inFace->AddNode(support,124,new TGeoCombiTrans("aco8_82",dx1,ms+dy11,dyz,idrotm232));
	inFace->AddNode(support,125,new TGeoCombiTrans("aco8_122",dyx11+des,ms+dyy11-des,dyz,idrotm232));
	inFace->AddNode(support,126,new TGeoCombiTrans("aco8_162",dyx21+des,ms+dyy21-des,dyz,idrotm232));

	// Put "support" on Out-Face

	outFace->AddNode(support,121,new TGeoCombiTrans("aco8_81",-dyx,dyy+ms,dyz,idrotm231));
	outFace->AddNode(support,122,new TGeoCombiTrans("aco8_121",-dyx1-des,ms+dyy1-des,dyz,idrotm231));
	outFace->AddNode(support,123,new TGeoCombiTrans("aco8_161",-dyx2-des,ms+dyy2-des,dyz,idrotm231));
	outFace->AddNode(support,124,new TGeoCombiTrans("aco8_82",-dx1,dy11+ms,dyz,idrotm231));
	outFace->AddNode(support,125,new TGeoCombiTrans("aco8_122",-dyx11-des,ms+dyy11-des,dyz,idrotm231));
	outFace->AddNode(support,126,new TGeoCombiTrans("aco8_162",-dyx21-des,ms+dyy21-des,dyz,idrotm231));
	
	aCORDE->AddNode(inFace,1);//---> volume of supports & bars in-face
	aCORDE->AddNode(upFace,2);//---> volume of supports & bars up-face
	aCORDE->AddNode(outFace,3);//---> volume of supports & bars out-face
	aCORDE->AddNode(modules,4);//---> volume of ALL ACORDE's Modules
	alice->AddNode(aCORDE,1);//---> put volume of ACORDE over ALICE's volume



}
//__________________________________________________________________________

void AliACORDEv1::Init()
{
  // Initialise L3 magnet after it has been built
  Int_t i;
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" ACORDEv1_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    // Here the ACORDEv initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
 // AliACORDE::Init();  
}
//____________________________________________________________________________
void AliACORDEv1::StepManager()
{

  //
  // Called for every step in the Cosmic Ray Trigger
  //


  // volume: 
  //  [0] = module number 1-60 (1==>(0-0), 60 (5-9)
  //  [1] = Plastic number: 0 (down) to 1 (up)
  static Int_t   vol[2]; 
  //
  // hit
  // [0] = PID
  // [1-3] = x, y, z 
  // [4] = time 
  // [5-7] = px, py, pz
  // [8] = energy 
  // [9] = energy loss
  // [10] = length of track through plastic
  static Float_t hits[11];

  // local static variables
  static Float_t eloss;
  static Float_t step;
  // scintillator volume
 static Int_t idScint = gMC->VolId("ACORDE2");
  // local variables
  Int_t copy;
  TLorentzVector pos;
  TLorentzVector mom;

  // only charged tracks
  if ( !gMC->TrackCharge() || !gMC->IsTrackAlive() ) return;

  // only in sensitive material
  if (gMC->CurrentVolID(copy) == idScint) {

    step  += gMC->TrackStep();
    eloss += gMC->Edep();
    // set all hit variables except eloss which is resetted
    // set volume variables
    if (gMC->IsTrackEntering()) {
      eloss = 0.0;
      step = 0.0;
      gMC->TrackPosition(pos);
      gMC->TrackMomentum(mom);
      // hit
      // [0] = PID
      // [1-3] = x, y, z 
      // [4] = time 
      // [5-7] = px, py, pz
      // [8] = energy 
      // [9] = energy loss
      hits[0]  = (Float_t ) gMC->TrackPid(); 


      hits[1] = pos[0]; 
      hits[2] = pos[1]; 
      hits[3] = pos[2]; 
      hits[4] = gMC->TrackTime();
      hits[5] = mom[0];
      hits[6] = mom[1];
      hits[7] = mom[2];
      hits[8] = gMC->Etot();
      // volume: 
      //  [0] = module number 1-60 (1==>(0-0), 60 (5-9)
      //  [1] = Plastic number: 0 (down) to 1 (up)
      Int_t copyPlastic; // plastic: down=1, up=2
      Int_t copyModule; // module: 1-60
      gMC->CurrentVolID(copyPlastic);
      gMC->CurrentVolOffID(1, copyModule);
      // module
      vol[0] = copyModule;
      // plastic: 0 = down, 1 = up
      vol[1] = copyPlastic - 4 ; // !!!!!!!
    // vol[1] = copyPlastic;
    } // end if gMC->IsTrackEntering()

    // set hit[9] = total energy loss and book hit
    if( gMC->IsTrackExiting() || 
	gMC->IsTrackStop() || 
	gMC->IsTrackDisappeared()){
      hits[9] = eloss;
      hits[10] = step;
      eloss = 0.0;
      step = 0.0;
      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol, hits);
     }
  } 



}

//_____________________________________________________________________________
void AliACORDEv1::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add an ACORDE hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliACORDEhit(fIshunt,track,vol,hits);
}

//_____________________________________________________________________________
void AliACORDEv1::AddDigits(Int_t* track, Int_t module, Float_t time)
{
  
  // Adds Digit
  
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliACORDEdigit(track,module,time);
}
//_____________________________________________________________________________

void AliACORDEv1::MakeBranch(Option_t *option)
{
// Creates new branches in the current Root Tree
    
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  AliDebug(2,Form("fBufferSize = %d",fBufferSize));
  const char *cH = strstr(option,"H");
  if (fHits   && fLoader->TreeH() && cH) {
    fLoader->TreeH()->Branch(branchname,&fHits, fBufferSize);
    AliDebug(2,Form("Making Branch %s for hits",branchname));
  }     
  const char *cD = strstr(option,"D");
  if (fDigits   && fLoader->TreeD() && cD) {
    fLoader->TreeD()->Branch(branchname,&fDigits, fBufferSize);
    AliDebug(2,Form("Making Branch %s for digits",branchname));
  }  
}

//_____________________________________________________________________________
void AliACORDEv1::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  // 

	// The alignable volumes are only the "ACORDE_MODULE_%d"
	//
	//	Structure of ACORDE's Geometry
	//
	//	ALIC_1
	//	    |---> ACORDE_1	
	//			|----> ACORDE_1_1 (in_face) ---
	//			|----> ACORDE_2_2 (up_face)    |--> BARS&SUPPORTS
	//			|----> ACORDE_3_3 (out_face)---
	//			|----> ACORDE_MODULES_4        |--> ACORDE'S MODULES
	//		
	//
	//     Send comments to: Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>

	TString vpstr1 = "ALIC_1/ACORDE_1/ALL_ACORDE_MODULES_4/ACORDE_MODULE_";
	TString snstr1 = "ACORDE/Array";
	TString volpath, symname;
	for(Int_t dy=1; dy<61 ; dy++)
	{
		volpath = vpstr1;
		volpath += dy;
		symname = snstr1;
		symname += dy;
		if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data()))
	        AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", symname.Data(),volpath.Data()));
	}
}
