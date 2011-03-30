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
//  with all the detectors. It includes the last survey of 2009.	     //
//  It include geometry and hits (position and momentum)                     //
//                                                                           //
//   Author: Mario Rodriguez Cahuantzi, FCFM-BUAP, Puebla, Pue. Mexico       //
//                                                                           //
//                  Send comments to:                                        //
//									     //
//      Arturo Fernandez Tellez   	<afernand@fcfm.buap.mx>              //
//      Eleazar Cuautle Flores    	<ecuautle@nucleares.unam.mx>         //
//	Mario Rodriguez	Cahuantzi 	<mrodrigu@mail.cern.ch>		     //	
//									     //
//			Puebla, Pue. Mexico December 2007                    //
//									     //
//	Last Update: Nov. 17th 2009				             //
//	Mario Rodriguez Cahuantzi 					     //
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

	// Call the global constants for the Modules
	
	AliACORDEConstants* constants = AliACORDEConstants::Instance();

	// Get the Alice Volume

	TGeoVolume *alice = gGeoManager->GetVolume("ALIC");

	// Define some materials & medium

	//*** Aluminium ***

	TGeoMedium* aluminium = gGeoManager->GetMedium("ACORDE_ALU_C0");

	//*** Scintillator ***

	TGeoMedium* scintillator = gGeoManager->GetMedium("ACORDE_CPV scint.1");
	
	//Define the mother volume for the ACORDE detector

	TGeoVolume *acorde = new TGeoVolumeAssembly("ACORDE");


	// Define 2 main-daughter volumes for ACORDE
	
	TGeoVolume *supportBars = new TGeoVolumeAssembly("ACORDE_SUPPORTS_BARS");
	TGeoVolume *acordeModules = new TGeoVolumeAssembly("ALL_ACORDE_MODULES");

  
	// Define rotation Matrix for Side's faces in Alice
	
	TGeoRotation *idrotm231 = new TGeoRotation("idrotm231",90, 45, 90, 135, 0, 0);
	TGeoRotation *idrotm232 = new TGeoRotation("idrotm232",90, 315, 90, 45, 0, 0);

	// Begin the Geometry for ACORDE

	// *** Definition of ACORDE's Modules ***

	// Define Measures of ACORDE's Modules

	Float_t acoFrameBox1[3],acoFrameBox2[3];

	acoFrameBox1[0] = constants->ModuleLength()/2;
	acoFrameBox1[1] = constants->ModuleHeight()/2;
	acoFrameBox1[2] = 2.50;

	acoFrameBox2[0] = 20.0;
	acoFrameBox2[1] = constants->ModuleHeight()/2;
	acoFrameBox2[2] = 10.0;

	// Define Measures of Scintillators
	
	Float_t acoScinBox[3];	
	acoScinBox[0] = constants->PlasticLength()/2;
	acoScinBox[1] = constants->PlasticHeight()/2;
	acoScinBox[2] = constants->PlasticWidth()/2;


	// Create the Modules of ACORDE, 1 aluminium frame and two scintillator plastics

	//*** Aluminium frame ***

	TGeoBBox *acordeModFrameBoxL = new TGeoBBox("acordeModFrameBoxL",acoFrameBox1[0],acoFrameBox1[1],acoFrameBox1[2]);
	TGeoBBox *acordeModFrameBoxH = new TGeoBBox("acordeModFrameBoxH",acoFrameBox2[0],acoFrameBox2[1],acoFrameBox2[2]);
	TGeoVolume *acordeModFrameVolumeL = new TGeoVolume("ACORDEMODFRAMEVOLUMEL",acordeModFrameBoxL,aluminium);
	TGeoVolume *acordeModFrameVolumeH = new TGeoVolume("ACORDEMODFRAMEVOLUMEH",acordeModFrameBoxH,aluminium);


	//*** Scintillators ***

	TGeoBBox *acordeScintillatorBox = new TGeoBBox("acordeScintillatorBox",acoScinBox[0],acoScinBox[1],acoScinBox[2]);
	TGeoVolume *acordeScintillatorVolume = new TGeoVolume("ACORDESCINTILLATORMODULE",acordeScintillatorBox,scintillator);


	// Here I create a single ACORDE module and then we make 60 copies of it

	TGeoVolume *acordeSingleModule = new TGeoVolumeAssembly("ACORDE_MODULE");
	acordeSingleModule->AddNode(acordeModFrameVolumeL,1,new TGeoTranslation("acordeFrame_01",0,0,12.5));
	acordeSingleModule->AddNode(acordeModFrameVolumeL,2,new TGeoTranslation("acordeFrame_02",0,0,-12.5));
	acordeSingleModule->AddNode(acordeModFrameVolumeH,3,new TGeoTranslation("acordeFrame_03",130,0,0));
	acordeSingleModule->AddNode(acordeModFrameVolumeH,4,new TGeoTranslation("acordeFrame_04",-130,0,0));
	acordeSingleModule->AddNode(acordeScintillatorVolume,5, new TGeoTranslation("acordeScintillator_01",0,1,0));
	acordeSingleModule->AddNode(acordeScintillatorVolume,6, new TGeoTranslation("acordeScintillator_01",0,-1,0));

	// Put the Modules of In-Face
	
	for(Int_t iAcordeModule=1;iAcordeModule<9;iAcordeModule++)
	{

		Float_t posx = constants->CenterModulePositionX(iAcordeModule);
		Float_t posy = constants->CenterModulePositionY(iAcordeModule);
		Float_t posz = constants->CenterModulePositionZ(iAcordeModule);	
		
		acordeModules->AddNode(acordeSingleModule,iAcordeModule,
			new TGeoCombiTrans("aco01",posx,posy,posz,idrotm232));
	}

	for(Int_t iAcordeModule=10;iAcordeModule<20;iAcordeModule++)
	{
		Float_t posx = constants->CenterModulePositionX(iAcordeModule);
		Float_t posy = constants->CenterModulePositionY(iAcordeModule);
		Float_t posz = constants->CenterModulePositionZ(iAcordeModule);

		acordeModules->AddNode(acordeSingleModule,iAcordeModule,
			new TGeoCombiTrans("aco01",posx,posy,posz,idrotm232));
	}

	// Put he Modules of Up-Face

	for(Int_t iAcordeModule=20;iAcordeModule<40;iAcordeModule++)
	{
		Float_t posx = constants->CenterModulePositionX(iAcordeModule);
		Float_t posy = constants->CenterModulePositionY(iAcordeModule);
		Float_t posz = constants->CenterModulePositionZ(iAcordeModule);	

		acordeModules->AddNode(acordeSingleModule,iAcordeModule,new TGeoTranslation("aco01",posx,posy,posz));
	}

	// Put the Modules of Out-Face

	for(Int_t iAcordeModule=40;iAcordeModule<50;iAcordeModule++)
	{
		Float_t posx = constants->CenterModulePositionX(iAcordeModule);
		Float_t posy = constants->CenterModulePositionY(iAcordeModule);
		Float_t posz = constants->CenterModulePositionZ(iAcordeModule);	

		acordeModules->AddNode(acordeSingleModule,iAcordeModule,
			new TGeoCombiTrans("aco01",posx,posy,posz,idrotm231));
	}

	// Put the Modules of Out-Face

	for(Int_t iAcordeModule=51;iAcordeModule<59;iAcordeModule++)
	{
		Float_t posx = constants->CenterModulePositionX(iAcordeModule);
		Float_t posy = constants->CenterModulePositionY(iAcordeModule);
		Float_t posz = constants->CenterModulePositionZ(iAcordeModule);	
		acordeModules->AddNode(acordeSingleModule,iAcordeModule,
					new TGeoCombiTrans("aco01",posx,posy,posz,idrotm231));
	}

	// Put the 4-central modules (Old-ITS modules)
	if (Get4CentralModulesGeometry()) 
	{ 
		acordeModules->AddNode(acordeSingleModule,0,
			new TGeoTranslation("Mod0_0",constants->CenterModulePositionX(0),constants->CenterModulePositionY(0),constants->CenterModulePositionZ(0)));
		acordeModules->AddNode(acordeSingleModule,9,
			new TGeoTranslation("Mod0_9",constants->CenterModulePositionX(9),constants->CenterModulePositionY(9),constants->CenterModulePositionZ(9)));
		acordeModules->AddNode(acordeSingleModule,50,
			new TGeoTranslation("Mod0_50",constants->CenterModulePositionX(50),constants->CenterModulePositionY(50),constants->CenterModulePositionZ(50)));
		acordeModules->AddNode(acordeSingleModule,59,
			new TGeoTranslation("Mod0_59",constants->CenterModulePositionX(59),constants->CenterModulePositionY(59),constants->CenterModulePositionZ(59)));
	}

	// Create a dummy support & bars of Aluminium (it doesn't exist a survey of this structure)

	Float_t boxLongSupport[3],boxThinSupport[3];

	boxLongSupport[0]=10.0;
	boxLongSupport[1]=0.5;
	boxLongSupport[2]=500.0;

	boxThinSupport[0]=1.0;
	boxThinSupport[1]=7.0;
	boxThinSupport[2]=500.0;

	TGeoBBox *acordeLongSupport = new TGeoBBox("ACORDELONGSUPPORT",boxLongSupport[0],boxLongSupport[1],boxLongSupport[2]);
	TGeoBBox *acordeThinSupport = new TGeoBBox("ACORDETHINSUPPORT",boxThinSupport[0],boxThinSupport[1],boxThinSupport[2]);

	TGeoVolume *acordeLSupport = new TGeoVolume("ACORDELS",acordeLongSupport,aluminium);
	TGeoVolume *acordeTSupport = new TGeoVolume("ACORDETS",acordeThinSupport,aluminium);
	TGeoVolume *acordeMainSupport = new TGeoVolumeAssembly("ACORDE_SUPPORT"); 
	acordeMainSupport->AddNode(acordeLSupport,1,new TGeoTranslation("ACOLSA",0,7.5,0));
	acordeMainSupport->AddNode(acordeLSupport,2,new TGeoTranslation("ACOLSB",0,-7.5,0));
	acordeMainSupport->AddNode(acordeTSupport,3);

	// Set the values for the bars support
	
	Float_t boxSingleBar[3];
	boxSingleBar[0]=10;
	boxSingleBar[1]=37;//36.722; // Correction to avoid overlaps with the L3 magnet
	boxSingleBar[2]=10;

	Float_t boxUnionUp[3];
	boxUnionUp[0]=10;
	boxUnionUp[1]=0.5;
	boxUnionUp[2]=15;

	Float_t boxUnionDown[3];
	boxUnionDown[0]=20;
	boxUnionDown[1]=1;
	boxUnionDown[2]=20;
	
	// Volume and Box for the bar
	
	TGeoBBox *acordeSingleBarSupport = new TGeoBBox("ACORDESBARS",boxSingleBar[0],boxSingleBar[1],boxSingleBar[2]);
	TGeoVolume *acordeBarSupport = new TGeoVolume("ACORDEBARSUPPORT",acordeSingleBarSupport,aluminium);

	// Volume and Box for the bar union with the long support (Up-with supports and Down with L3 magnet)
	
	TGeoBBox *acordeSingleBoxUnionUp = new TGeoBBox("ACORDEBUP",boxUnionUp[0],boxUnionUp[1],boxUnionUp[2]);
	TGeoVolume *acordeBoxUnionUp = new TGeoVolume("ACORDEBOXUNIONUP",acordeSingleBoxUnionUp,aluminium);
	
	TGeoBBox *acordeSingleBoxUnionDown = new TGeoBBox("ACORDEBDOWN",boxUnionDown[0],boxUnionDown[1],boxUnionDown[2]);
	TGeoVolume *acordeBoxUnionDown = new TGeoVolume("ACORDEBOXUNIONDOWN",acordeSingleBoxUnionDown,aluminium);

	TGeoVolume *acordeMainBar = new TGeoVolumeAssembly("ACORDE_BAR");
	acordeMainBar->AddNode(acordeBoxUnionUp,1,new TGeoTranslation("ACOBAR01",0,boxSingleBar[1]+boxUnionUp[1],0));
	acordeMainBar->AddNode(acordeBarSupport,2);
	acordeMainBar->AddNode(acordeBoxUnionDown,3,new TGeoTranslation("ACOBAR01",0,-boxSingleBar[1]-boxUnionDown[1],0));

	// Volume for the Full support (supports and bars) UP face of L3 Magnet
	Float_t supportPosXIn = constants->CenterModulePositionX(20);
	Float_t supportPosY = 859.044-7.5-5.5; // Minimum module position Y less the heigh of the support
	Float_t supportPosZ = 0;
	Float_t supportPosXOut = constants->CenterModulePositionX(30);
	Float_t deltaXA = 120.;
	Float_t deltaXB = 60.;

	TGeoVolume *acordeFullSupportUpFace = new TGeoVolumeAssembly("ACORDE_FULL_SUPPORT_UPFACE");
	acordeFullSupportUpFace->AddNode(acordeMainSupport,1,new TGeoTranslation("ACOFSB01",supportPosXIn+deltaXA,supportPosY,supportPosZ));
	acordeFullSupportUpFace->AddNode(acordeMainSupport,2,new TGeoTranslation("ACOFSB02",supportPosXIn+deltaXB,supportPosY,supportPosZ));
	acordeFullSupportUpFace->AddNode(acordeMainSupport,3,new TGeoTranslation("ACOFSB03",supportPosXIn-deltaXA,supportPosY,supportPosZ));
	acordeFullSupportUpFace->AddNode(acordeMainSupport,4,new TGeoTranslation("ACOFSB04",supportPosXIn-deltaXB,supportPosY,supportPosZ));
	acordeFullSupportUpFace->AddNode(acordeMainSupport,5,new TGeoTranslation("ACOFSB05",supportPosXOut+deltaXA,supportPosY,supportPosZ));
	acordeFullSupportUpFace->AddNode(acordeMainSupport,6,new TGeoTranslation("ACOFSB06",supportPosXOut+deltaXB,supportPosY,supportPosZ));
	acordeFullSupportUpFace->AddNode(acordeMainSupport,7,new TGeoTranslation("ACOFSB07",supportPosXOut-deltaXA,supportPosY,supportPosZ));
	acordeFullSupportUpFace->AddNode(acordeMainSupport,8,new TGeoTranslation("ACOFSB08",supportPosXOut-deltaXB,supportPosY,supportPosZ));
	
	// Put the bars in the main volume acordeFullSupportUpFace

	Float_t barPosXIn = constants->CenterModulePositionX(20);
	Float_t barPosY = supportPosY-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
	Float_t barPosZ = 0;
	Int_t barIndex = 9;
	for(Int_t iBarModule = 20; iBarModule<30; iBarModule+=2)
	{
		barPosZ = constants->CenterModulePositionZ(iBarModule);
		acordeFullSupportUpFace->AddNode(acordeMainBar,barIndex,new TGeoTranslation("ACOFSB09",barPosXIn+deltaXA,barPosY,barPosZ));	
		acordeFullSupportUpFace->AddNode(acordeMainBar,barIndex+1,new TGeoTranslation("ACOFSB09",barPosXIn+deltaXB,barPosY,barPosZ));	
	        acordeFullSupportUpFace->AddNode(acordeMainBar,barIndex+2,new TGeoTranslation("ACOFSB09",barPosXIn-deltaXA,barPosY,barPosZ));
                acordeFullSupportUpFace->AddNode(acordeMainBar,barIndex+3,new TGeoTranslation("ACOFSB09",barPosXIn-deltaXB,barPosY,barPosZ));

		barIndex+=4;
	}

	Float_t barPosXOut = constants->CenterModulePositionX(30);
        for(Int_t iBarModule = 30; iBarModule<40; iBarModule+=2)
        {
                barPosZ = constants->CenterModulePositionZ(iBarModule);
                acordeFullSupportUpFace->AddNode(acordeMainBar,barIndex,new TGeoTranslation("ACOFSB09",barPosXOut+deltaXA,barPosY,barPosZ));
                acordeFullSupportUpFace->AddNode(acordeMainBar,barIndex+1,new TGeoTranslation("ACOFSB09",barPosXOut+deltaXB,barPosY,barPosZ));
	        acordeFullSupportUpFace->AddNode(acordeMainBar,barIndex+2,new TGeoTranslation("ACOFSB09",barPosXOut-deltaXA,barPosY,barPosZ));
              	acordeFullSupportUpFace->AddNode(acordeMainBar,barIndex+3,new TGeoTranslation("ACOFSB09",barPosXOut-deltaXB,barPosY,barPosZ));

                barIndex+=4;
        }

	// Supports and bars for the InSide of L3 Magnet

	TGeoVolume *acordeFullSupportInFace = new TGeoVolumeAssembly("ACORDE_FULL_SUPPORT_INFACE");
	Float_t supportPosXInA = constants->CenterModulePositionX(1)+deltaXA;
	Float_t supportPosYInA = 592.017-7.5-5.5;
	Float_t x0 = constants->CenterModulePositionX(1);
	Float_t y0 = constants->CenterModulePositionY(1);
	Float_t theta = -1*TMath::Pi()/4;
	Float_t supportPosXInPA = x0+(supportPosXInA-x0)*TMath::Cos(theta)-(supportPosYInA-y0)*TMath::Sin(theta);
	Float_t supportPosYInPA = y0+(supportPosYInA-y0)*TMath::Cos(theta)+(supportPosXInA-x0)*TMath::Sin(theta);
	Float_t supportPosZInIn = 0;

	Float_t supportPosXInB = constants->CenterModulePositionX(1)+deltaXB;
	Float_t supportPosXInPB = x0+(supportPosXInB-x0)*TMath::Cos(theta)-(supportPosYInA-y0)*TMath::Sin(theta);
	Float_t supportPosYInPB = y0+(supportPosYInA-y0)*TMath::Cos(theta)+(supportPosXInB-x0)*TMath::Sin(theta);

        Float_t supportPosXInC = constants->CenterModulePositionX(1)-deltaXA;
        Float_t supportPosXInPC = x0+(supportPosXInC-x0)*TMath::Cos(theta)-(supportPosYInA-y0)*TMath::Sin(theta);
        Float_t supportPosYInPC = y0+(supportPosYInA-y0)*TMath::Cos(theta)+(supportPosXInC-x0)*TMath::Sin(theta);

        Float_t supportPosXInD = constants->CenterModulePositionX(1)-deltaXB;
        Float_t supportPosXInPD = x0+(supportPosXInD-x0)*TMath::Cos(theta)-(supportPosYInA-y0)*TMath::Sin(theta);
        Float_t supportPosYInPD = y0+(supportPosYInA-y0)*TMath::Cos(theta)+(supportPosXInD-x0)*TMath::Sin(theta);


	acordeFullSupportInFace->AddNode(acordeMainSupport,1,new TGeoCombiTrans("ACOFSBIN01",supportPosXInPA,supportPosYInPA,supportPosZInIn,idrotm232));
	acordeFullSupportInFace->AddNode(acordeMainSupport,2,new TGeoCombiTrans("ACOFSBIN02",supportPosXInPB,supportPosYInPB,supportPosZInIn,idrotm232));
	acordeFullSupportInFace->AddNode(acordeMainSupport,3,new TGeoCombiTrans("ACOFSBIN03",supportPosXInPC,supportPosYInPC,supportPosZInIn,idrotm232));
	acordeFullSupportInFace->AddNode(acordeMainSupport,4,new TGeoCombiTrans("ACOFSBIN04",supportPosXInPD,supportPosYInPD,supportPosZInIn,idrotm232));

        Float_t supportPosXInE = constants->CenterModulePositionX(10)+deltaXA;
        Float_t supportPosYInE = 806.312-7.5-5.5;
        Float_t x00 = constants->CenterModulePositionX(10);
        Float_t y00 = constants->CenterModulePositionY(10);
        Float_t supportPosXInPE = x00+(supportPosXInE-x00)*TMath::Cos(theta)-(supportPosYInE-y00)*TMath::Sin(theta);
        Float_t supportPosYInPE = y00+(supportPosYInE-y00)*TMath::Cos(theta)+(supportPosXInE-x00)*TMath::Sin(theta);

        Float_t supportPosXInF = constants->CenterModulePositionX(10)+deltaXB;
        Float_t supportPosXInPF = x00+(supportPosXInF-x00)*TMath::Cos(theta)-(supportPosYInE-y00)*TMath::Sin(theta);
        Float_t supportPosYInPF = y00+(supportPosYInE-y00)*TMath::Cos(theta)+(supportPosXInF-x00)*TMath::Sin(theta);

        Float_t supportPosXInG = constants->CenterModulePositionX(10)-deltaXA+100;
        Float_t supportPosXInPG = x00+(supportPosXInG-x00)*TMath::Cos(theta)-(supportPosYInE-y00)*TMath::Sin(theta);
        Float_t supportPosYInPG = y00+(supportPosYInE-y00)*TMath::Cos(theta)+(supportPosXInG-x00)*TMath::Sin(theta);

        acordeFullSupportInFace->AddNode(acordeMainSupport,5,new TGeoCombiTrans("ACOFSBIN05",supportPosXInPE,supportPosYInPE,supportPosZInIn,idrotm232));
        acordeFullSupportInFace->AddNode(acordeMainSupport,6,new TGeoCombiTrans("ACOFSBIN06",supportPosXInPF,supportPosYInPF,supportPosZInIn,idrotm232));
        acordeFullSupportInFace->AddNode(acordeMainSupport,7,new TGeoCombiTrans("ACOFSBIN07",supportPosXInPG,supportPosYInPG,supportPosZInIn,idrotm232));

        // Put the bars in the main volume acordeFullSupportInFace

	Int_t barIndexIn = 8;

	for (Int_t iBarModule=1; iBarModule<9; iBarModule+=2)
	{
        	Float_t barPosXInIn = constants->CenterModulePositionX(iBarModule)+deltaXA;
		Float_t barPosYInIn = supportPosYInA-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZInIn = constants->CenterModulePositionZ(iBarModule);
		Float_t barPosXInP = x0+(barPosXInIn-x0)*TMath::Cos(theta)-(barPosYInIn-y0)*TMath::Sin(theta);
		Float_t barPosYInP = y0+(barPosYInIn-y0)*TMath::Cos(theta)+(barPosXInIn-x0)*TMath::Sin(theta);

                Float_t barPosXInInA = constants->CenterModulePositionX(iBarModule)+deltaXB;
                Float_t barPosYInInA = supportPosYInA-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
                Float_t barPosZInInA = constants->CenterModulePositionZ(iBarModule);
                Float_t barPosXInPA = x0+(barPosXInInA-x0)*TMath::Cos(theta)-(barPosYInInA-y0)*TMath::Sin(theta);
                Float_t barPosYInPA = y0+(barPosYInInA-y0)*TMath::Cos(theta)+(barPosXInInA-x0)*TMath::Sin(theta);

                Float_t barPosXInInB = constants->CenterModulePositionX(iBarModule)-deltaXA;
                Float_t barPosYInInB = supportPosYInA-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
                Float_t barPosZInInB = constants->CenterModulePositionZ(iBarModule);
                Float_t barPosXInPB = x0+(barPosXInInB-x0)*TMath::Cos(theta)-(barPosYInInB-y0)*TMath::Sin(theta);
                Float_t barPosYInPB = y0+(barPosYInInB-y0)*TMath::Cos(theta)+(barPosXInInB-x0)*TMath::Sin(theta);

                Float_t barPosXInInC = constants->CenterModulePositionX(iBarModule)-deltaXB;
                Float_t barPosYInInC = supportPosYInA-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
                Float_t barPosZInInC = constants->CenterModulePositionZ(iBarModule);
                Float_t barPosXInPC = x0+(barPosXInInC-x0)*TMath::Cos(theta)-(barPosYInInC-y0)*TMath::Sin(theta);
                Float_t barPosYInPC = y0+(barPosYInInC-y0)*TMath::Cos(theta)+(barPosXInInC-x0)*TMath::Sin(theta);

		acordeFullSupportInFace->AddNode(acordeMainBar,barIndexIn,new TGeoCombiTrans("ACOFSBIN08",barPosXInP,barPosYInP,barPosZInIn,idrotm232));
		acordeFullSupportInFace->AddNode(acordeMainBar,barIndexIn+1,new TGeoCombiTrans("ACOFSBIN09",barPosXInPA,barPosYInPA,barPosZInInA,idrotm232));
                acordeFullSupportInFace->AddNode(acordeMainBar,barIndexIn+2,new TGeoCombiTrans("ACOFSBIN10",barPosXInPB,barPosYInPB,barPosZInInB,idrotm232));
                acordeFullSupportInFace->AddNode(acordeMainBar,barIndexIn+3,new TGeoCombiTrans("ACOFSBIN11",barPosXInPC,barPosYInPC,barPosZInInC,idrotm232));

		barIndexIn+=4;
	}

        for (Int_t iBarModule=10; iBarModule<20; iBarModule+=2)
        {
                Float_t barPosXInIn = constants->CenterModulePositionX(iBarModule)+deltaXA;
                Float_t barPosYInIn = supportPosYInE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
                Float_t barPosZInIn = constants->CenterModulePositionZ(iBarModule);
                Float_t barPosXInP = x00+(barPosXInIn-x00)*TMath::Cos(theta)-(barPosYInIn-y00)*TMath::Sin(theta);
                Float_t barPosYInP = y00+(barPosYInIn-y00)*TMath::Cos(theta)+(barPosXInIn-x00)*TMath::Sin(theta);

                Float_t barPosXInInA = constants->CenterModulePositionX(iBarModule)+deltaXB;
                Float_t barPosYInInA = supportPosYInE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
                Float_t barPosZInInA = constants->CenterModulePositionZ(iBarModule);
                Float_t barPosXInPA = x00+(barPosXInInA-x00)*TMath::Cos(theta)-(barPosYInInA-y00)*TMath::Sin(theta);
                Float_t barPosYInPA = y00+(barPosYInInA-y00)*TMath::Cos(theta)+(barPosXInInA-x00)*TMath::Sin(theta);

                Float_t barPosXInInB = constants->CenterModulePositionX(iBarModule)-deltaXA+100;
                Float_t barPosYInInB = supportPosYInE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
                Float_t barPosZInInB = constants->CenterModulePositionZ(iBarModule);
                Float_t barPosXInPB = x00+(barPosXInInB-x00)*TMath::Cos(theta)-(barPosYInInB-y00)*TMath::Sin(theta);
                Float_t barPosYInPB = y00+(barPosYInInB-y00)*TMath::Cos(theta)+(barPosXInInB-x00)*TMath::Sin(theta);

                acordeFullSupportInFace->AddNode(acordeMainBar,barIndexIn,new TGeoCombiTrans("ACOFSBIN08",barPosXInP,barPosYInP,barPosZInIn,idrotm232));
                acordeFullSupportInFace->AddNode(acordeMainBar,barIndexIn+1,new TGeoCombiTrans("ACOFSBIN09",barPosXInPA,barPosYInPA,barPosZInInA,idrotm232));
                acordeFullSupportInFace->AddNode(acordeMainBar,barIndexIn+2,new TGeoCombiTrans("ACOFSBIN10",barPosXInPB,barPosYInPB,barPosZInInB,idrotm232));

                barIndexIn+=4;
        }


	// Construction of the support and bars for the OutSide


        TGeoVolume *acordeFullSupportOutFace = new TGeoVolumeAssembly("ACORDE_FULL_SUPPORT_OUTFACE");
        Float_t supportPosXOutA = constants->CenterModulePositionX(45)+deltaXA-90;
        Float_t supportPosYOutA = 807.2915-7.5-5.5;
        Float_t x000 = constants->CenterModulePositionX(45);
        Float_t y000 = constants->CenterModulePositionY(45);
        Float_t theta1 = 1*TMath::Pi()/4;

        Float_t supportPosXOutPA = x000+(supportPosXOutA-x000)*TMath::Cos(theta1)-(supportPosYOutA-y000)*TMath::Sin(theta1);
        Float_t supportPosYOutPA = y000+(supportPosYOutA-y000)*TMath::Cos(theta1)+(supportPosXOutA-x000)*TMath::Sin(theta1);
        Float_t supportPosZOutIn = 0;


        Float_t supportPosXOutC = constants->CenterModulePositionX(45)-deltaXA;
        Float_t supportPosXOutPC = x000+(supportPosXOutC-x000)*TMath::Cos(theta1)-(supportPosYOutA-1-y000)*TMath::Sin(theta1);
        Float_t supportPosYOutPC = y000+(supportPosYOutA-1-y000)*TMath::Cos(theta1)+(supportPosXOutC-x000)*TMath::Sin(theta1);

        Float_t supportPosXOutD = constants->CenterModulePositionX(45)-deltaXB;
        Float_t supportPosXOutPD = x000+(supportPosXOutD-x000)*TMath::Cos(theta1)-(supportPosYOutA-1-y000)*TMath::Sin(theta1);
        Float_t supportPosYOutPD = y000+(supportPosYOutA-1-y000)*TMath::Cos(theta1)+(supportPosXOutD-x000)*TMath::Sin(theta1);


        acordeFullSupportOutFace->AddNode(acordeMainSupport,1,new TGeoCombiTrans("ACOFSBOUT01",supportPosXOutPA,supportPosYOutPA,supportPosZOutIn,idrotm231));
	acordeFullSupportOutFace->AddNode(acordeMainSupport,2,new TGeoCombiTrans("ACOFSBOUT02",supportPosXOutPC,supportPosYOutPC,supportPosZOutIn,idrotm231));
	acordeFullSupportOutFace->AddNode(acordeMainSupport,3,new TGeoCombiTrans("ACOFSBOUT03",supportPosXOutPD,supportPosYOutPD,supportPosZOutIn,idrotm231));
	

	Float_t supportPosXOutE = constants->CenterModulePositionX(52)+deltaXA;
	Float_t supportPosYOutE = 585.616-7.5-5.5-0.1;
	Float_t x0000 = constants->CenterModulePositionX(52);
	Float_t y0000 = constants->CenterModulePositionY(52);

	Float_t supportPosXOutEP = x0000+(supportPosXOutE-x0000)*TMath::Cos(theta1)-(supportPosYOutE-y0000)*TMath::Sin(theta1);
	Float_t supportPosYOutEP = y0000+(supportPosYOutE-y0000)*TMath::Cos(theta1)+(supportPosXOutE-x0000)*TMath::Sin(theta1);
	
	Float_t supportPosXOutF = constants->CenterModulePositionX(52)+deltaXB;

	Float_t supportPosXOutFP = x0000+(supportPosXOutF-x0000)*TMath::Cos(theta1)-(supportPosYOutE-y0000)*TMath::Sin(theta1);
	Float_t supportPosYOutFP = y0000+(supportPosYOutE-y0000)*TMath::Cos(theta1)+(supportPosXOutF-x0000)*TMath::Sin(theta1);

	Float_t supportPosXOutG = constants->CenterModulePositionX(52)-deltaXA;
	Float_t supportPosXOutGP = x0000+(supportPosXOutG-x0000)*TMath::Cos(theta1)-(supportPosYOutE-y0000)*TMath::Sin(theta1);
	Float_t supportPosYOutGP = y0000+(supportPosYOutE-y0000)*TMath::Cos(theta1)+(supportPosXOutG-x0000)*TMath::Sin(theta1);

	Float_t supportPosXOutH = constants->CenterModulePositionX(52)-deltaXB;
	Float_t supportPosXOutHP = x0000+(supportPosXOutH-x0000)*TMath::Cos(theta1)-(supportPosYOutE-0.4-y0000)*TMath::Sin(theta1);
	Float_t supportPosYOutHP = y0000+(supportPosYOutE-0.4-y0000)*TMath::Cos(theta1)+(supportPosXOutH-x0000)*TMath::Sin(theta1);


	acordeFullSupportOutFace->AddNode(acordeMainSupport,4,new TGeoCombiTrans("ACOFSBOUT04",supportPosXOutEP,supportPosYOutEP,supportPosZOutIn,idrotm231));
	acordeFullSupportOutFace->AddNode(acordeMainSupport,5,new TGeoCombiTrans("ACOFSBOUT05",supportPosXOutFP,supportPosYOutFP,supportPosZOutIn,idrotm231));
	acordeFullSupportOutFace->AddNode(acordeMainSupport,6,new TGeoCombiTrans("ACOFSBOUT06",supportPosXOutGP,supportPosYOutGP,supportPosZOutIn,idrotm231));
	acordeFullSupportOutFace->AddNode(acordeMainSupport,7,new TGeoCombiTrans("ACOFSBOUT07",supportPosXOutHP,supportPosYOutHP,supportPosZOutIn,idrotm231));

	// Put the bars of the PutFace Side 2 L3-Magnet


	Int_t indexBar0=8;
	for(Int_t iAcoBar = 40;iAcoBar < 50 ; iAcoBar+=2)
	{
		Float_t barPosXOutIn = constants->CenterModulePositionX(iAcoBar)+deltaXA-90;
		Float_t barPosYOutIn = supportPosYOutA-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutIn = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutP = x000+(barPosXOutIn-x000)*TMath::Cos(theta1)-(barPosYOutIn-y000)*TMath::Sin(theta1);
		Float_t barPosYOutP = y000+(barPosYOutIn-y000)*TMath::Cos(theta1)+(barPosXOutIn-x000)*TMath::Sin(theta1);

		Float_t barPosXOutInA = constants->CenterModulePositionX(iAcoBar)-deltaXA;
		Float_t barPosYOutInA = supportPosYOutA-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5-1.;
        	Float_t barPosZOutInA = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPA = x000+(barPosXOutInA-x000)*TMath::Cos(theta1)-(barPosYOutInA-y000)*TMath::Sin(theta1);
		Float_t barPosYOutPA = y000+(barPosYOutInA-y000)*TMath::Cos(theta1)+(barPosXOutInA-x000)*TMath::Sin(theta1);


		Float_t barPosXOutInB = constants->CenterModulePositionX(iAcoBar)-deltaXB;
		Float_t barPosYOutInB = supportPosYOutA-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5-1.;
        	Float_t barPosZOutInB = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPB = x000+(barPosXOutInB-x000)*TMath::Cos(theta1)-(barPosYOutInB-y000)*TMath::Sin(theta1);
		Float_t barPosYOutPB = y000+(barPosYOutInB-y000)*TMath::Cos(theta1)+(barPosXOutInB-x000)*TMath::Sin(theta1);
	
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutP,barPosYOutP,barPosZOutIn,idrotm231));
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+1,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPA,barPosYOutPA,barPosZOutInA,idrotm231));
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+2,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPB,barPosYOutPB,barPosZOutInB,idrotm231));

		indexBar0+=4;
	}
	for(Int_t iAcoBar = 51;iAcoBar < 54 ; iAcoBar+=2)
	{

		Float_t barPosXOutInC = constants->CenterModulePositionX(iAcoBar)+deltaXA;
		Float_t barPosYOutInC = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutInC = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPC = x0000+(barPosXOutInC-x0000)*TMath::Cos(theta1)-(barPosYOutInC-y0000)*TMath::Sin(theta1);
		Float_t barPosYOutPC = y0000+(barPosYOutInC-y0000)*TMath::Cos(theta1)+(barPosXOutInC-x0000)*TMath::Sin(theta1);


		Float_t barPosXOutInD = constants->CenterModulePositionX(iAcoBar)+deltaXB;
		Float_t barPosYOutInD = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutInD = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPD = x0000+(barPosXOutInD-x0000)*TMath::Cos(theta1)-(barPosYOutInD-y0000)*TMath::Sin(theta1);
		Float_t barPosYOutPD = y0000+(barPosYOutInD-y0000)*TMath::Cos(theta1)+(barPosXOutInD-x0000)*TMath::Sin(theta1);


		Float_t barPosXOutInE = constants->CenterModulePositionX(iAcoBar)-deltaXA;
		Float_t barPosYOutInE = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutInE = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPE = x0000+(barPosXOutInE-x0000)*TMath::Cos(theta1)-(barPosYOutInE-y0000)*TMath::Sin(theta1);
		Float_t barPosYOutPE = y0000+(barPosYOutInE-y0000)*TMath::Cos(theta1)+(barPosXOutInE-x0000)*TMath::Sin(theta1);
		Float_t barPosXOutInF = constants->CenterModulePositionX(iAcoBar)-deltaXB;
		Float_t barPosYOutInF = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutInF = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPF = x0000+(barPosXOutInF-x0000)*TMath::Cos(theta1)-(barPosYOutInF-0.4-y0000)*TMath::Sin(theta1);
		Float_t barPosYOutPF = y0000+(barPosYOutInF-0.4-y0000)*TMath::Cos(theta1)+(barPosXOutInF-x0000)*TMath::Sin(theta1);

		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPC,barPosYOutPC,barPosZOutInC,idrotm231));
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+1,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPD,barPosYOutPD,barPosZOutInD,idrotm231));
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+2,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPE,barPosYOutPE,barPosZOutInE,idrotm231));
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+3,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPF,barPosYOutPF,barPosZOutInF,idrotm231));
		indexBar0+=4;
	}


	for(Int_t iAcoBar = 57;iAcoBar < 58 ; iAcoBar+=2)
	{

		Float_t barPosXOutInC = constants->CenterModulePositionX(iAcoBar)+deltaXA;
		Float_t barPosYOutInC = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutInC = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPC = x0000+(barPosXOutInC-x0000)*TMath::Cos(theta1)-(barPosYOutInC-y0000)*TMath::Sin(theta1);
		Float_t barPosYOutPC = y0000+(barPosYOutInC-y0000)*TMath::Cos(theta1)+(barPosXOutInC-x0000)*TMath::Sin(theta1);


		Float_t barPosXOutInD = constants->CenterModulePositionX(iAcoBar)+deltaXB;
		Float_t barPosYOutInD = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutInD = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPD = x0000+(barPosXOutInD-x0000)*TMath::Cos(theta1)-(barPosYOutInD-y0000)*TMath::Sin(theta1);
		Float_t barPosYOutPD = y0000+(barPosYOutInD-y0000)*TMath::Cos(theta1)+(barPosXOutInD-x0000)*TMath::Sin(theta1);


		Float_t barPosXOutInE = constants->CenterModulePositionX(iAcoBar)-deltaXA;
		Float_t barPosYOutInE = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutInE = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPE = x0000+(barPosXOutInE-x0000)*TMath::Cos(theta1)-(barPosYOutInE-y0000)*TMath::Sin(theta1);
		Float_t barPosYOutPE = y0000+(barPosYOutInE-y0000)*TMath::Cos(theta1)+(barPosXOutInE-x0000)*TMath::Sin(theta1);
		Float_t barPosXOutInF = constants->CenterModulePositionX(iAcoBar)-deltaXB;
		Float_t barPosYOutInF = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        	Float_t barPosZOutInF = constants->CenterModulePositionZ(iAcoBar);
		Float_t barPosXOutPF = x0000+(barPosXOutInF-x0000)*TMath::Cos(theta1)-(barPosYOutInF-0.4-y0000)*TMath::Sin(theta1);
		Float_t barPosYOutPF = y0000+(barPosYOutInF-0.4-y0000)*TMath::Cos(theta1)+(barPosXOutInF-x0000)*TMath::Sin(theta1);

		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPC,barPosYOutPC,barPosZOutInC,idrotm231));
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+1,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPD,barPosYOutPD,barPosZOutInD,idrotm231));
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+2,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPE,barPosYOutPE,barPosZOutInE,idrotm231));
		acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+3,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPF,barPosYOutPF,barPosZOutInF,idrotm231));
		indexBar0+=4;
	}


	Float_t barPosXOutInC = constants->CenterModulePositionX(52)+deltaXA;
	Float_t barPosYOutInC = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
        Float_t barPosZOutInC = constants->CenterModulePositionZ(55);
	Float_t barPosXOutPC = x0000+(barPosXOutInC-x0000)*TMath::Cos(theta1)-(barPosYOutInC-y0000)*TMath::Sin(theta1);
	Float_t barPosYOutPC = y0000+(barPosYOutInC-y0000)*TMath::Cos(theta1)+(barPosXOutInC-x0000)*TMath::Sin(theta1);


	Float_t barPosXOutInD = constants->CenterModulePositionX(52)+deltaXB;
	Float_t barPosYOutInD = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
	Float_t barPosXOutPD = x0000+(barPosXOutInD-x0000)*TMath::Cos(theta1)-(barPosYOutInD-y0000)*TMath::Sin(theta1);
	Float_t barPosYOutPD = y0000+(barPosYOutInD-y0000)*TMath::Cos(theta1)+(barPosXOutInD-x0000)*TMath::Sin(theta1);


	Float_t barPosXOutInE = constants->CenterModulePositionX(52)-deltaXA;
	Float_t barPosYOutInE = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
	Float_t barPosXOutPE = x0000+(barPosXOutInE-x0000)*TMath::Cos(theta1)-(barPosYOutInE-y0000)*TMath::Sin(theta1);
	Float_t barPosYOutPE = y0000+(barPosYOutInE-y0000)*TMath::Cos(theta1)+(barPosXOutInE-x0000)*TMath::Sin(theta1);
	Float_t barPosXOutInF = constants->CenterModulePositionX(52)-deltaXB;
	Float_t barPosYOutInF = supportPosYOutE-boxSingleBar[1]-boxUnionUp[1]-boxUnionDown[1]-7.5;
	Float_t barPosXOutPF = x0000+(barPosXOutInF-x0000)*TMath::Cos(theta1)-(barPosYOutInF-0.4-y0000)*TMath::Sin(theta1);
	Float_t barPosYOutPF = y0000+(barPosYOutInF-0.4-y0000)*TMath::Cos(theta1)+(barPosXOutInF-x0000)*TMath::Sin(theta1);

	acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPC,barPosYOutPC,barPosZOutInC,idrotm231));
	acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+1,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPD,barPosYOutPD,barPosZOutInC,idrotm231));
	acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+2,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPE,barPosYOutPE,barPosZOutInC,idrotm231));
	acordeFullSupportOutFace->AddNode(acordeMainBar,indexBar0+3,new TGeoCombiTrans("ACOFSBOUT08",barPosXOutPF,barPosYOutPF,barPosZOutInC,idrotm231));
	




	supportBars->AddNode(acordeFullSupportUpFace,1);
	supportBars->AddNode(acordeFullSupportInFace,2);
	supportBars->AddNode(acordeFullSupportOutFace,3);
	acorde->AddNode(acordeModules,1);
	acorde->AddNode(supportBars,2);
	alice->AddNode(acorde,1);//---> put volume of ACORDE over ALICE's volume



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
 static Int_t idScint = gMC->VolId("ACORDESCINTILLATORMODULE");
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
  snprintf(branchname,9,"%s",GetName());
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
	//			|----> ALL_ACORDE_MODULES_1/ACORDE_MODULE_%d (d:0->to->59)
	//			|----> ACORDE_SUPPORT_BARS_2    |--> BARS&SUPPORTS
	//
	//     Send comments to: Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>

	TString vpstr1 = "ALIC_1/ACORDE_1/ALL_ACORDE_MODULES_1/ACORDE_MODULE_";
	TString snstr1 = "ACORDE/Array";
	TString volpath, symname;
	for(Int_t dy=0; dy<60 ; dy++)
	{
		volpath = vpstr1;
		volpath += dy;
		symname = snstr1;
		symname += dy;
		if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data()))
	        AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", symname.Data(),volpath.Data()));
	}
}
