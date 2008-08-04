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
#include <TGeometry.h>
#include <TMath.h>
#include <TTUBE.h>
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
#include <TBRIK.h>
#include <TNode.h>
 

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
void AliACORDEv1::BuildGeometry()
{

  // not needed anymore

}

//_____________________________________________________________________________
void AliACORDEv1::CreateGeometry()
{
  CreateAcorde();
  if (GetCreateCavern()) CreateCavern();
}


void AliACORDEv1::CreateCavern()
{

	// Create the mother volume, the one which contain all the material
	//above the hall

	TGeoManager *acorde = new TGeoManager("ACORDE", "Geometry of ACORDE");	

	//---> define some materials

	TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);

	//---> define some media

	TGeoMedium *vacuum = new TGeoMedium("Vacuum",1, matVacuum);

	//---> define the measures

	Double_t dx1 = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	Double_t dy1 = AliACORDEConstants::Instance()->Depth();
	Double_t dz1 = dx1;

	//---> define the box for the mother volume

	TGeoVolume *aCORDE = acorde->MakeBox("ACORDE", vacuum, dx1, dy1, dz1);
	acorde->SetTopVolume(aCORDE);

	//---> create shafts&molasse

	CreateShafts();
	CreateMolasse();
}


void AliACORDEv1::CreateShafts()
{

	//---> This shaft is composes by an open tube down in the hall
	//---> and a cilinder above the level of the celling
	//---> Every structure relative to the shaft will be put into this volume


	TGeoManager *acorde = new TGeoManager("ACORDE2007", "Geometry of ACORDE");	

	//---> define some materials

	TGeoMaterial *matVacuum = new TGeoMaterial("Al", 0,0,0);

	//---> define some media

	TGeoMedium *vacuum = new TGeoMedium("Vacuum",1, matVacuum);
	TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
	TGeoMedium *al = new TGeoMedium("Root Material",2, matAl);
	
	
	//---> Access to shafts
	//---> define the Matrix Rotation&other variables

	TGeoRotation *rot1 = new TGeoRotation("rot1", 0.0, 0.0, 90.0, 0.0, 90.0, 90.0);
	
	Float_t ptube[5];
	ptube[0]=0;
	ptube[1]=1250;
	ptube[2]=5150/2;
	ptube[3]=360;
	ptube[4]=360;
	
 	Float_t ptubs[5];

	//---> The open section of the PX24
	ptubs[0] = 1150; //---> Inner radius
	ptubs[1] = 1250; //---> Outer radius
	ptubs[2] = 1300; //---> Half length
	ptubs[3] = 180 + kRaddeg*TMath::ASin(1070/ptubs[0]); //---> starting angle
	ptubs[4] = 180 -  kRaddeg*TMath::ASin(1070/ptubs[0]);
		
	//---> Set position for the tubes 

	TGeoTranslation *tr2 = new TGeoTranslation(0,0,-ptube[2]+ptubs[2]);

	//---> define the cilinders to hold the main structure in the shaft

	TGeoVolume *o = acorde->MakeBox("O", vacuum, 25., 25., 5.);
	TGeoVolume *cSF1 = acorde->MakeTubs("CSF1",al,ptube[0],ptube[1],ptube[2],ptube[3],ptube[4]);
	TGeoVolume *cSF2 = acorde->MakeTubs("CSF2",al,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	o->AddNode(cSF1, 1);
	cSF1->AddNode(cSF2,1,tr2);

	//---> definition of the other part of the shaft

  	ptube[0] = ptubs[0]; // Inner radius
	ptube[1] = ptubs[1]; // Outer radius
	ptube[2] = 5150/2 - ptubs[2]; // Half lenght
	TGeoVolume *cSF3 = acorde->MakeTubs("CSF3",al,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	TGeoTranslation *tr3 = new TGeoTranslation(0,0,5150/2-ptube[2]);
	cSF1->AddNode(cSF3,1,tr3);

	//---> define concrete walls along the shaft (next to the elevator)
	
	Float_t pbox[3];
	pbox[0]=480/2;
	pbox[1]=120/2;
	pbox[2]=5150/2;
	TGeoVolume *cSW1 = acorde->MakeBox("CSW1",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br1 = new TGeoTranslation(820+pbox[0],150+pbox[1],0);
	TGeoTranslation *br1a = new TGeoTranslation(820+pbox[0],-300-pbox[1],0);
	cSF1->AddNode(cSW1,1,br1);
	cSF1->AddNode(cSW1,1,br1a);

	pbox[0] = 120/2;  // Half length in X
	pbox[1] = 750/2;  // Half length in Y
	pbox[2] = 5150/2; // Half length in Z
	TGeoVolume *cSW2 = acorde->MakeBox("CSW2",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br2 = new TGeoTranslation(820-60,150+pbox[1],0);
	cSF1->AddNode(cSW2,1,br2);


	pbox[0] = 120/2;  // Half length in X
	pbox[1] = 600/2;  // Half lenght in Y
	pbox[2] = 5150/2; // Half length in Z
	TGeoVolume *cSW3 = acorde->MakeBox("CSW3",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br3 = new TGeoTranslation(820-60,-300-pbox[1],0);
	cSF1->AddNode(cSW3,1,br3);

	pbox[0] = 400/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 300/2; // Half length in Z
	TGeoVolume *cSW4 = acorde->MakeBox("CSW4",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br4 = new TGeoTranslation(pbox[1]-pbox[0],0,3000-5150/2-pbox[2]);
	cSF1->AddNode(cSW4,1,br4);


	pbox[0] = 1400/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 170/2; // Half length in Z
	TGeoVolume *cSW5 = acorde->MakeBox("CSW5",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br5 = new TGeoTranslation(0,0,3000-5150/2-130);
	cSF1->AddNode(cSW5,1,br5);


	pbox[0] = 170/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 300/2; // Half length in Z
	TGeoVolume *cSW6 = acorde->MakeBox("CSW6",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br6 = new TGeoTranslation(-1400/2-pbox[0],0,3000-5150/2-pbox[2]);
	cSF1->AddNode(cSW6,1,br6);


	pbox[0] = 100/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 450/2; // Half length in Z
	TGeoVolume *cSW7 = acorde->MakeBox("CSW7",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br7 = new TGeoTranslation(-1400/2-170-pbox[0],0,3000-5150/2+pbox[2]);
	cSF1->AddNode(cSW7,1,br7);


	pbox[0] = 300/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 170/2; // Half length in Z
	TGeoVolume *cSW8 = acorde->MakeBox("CSW8",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br8 = new TGeoTranslation(-2300/2+pbox[0],0,2500-5150/2);
	cSF1->AddNode(cSW8,1,br8);

	//---> put the shaft into the mother volume

	TGeoCombiTrans *br = new TGeoCombiTrans(0,AliACORDEConstants::Instance()->Depth()-5150/2,2300,rot1);
	cSF1->AddNode(cSF1,1,br);


	//---> PM25 Access Shafts

	ptube[0]=910/2;
	ptube[1]=ptube[0]+100;
	ptube[2]=(5150-1166)/2;
	TGeoVolume *cSF4 = acorde->MakeTubs("CSF4",vacuum,pbox[0],pbox[1],pbox[2],360,360);
	TGeoCombiTrans *tr4 = new TGeoCombiTrans(2100,AliACORDEConstants::Instance()->Depth()-ptube[2],0,rot1);
	cSF4->AddNode(cSF4,1,tr4);


	//---> PGC2 Access shaft

	ptube[0]=1100/2;
	ptube[1]=ptube[0]+100;
	ptube[2]=(5150-690)/2;
	TGeoVolume *cSF5 = acorde->MakeTubs("CSF5",vacuum,pbox[0],pbox[1],pbox[2],360,360);
	TGeoCombiTrans *tr5 = new TGeoCombiTrans(-375,AliACORDEConstants::Instance()->Depth()-ptube[2],-1900-2987.7,rot1);
	cSF5->AddNode(cSF5,1,tr5);

}


void AliACORDEv1::CreateMolasse()

{
	// create a big molasse for ACORDE detector
	TGeoManager *acorde = new TGeoManager("ACORDE2007", "Geometry of ACORDE");	

	//---> define some media
	

	TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
	TGeoMedium *vacuum = new TGeoMedium("Root Material",2, matAl);

	Float_t px24radius = 2300/2;
	Float_t px24X = 0;
	Float_t px24Z = 2300;
	Float_t pm25radius = 910/2;
	Float_t pm25X = 2100;
	Float_t pm25Z = 0;
	Float_t pgc2radius = 1100/2;
	Float_t pgc2X = -375;
	Float_t pgc2Z = -(1900 + 2987.7);
	Float_t concreteWidth = 100; //---> Standard width of the hall walls.


	//---> Create a local mother volume.
	Float_t pbox[3];
	pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = pbox[0];
	TGeoVolume *cM01 = acorde->MakeBox("CM01", vacuum, pbox[0],pbox[1],pbox[2]);

	//---> Now put the molasse exactly above the hall. OK
	//---> Above the ceiling
	
	Float_t ptubs[5];
	ptubs[0] = 1170;
	ptubs[1] = 2100 - pm25radius;
	ptubs[2] = 1900/2 + px24radius;
	ptubs[3] = 0;
	ptubs[4] = 180;
	TGeoVolume *cM02 = acorde->MakeTubs("CM02",vacuum,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	TGeoTranslation *tr2 = new TGeoTranslation(0,500-AliACORDEConstants::Instance()->Depth()/2,ptubs[2]-1900);
	cM01->AddNode(cM02,1,tr2);


	//---> Molasse around the RB24/26 Wall. OK

	ptubs[0] = 220 + 1600;
	ptubs[1] = AliACORDEConstants::Instance()->Depth() - ptubs[0];
	ptubs[2] = 2987.7/2 - 1100/4 - concreteWidth/2;
	ptubs[3] = 0;
	ptubs[4] = 180;
	TGeoVolume *cM03 = acorde->MakeTubs("CM03",vacuum,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	TGeoTranslation *tr3 = new TGeoTranslation(70,40-AliACORDEConstants::Instance()->Depth()/2,-ptubs[2]-1900);
	cM01->AddNode(cM03,1,tr3);


	//---> A big block above the RB24/26 wall. OK

	pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	pbox[1] = (AliACORDEConstants::Instance()->Depth() - 220 - 1600)/2;
	pbox[2] = 2987.7/2 - 1100/4 - concreteWidth/2;
	TGeoVolume *cM04 = acorde->MakeBox("CM04", vacuum, pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr4 = new TGeoTranslation(0,AliACORDEConstants::Instance()->Depth()/2-pbox[1],-1900-pbox[2]);
	cM01->AddNode(cM04,1,tr4);




	//---> Small blocks below the volume CMO4 on both sides of the wall RB24/26. OK

	pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)-ptubs[0])/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2 - pbox[1];
	TGeoVolume *cM17 = acorde->MakeBox("CM17", vacuum, pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr17 = new TGeoTranslation(AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - pbox[0],-AliACORDEConstants::Instance()->Depth()/2 + pbox[1],-1900 - pbox[2]);
	TGeoTranslation *tr17a = new TGeoTranslation(-AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)+ pbox[0],-AliACORDEConstants::Instance()->Depth()/2 + pbox[1], -1900 - pbox[2]);
	cM01->AddNode(cM17,1,tr17);
	cM01->AddNode(cM17,2,tr17a);


	//---> And a big block of molasse above the hall up to the surface. OK

	pbox[0] = pm25X - pm25radius;
	pbox[1] = (AliACORDEConstants::Instance()->Depth()-500-1170)/2;
	pbox[2] = (1900 + 1150)/2;
	TGeoVolume *cM05 = acorde->MakeBox("CM05", vacuum, pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr5 = new TGeoTranslation(0,AliACORDEConstants::Instance()->Depth()/2-pbox[1], pbox[2]-1900);
	cM01->AddNode(cM05,1,tr5);


	//---> Small blocks of molasse betwen the blocks CMO2, CMO5 and PM25. Ok

	pbox[0] = (pm25X - pm25radius - 1170)/2;
	pbox[1] = 1000;
	TGeoVolume *cM16 = acorde->MakeBox("CM16", vacuum, pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr16 = new TGeoTranslation(1170 + pbox[0], -AliACORDEConstants::Instance()->Depth()/2+pbox[1], pbox[2] - 1900);
	cM01->AddNode(cM16,1,tr16);


	//---> Molasse around the shafts.

	TGeoRotation *rot2 = new TGeoRotation("rot1",0, 0, 90, 0, 90, 90 );

	//---> Around the PX24, the open section. OK

	ptubs[0] = px24radius + concreteWidth;
	ptubs[1] = ptubs[0] + 1000;
	ptubs[2] = (2300 - (5150 - AliACORDEConstants::Instance()->Depth()))/2;
	ptubs[3] = 180 + kRaddeg*TMath::ASin(1070/ptubs[0]);
	ptubs[4] = 180 -  kRaddeg*TMath::ASin(1070/ptubs[0]);
	TGeoVolume *cM06 = acorde->MakeTubs("CM06", vacuum,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	TGeoTranslation *tr6 = new TGeoTranslation(px24X, ptubs[2] - AliACORDEConstants::Instance()->Depth()/2, px24Z);
	cM01->AddNode(cM06,1,tr6);


	//---> Around the PX24, the closed section. OK

	Float_t ptube[3];
	ptube[0] = px24radius + concreteWidth;
	ptube[1] = ptube[0] + 1000;
	ptube[2] = (5150 - 2300)/2;
	TGeoVolume *cM07 = acorde->MakeTubs("CM07", vacuum,ptube[0],ptube[1],ptubs[2],ptube[3],ptube[4]);
	TGeoTranslation *tr7 = new TGeoTranslation(px24X, AliACORDEConstants::Instance()->Depth()/2-ptube[2], px24Z);
	cM01->AddNode(cM07,1,tr7);


	//---> Around PM25. OK

	ptube[0] = pm25radius + concreteWidth;
	ptube[1] = ptube[0] + 400;
	ptube[2] = AliACORDEConstants::Instance()->Depth()/2;
	TGeoVolume *cM08 = acorde->MakeTubs("CM08", vacuum,ptube[0],ptube[1],ptube[2],ptube[3],ptube[4]);
	TGeoCombiTrans *tr8 = new TGeoCombiTrans(pm25X, 0, pm25Z,rot2);
	cM01->AddNode(cM08,1,tr8);


	//---> On both sides of the PM25 along the HALL.

	pbox[0] = (2100 + pm25radius - 1170)/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (3*px24radius - pm25radius)/2;
	TGeoVolume *cM18 = acorde->MakeBox("CM18",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr18 = new TGeoTranslation(2100, 0, pbox[2] + pm25radius);
	cM01->AddNode(cM18,1,tr18);
  
  	pbox[2] = (1900 - pm25radius)/2;
	TGeoVolume *cM19 = acorde->MakeBox("CM19",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr19 = new TGeoTranslation(2100, 0, -pbox[2] - pm25radius);
  	cM01->AddNode(cM19,1,tr19);


	//---> Around the PGC2. OK

	ptube[0] = pgc2radius + concreteWidth;
	ptube[1] = 2987.7 - 740;
	ptube[2] = AliACORDEConstants::Instance()->Depth()/2;
	TGeoVolume *cM09 = acorde->MakeTubs("CM09",vacuum,ptube[0],ptube[1],ptube[2],ptube[3],ptube[4]);
	TGeoCombiTrans *tr09 = new TGeoCombiTrans(pgc2X, 0, pgc2Z,rot2);
	cM01->AddNode(cM09,1,tr09);

	//---> On both sides of the PGC2.OK

	pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)-1100 - 375)/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = pgc2radius + concreteWidth;
	TGeoVolume *cM10 = acorde->MakeBox("CM10",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr10 = new TGeoTranslation(AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - pbox[0], 0, pgc2Z);
	TGeoTranslation *tr10a = new TGeoTranslation(-AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) + pbox[0], 0, pgc2Z);
	cM01->AddNode(cM10,1,tr10);
	cM01->AddNode(cM10,2,tr10a);


	//---> big block of molasse behind the PX24. OK

	pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (pbox[0] - (2300 + 1150 + 100))/2;
	TGeoVolume *cM12 = acorde->MakeBox("CM12",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr12 = new TGeoTranslation(px24X, 0, px24Z + px24radius + concreteWidth + pbox[2]);	
	cM01->AddNode(cM12,1,tr12);


	//---> big block of molasse in the opposite side of the PM25. OK

	pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)-1150)/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (1900 + 2300 + 1150)/2;
	TGeoVolume *cM13 = acorde->MakeBox("CM13",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr13 = new TGeoTranslation(-1150 - pbox[0], 0, pbox[2] - 1900);	
	cM01->AddNode(cM13,1,tr13);
 

	//---> big block of molasse behind the PM25. OK

	pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)-(2100 + 910/2 + 100))/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (1900 + 2300 + 1150)/2;
	TGeoVolume *cM14 = acorde->MakeBox("CM14",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr14 = new TGeoTranslation(pm25X + pm25radius + concreteWidth + pbox[0], 0, pbox[2] - 1900);
	cM01->AddNode(cM14,1,tr14);


	//---> big block of molasse behind the PGC2. OK

	pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (pbox[0] - (2987.7 + 1900 + 1100/2 + 100))/2;
	TGeoVolume *cM15 = acorde->MakeBox("CM15",vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr15 = new TGeoTranslation(0, 0, -pbox[0] + pbox[2]);
	TGeoTranslation *tr15a = new TGeoTranslation(0,AliACORDEConstants::Instance()->Depth()/2,0);
	cM01->AddNode(cM15,1,tr15);
	cM01->AddNode(cM01,1,tr15a);

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
	Float_t w1 = 8;

	


	// Put the Modules of In-Face
	
	count=1;
	for(Int_t i=1;i<9;i++){

		Float_t posx = constants->ModulePositionX(i);
		Float_t posy = constants->ModulePositionY(i);
		Float_t posz = constants->ModulePositionZ(i);	
                Int_t moduleElectronicID = constants->ModuleElectronicChannel(i); 
		
		modules->AddNode(acomodule,moduleElectronicID,
			new TGeoCombiTrans("aco01",posx,posy-w1,posz,idrotm232));
		count++;

	}

	count=9;
	for(Int_t i=10;i<20;i++){
		Float_t posx = constants->ModulePositionX(i);
		Float_t posy = constants->ModulePositionY(i);
		Float_t posz = constants->ModulePositionZ(i);
		Int_t moduleElectronicID = constants->ModuleElectronicChannel(i); 	

		modules->AddNode(acomodule,moduleElectronicID,
			new TGeoCombiTrans("aco01",posx,posy-w1,posz,idrotm232));
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
			new TGeoCombiTrans("aco01",posx,posy-w1,posz,idrotm231));
		count++;
	}

	// Put the Modules of Out-Face

	count=11;
	for(Int_t i=51;i<59;i++){
		Float_t posx = constants->ModulePositionX(i);
		Float_t posy = constants->ModulePositionY(i);
		Float_t posz = constants->ModulePositionZ(i);	
		Int_t moduleElectronicID = constants->ModuleElectronicChannel(i); 

	if ((i==57) || (i==56))
		 modules->AddNode(acomodule,moduleElectronicID,
					new TGeoCombiTrans("aco01",posx,posy-w1,posz-w1,idrotm231));
	else
		modules->AddNode(acomodule,moduleElectronicID,
			new TGeoCombiTrans("aco01",posx,posy-w1,posz,idrotm231));
		count++;
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
					constants->ModulePositionX(ma)-0.5*293+dy2,
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));

		upFace->AddNode(aCORDE00,count,aco00);

		TGeoTranslation *aco00q1=new TGeoTranslation("aco00q1",
					-(constants->ModulePositionX(ma)-0.5*293+dy2),
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+1,aco00q1);

		TGeoTranslation *aco00q2=new TGeoTranslation("aco00q2",
					constants->ModulePositionX(ma)+0.5*293-dy2,
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+2,aco00q2);

		TGeoTranslation *aco00q3=new TGeoTranslation("aco00q3",
					-(constants->ModulePositionX(ma)+0.5*293-dy2),
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+3,aco00q3);
		count=count+4;

		ma++;
	}

	count=41;
	for(Int_t ma=25;ma<=29;ma++)
	{
		TGeoTranslation *aco00=new TGeoTranslation("aco00",
					constants->ModulePositionX(ma)-0.5*293+dy2,
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count,aco00);

		TGeoTranslation *aco00q1=new TGeoTranslation("aco00q1",
					-(constants->ModulePositionX(ma)-0.5*293+dy2),
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+1,aco00q1);

		TGeoTranslation *aco00q2=new TGeoTranslation("aco00q2",
					constants->ModulePositionX(ma)+0.5*293-dy2,
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+2,aco00q2);

		TGeoTranslation *aco00q3=new TGeoTranslation("aco00q3",
					-(constants->ModulePositionX(ma)+0.5*293-dy2),
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		upFace->AddNode(aCORDE00,count+3,aco00q3);
		count=count+4;
		ma++;
	}

	TGeoTranslation *c1 = new TGeoTranslation ("c1",
				constants->ModulePositionX(20)-0.5*293,
				constants->ModulePositionY(20)-box[1]-z1,
				constants->ModulePositionZ(20)-40);
	TGeoTranslation *c2 = new TGeoTranslation ("c2",
				constants->ModulePositionX(23)-0.5*293,
				constants->ModulePositionY(23)-box[1]-z1,
				constants->ModulePositionZ(23)-40);
	TGeoTranslation *c3 = new TGeoTranslation ("c3",
				constants->ModulePositionX(24)-0.5*293,
				constants->ModulePositionY(24)-box[1]-z1,
				constants->ModulePositionZ(25)-40);
	TGeoTranslation *c4 = new TGeoTranslation ("c4",
				constants->ModulePositionX(27)-0.5*293,
				constants->ModulePositionY(27)-box[1]-z1,
				constants->ModulePositionZ(28)-40);
	upFace->AddNode(aCORDE00,57,c1);
	upFace->AddNode(aCORDE00,58,c2);
	upFace->AddNode(aCORDE00,59,c3);
	upFace->AddNode(aCORDE00,60,c4);


	// Construct Bars for lateral supports (up-face)

	TGeoTranslation *aco00=new TGeoTranslation("aco00",
				constants->ModulePositionX(20)+0.5*293-dy,
				constants->ModulePositionY(20)-box[1]-z1,
				constants->ModulePositionZ(20)-40);
	upFace->AddNode(aCORDE00,61,aco00);

	TGeoTranslation *aco00q1=new TGeoTranslation("aco00q1",
				constants->ModulePositionX(23)+0.5*293-dy,
				constants->ModulePositionY(23)-box[1]-z1,
				constants->ModulePositionZ(23)-40);
	upFace->AddNode(aCORDE00,62,aco00q1);

	TGeoTranslation *aco00q2=new TGeoTranslation("aco00q2",
				constants->ModulePositionX(24)+0.5*293-dy,
				constants->ModulePositionY(24)-box[1]-z1,
				constants->ModulePositionZ(25)-40);
	upFace->AddNode(aCORDE00,63,aco00q2);

	TGeoTranslation *aco00q3=new TGeoTranslation("aco00q3",
				constants->ModulePositionX(27)+0.5*293-dy,
				constants->ModulePositionY(27)-box[1]-z1,
				constants->ModulePositionZ(28)-40);
	upFace->AddNode(aCORDE00,64,aco00q3);


	TGeoTranslation *aco01=new TGeoTranslation("aco01",
				constants->ModulePositionX(30)-0.5*293+dy,
				constants->ModulePositionY(30)-box[1]-z1,
				constants->ModulePositionZ(30)-40);
	upFace->AddNode(aCORDE00,65,aco01);

	TGeoTranslation *aco01q1=new TGeoTranslation("aco01q1",
				constants->ModulePositionX(33)-0.5*293+dy,
				constants->ModulePositionY(33)-box[1]-z1,
				constants->ModulePositionZ(33)-40);
	upFace->AddNode(aCORDE00,66,aco01q1);

	TGeoTranslation *aco01q2=new TGeoTranslation("aco01q2",
				constants->ModulePositionX(34)-0.5*293+dy,
				constants->ModulePositionY(34)-box[1]-z1,
				constants->ModulePositionZ(35)-40);
	upFace->AddNode(aCORDE00,67,aco01q2);

	TGeoTranslation *aco01q3=new TGeoTranslation("aco01q3",
				constants->ModulePositionX(37)-0.5*293+dy,
				constants->ModulePositionY(37)-box[1]-z1,
				constants->ModulePositionZ(38)-40);
	upFace->AddNode(aCORDE00,68,aco01q3);



	// Acorde's support bars (side's faces)

	//*** In Face ***

//	box[0]=39;
	box[0]=27;
	box[1]=5;
	box[2]=5;
	Float_t kro=3;
	Float_t q1=0;
	Float_t posx=constants->ModulePositionX(0)+0.5*293*0.7071-56*0.7071-18;
	Float_t posy=constants->ModulePositionY(0)-0.5*293*0.7071-56*0.7071+3-q1+kro;
	Float_t posz=constants->ModulePositionZ(0);

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

		Float_t posx1=constants->ModulePositionX(dyA)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->ModulePositionY(dyA)-0.1*293*0.7071-56*0.7071+3-des-q1+kro;
		Float_t posza=constants->ModulePositionZ(dyA);
		Float_t posx2=constants->ModulePositionX(dyA)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->ModulePositionY(dyA)+0.27*293*0.7071-56*0.7071+3-des-q1+kro;
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

		Float_t posx1=constants->ModulePositionX(dyb)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->ModulePositionY(dyb)-0.1*293*0.7071-56*0.7071+3-des-q1+kro;
		Float_t poszb=constants->ModulePositionZ(dyb+10);
		Float_t posx2=constants->ModulePositionX(dyb)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->ModulePositionY(dyb)+0.27*293*0.7071-56*0.7071+3-des-q1+kro;
		TGeoCombiTrans *aco7q6 = new TGeoCombiTrans("aco7q6",posx1,posy1,poszb,idrotm231);
		TGeoCombiTrans *aco7q7 = new TGeoCombiTrans("aco7q7",posx2,posy2,poszb,idrotm231);
		inFace->AddNode(aCORDE7,count,aco7q6);
		inFace->AddNode(aCORDE7,count+1,aco7q7);
		count=count+2;
		dyb++;
	}	



	Float_t posxq1=constants->ModulePositionX(10)+0.5*293*0.7071-56*0.7071-18;
	Float_t posyq1=constants->ModulePositionY(10)-0.5*293*0.7071-56*0.7071+3-q1+kro;
	Float_t poszq1=constants->ModulePositionZ(10);
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

		Float_t posx1=constants->ModulePositionX(dyc)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->ModulePositionY(dyc)-0.1*293*0.7071-56*0.7071+3-des-0.8+kro;
		Float_t poszc=constants->ModulePositionZ(dyc);
		Float_t posx2=constants->ModulePositionX(dyc)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->ModulePositionY(dyc)+0.27*293*0.7071-56*0.7071+3-des-1.5-0.8+kro;
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

		Float_t posx1=constants->ModulePositionX(dyd)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->ModulePositionY(dyd)-0.1*293*0.7071-56*0.7071+3-des-q1-0.8+kro;
		Float_t poszd=constants->ModulePositionZ(dyd);
		Float_t posx2=constants->ModulePositionX(dyd)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->ModulePositionY(dyd)+0.27*293*0.7071-56*0.7071+3-des-1.5-q1-0.8+kro;
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
	Float_t posxqa=constants->ModulePositionX(50)-0.5*293*0.7071+56*0.7071+18;
	Float_t posyqa=constants->ModulePositionY(50)-0.5*293*0.7071-56*0.7071+3-s1+kro;
	Float_t poszqa=constants->ModulePositionZ(50);
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
					constants->ModulePositionX(50)-0.1*293*0.7071
					+56*0.7071+18-des,
					constants->ModulePositionY
					(50)-0.1*293*0.7071-56*0.7071+3-des-s1,
					constants->ModulePositionZ(45),idrotm232);
	TGeoCombiTrans *aco7q21 = new TGeoCombiTrans("aco7q21",
					constants->ModulePositionX(50)+0.27*293*0.7071
					+56*0.7071+18-des,
					constants->ModulePositionY(50)
					+0.27*293*0.7071-56*0.7071+3-des-s1,
					constants->ModulePositionZ(45),idrotm232);
	outFace->AddNode(aCORDE7,19,aco7q16);
	outFace->AddNode(aCORDE7,20,aco7q17);
	outFace->AddNode(aCORDE7,21,aco7q18);
	outFace->AddNode(aCORDE7,22,aco7q19);
	outFace->AddNode(aCORDE7,23,aco7q20);
	outFace->AddNode(aCORDE7,24,aco7q21);


	count=25;
	for(Int_t dye=50;dye<=54;dye++)
	{

		Float_t posx1=constants->ModulePositionX(dye)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->ModulePositionY(dye)-0.1*293*0.7071-56*0.7071+3-des-s1+kro;
		Float_t posze=constants->ModulePositionZ(dye);
		Float_t posx2=constants->ModulePositionX(dye)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->ModulePositionY(dye)+0.27*293*0.7071-56*0.7071+3-des-s1+kro;
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

		Float_t posx1=constants->ModulePositionX(dyf)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->ModulePositionY(dyf)-0.1*293*0.7071-56*0.7071+3-des-s1+kro;
		Float_t poszf=constants->ModulePositionZ(dyf-10);
		Float_t posx2=constants->ModulePositionX(dyf)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->ModulePositionY(dyf)+0.27*293*0.7071-56*0.7071+3-des-s1+kro;
		TGeoCombiTrans *aco7q24 = new TGeoCombiTrans("aco7q24",posx1,posy1,poszf,idrotm232);
		TGeoCombiTrans *aco7q25 = new TGeoCombiTrans("aco7q25",posx2,posy2,poszf,idrotm232);
		outFace->AddNode(aCORDE7,count,aco7q24);
		outFace->AddNode(aCORDE7,count+1,aco7q25);
		count=count+2;
		dyf++;
	}


	Float_t posxqb=constants->ModulePositionX(40)-0.5*293*0.7071+56*0.7071+18;
	Float_t posyqb=constants->ModulePositionY(40)-0.5*293*0.7071-56*0.7071+3-s1+kro;
	Float_t poszqb=constants->ModulePositionZ(40);
	TGeoCombiTrans *aco7q26 = new TGeoCombiTrans("aco7q26",
					posxqb,posyqb,poszqb-4*dy,idrotm232);
	TGeoCombiTrans *aco7q27 = new TGeoCombiTrans("aco7q27",
					posxqb,posyqb,
					constants->ModulePositionZ(43)-4*dy,idrotm232);
	TGeoCombiTrans *aco7q28 = new TGeoCombiTrans("aco7q28",
					posxqb,posyqb,
					constants->ModulePositionZ(45)-4*dy,idrotm232);
	TGeoCombiTrans *aco7q29 = new TGeoCombiTrans("aco7q29",posxqb,posyqb,
					constants->ModulePositionZ(48)-4*dy,idrotm232);
	outFace->AddNode(aCORDE7,41,aco7q26);
	outFace->AddNode(aCORDE7,42,aco7q27);
	outFace->AddNode(aCORDE7,43,aco7q28);
	outFace->AddNode(aCORDE7,44,aco7q29);

	count=45;
	for(Int_t dyg=40;dyg<=44;dyg++)
	{

		Float_t posx1=constants->ModulePositionX(dyg)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->ModulePositionY(dyg)-0.1*293*0.7071-56*0.7071+3-des-s1+kro;
		Float_t poszg=constants->ModulePositionZ(dyg);
		Float_t posx2=constants->ModulePositionX(dyg)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->ModulePositionY(dyg)+0.27*293*0.7071-56*0.7071+3-des-s1+kro;
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

		Float_t posx1=constants->ModulePositionX(dyh)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->ModulePositionY(dyh)-0.1*293*0.7071-56*0.7071+3-des-s1+kro;
		Float_t poszh=constants->ModulePositionZ(dyh);
		Float_t posx2=constants->ModulePositionX(dyh)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->ModulePositionY(dyh)+0.27*293*0.7071-56*0.7071+3-des-s1+kro;
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
	Float_t posx1=constants->ModulePositionX(0)+0.5*293*0.7071-4*box[0]-8+re;
	Float_t posy1=constants->ModulePositionY(0)-0.5*293*0.7071-box[1]-18-2+sm;
	Float_t posz1=constants->ModulePositionZ(0);

	TGeoBBox *acorde7q1 = new TGeoBBox("acorde7q1",box[0],box[1],box[2]);

	TGeoVolume *aCORDE7q1 = new TGeoVolume("ACORDE7_1",acorde7q1,al);
	TGeoTranslation *aco71 = new TGeoTranslation("aco71",posx1,posy1,posz1-4*dy);
	TGeoTranslation *aco72 = new TGeoTranslation("aco72",posx1,posy1,
					constants->ModulePositionZ(3)-4*dy);
	TGeoTranslation *aco73 = new TGeoTranslation("aco73",posx1,posy1,
					constants->ModulePositionZ(5)-4*dy);
	TGeoTranslation *aco74 = new TGeoTranslation("aco74",posx1,posy1,
					constants->ModulePositionZ(8)-4*dy);
	inFace->AddNode(aCORDE7q1,67,aco71);
	inFace->AddNode(aCORDE7q1,68,aco72);
	inFace->AddNode(aCORDE7q1,69,aco73);
	inFace->AddNode(aCORDE7q1,70,aco74);


	count=71;
	for(Int_t dyi=0;dyi<=4;dyi++)
	{

		Float_t posx1a=constants->ModulePositionX(dyi)+0.1*293*0.7071-4*box[0]-8+des+re;
		Float_t posy1a=constants->ModulePositionY(dyi)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1a=constants->ModulePositionZ(dyi);
		Float_t dyx2=constants->ModulePositionX(dyi)-0.27*293*0.7071-4*box[0]-8+des+re;
		Float_t dyy2=constants->ModulePositionY(dyi)+0.27*293*0.7071-box[1]-18-2-des+sm;
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

		Float_t posx1b=constants->ModulePositionX(dyj)+0.1*293*0.7071-4*box[0]-8+des+re;
		Float_t posy1b=constants->ModulePositionY(dyj)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1b=constants->ModulePositionZ(dyj+10);
		Float_t dyx2=constants->ModulePositionX(dyj)-0.27*293*0.7071-4*box[0]-8+des+re;
		Float_t dyy2=constants->ModulePositionY(dyj)+0.27*293*0.7071-box[1]-18-2-des+sm;
		TGeoTranslation *aco75=new TGeoTranslation("aco75",posx1b,posy1b,posz1b);
		TGeoTranslation *aco76=new TGeoTranslation("aco76",dyx2,dyy2,posz1b);
		inFace->AddNode(aCORDE7q1,count,aco75);
		inFace->AddNode(aCORDE7q1,count+1,aco76);
		count=count+2;
		dyj++;
	}


	Float_t posx1q1=constants->ModulePositionX(10)+0.5*293*0.7071-4*box[0]-8+re;
	Float_t posy1q1=constants->ModulePositionY(10)-0.5*293*0.7071-box[1]-18-2+sm;
	Float_t posz1q1=constants->ModulePositionZ(10);
	TGeoTranslation *aco77=new TGeoTranslation("aco77",posx1q1,posy1q1,posz1q1-4*dy);
	TGeoTranslation *aco78=new TGeoTranslation("aco78",posx1q1,posy1q1,
					constants->ModulePositionZ(13)-4*dy);

	TGeoTranslation *aco79=new TGeoTranslation("aco79",posx1q1,posy1q1,
					constants->ModulePositionZ(15)-4*dy);
	TGeoTranslation *aco710=new TGeoTranslation("aco710",posx1q1,posy1q1,
					constants->ModulePositionZ(18)-4*dy);
	inFace->AddNode(aCORDE7q1,91,aco77);
	inFace->AddNode(aCORDE7q1,92,aco78);
	inFace->AddNode(aCORDE7q1,93,aco79);
	inFace->AddNode(aCORDE7q1,94,aco710);

	count=95;
	for(Int_t dyk=10;dyk<=14;dyk++)
	{

		Float_t posx1c=constants->ModulePositionX(dyk)+0.1*293*0.7071-4*box[0]-8+des+re+.83;
		Float_t posy1c=constants->ModulePositionY(dyk)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1c=constants->ModulePositionZ(dyk);
		Float_t dyx2=constants->ModulePositionX(dyk)-0.27*293*0.7071-4*box[0]-4+des+re+0.83;
		Float_t dyy2=constants->ModulePositionY(dyk)+0.27*293*0.7071-box[1]-18-5-des+sm;
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

		Float_t posx1d=constants->ModulePositionX(dyl)+0.1*293*0.7071-4*box[0]-8+des+re+0.83;
		Float_t posy1d=constants->ModulePositionY(dyl)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1d=constants->ModulePositionZ(dyl);
		Float_t dyx2=constants->ModulePositionX(dyl)-0.27*293*0.7071-4*box[0]-4+des+re+0.83;
		Float_t dyy2=constants->ModulePositionY(dyl)+0.27*293*0.7071-box[1]-18-5-des;
		TGeoTranslation *aco713=new TGeoTranslation("aco713",posx1d,posy1d,posz1d);
		TGeoTranslation *aco714=new TGeoTranslation("aco714",dyx2,dyy2,posz1d);
		inFace->AddNode(aCORDE7q1,count,aco713);
		inFace->AddNode(aCORDE7q1,count+1,aco714);
		count=count+2;
		dyl++;
	}

		//*** Out-Face ***

	Float_t posx1qa=constants->ModulePositionX(50)-0.5*293*0.7071+4*box[0]+8-re-1;
	Float_t posy1qa=constants->ModulePositionY(50)-0.5*293*0.7071-box[1]-18-2+sm-2.5;
	Float_t posz1qa=constants->ModulePositionZ(50);
	TGeoTranslation *aco715=new TGeoTranslation("aco715",posx1qa,posy1qa,posz1qa-4*dy);
	TGeoTranslation *aco716=new TGeoTranslation("aco716",posx1qa,posy1qa,
				constants->ModulePositionZ(43)-4*dy);
	TGeoTranslation *aco717=new TGeoTranslation("aco717",posx1qa,posy1qa,
				constants->ModulePositionZ(55)-4*dy);
	TGeoTranslation *aco718=new TGeoTranslation("aco718",posx1qa,posy1qa,
				constants->ModulePositionZ(58)-4*dy);
	TGeoTranslation *aco719=new TGeoTranslation("aco719",
				constants->ModulePositionX(50)-0.1*293*0.7071+4*box[0]+8-des-re-1,		
				constants->ModulePositionY(50)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5,
				constants->ModulePositionZ(45));
	TGeoTranslation *aco720=new TGeoTranslation("aco720",
				constants->ModulePositionX(50)+0.27*293*0.7071+4*box[0]+8-des-re-1,
				constants->ModulePositionY(50)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5,
				constants->ModulePositionZ(45));


	outFace->AddNode(aCORDE7q1,115,aco715);
	outFace->AddNode(aCORDE7q1,116,aco716);
	outFace->AddNode(aCORDE7q1,117,aco717);
	outFace->AddNode(aCORDE7q1,118,aco718);
	outFace->AddNode(aCORDE7q1,119,aco719);
	outFace->AddNode(aCORDE7q1,120,aco720);




	count=65;
	for(Int_t dym=50;dym<=54;dym++)
	{

		Float_t posx1e=constants->ModulePositionX(dym)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1e=constants->ModulePositionY(dym)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1e=constants->ModulePositionZ(dym);
		Float_t dyx2=constants->ModulePositionX(dym)+0.27*293*0.7071+4*box[0]+8-des-re-1;
		Float_t dyy2=constants->ModulePositionY(dym)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5;
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

		Float_t posx1f=constants->ModulePositionX(dyn)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1f=constants->ModulePositionY(dyn)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1f=constants->ModulePositionZ(dyn-10);
		Float_t dyx2=constants->ModulePositionX(dyn)+0.27*293*0.7071+4*box[0]+8-des-re-1;
		Float_t dyy2=constants->ModulePositionY(dyn)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5;
		TGeoTranslation *aco723=new TGeoTranslation("aco723",posx1f,posy1f,posz1f);
		TGeoTranslation *aco724=new TGeoTranslation("aco724",dyx2,dyy2,posz1f);
		outFace->AddNode(aCORDE7q1,count,aco723);
		outFace->AddNode(aCORDE7q1,count+1,aco724);
		count=count+2;
		dyn++;
	}


	Float_t posx1qb=constants->ModulePositionX(40)-0.5*293*0.7071+4*box[0]+5;
	Float_t posy1qb=constants->ModulePositionY(40)-0.5*293*0.7071-box[1]-18-2;
	Float_t posz1qb=constants->ModulePositionZ(40);
	TGeoTranslation *aco725=new TGeoTranslation("aco725",posx1qb,posy1qb,posz1qb-4*dy);
	TGeoTranslation *aco726=new TGeoTranslation("aco726",posx1qb,posy1qb,
				constants->ModulePositionZ(43)-4*dy);
	TGeoTranslation *aco727=new TGeoTranslation("aco727",posx1qb,posy1qb,
				constants->ModulePositionZ(45)-4*dy);
	TGeoTranslation *aco728=new TGeoTranslation("aco728",posx1qb,posy1qb,
				constants->ModulePositionZ(48)-4*dy);
	outFace->AddNode(aCORDE7q1,85,aco725);
	outFace->AddNode(aCORDE7q1,86,aco726);
	outFace->AddNode(aCORDE7q1,87,aco727);
	outFace->AddNode(aCORDE7q1,88,aco728);



	count=89;
	for(Int_t dyo=40;dyo<=44;dyo++)
	{

		Float_t posx1g=constants->ModulePositionX(dyo)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1g=constants->ModulePositionY(dyo)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1g=constants->ModulePositionZ(dyo);
		Float_t dyx2=constants->ModulePositionX(dyo)+0.27*293*0.7071+4*box[0]+4-des-re-1+2.8;
		Float_t dyy2=constants->ModulePositionY(dyo)+0.27*293*0.7071-box[1]-18-5-des+sm-2.5+3;
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

		Float_t posx1h=constants->ModulePositionX(dyp)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1h=constants->ModulePositionY(dyp)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1h=constants->ModulePositionZ(dyp);
		Float_t dyx2=constants->ModulePositionX(dyp)+0.27*293*0.7071+4*box[0]+4-des-re-1+2.8;
		Float_t dyy2=constants->ModulePositionY(dyp)+0.27*293*0.7071-box[1]-18-5-des+sm-2.5+3;
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

	Float_t sx=constants->ModulePositionX(24)-0.5*293;
	Float_t sy=constants->ModulePositionY(24)-box[1];
	Float_t sz=0;
	Float_t sx2=constants->ModulePositionX(24)+0.5*293-dy;
	Float_t sy2=constants->ModulePositionY(24)-box[1];
	Float_t sx4=constants->ModulePositionX(24)-0.5*293+dy2;
	Float_t sy4=constants->ModulePositionY(24)-box[1];
	Float_t sx5=constants->ModulePositionX(24)+0.5*293-dy2;
	Float_t sy5=constants->ModulePositionY(24)-box[1];

	Float_t dyx=constants->ModulePositionX(4)+0.5*293*0.7071-box[0];
	Float_t dyy=constants->ModulePositionY(4)-0.5*293*0.7071-box[1];
	Float_t dyz=0;
	Float_t dyx1=constants->ModulePositionX(4)+0.1*293*0.7071-box[0];
	Float_t dyy1=constants->ModulePositionY(4)-0.1*293*0.7071-box[1];
	Float_t dyx2=constants->ModulePositionX(4)-0.27*293*0.7071-box[0];
	Float_t dyy2=constants->ModulePositionY(4)+0.27*293*0.7071-box[1];


	Float_t dx1=constants->ModulePositionX(14)+0.5*293*0.7071-box[0];
	Float_t dy11=constants->ModulePositionY(14)-0.5*293*0.7071-box[1];
	Float_t dyx11=constants->ModulePositionX(14)+0.1*293*0.7071-box[0];
	Float_t dyy11=constants->ModulePositionY(14)-0.1*293*0.7071-box[1];
	Float_t dyx21=constants->ModulePositionX(14)-0.27*293*0.7071-box[0];
	Float_t dyy21=constants->ModulePositionY(14)+0.27*293*0.7071-box[1];


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
//	aCORDE->AddNode(inFacem,4);//---> volume of modules in-face
//	aCORDE->AddNode(upFacem,5);//---> volume of modules up-face
//	aCORDE->AddNode(outFacem,6);//---> volume of modules out-face
        aCORDE->AddNode(modules,4);//---> volume of ALL ACORDE's Modules
	alice->AddNode(aCORDE,1);//---> put volume of ACORDE over ALICE's volume



}


//_____________________________________________________________________________
void AliACORDEv1::DrawDetector() const
{

  // not needed anymore

}

//____________________________________________________________________________

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



//_____________________________________________________________________________
void AliACORDEv1::MakeBranch(Option_t *option)
{
// Creates new branches in the current Root Tree
    
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  AliDebug(2,Form("fBufferSize = %d",fBufferSize));
  const char *cH = strstr(option,"H");
  if (fHits   && TreeH() && cH) {
    TreeH()->Branch(branchname,&fHits, fBufferSize);
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
