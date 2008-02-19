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
//  This class contains the functions for version 0 of the ALICE Cosmic Ray  //
//  Trigger. This version will be used to simulation comic rays in alice with//
//  all the detectors. It include geometry and hits (position and momentum)  //
//                                                                           //
//                  Send comments to:                                        //
//									     //
//      Arturo Fernandez Tellez   	<afernand@fcfm.buap.mx>              //
//      Enrique Gamez             	<egamez@fcfm.buap.mx>                //
//      Eleazar Cuautle Flores    	<ecuautle@nucleares.unam.mx>         //
//	Mario Rodriguez	Cahuantzi 	<mrodrigu@mail.cern.ch>		     //	
//									     //
//			Puebla, Pue. Mexico December 2007                    //
///////////////////////////////////////////////////////////////////////////////


//
#include <Riostream.h>
#include <TGeoMatrix.h>
#include <TGeometry.h>
#include <TMath.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TVirtualMC.h>
#include <TString.h>
#include <TSystem.h>

#include "AliConst.h"
#include "AliRun.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoPcon.h"
#include "TGeoTube.h"
#include "TGeoPgon.h"
#include "TGeoTrd1.h"
#include "TGeoCompositeShape.h"
#include "TGeoPara.h"
//

#include "AliACORDEv1.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVirtualMC.h>
#include <TPDGCode.h>
#include <TGeometry.h>
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

	TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);

	//---> define the measures

	Double_t dx1 = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	Double_t dy1 = AliACORDEConstants::Instance()->Depth();
	Double_t dz1 = dx1;

	//---> define the box for the mother volume

	TGeoVolume *ACORDE = acorde->MakeBox("ACORDE", Vacuum, dx1, dy1, dz1);
	acorde->SetTopVolume(ACORDE);

	//---> create shafts&molasse

	CreateShafts();
	CreateMolasse();
}

void AliACORDEv1::CreateShafts()
{


	TGeoManager *acorde = new TGeoManager("ACORDE2007", "Geometry of ACORDE");	

	//---> define some materials

	TGeoMaterial *matVacuum = new TGeoMaterial("Al", 0,0,0);

	//---> define some media

	TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
	TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
	TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
	
	//---> This shaft is composes by an open tube down in the hall
	//---> and a cilinder above the level of the celling
	//---> Every structure relative to the shaft will be put into this volume

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

	TGeoVolume *O = acorde->MakeBox("O", Vacuum, 25., 25., 5.);
	TGeoVolume *CSF1 = acorde->MakeTubs("CSF1",Al,ptube[0],ptube[1],ptube[2],ptube[3],ptube[4]);
	TGeoVolume *CSF2 = acorde->MakeTubs("CSF2",Al,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	O->AddNode(CSF1, 1);
	CSF1->AddNode(CSF2,1,tr2);

	//---> definition of the other part of the shaft

  	ptube[0] = ptubs[0]; // Inner radius
	ptube[1] = ptubs[1]; // Outer radius
	ptube[2] = 5150/2 - ptubs[2]; // Half lenght
	TGeoVolume *CSF3 = acorde->MakeTubs("CSF3",Al,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	TGeoTranslation *tr3 = new TGeoTranslation(0,0,5150/2-ptube[2]);
	CSF1->AddNode(CSF3,1,tr3);

	//---> define concrete walls along the shaft (next to the elevator)
	
	Float_t pbox[3];
	pbox[0]=480/2;
	pbox[1]=120/2;
	pbox[2]=5150/2;
	TGeoVolume *CSW1 = acorde->MakeBox("CSW1",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br1 = new TGeoTranslation(820+pbox[0],150+pbox[1],0);
	TGeoTranslation *br1_a = new TGeoTranslation(820+pbox[0],-300-pbox[1],0);
	CSF1->AddNode(CSW1,1,br1);
	CSF1->AddNode(CSW1,1,br1_a);

	pbox[0] = 120/2;  // Half length in X
	pbox[1] = 750/2;  // Half length in Y
	pbox[2] = 5150/2; // Half length in Z
	TGeoVolume *CSW2 = acorde->MakeBox("CSW2",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br2 = new TGeoTranslation(820-60,150+pbox[1],0);
	CSF1->AddNode(CSW2,1,br2);


	pbox[0] = 120/2;  // Half length in X
	pbox[1] = 600/2;  // Half lenght in Y
	pbox[2] = 5150/2; // Half length in Z
	TGeoVolume *CSW3 = acorde->MakeBox("CSW3",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br3 = new TGeoTranslation(820-60,-300-pbox[1],0);
	CSF1->AddNode(CSW3,1,br3);

	pbox[0] = 400/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 300/2; // Half length in Z
	TGeoVolume *CSW4 = acorde->MakeBox("CSW4",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br4 = new TGeoTranslation(pbox[1]-pbox[0],0,3000-5150/2-pbox[2]);
	CSF1->AddNode(CSW4,1,br4);


	pbox[0] = 1400/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 170/2; // Half length in Z
	TGeoVolume *CSW5 = acorde->MakeBox("CSW5",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br5 = new TGeoTranslation(0,0,3000-5150/2-130);
	CSF1->AddNode(CSW5,1,br5);


	pbox[0] = 170/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 300/2; // Half length in Z
	TGeoVolume *CSW6 = acorde->MakeBox("CSW6",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br6 = new TGeoTranslation(-1400/2-pbox[0],0,3000-5150/2-pbox[2]);
	CSF1->AddNode(CSW6,1,br6);


	pbox[0] = 100/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 450/2; // Half length in Z
	TGeoVolume *CSW7 = acorde->MakeBox("CSW7",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br7 = new TGeoTranslation(-1400/2-170-pbox[0],0,3000-5150/2+pbox[2]);
	CSF1->AddNode(CSW7,1,br7);


	pbox[0] = 300/2;  // Half length in X
	pbox[1] = 2300/2;  // Half lenght in Y
	pbox[2] = 170/2; // Half length in Z
	TGeoVolume *CSW8 = acorde->MakeBox("CSW8",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *br8 = new TGeoTranslation(-2300/2+pbox[0],0,2500-5150/2);
	CSF1->AddNode(CSW8,1,br8);

	//---> put the shaft into the mother volume

	TGeoCombiTrans *br = new TGeoCombiTrans(0,AliACORDEConstants::Instance()->Depth()-5150/2,2300,rot1);
	CSF1->AddNode(CSF1,1,br);


	//---> PM25 Access Shafts

	ptube[0]=910/2;
	ptube[1]=ptube[0]+100;
	ptube[2]=(5150-1166)/2;
	TGeoVolume *CSF4 = acorde->MakeTubs("CSF4",Vacuum,pbox[0],pbox[1],pbox[2],360,360);
	TGeoCombiTrans *tr4 = new TGeoCombiTrans(2100,AliACORDEConstants::Instance()->Depth()-ptube[2],0,rot1);
	CSF4->AddNode(CSF4,1,tr4);


	//---> PGC2 Access shaft

	ptube[0]=1100/2;
	ptube[1]=ptube[0]+100;
	ptube[2]=(5150-690)/2;
	TGeoVolume *CSF5 = acorde->MakeTubs("CSF5",Vacuum,pbox[0],pbox[1],pbox[2],360,360);
	TGeoCombiTrans *tr5 = new TGeoCombiTrans(-375,AliACORDEConstants::Instance()->Depth()-ptube[2],-1900-2987.7,rot1);
	CSF5->AddNode(CSF5,1,tr5);

}


void AliACORDEv1::CreateMolasse()

{

	TGeoManager *acorde = new TGeoManager("ACORDE2007", "Geometry of ACORDE");	

	//---> define some media


	TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
	TGeoMedium *Vacuum = new TGeoMedium("Root Material",2, matAl);

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
	TGeoVolume *CM01 = acorde->MakeBox("CM01", Vacuum, pbox[0],pbox[1],pbox[2]);

	//---> Now put the molasse exactly above the hall. OK
	//---> Above the ceiling
	
	Float_t ptubs[5];
	ptubs[0] = 1170;
	ptubs[1] = 2100 - pm25radius;
	ptubs[2] = 1900/2 + px24radius;
	ptubs[3] = 0;
	ptubs[4] = 180;
	TGeoVolume *CM02 = acorde->MakeTubs("CM02",Vacuum,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	TGeoTranslation *tr2 = new TGeoTranslation(0,500-AliACORDEConstants::Instance()->Depth()/2,ptubs[2]-1900);
	CM01->AddNode(CM02,1,tr2);


	//---> Molasse around the RB24/26 Wall. OK

	ptubs[0] = 220 + 1600;
	ptubs[1] = AliACORDEConstants::Instance()->Depth() - ptubs[0];
	ptubs[2] = 2987.7/2 - 1100/4 - concreteWidth/2;
	ptubs[3] = 0;
	ptubs[4] = 180;
	TGeoVolume *CM03 = acorde->MakeTubs("CM03",Vacuum,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	TGeoTranslation *tr3 = new TGeoTranslation(70,40-AliACORDEConstants::Instance()->Depth()/2,-ptubs[2]-1900);
	CM01->AddNode(CM03,1,tr3);


	//---> A big block above the RB24/26 wall. OK

	pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	pbox[1] = (AliACORDEConstants::Instance()->Depth() - 220 - 1600)/2;
	pbox[2] = 2987.7/2 - 1100/4 - concreteWidth/2;
	TGeoVolume *CM04 = acorde->MakeBox("CM04", Vacuum, pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr4 = new TGeoTranslation(0,AliACORDEConstants::Instance()->Depth()/2-pbox[1],-1900-pbox[2]);
	CM01->AddNode(CM04,1,tr4);



	//---> Small blocks below the volume CMO4 on both sides of the wall RB24/26. OK

	pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)-ptubs[0])/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2 - pbox[1];
	TGeoVolume *CM17 = acorde->MakeBox("CM17", Vacuum, pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr17 = new TGeoTranslation(AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - pbox[0],-AliACORDEConstants::Instance()->Depth()/2 + pbox[1],-1900 - pbox[2]);
	TGeoTranslation *tr17_a = new TGeoTranslation(-AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)+ pbox[0],-AliACORDEConstants::Instance()->Depth()/2 + pbox[1], -1900 - pbox[2]);
	CM01->AddNode(CM17,1,tr17);
	CM01->AddNode(CM17,2,tr17_a);


	//---> And a big block of molasse above the hall up to the surface. OK

	pbox[0] = pm25X - pm25radius;
	pbox[1] = (AliACORDEConstants::Instance()->Depth()-500-1170)/2;
	pbox[2] = (1900 + 1150)/2;
	TGeoVolume *CM05 = acorde->MakeBox("CM05", Vacuum, pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr5 = new TGeoTranslation(0,AliACORDEConstants::Instance()->Depth()/2-pbox[1], pbox[2]-1900);
	CM01->AddNode(CM05,1,tr5);


	//---> Small blocks of molasse betwen the blocks CMO2, CMO5 and PM25. Ok

	pbox[0] = (pm25X - pm25radius - 1170)/2;
	pbox[1] = 1000;
	TGeoVolume *CM16 = acorde->MakeBox("CM16", Vacuum, pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr16 = new TGeoTranslation(1170 + pbox[0], -AliACORDEConstants::Instance()->Depth()/2+pbox[1], pbox[2] - 1900);
	CM01->AddNode(CM16,1,tr16);


	//---> Molasse around the shafts.

	TGeoRotation *rot2 = new TGeoRotation("rot1",0, 0, 90, 0, 90, 90 );

	//---> Around the PX24, the open section. OK

	ptubs[0] = px24radius + concreteWidth;
	ptubs[1] = ptubs[0] + 1000;
	ptubs[2] = (2300 - (5150 - AliACORDEConstants::Instance()->Depth()))/2;
	ptubs[3] = 180 + kRaddeg*TMath::ASin(1070/ptubs[0]);
	ptubs[4] = 180 -  kRaddeg*TMath::ASin(1070/ptubs[0]);
	TGeoVolume *CM06 = acorde->MakeTubs("CM06", Vacuum,ptubs[0],ptubs[1],ptubs[2],ptubs[3],ptubs[4]);
	TGeoTranslation *tr6 = new TGeoTranslation(px24X, ptubs[2] - AliACORDEConstants::Instance()->Depth()/2, px24Z);
	CM01->AddNode(CM06,1,tr6);


	//---> Around the PX24, the closed section. OK

	Float_t ptube[3];
	ptube[0] = px24radius + concreteWidth;
	ptube[1] = ptube[0] + 1000;
	ptube[2] = (5150 - 2300)/2;
	TGeoVolume *CM07 = acorde->MakeTubs("CM07", Vacuum,ptube[0],ptube[1],ptubs[2],ptube[3],ptube[4]);
	TGeoTranslation *tr7 = new TGeoTranslation(px24X, AliACORDEConstants::Instance()->Depth()/2-ptube[2], px24Z);
	CM01->AddNode(CM07,1,tr7);


	//---> Around PM25. OK

	ptube[0] = pm25radius + concreteWidth;
	ptube[1] = ptube[0] + 400;
	ptube[2] = AliACORDEConstants::Instance()->Depth()/2;
	TGeoVolume *CM08 = acorde->MakeTubs("CM08", Vacuum,ptube[0],ptube[1],ptube[2],ptube[3],ptube[4]);
	TGeoCombiTrans *tr8 = new TGeoCombiTrans(pm25X, 0, pm25Z,rot2);
	CM01->AddNode(CM08,1,tr8);


	//---> On both sides of the PM25 along the HALL.

	pbox[0] = (2100 + pm25radius - 1170)/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (3*px24radius - pm25radius)/2;
	TGeoVolume *CM18 = acorde->MakeBox("CM18",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr18 = new TGeoTranslation(2100, 0, pbox[2] + pm25radius);
	CM01->AddNode(CM18,1,tr18);
  
  	pbox[2] = (1900 - pm25radius)/2;
	TGeoVolume *CM19 = acorde->MakeBox("CM19",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr19 = new TGeoTranslation(2100, 0, -pbox[2] - pm25radius);
  	CM01->AddNode(CM19,1,tr19);


	//---> Around the PGC2. OK

	ptube[0] = pgc2radius + concreteWidth;
	ptube[1] = 2987.7 - 740;
	ptube[2] = AliACORDEConstants::Instance()->Depth()/2;
	TGeoVolume *CM09 = acorde->MakeTubs("CM09",Vacuum,ptube[0],ptube[1],ptube[2],ptube[3],ptube[4]);
	TGeoCombiTrans *tr09 = new TGeoCombiTrans(pgc2X, 0, pgc2Z,rot2);
	CM01->AddNode(CM09,1,tr09);

	//---> On both sides of the PGC2.OK

	pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)-1100 - 375)/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = pgc2radius + concreteWidth;
	TGeoVolume *CM10 = acorde->MakeBox("CM10",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr10 = new TGeoTranslation(AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - pbox[0], 0, pgc2Z);
	TGeoTranslation *tr10_a = new TGeoTranslation(-AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) + pbox[0], 0, pgc2Z);
	CM01->AddNode(CM10,1,tr10);
	CM01->AddNode(CM10,2,tr10_a);


	//---> big block of molasse behind the PX24. OK

	pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (pbox[0] - (2300 + 1150 + 100))/2;
	TGeoVolume *CM12 = acorde->MakeBox("CM12",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr12 = new TGeoTranslation(px24X, 0, px24Z + px24radius + concreteWidth + pbox[2]);	
	CM01->AddNode(CM12,1,tr12);


	//---> big block of molasse in the opposite side of the PM25. OK

	pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)-1150)/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (1900 + 2300 + 1150)/2;
	TGeoVolume *CM13 = acorde->MakeBox("CM13",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr13 = new TGeoTranslation(-1150 - pbox[0], 0, pbox[2] - 1900);	
	CM01->AddNode(CM13,1,tr13);
 

	//---> big block of molasse behind the PM25. OK

	pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)-(2100 + 910/2 + 100))/2;
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (1900 + 2300 + 1150)/2;
	TGeoVolume *CM14 = acorde->MakeBox("CM14",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr14 = new TGeoTranslation(pm25X + pm25radius + concreteWidth + pbox[0], 0, pbox[2] - 1900);
	CM01->AddNode(CM14,1,tr14);


	//---> big block of molasse behind the PGC2. OK

	pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
	pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
	pbox[2] = (pbox[0] - (2987.7 + 1900 + 1100/2 + 100))/2;
	TGeoVolume *CM15 = acorde->MakeBox("CM15",Vacuum,pbox[0],pbox[1],pbox[2]);
	TGeoTranslation *tr15 = new TGeoTranslation(0, 0, -pbox[0] + pbox[2]);
	TGeoTranslation *tr15_a = new TGeoTranslation(0,AliACORDEConstants::Instance()->Depth()/2,0);
	CM01->AddNode(CM15,1,tr15);
	CM01->AddNode(CM01,1,tr15_a);

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
  // |		ACORDE_MODULE--> Volume that contains a full Acorde's module	     |	
  // |		ACORDE_1--> Volume that contains the bars&supports in-face	     |
  // |		ACORDE_2--> Volume that contains the bars&supports up-face	     | 	
  // |		ACORDE_3--> Volume that contains the bars&supports out-face	     |
  // |		ACORDE_4--> Volume that contains the modules of ACORDE in-face	     |
  // |		ACORDE_5--> Volume that contains the modules of ACORDE up-face	     |
  // |		ACORDE_6--> Volume that contains the modules of ACORDE out-face	     |
  // |										     |	
  // |_______________________________________________________________________________|


	// Call the global constants for the Modules
	
	AliACORDEConstants* constants = AliACORDEConstants::Instance();

	// Get the Alice Volume

	TGeoVolume *alice = gGeoManager->GetVolume("ALIC");

	// Define some materials & medium

	//*** Support & Bars***

	TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
	TGeoMedium *Al = new TGeoMedium("Al",2, matAl);

	//*** Scintillators ***

	TGeoMixture *mat17 = new TGeoMixture("SCINTILLATOR",2,1.03200);
    	mat17->SetUniqueID(17);
    	mat17->DefineElement(0,12.01,6,0.9153266);
    	mat17->DefineElement(1,1.01,1,0.8467343E-01);
	TGeoMedium *med6 = new TGeoMedium("SCINTILLATOR",6,17,1,0,0,1,
					  0.1000000E+11,0.3081371E-01,
					  0.1000000E-01,0.8135138E-02);
	
	//Define a FULL-ACORDE-VOLUME

	TGeoVolume *ACORDE = new TGeoVolumeAssembly("ACORDE");


	// Define 6 master volumes for ACORDE
	
	TGeoVolume *in_face = new TGeoVolumeAssembly("ACORDE_1");
	TGeoVolume *up_face = new TGeoVolumeAssembly("ACORDE_2");
	TGeoVolume *out_face = new TGeoVolumeAssembly("ACORDE_3");

	TGeoVolume *in_facem = new TGeoVolumeAssembly("ACORDE_4");
	TGeoVolume *up_facem = new TGeoVolumeAssembly("ACORDE_5");
	TGeoVolume *out_facem = new TGeoVolumeAssembly("ACORDE_6");


	// Define global variables

	Float_t box[3];
	Int_t count;
	Float_t dy=10;//-->displacement of the support and bars of ACORDE
	Float_t dy2=66.5;//-->displacement of the support and bars of ACORDE
	Float_t placed_at;
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
	TGeoVolume *ACORDE1_a = new TGeoVolume("ACORDE1_a",acorde1,Al);
	TGeoVolume *ACORDE10 = new TGeoVolume("ACORDE10",acorde10,Al);

	//*** Scintillators ***

	TGeoBBox *acorde2 = new TGeoBBox("acorde2",pbox[0],pbox[1],pbox[2]);
	TGeoVolume *ACORDE2 = new TGeoVolume("ACORDE2",acorde2,med6);


	// Here I define & construct a Master Volume ("ACORDE_MODULE") for one Module in ACORDE

	TGeoVolume *module = new TGeoVolumeAssembly("ACORDE_MODULE");
	module->AddNode(ACORDE1_a,1,new TGeoTranslation("aco1",0,0,13));
	module->AddNode(ACORDE1_a,2,new TGeoTranslation("aco10",0,0,-13));
	module->AddNode(ACORDE10,3,new TGeoTranslation("aco10",293/2+5,0,0));
	module->AddNode(ACORDE10,4,new TGeoTranslation("aco10",-293/2-5,0,0));
        placed_at = pbox[1]+constants->ProfileThickness()-constants->ModuleHeight()/2+small;
	module->AddNode(ACORDE2,5,new TGeoTranslation("aco2",placed_at,0,0));
        placed_at = placed_at + 2.0*pbox[1]+small;
	module->AddNode(ACORDE2,6,new TGeoTranslation("aco2",placed_at,-1,0));
	Float_t w1 = 8;

	


	// Put the Modules of In-Face
	
	count=1;
	for(Int_t i=2;i<10;i++){

		Float_t posx = constants->ModulePositionX(i-1);
		Float_t posy = constants->ModulePositionY(i-1);
		Float_t posz = constants->ModulePositionZ(i-1);	

		in_facem->AddNode(module,i-1,
			new TGeoCombiTrans("aco01",posx,posy-w1,posz,idrotm232));
		count++;

	}

	count=9;
	for(Int_t i=11;i<21;i++){
		Float_t posx = constants->ModulePositionX(i-1);
		Float_t posy = constants->ModulePositionY(i-1);
		Float_t posz = constants->ModulePositionZ(i-1);	

		in_facem->AddNode(module,i-1,
			new TGeoCombiTrans("aco01",posx,posy-w1,posz,idrotm232));
	}

	// Put he Modules of Up-Face

	count=1;
	for(Int_t i=21;i<41;i++){
		Float_t posx = constants->ModulePositionX(i-1);
		Float_t posy = constants->ModulePositionY(i-1);
		Float_t posz = constants->ModulePositionZ(i-1);	

		up_facem->AddNode(module,i-1,new TGeoTranslation("aco01",posx,posy,posz));
		count++;
	}

	// Put the Modules of Out-Face

	count=1;
	for(Int_t i=41;i<51;i++){
		Float_t posx = constants->ModulePositionX(i-1);
		Float_t posy = constants->ModulePositionY(i-1);
		Float_t posz = constants->ModulePositionZ(i-1);	

		out_facem->AddNode(module,i-1,
			new TGeoCombiTrans("aco01",posx,posy-w1,posz,idrotm231));
		count++;
	}

	// Put the Modules of Out-Face

	count=11;
	for(Int_t i=52;i<60;i++){
		Float_t posx = constants->ModulePositionX(i-1);
		Float_t posy = constants->ModulePositionY(i-1);
		Float_t posz = constants->ModulePositionZ(i-1);	
	if ((i==57) || (i==56))
		 out_facem->AddNode(module,i-1,
					new TGeoCombiTrans("aco01",posx,posy-w1,posz-w1,idrotm231));
	else
		out_facem->AddNode(module,i-1,
			new TGeoCombiTrans("aco01",posx,posy-w1,posz,idrotm231));
		count++;
	}


	// Put th Modules ITS-ACORDE

	if (GetITSGeometry()) {

		up_facem->AddNode(module,71,new TGeoTranslation("its1",
				constants->ExtraModulePositionX(),
				constants->ExtraModulePositionY(),
				constants->ExtraModulePositionZ(0)));

		up_facem->AddNode(module,72,new TGeoTranslation("its2",
				constants->ExtraModulePositionX(),
				constants->ExtraModulePositionY(),
				constants->ExtraModulePositionZ(1)));

		up_facem->AddNode(module,73,new TGeoTranslation("its3",
				constants->ExtraModulePositionX(),
				constants->ExtraModulePositionY(),
				constants->ExtraModulePositionZ(2)));

		up_facem->AddNode(module,74,new TGeoTranslation("its4",
				constants->ExtraModulePositionX(),
				constants->ExtraModulePositionY(),
				constants->ExtraModulePositionZ(3)));


		} 
	else {


		up_facem->AddNode(module,61,new TGeoTranslation("its1",
				constants->ModulePositionX(0),
				constants->ModulePositionY(0),
				constants->ModulePositionZ(0)));

		up_facem->AddNode(module,62,new TGeoTranslation("its2",
				constants->ModulePositionX(9),
				constants->ModulePositionY(9),
				constants->ModulePositionZ(9)));

		up_facem->AddNode(module,63,new TGeoTranslation("its3",
				constants->ModulePositionX(50),
				constants->ModulePositionY(50),
				constants->ModulePositionZ(50)));

		up_facem->AddNode(module,64,new TGeoTranslation("its4",
				constants->ModulePositionX(59),
				constants->ModulePositionY(59),
				constants->ModulePositionZ(59)));

		} // end if (fITSGeometry)



	//*** Begin the structure of support & bars for ACORDE ***

	// Define a volume for the bars (up-face)

	box[0]=5;
	box[1]=40;
//	box[1]=29;
	box[2]=5;
	Float_t z1 = 21 ;
	TGeoBBox *acorde00 = new TGeoBBox("acorde00",box[0],box[1],box[2]);

	TGeoVolume *ACORDE00 = new TGeoVolume("ACORDE00",acorde00,Al);

	count=25;
	for (Int_t ma=20;ma<=24;ma++)
	{
		TGeoTranslation *aco00=new TGeoTranslation("aco00",
					constants->ModulePositionX(ma)-0.5*293+dy2,
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));

		up_face->AddNode(ACORDE00,count,aco00);

		TGeoTranslation *aco00_1=new TGeoTranslation("aco00_1",
					-(constants->ModulePositionX(ma)-0.5*293+dy2),
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		up_face->AddNode(ACORDE00,count+1,aco00_1);

		TGeoTranslation *aco00_2=new TGeoTranslation("aco00_2",
					constants->ModulePositionX(ma)+0.5*293-dy2,
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		up_face->AddNode(ACORDE00,count+2,aco00_2);

		TGeoTranslation *aco00_3=new TGeoTranslation("aco00_3",
					-(constants->ModulePositionX(ma)+0.5*293-dy2),
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		up_face->AddNode(ACORDE00,count+3,aco00_3);
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
		up_face->AddNode(ACORDE00,count,aco00);

		TGeoTranslation *aco00_1=new TGeoTranslation("aco00_1",
					-(constants->ModulePositionX(ma)-0.5*293+dy2),
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		up_face->AddNode(ACORDE00,count+1,aco00_1);

		TGeoTranslation *aco00_2=new TGeoTranslation("aco00_2",
					constants->ModulePositionX(ma)+0.5*293-dy2,
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		up_face->AddNode(ACORDE00,count+2,aco00_2);

		TGeoTranslation *aco00_3=new TGeoTranslation("aco00_3",
					-(constants->ModulePositionX(ma)+0.5*293-dy2),
					constants->ModulePositionY(ma)-box[1]-z1,
					constants->ModulePositionZ(ma));
		up_face->AddNode(ACORDE00,count+3,aco00_3);
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
	up_face->AddNode(ACORDE00,57,c1);
	up_face->AddNode(ACORDE00,58,c2);
	up_face->AddNode(ACORDE00,59,c3);
	up_face->AddNode(ACORDE00,60,c4);


	// Construct Bars for lateral supports (up-face)

	TGeoTranslation *aco00=new TGeoTranslation("aco00",
				constants->ModulePositionX(20)+0.5*293-dy,
				constants->ModulePositionY(20)-box[1]-z1,
				constants->ModulePositionZ(20)-40);
	up_face->AddNode(ACORDE00,61,aco00);

	TGeoTranslation *aco00_1=new TGeoTranslation("aco00_1",
				constants->ModulePositionX(23)+0.5*293-dy,
				constants->ModulePositionY(23)-box[1]-z1,
				constants->ModulePositionZ(23)-40);
	up_face->AddNode(ACORDE00,62,aco00_1);

	TGeoTranslation *aco00_2=new TGeoTranslation("aco00_2",
				constants->ModulePositionX(24)+0.5*293-dy,
				constants->ModulePositionY(24)-box[1]-z1,
				constants->ModulePositionZ(25)-40);
	up_face->AddNode(ACORDE00,63,aco00_2);

	TGeoTranslation *aco00_3=new TGeoTranslation("aco00_3",
				constants->ModulePositionX(27)+0.5*293-dy,
				constants->ModulePositionY(27)-box[1]-z1,
				constants->ModulePositionZ(28)-40);
	up_face->AddNode(ACORDE00,64,aco00_3);


	TGeoTranslation *aco01=new TGeoTranslation("aco01",
				constants->ModulePositionX(30)-0.5*293+dy,
				constants->ModulePositionY(30)-box[1]-z1,
				constants->ModulePositionZ(30)-40);
	up_face->AddNode(ACORDE00,65,aco01);

	TGeoTranslation *aco01_1=new TGeoTranslation("aco01_1",
				constants->ModulePositionX(33)-0.5*293+dy,
				constants->ModulePositionY(33)-box[1]-z1,
				constants->ModulePositionZ(33)-40);
	up_face->AddNode(ACORDE00,66,aco01_1);

	TGeoTranslation *aco01_2=new TGeoTranslation("aco01_2",
				constants->ModulePositionX(34)-0.5*293+dy,
				constants->ModulePositionY(34)-box[1]-z1,
				constants->ModulePositionZ(35)-40);
	up_face->AddNode(ACORDE00,67,aco01_2);

	TGeoTranslation *aco01_3=new TGeoTranslation("aco01_3",
				constants->ModulePositionX(37)-0.5*293+dy,
				constants->ModulePositionY(37)-box[1]-z1,
				constants->ModulePositionZ(38)-40);
	up_face->AddNode(ACORDE00,68,aco01_3);



	// Acorde's support bars (side's faces)

	//*** In Face ***

	box[0]=39;
//	box[0]=29;
	box[1]=5;
	box[2]=5;
	Float_t q1=0;
	Float_t posx=constants->ModulePositionX(0)+0.5*293*0.7071-56*0.7071-18;
	Float_t posy=constants->ModulePositionY(0)-0.5*293*0.7071-56*0.7071+3-q1;
	Float_t posz=constants->ModulePositionZ(0);

	TGeoBBox *acorde7 = new TGeoBBox("acorde7",box[0],box[1],box[2]);

	TGeoVolume *ACORDE7 = new TGeoVolume("ACORDE7",acorde7,Al);

	TGeoCombiTrans *aco7 = new TGeoCombiTrans("aco7",posx,posy,posz-4*dy,idrotm231);
	TGeoCombiTrans *aco7_1 = new TGeoCombiTrans("aco7_1",posx,posy,
					constants->ModulePositionZ(3)-4*dy,idrotm231);
	TGeoCombiTrans *aco7_2 = new TGeoCombiTrans("aco7_2",posx,posy,
					constants->ModulePositionZ(5)-4*dy,idrotm231);
	TGeoCombiTrans *aco7_3 = new TGeoCombiTrans("aco7_3",posx,posy,
					constants->ModulePositionZ(8)-4*dy,idrotm231);

	in_face->AddNode(ACORDE7,20,aco7);
	in_face->AddNode(ACORDE7,21,aco7_1);
	in_face->AddNode(ACORDE7,22,aco7_2);
	in_face->AddNode(ACORDE7,23,aco7_3);


	count=24;
	for(Int_t dy=0;dy<=4;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-56*0.7071+3-des-q1;
		Float_t posz=constants->ModulePositionZ(dy);
		Float_t posx2=constants->ModulePositionX(dy)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->ModulePositionY(dy)+0.27*293*0.7071-56*0.7071+3-des-q1;
		TGeoCombiTrans *aco7_4 = new TGeoCombiTrans("aco7_4",posx1,posy1,posz,idrotm231);
		TGeoCombiTrans *aco7_5 = new TGeoCombiTrans("aco7_5",posx2,posy2,posz,idrotm231);
		in_face->AddNode(ACORDE7,count,aco7_4);
		in_face->AddNode(ACORDE7,count+1,aco7_5);
		count=count+2;
		dy++;
	}	


	count=34;
	for(Int_t dy=5;dy<=9;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-56*0.7071+3-des-q1;
		Float_t posz=constants->ModulePositionZ(dy+10);
		Float_t posx2=constants->ModulePositionX(dy)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->ModulePositionY(dy)+0.27*293*0.7071-56*0.7071+3-des-q1;
		TGeoCombiTrans *aco7_6 = new TGeoCombiTrans("aco7_6",posx1,posy1,posz,idrotm231);
		TGeoCombiTrans *aco7_7 = new TGeoCombiTrans("aco7_7",posx2,posy2,posz,idrotm231);
		in_face->AddNode(ACORDE7,count,aco7_6);
		in_face->AddNode(ACORDE7,count+1,aco7_7);
		count=count+2;
		dy++;
	}	



	Float_t posx_1=constants->ModulePositionX(10)+0.5*293*0.7071-56*0.7071-18;
	Float_t posy_1=constants->ModulePositionY(10)-0.5*293*0.7071-56*0.7071+3-q1;
	Float_t posz_1=constants->ModulePositionZ(10);
	TGeoCombiTrans *aco7_8 = new TGeoCombiTrans("aco7_8",posx_1,posy_1,posz_1-4*dy,idrotm231);
	TGeoCombiTrans *aco7_9 = new TGeoCombiTrans("aco7_9",posx_1,posy_1,
					constants->ModulePositionZ(13)-4*dy,idrotm231);
	TGeoCombiTrans *aco7_10 = new TGeoCombiTrans("aco7_10",posx_1,posy_1,
					constants->ModulePositionZ(15)-4*dy,idrotm231);
	TGeoCombiTrans *aco7_11 = new TGeoCombiTrans("aco7_11",posx_1,posy_1,
					constants->ModulePositionZ(18)-4*dy,idrotm231);
	in_face->AddNode(ACORDE7,44,aco7_8);
	in_face->AddNode(ACORDE7,45,aco7_9);
	in_face->AddNode(ACORDE7,46,aco7_10);
	in_face->AddNode(ACORDE7,47,aco7_11);


	count=48;
	for(Int_t dy=10;dy<=14;dy++)

	{

		Float_t posx1=constants->ModulePositionX(dy)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-56*0.7071+3-des-0.8;
		Float_t posz=constants->ModulePositionZ(dy);
		Float_t posx2=constants->ModulePositionX(dy)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->ModulePositionY(dy)+0.27*293*0.7071-56*0.7071+3-des-1.5-0.8;
		TGeoRotation *rot1 = new TGeoRotation();
		rot1->RotateZ(70);
		TGeoCombiTrans *aco7_12 = new TGeoCombiTrans("aco7_12",posx1,posy1,posz,idrotm231);
		TGeoCombiTrans *aco7_13 = new TGeoCombiTrans("aco7_13",posx2+15,posy2-10,posz,rot1);
		in_face->AddNode(ACORDE7,count,aco7_12);
		in_face->AddNode(ACORDE7,count+1,aco7_13);// bars 25 grades
		count=count+2;
		dy++;
	}


	count=57;
	for(Int_t dy=15;dy<=19;dy++)

	{

		Float_t posx1=constants->ModulePositionX(dy)+0.1*293*0.7071-56*0.7071-18+des;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-56*0.7071+3-des-q1-0.8;
		Float_t posz=constants->ModulePositionZ(dy);
		Float_t posx2=constants->ModulePositionX(dy)-0.27*293*0.7071-56*0.7071-18+des;
		Float_t posy2=constants->ModulePositionY(dy)+0.27*293*0.7071-56*0.7071+3-des-1.5-q1-0.8;
		TGeoRotation *rot1 = new TGeoRotation();
		rot1->RotateZ(70);
		TGeoCombiTrans *aco7_14 = new TGeoCombiTrans("aco7_14",posx1,posy1,posz,idrotm231);
		TGeoCombiTrans *aco7_15 = new TGeoCombiTrans("aco7_15",posx2+15,posy2-10,posz,rot1);
		in_face->AddNode(ACORDE7,count,aco7_14);
		in_face->AddNode(ACORDE7,count+1,aco7_15);// bars 25 grades
		count=count+2;
		dy++;
	}


	//*** Out Face ***

	box[0]=39;
//	box[0]=29;
	box[1]=5;
	box[2]=5;
	Float_t s1=2.5;
	Float_t posx_a=constants->ModulePositionX(50)-0.5*293*0.7071+56*0.7071+18;
	Float_t posy_a=constants->ModulePositionY(50)-0.5*293*0.7071-56*0.7071+3-s1;
	Float_t posz_a=constants->ModulePositionZ(50);
	TGeoCombiTrans *aco7_16 = new TGeoCombiTrans("aco7_16",
					posx_a,posy_a,posz_a-4*dy,idrotm232);
	TGeoCombiTrans *aco7_17 = new TGeoCombiTrans("aco7_17",
					posx_a,posy_a,
					constants->ModulePositionZ(43)-4*dy,idrotm232);
	TGeoCombiTrans *aco7_18 = new TGeoCombiTrans("aco7_18",posx_a,posy_a,
					constants->ModulePositionZ(55)-4*dy,idrotm232);
	TGeoCombiTrans *aco7_19 = new TGeoCombiTrans("aco7_19",posx_a,posy_a,
					constants->ModulePositionZ(58)-4*dy,idrotm232);
	TGeoCombiTrans *aco7_20 = new TGeoCombiTrans("aco7_20",
					constants->ModulePositionX(50)-0.1*293*0.7071
					+56*0.7071+18-des,
					constants->ModulePositionY
					(50)-0.1*293*0.7071-56*0.7071+3-des-s1,
					constants->ModulePositionZ(45),idrotm232);
	TGeoCombiTrans *aco7_21 = new TGeoCombiTrans("aco7_21",
					constants->ModulePositionX(50)+0.27*293*0.7071
					+56*0.7071+18-des,
					constants->ModulePositionY(50)
					+0.27*293*0.7071-56*0.7071+3-des-s1,
					constants->ModulePositionZ(45),idrotm232);
	out_face->AddNode(ACORDE7,19,aco7_16);
	out_face->AddNode(ACORDE7,20,aco7_17);
	out_face->AddNode(ACORDE7,21,aco7_18);
	out_face->AddNode(ACORDE7,22,aco7_19);
	out_face->AddNode(ACORDE7,23,aco7_20);
	out_face->AddNode(ACORDE7,24,aco7_21);


	count=25;
	for(Int_t dy=50;dy<=54;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-56*0.7071+3-des-s1;
		Float_t posz=constants->ModulePositionZ(dy);
		Float_t posx2=constants->ModulePositionX(dy)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->ModulePositionY(dy)+0.27*293*0.7071-56*0.7071+3-des-s1;
		TGeoCombiTrans *aco7_22 = new TGeoCombiTrans("aco7_22",posx1,posy1,posz,idrotm232);
		TGeoCombiTrans *aco7_23 = new TGeoCombiTrans("aco7_23",posx2,posy2,posz,idrotm232);
		out_face->AddNode(ACORDE7,count,aco7_22);
		out_face->AddNode(ACORDE7,count+1,aco7_23);
		count=count+2;
		dy++;
	}


	count=35;
	for(Int_t dy=57;dy<=59;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-56*0.7071+3-des-s1;
		Float_t posz=constants->ModulePositionZ(dy-10);
		Float_t posx2=constants->ModulePositionX(dy)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->ModulePositionY(dy)+0.27*293*0.7071-56*0.7071+3-des-s1;
		TGeoCombiTrans *aco7_24 = new TGeoCombiTrans("aco7_24",posx1,posy1,posz,idrotm232);
		TGeoCombiTrans *aco7_25 = new TGeoCombiTrans("aco7_25",posx2,posy2,posz,idrotm232);
		out_face->AddNode(ACORDE7,count,aco7_24);
		out_face->AddNode(ACORDE7,count+1,aco7_25);
		count=count+2;
		dy++;
	}


	Float_t posx_b=constants->ModulePositionX(40)-0.5*293*0.7071+56*0.7071+18;
	Float_t posy_b=constants->ModulePositionY(40)-0.5*293*0.7071-56*0.7071+3-s1;
	Float_t posz_b=constants->ModulePositionZ(40);
	TGeoCombiTrans *aco7_26 = new TGeoCombiTrans("aco7_26",
					posx_b,posy_b,posz_b-4*dy,idrotm232);
	TGeoCombiTrans *aco7_27 = new TGeoCombiTrans("aco7_27",
					posx_b,posy_b,
					constants->ModulePositionZ(43)-4*dy,idrotm232);
	TGeoCombiTrans *aco7_28 = new TGeoCombiTrans("aco7_28",
					posx_b,posy_b,
					constants->ModulePositionZ(45)-4*dy,idrotm232);
	TGeoCombiTrans *aco7_29 = new TGeoCombiTrans("aco7_29",posx_b,posy_b,
					constants->ModulePositionZ(48)-4*dy,idrotm232);
	out_face->AddNode(ACORDE7,41,aco7_26);
	out_face->AddNode(ACORDE7,42,aco7_27);
	out_face->AddNode(ACORDE7,43,aco7_28);
	out_face->AddNode(ACORDE7,44,aco7_29);

	count=45;
	for(Int_t dy=40;dy<=44;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-56*0.7071+3-des-s1;
		Float_t posz=constants->ModulePositionZ(dy);
		Float_t posx2=constants->ModulePositionX(dy)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->ModulePositionY(dy)+0.27*293*0.7071-56*0.7071+3-des-s1;
		TGeoRotation *rot1 = new TGeoRotation();
		rot1->RotateZ(105);
		TGeoCombiTrans *aco7_30 = new TGeoCombiTrans("aco7_30",posx1,posy1,posz,idrotm232);
		TGeoCombiTrans *aco7_31 = new TGeoCombiTrans("aco7_31",posx2-15,posy2-10,posz,rot1);
		out_face->AddNode(ACORDE7,count,aco7_30);
		out_face->AddNode(ACORDE7,count+1,aco7_31);// bars 25 grades
		count=count+2;
		dy++;
	}


	count=55;
	for(Int_t dy=45;dy<=49;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)-0.1*293*0.7071+56*0.7071+18-des;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-56*0.7071+3-des-s1;
		Float_t posz=constants->ModulePositionZ(dy);
		Float_t posx2=constants->ModulePositionX(dy)+0.27*293*0.7071+56*0.7071+18-des;
		Float_t posy2=constants->ModulePositionY(dy)+0.27*293*0.7071-56*0.7071+3-des-s1;
		TGeoRotation *rot1 = new TGeoRotation();
		rot1->RotateZ(105);
		TGeoCombiTrans *aco7_32 = new TGeoCombiTrans("aco7_32",posx1,posy1,posz,idrotm232);
		TGeoCombiTrans *aco7_33 = new TGeoCombiTrans("aco7_33",posx2-15,posy2-10,posz,rot1);
		out_face->AddNode(ACORDE7,count,aco7_32);
		out_face->AddNode(ACORDE7,count+1,aco7_33);// bars 25 grades
		count=count+2;
		dy++;
	}



	// Set the bars non perpendicular at side faces

	//*** In-Face ***

	box[0]=5;
	box[1]=55.15;
//	box[1]=41.15;
	box[2]=5;
	Float_t sm=2;
	Float_t re=1;
	Float_t posx1=constants->ModulePositionX(0)+0.5*293*0.7071-4*box[0]-8+re;
	Float_t posy1=constants->ModulePositionY(0)-0.5*293*0.7071-box[1]-18-2+sm;
	Float_t posz1=constants->ModulePositionZ(0);

	TGeoBBox *acorde7_1 = new TGeoBBox("acorde7_1",box[0],box[1],box[2]);

	TGeoVolume *ACORDE7_1 = new TGeoVolume("ACORDE7_1",acorde7_1,Al);
	TGeoTranslation *aco71 = new TGeoTranslation("aco71",posx1,posy1,posz1-4*dy);
	TGeoTranslation *aco72 = new TGeoTranslation("aco72",posx1,posy1,
					constants->ModulePositionZ(3)-4*dy);
	TGeoTranslation *aco73 = new TGeoTranslation("aco73",posx1,posy1,
					constants->ModulePositionZ(5)-4*dy);
	TGeoTranslation *aco74 = new TGeoTranslation("aco74",posx1,posy1,
					constants->ModulePositionZ(8)-4*dy);
	in_face->AddNode(ACORDE7_1,67,aco71);
	in_face->AddNode(ACORDE7_1,68,aco72);
	in_face->AddNode(ACORDE7_1,69,aco73);
	in_face->AddNode(ACORDE7_1,70,aco74);


	count=71;
	for(Int_t dy=0;dy<=4;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)+0.1*293*0.7071-4*box[0]-8+des+re;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1=constants->ModulePositionZ(dy);
		Float_t dyx2=constants->ModulePositionX(dy)-0.27*293*0.7071-4*box[0]-8+des+re;
		Float_t dyy2=constants->ModulePositionY(dy)+0.27*293*0.7071-box[1]-18-2-des+sm;
		TGeoTranslation *aco75=new TGeoTranslation("aco75",posx1,posy1,posz1);
		TGeoTranslation *aco76=new TGeoTranslation("aco76",dyx2,dyy2,posz1);
		in_face->AddNode(ACORDE7_1,count,aco75);
		in_face->AddNode(ACORDE7_1,count+1,aco76);
		count=count+2;
		dy++;
	}


	count=81;
	for(Int_t dy=5;dy<=9;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)+0.1*293*0.7071-4*box[0]-8+des+re;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1=constants->ModulePositionZ(dy+10);
		Float_t dyx2=constants->ModulePositionX(dy)-0.27*293*0.7071-4*box[0]-8+des+re;
		Float_t dyy2=constants->ModulePositionY(dy)+0.27*293*0.7071-box[1]-18-2-des+sm;
		TGeoTranslation *aco75=new TGeoTranslation("aco75",posx1,posy1,posz1);
		TGeoTranslation *aco76=new TGeoTranslation("aco76",dyx2,dyy2,posz1);
		in_face->AddNode(ACORDE7_1,count,aco75);
		in_face->AddNode(ACORDE7_1,count+1,aco76);
		count=count+2;
		dy++;
	}


	Float_t posx1_1=constants->ModulePositionX(10)+0.5*293*0.7071-4*box[0]-8+re;
	Float_t posy1_1=constants->ModulePositionY(10)-0.5*293*0.7071-box[1]-18-2+sm;
	Float_t posz1_1=constants->ModulePositionZ(10);
	TGeoTranslation *aco77=new TGeoTranslation("aco77",posx1_1,posy1_1,posz1_1-4*dy);
	TGeoTranslation *aco78=new TGeoTranslation("aco78",posx1_1,posy1_1,
					constants->ModulePositionZ(13)-4*dy);

	TGeoTranslation *aco79=new TGeoTranslation("aco79",posx1_1,posy1_1,
					constants->ModulePositionZ(15)-4*dy);
	TGeoTranslation *aco710=new TGeoTranslation("aco710",posx1_1,posy1_1,
					constants->ModulePositionZ(18)-4*dy);
	in_face->AddNode(ACORDE7_1,91,aco77);
	in_face->AddNode(ACORDE7_1,92,aco78);
	in_face->AddNode(ACORDE7_1,93,aco79);
	in_face->AddNode(ACORDE7_1,94,aco710);

	count=95;
	for(Int_t dy=10;dy<=14;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)+0.1*293*0.7071-4*box[0]-8+des+re+.83;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1=constants->ModulePositionZ(dy);
		Float_t dyx2=constants->ModulePositionX(dy)-0.27*293*0.7071-4*box[0]-4+des+re+0.83;
		Float_t dyy2=constants->ModulePositionY(dy)+0.27*293*0.7071-box[1]-18-5-des+sm;
		TGeoTranslation *aco711=new TGeoTranslation("aco711",posx1,posy1,posz1);
		TGeoTranslation *aco712=new TGeoTranslation("aco712",dyx2,dyy2,posz1);
		in_face->AddNode(ACORDE7_1,count,aco711);
		in_face->AddNode(ACORDE7_1,count+1,aco712);
		count=count+2;
		dy++;
	}



	count=105;
	for(Int_t dy=15;dy<=19;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)+0.1*293*0.7071-4*box[0]-8+des+re+0.83;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-box[1]-18-2-des+sm;
		Float_t posz1=constants->ModulePositionZ(dy);
		Float_t dyx2=constants->ModulePositionX(dy)-0.27*293*0.7071-4*box[0]-4+des+re+0.83;
		Float_t dyy2=constants->ModulePositionY(dy)+0.27*293*0.7071-box[1]-18-5-des;
		TGeoTranslation *aco713=new TGeoTranslation("aco713",posx1,posy1,posz1);
		TGeoTranslation *aco714=new TGeoTranslation("aco714",dyx2,dyy2,posz1);
		in_face->AddNode(ACORDE7_1,count,aco713);
		in_face->AddNode(ACORDE7_1,count+1,aco714);
		count=count+2;
		dy++;
	}

		//*** Out-Face ***

	Float_t posx1_a=constants->ModulePositionX(50)-0.5*293*0.7071+4*box[0]+8-re-1;
	Float_t posy1_a=constants->ModulePositionY(50)-0.5*293*0.7071-box[1]-18-2+sm-2.5;
	Float_t posz1_a=constants->ModulePositionZ(50);
	TGeoTranslation *aco715=new TGeoTranslation("aco715",posx1_a,posy1_a,posz1_a-4*dy);
	TGeoTranslation *aco716=new TGeoTranslation("aco716",posx1_a,posy1_a,
				constants->ModulePositionZ(43)-4*dy);
	TGeoTranslation *aco717=new TGeoTranslation("aco717",posx1_a,posy1_a,
				constants->ModulePositionZ(55)-4*dy);
	TGeoTranslation *aco718=new TGeoTranslation("aco718",posx1_a,posy1_a,
				constants->ModulePositionZ(58)-4*dy);
	TGeoTranslation *aco719=new TGeoTranslation("aco719",
				constants->ModulePositionX(50)-0.1*293*0.7071+4*box[0]+8-des-re-1,		
				constants->ModulePositionY(50)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5,
				constants->ModulePositionZ(45));
	TGeoTranslation *aco720=new TGeoTranslation("aco720",
				constants->ModulePositionX(50)+0.27*293*0.7071+4*box[0]+8-des-re-1,
				constants->ModulePositionY(50)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5,
				constants->ModulePositionZ(45));


	out_face->AddNode(ACORDE7_1,115,aco715);
	out_face->AddNode(ACORDE7_1,116,aco716);
	out_face->AddNode(ACORDE7_1,117,aco717);
	out_face->AddNode(ACORDE7_1,118,aco718);
	out_face->AddNode(ACORDE7_1,119,aco719);
	out_face->AddNode(ACORDE7_1,120,aco720);




	count=65;
	for(Int_t dy=50;dy<=54;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1=constants->ModulePositionZ(dy);
		Float_t dyx2=constants->ModulePositionX(dy)+0.27*293*0.7071+4*box[0]+8-des-re-1;
		Float_t dyy2=constants->ModulePositionY(dy)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5;
		TGeoTranslation *aco721=new TGeoTranslation("aco721",posx1,posy1,posz1);
		TGeoTranslation *aco722=new TGeoTranslation("aco722",dyx2,dyy2,posz1);
		out_face->AddNode(ACORDE7_1,count,aco721);
		out_face->AddNode(ACORDE7_1,count+1,aco722);
		count=count+2;
		dy++;
	}



	count=75;
	for(Int_t dy=57;dy<=59;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1=constants->ModulePositionZ(dy-10);
		Float_t dyx2=constants->ModulePositionX(dy)+0.27*293*0.7071+4*box[0]+8-des-re-1;
		Float_t dyy2=constants->ModulePositionY(dy)+0.27*293*0.7071-box[1]-18-2-des+sm-2.5;
		TGeoTranslation *aco723=new TGeoTranslation("aco723",posx1,posy1,posz1);
		TGeoTranslation *aco724=new TGeoTranslation("aco724",dyx2,dyy2,posz1);
		out_face->AddNode(ACORDE7_1,count,aco723);
		out_face->AddNode(ACORDE7_1,count+1,aco724);
		count=count+2;
		dy++;
	}


	Float_t posx1_b=constants->ModulePositionX(40)-0.5*293*0.7071+4*box[0]+5;
	Float_t posy1_b=constants->ModulePositionY(40)-0.5*293*0.7071-box[1]-18-2;
	Float_t posz1_b=constants->ModulePositionZ(40);
	TGeoTranslation *aco725=new TGeoTranslation("aco725",posx1_b,posy1_b,posz1_b-4*dy);
	TGeoTranslation *aco726=new TGeoTranslation("aco726",posx1_b,posy1_b,
				constants->ModulePositionZ(43)-4*dy);
	TGeoTranslation *aco727=new TGeoTranslation("aco727",posx1_b,posy1_b,
				constants->ModulePositionZ(45)-4*dy);
	TGeoTranslation *aco728=new TGeoTranslation("aco728",posx1_b,posy1_b,
				constants->ModulePositionZ(48)-4*dy);
	out_face->AddNode(ACORDE7_1,85,aco725);
	out_face->AddNode(ACORDE7_1,86,aco726);
	out_face->AddNode(ACORDE7_1,87,aco727);
	out_face->AddNode(ACORDE7_1,88,aco728);



	count=89;
	for(Int_t dy=40;dy<=44;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1=constants->ModulePositionZ(dy);
		Float_t dyx2=constants->ModulePositionX(dy)+0.27*293*0.7071+4*box[0]+4-des-re-1+2.8;
		Float_t dyy2=constants->ModulePositionY(dy)+0.27*293*0.7071-box[1]-18-5-des+sm-2.5+3;
		TGeoTranslation *aco729=new TGeoTranslation("aco729",posx1,posy1,posz1);
		TGeoTranslation *aco730=new TGeoTranslation("aco730",dyx2,dyy2,posz1);
		out_face->AddNode(ACORDE7_1,count,aco729);
		out_face->AddNode(ACORDE7_1,count+1,aco730);
		count=count+2;
		dy++;
	}



	count=99;
	for(Int_t dy=45;dy<=49;dy++)
	{

		Float_t posx1=constants->ModulePositionX(dy)-0.1*293*0.7071+4*box[0]+8-des-re-1;
		Float_t posy1=constants->ModulePositionY(dy)-0.1*293*0.7071-box[1]-18-2-des+sm-2.5;
		Float_t posz1=constants->ModulePositionZ(dy);
		Float_t dyx2=constants->ModulePositionX(dy)+0.27*293*0.7071+4*box[0]+4-des-re-1+2.8;
		Float_t dyy2=constants->ModulePositionY(dy)+0.27*293*0.7071-box[1]-18-5-des+sm-2.5+3;
		TGeoTranslation *aco729=new TGeoTranslation("aco729",posx1,posy1,posz1);
		TGeoTranslation *aco730=new TGeoTranslation("aco730",dyx2,dyy2,posz1);
		out_face->AddNode(ACORDE7_1,count,aco729);
		out_face->AddNode(ACORDE7_1,count+1,aco730);
		count=count+2;
		dy++;
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
	TGeoVolume *ACORDE8 = new TGeoVolume("ACORDE8",acorde8,Al);

	TGeoBBox *acorde9 = new TGeoBBox("acorde9",tbox[0],tbox[1],tbox[2]);
	TGeoVolume *ACORDE9 = new TGeoVolume("ACORDE9",acorde9,Al);

	support->AddNode(ACORDE8,1,new TGeoTranslation(0,-5,0));
	support->AddNode(ACORDE8,2,new TGeoTranslation(0,-dy1,0));
	support->AddNode(ACORDE9,3,new TGeoTranslation(0,-tbox[1]-5.5,0));


	// Put "support" on Up-Face

	up_face->AddNode(support,69,new TGeoTranslation("aco8",sx,sy,sz));
	up_face->AddNode(support,70,new TGeoTranslation("aco8_2",sx2,sy2,sz));
	up_face->AddNode(support,71,new TGeoTranslation("aco8_4",sx4,sy4,sz));
	up_face->AddNode(support,72,new TGeoTranslation("aco8_6",sx5,sy5,sz));
	up_face->AddNode(support,73,new TGeoTranslation("aco8_2",-sx2,sy2,sz));
	up_face->AddNode(support,74,new TGeoTranslation("aco8_4",-sx4,sy4,sz));
	up_face->AddNode(support,75,new TGeoTranslation("aco8_6",-sx5,sy5,sz));

	// Put "support" on In-Face
	Float_t ms = 1.3;
	in_face->AddNode(support,121,new TGeoCombiTrans("aco8_81",dyx,dyy+ms,dyz,idrotm232));
	in_face->AddNode(support,122,new TGeoCombiTrans("aco8_121",dyx1+des,ms+dyy1-des,dyz,idrotm232));
	in_face->AddNode(support,123,new TGeoCombiTrans("aco8_161",dyx2+des,ms+dyy2-des,dyz,idrotm232));
	in_face->AddNode(support,124,new TGeoCombiTrans("aco8_82",dx1,ms+dy11,dyz,idrotm232));
	in_face->AddNode(support,125,new TGeoCombiTrans("aco8_122",dyx11+des,ms+dyy11-des,dyz,idrotm232));
	in_face->AddNode(support,126,new TGeoCombiTrans("aco8_162",dyx21+des,ms+dyy21-des,dyz,idrotm232));

	// Put "support" on Out-Face

	out_face->AddNode(support,121,new TGeoCombiTrans("aco8_81",-dyx,dyy+ms,dyz,idrotm231));
	out_face->AddNode(support,122,new TGeoCombiTrans("aco8_121",-dyx1-des,ms+dyy1-des,dyz,idrotm231));
	out_face->AddNode(support,123,new TGeoCombiTrans("aco8_161",-dyx2-des,ms+dyy2-des,dyz,idrotm231));
	out_face->AddNode(support,124,new TGeoCombiTrans("aco8_82",-dx1,dy11+ms,dyz,idrotm231));
	out_face->AddNode(support,125,new TGeoCombiTrans("aco8_122",-dyx11-des,ms+dyy11-des,dyz,idrotm231));
	out_face->AddNode(support,126,new TGeoCombiTrans("aco8_162",-dyx21-des,ms+dyy21-des,dyz,idrotm231));
	
	ACORDE->AddNode(in_face,1);//---> volume of supports & bars in-face
	ACORDE->AddNode(up_face,2);//---> volume of supports & bars up-face
	ACORDE->AddNode(out_face,3);//---> volume of supports & bars out-face
	ACORDE->AddNode(in_facem,4);//---> volume of modules in-face
	ACORDE->AddNode(up_facem,5);//---> volume of modules up-face
	ACORDE->AddNode(out_facem,6);//---> volume of modules out-face
	alice->AddNode(ACORDE,1);//---> put volume of ACORDE over ALICE's volume



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

    printf("Hay particula cargada en el volumen %d \n",idScint);
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
      hits[1] = pos[0]; printf("PosX %f \n",hits[1]);
      hits[2] = pos[1]; printf("PosY %f \n",hits[2]);
      hits[3] = pos[2]; printf("PosZ %f \n",hits[3]);
      hits[4] = gMC->TrackTime();printf("TimeTracking %f \n",hits[4]);
      hits[5] = mom[0]; printf("MomentoX %f \n",hits[5]);
      hits[6] = mom[1]; printf("MomentoY %f \n",hits[6]);
      hits[7] = mom[2]; printf("MomentoZ %f \n",hits[7]);
      hits[8] = gMC->Etot();printf("EnergiaTotal %f \n",hits[8]);
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
      vol[1] = copyPlastic;
    } // end if gMC->IsTrackEntering()

    // set hit[9] = total energy loss and book hit
    if( gMC->IsTrackExiting() || 
	gMC->IsTrackStop() || 
	gMC->IsTrackDisappeared()){
      hits[9] = eloss;printf("Energia Perdida %f \n",hits[9]);
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
  // Add a ACORDE hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliACORDEhit(fIshunt,track,vol,hits);
}

