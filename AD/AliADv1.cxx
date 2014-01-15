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

/* $Id: AliAD.cxx  $ */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                  AD (ALICE Diffractive)  Detector                     //
//                                                                       //
//  This class contains the base procedures for the AD  detector         //
//  Default geometry of 2013: 16 modules                                 //
//  All comments should be sent to :                                     //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// --- Standard libraries ---
#include <Riostream.h>

// --- ROOT libraries ---
#include <TMath.h>
#include <TString.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoShape.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGeoCompositeShape.h>
#include <TGeoGlobalMagField.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoArb8.h>
#include <TClonesArray.h>
#include <TGeoTrd2.h>
#include <TParticle.h>

#include <TH2F.h>
#include <TCanvas.h>

// --- AliRoot header files ---


#include "AliADhit.h"
#include "AliADdigit.h"
#include "AliADv1.h"
#include "AliLog.h"
#include "AliConst.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"


ClassImp(AliADv1)
//__________________________________________________________________
AliADv1::AliADv1()
  : AliAD(),
  fADCLightYield(93.75),
  fADCPhotoCathodeEfficiency(0.18),
  fADALightYield(93.75),
  fADAPhotoCathodeEfficiency(0.18)

{
   // Default Constructor
    fHits = 0;
}

//_____________________________________________________________________________
AliADv1::AliADv1(const char *name, const char *title) : 
  AliAD(name,title),  
  fADCLightYield(93.75),
  fADCPhotoCathodeEfficiency(0.18),
  fADALightYield(93.75),
  fADAPhotoCathodeEfficiency(0.18)
{
   // Standard constructor for AD Detector
  
   AliModule* pipe = gAlice->GetModule("PIPE");
   if( (!pipe) ) {
      Error("Constructor","AD needs PIPE!!!\n");
      exit(1);
   } 
   fHits = new TClonesArray("AliADhit",400);
   gAlice->GetMCApp()->AddHitList(fHits);
}

//_____________________________________________________________________________
AliADv1::~AliADv1()
{
	// default destructor
}
//_____________________________________________________________________________
void AliADv1::Init()
{
  // Initialise L3 magnet after it has been built
  Int_t i;
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" ADv1_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

//_____________________________________________________________________________
void AliADv1::CreateGeometry()
{
  //
  // Create the geometry for the AD arrays
  //
  
  CreateAD();
  
}


//_____________________________________________________________________________
void AliADv1::CreateAD()
{

	// here we create AD: ADA & ADC

	// Get ALICE volume

	TGeoVolume *alice = gGeoManager->GetVolume("ALIC");

	// Define the mother volume for AD

	TGeoVolume *ad = new TGeoVolumeAssembly("AD");

	// Get medium

	TGeoMedium *medADASci		= gGeoManager->GetMedium("AD_NE102"); // AD Scin. 
	TGeoMedium *medADALG		= gGeoManager->GetMedium("AD_PMMA");  // lightGuide
	TGeoMedium *medADAPMGlass	= gGeoManager->GetMedium("AD_Glass"); // Glass for Aluminium simulation
	TGeoMedium *medADAPMAlum	= gGeoManager->GetMedium("AD_Alum");  // Aluminium 


   	///  ADA Scintillator Pad Measures
   	const Double_t kADATriangleSide   = 7.8;   //
   	const Double_t kADACellSide       = 20.0;
   	const Double_t kADACellThickness  = 2.0;  // Half thickness
   	const int kColorADA = kGreen;
 

	// Creation of the Box's pad

   	new TGeoBBox( "ADAbox", kADACellSide/2.0-kADATriangleSide/2., kADACellSide/2.0, kADACellThickness );
   	const Double_t boxSide2 = kADACellSide/2.0-kADATriangleSide/2.;
   	new TGeoBBox( "ADAbox1", kADATriangleSide/2., boxSide2, kADACellThickness );
	
	// translation

   	TGeoTranslation *trada2 = new TGeoTranslation( -kADACellSide/2.0,  kADACellSide/2.0 - boxSide2, 0. );
   	trada2->SetName( "trada2" );
   	trada2->RegisterYourself();

   	TGeoArb8* sADAtriang = new TGeoArb8( "ADAtriang", kADACellThickness );
   	for ( int iz = 0; iz < 2; iz++ ) {
      		sADAtriang->SetVertex( 0+iz*4, kADACellSide/2.0, kADACellSide/2.0 );
      		sADAtriang->SetVertex( 1+iz*4, kADACellSide/2.0, (kADACellSide/2.0)-kADATriangleSide );
      		sADAtriang->SetVertex( 2+iz*4, kADACellSide/2.0, (kADACellSide/2.0)-kADATriangleSide );
      		sADAtriang->SetVertex( 3+iz*4, kADACellSide/2.0-kADATriangleSide, kADACellSide/2.0 );
   	}    
   	TGeoTranslation *trada1 = new TGeoTranslation( -kADACellSide+kADATriangleSide/2. , -kADACellSide+kADATriangleSide, 0. );
   	trada1->SetName( "trada1" );
   	trada1->RegisterYourself();

   	TGeoCompositeShape *sADA1 = new TGeoCompositeShape ( "sADA1s1", "ADAbox+(ADAbox1:trada2)+(ADAtriang:trada1)" );
   	TGeoVolume *vADA1 = new TGeoVolume( "ADApad", sADA1, medADASci );
   	vADA1->SetLineColor( kColorADA ); 

   	/// Light guide
   	Double_t kADALGThickness   = 2.0; // Half thickness
   	Double_t kADALGSideScint   = 20.0;
   	Double_t kADALGHeigth      = TMath::Sqrt( kADALGSideScint*kADALGSideScint - (10.5 * 10.5) );
   	Double_t kADALGSideCoupling = 3.0; 
   	const int kColorADALG = kYellow;
 
 	  // Triangle
   	TGeoTrd2* sADALGtriang = new TGeoTrd2( kADALGThickness, kADALGThickness, kADALGSideScint/2., kADALGSideCoupling/2., kADALGHeigth/2.0);
   	TGeoVolume *vADALGtriang = new TGeoVolume( "ADALG", sADALGtriang, medADALG );
   	vADALGtriang->SetLineColor( kColorADALG ); 

   	// Coupling
   	Double_t kADALGCoupling = 5.0; 
   	TGeoTube* sADACouplTube = new TGeoTube( "ADACouplTube", 0, 1.4, kADALGCoupling/2. );
   	TGeoVolume * vADACoupling = new TGeoVolume( "ADACoupling", sADACouplTube, medADALG );
   	vADACoupling->SetLineColor( kColorADALG ); 
   
   	TGeoVolume * vADALG1  = new TGeoVolumeAssembly( "ADALG" );
   	vADALG1->AddNode( vADALGtriang, 1 );   
   	vADALG1->AddNode( vADACoupling, 1, new TGeoTranslation(0., 0., kADALGHeigth/2.+kADALGCoupling/2.) );   
   	vADALG1->SetLineColor( kColorADALG ); 
      
   	/// PMT  Hamamatsu R5946
   	Double_t kADAPMR1 = 1.95;          // 3.9 cm diameter
   	Double_t kADAPMR2 = 2.15;          // + 2 mm?? aluminium case 
   	Double_t kADAPMlength = 6.4;      // 5 cm PMT + 1.4 socket 
   	const int kColorPMG   = kWhite;
   	const int kColorPMA   = kGray;
 
   	TGeoTube *sADAPMg1   = new TGeoTube( "sADAPMg", 0., kADAPMR1, kADAPMlength/2. );
   	TGeoVolume *vADAPMg1 = new TGeoVolume( "ADAPMg", sADAPMg1, medADAPMGlass );
   	vADAPMg1->SetLineColor( kColorPMG );
   	TGeoTube *sADAPMa1   = new TGeoTube( "ADAPMa", kADAPMR1, kADAPMR2, kADAPMlength/2. );
   	TGeoVolume *vADAPMa1 = new TGeoVolume( "ADAPMa", sADAPMa1, medADAPMAlum );
   	vADAPMa1->SetLineColor( kColorPMA );
   	TGeoVolume *vADAPM1  = new TGeoVolumeAssembly("ADAPM");
   	vADAPM1->AddNode( vADAPMg1, 1 );
   	vADAPM1->AddNode( vADAPMa1, 1 );
 
    	/// Sector (Assembly:  Scintillator Pad + Light guide + PM )
   	TGeoVolume *secADA  = new TGeoVolumeAssembly( "ADAsec" ); 
   	// Add PAD
   	secADA->AddNode( vADA1, 1, new TGeoTranslation( kADACellSide/2.0+kADATriangleSide/2. + 0.05, kADACellSide/2.0 + 0.05, 0. ) );
   	// Add Light Guide
   	TGeoCombiTrans *transrot = new TGeoCombiTrans( kADACellSide + kADALGHeigth/2.0 + 0.05, kADALGSideScint/2.0 + 0.05, 0., 
                                                  new TGeoRotation("rot",90.,90.,90.) );
   	secADA->AddNode( vADALG1, 1, transrot );   
   	// Add PM
   	transrot = new TGeoCombiTrans( kADACellSide + kADALGHeigth  + kADALGCoupling + kADAPMlength/2. + 0.05, kADACellSide/2.0 + 0.05, 0, 
                                  new TGeoRotation("rot",90.,90.,0.) );
   	secADA->AddNode(vADAPM1, 1, transrot);

   	// TODO: Add mechanical support
   
   	/// Assembling ADA adding 4 sectors                                       //  Sectors
   	TGeoVolume *vADAarray = new TGeoVolumeAssembly( "ADA" );                  //        ^ y
   	vADAarray->AddNode( secADA, 1 );                                          //        |   
   	vADAarray->AddNode( secADA, 2, new TGeoRotation("rot",0.  , 180.,0.) );   //   4    |   1
   	vADAarray->AddNode( secADA, 3, new TGeoRotation("rot",180., 0.,  0.) );   // --------------->  x     
   	vADAarray->AddNode( secADA, 4, new TGeoRotation("rot",180., 180.,0.) );   //   3    |   2
        TGeoRotation *rotADA90 = new TGeoRotation("adarot",90,0,0);
	// here I add ADA to AD volume
	const Float_t kPosADA = 1700.0;
	ad->AddNode(vADAarray,1, new TGeoCombiTrans("ada",0,0,kPosADA,rotADA90));                                                                     //        |
 
	if (GetADAToInstalled())
	{
		const Float_t kPosADA2 = 1695.0;
		ad->AddNode(vADAarray,2, new TGeoCombiTrans("ada",0,0,kPosADA2,rotADA90));
	}

	// Creation of ADC

	// Get Medium for ADC (in principle is the same as ADA, but I keep the previous variables)

   	TGeoMedium *medADCSci     = gGeoManager->GetMedium("AD_NE102");
   	TGeoMedium *medADCLG      = gGeoManager->GetMedium("AD_PMMA");
   	TGeoMedium *medADCPMGlass = gGeoManager->GetMedium("AD_Glass");
   	TGeoMedium *medADCPMAlum  = gGeoManager->GetMedium("AD_Alum");
 
   	/// Creation of assembly of one ADC sector
 
   	/// ADC Scintillator Pad 
   	const Double_t kADCCellSide = 30.;
   	const Double_t kADCCellThickness = 4.0;
   	const int kColorADC = kGreen;
   
   	new TGeoBBox( "ADCbox0", kADCCellSide/4.0, kADCCellSide/2.0, kADCCellThickness/2.0 );
   	new TGeoBBox( "ADCbox1", kADCCellSide/4.0, kADCCellSide/4.0, kADCCellThickness/2.0 );
   	new TGeoBBox( "ADCbox2", 2.5, 5.5, kADCCellThickness/2.0 );
   	TGeoTranslation *tradd1 = new TGeoTranslation( -kADCCellSide/2.0, kADCCellSide/4.0, 0. );
   	TGeoTranslation *tradd2 = new TGeoTranslation( -kADCCellSide/4.0 - 2.5, -kADCCellSide/2.0 + 5.5 , 0. );
   	tradd1->SetName( "tradd1" );
   	tradd2->SetName( "tradd2" );
   	tradd1->RegisterYourself();
   	tradd2->RegisterYourself();
   	TGeoCompositeShape *sADC1 = new TGeoCompositeShape ( "sADCpad", "ADCbox0+(ADCbox1:tradd1)+(ADCbox2:tradd2)" );
   	TGeoVolume *vADC = new TGeoVolume( "ADCpad", sADC1, medADCSci );      
   	vADC->SetLineColor( kColorADC );

   	/// Light guide
   	const Double_t kADCLGThickness    = 4.0;
   	const Double_t kADCLGHeigth       = 28.95;        // Dist from scint to coupling
   	const Double_t kADCLGSideScint    = kADCCellSide; // 30.0 
   	const Double_t kADCLGSideCoupling = 4.0; 
   	const int kColorADCLG = kYellow;
  
   	// Triangle
   	TGeoTrd2* sADCLGtriang = new TGeoTrd2( kADCLGThickness/2., kADCLGThickness/2., kADCLGSideScint/2., kADCLGSideCoupling/2., kADCLGHeigth/2.0);
   	TGeoVolume *vADCLGtriang = new TGeoVolume( "ADCLG", sADCLGtriang, medADCLG );
   	vADCLGtriang->SetLineColor( kColorADCLG ); 

   	// Coupling
   	Double_t kADCLGCoupling = 10.0; // Total lenght
   	new TGeoCone( "ADCCouplCone", kADCLGCoupling/4., 0, kADCLGSideCoupling/TMath::Sqrt(2.), 0, 1.4 );
   	new TGeoBBox( "ADCCouplBox", kADCLGSideCoupling/2., kADCLGSideCoupling/2., kADCLGCoupling/4. );
   	new TGeoTube( "ADCCouplTube", 0, 1.4, kADCLGCoupling/4. );
   	TGeoTranslation *tradd3 = new TGeoTranslation(0, 0, kADCLGCoupling/2. );
   	tradd3->SetName( "tradd3" );
   	tradd3->RegisterYourself();
   
   	TGeoCompositeShape * sADCCoupling = new TGeoCompositeShape ( "sADCCoupling", "ADCCouplBox * ADCCouplCone + (ADCCouplTube:tradd3)" );
   	TGeoVolume * vADCCoupling = new TGeoVolume( "ADCCoupling", sADCCoupling, medADCLG );
   	vADCCoupling->SetLineColor( kColorADCLG ); 
   
   	TGeoVolume * vADCLG  = new TGeoVolumeAssembly( "ADCLG" );
   	vADCLG->AddNode( vADCLGtriang, 1 );   
   	vADCLG->AddNode( vADCCoupling, 1, new TGeoTranslation(0., 0., kADCLGHeigth/2.+kADCLGCoupling/4.) );   
   	vADCLG->SetLineColor( kColorADCLG ); 
   
   	/// PM  Hamamatsu R5946  
   	const Double_t kADCPMR1 = 1.95;          // 3.9 cm diameter
   	const Double_t kADCPMR2 = 2.15;          // + 2 mm?? aluminium case 
   	const Double_t kADCPMlength = 6.4;       // 5 cm PMT + 1.4 socket 
   	//const int kColorPMG   = kWhite;
   	////const int kColorPMA   = kGray;
 
   	TGeoTube *sADCPMg   = new TGeoTube( "sADCPMg", 0., kADCPMR1, kADCPMlength/2. );
   	TGeoVolume *vADCPMg = new TGeoVolume( "ADCPMg", sADCPMg, medADCPMGlass );
   	vADCPMg->SetLineColor(kColorPMG);
   	TGeoTube *sADCPMa   = new TGeoTube( "ADCPMa", kADCPMR1, kADCPMR2, kADCPMlength/2. );
   	TGeoVolume *vADCPMa = new TGeoVolume( "ADCPMa", sADCPMa, medADCPMAlum );
   	vADCPMa->SetLineColor( kColorPMA );
   	TGeoVolume *vADCPM  = new TGeoVolumeAssembly( "ADCPM" );
   	vADCPM->AddNode( vADCPMg, 1 );
   	vADCPM->AddNode( vADCPMa, 1 );

   	/// Sector (Asembly:  Scintillator Pad + Light guide + PM )
   	TGeoVolume *secADC  = new TGeoVolumeAssembly("ADCsec");
   	// Add PAD
   	TGeoCombiTrans *transrot1 = new TGeoCombiTrans( 3*kADCCellSide/4.0, kADCCellSide/2.0, 0., new TGeoRotation("rot",0.,0.,0.) );
   	secADC->AddNode( vADC, 1, transrot1 );
   	// Add Light Guide
   	transrot1 = new TGeoCombiTrans( kADCCellSide + kADCLGHeigth/2.0, kADCLGSideScint/2.0, 0., new TGeoRotation("rot",90.,90.,90.) );
   	secADC->AddNode( vADCLG, 1, transrot1 );
   	// Add PM
   	transrot1 = new TGeoCombiTrans( kADCCellSide + kADCLGHeigth  + kADCLGCoupling + kADCPMlength/2., kADCCellSide/2.0, 0, new TGeoRotation("rot",90.,90.,0.) );
   	secADC->AddNode( vADCPM, 1, transrot1 );

   	// TODO: Add mechanical support

   	/// Assembling ADC adding the 4 sectors                                   //  Sectors
   	TGeoVolume *vADCarray = new TGeoVolumeAssembly("ADC");                    //        ^ y
   	vADCarray->AddNode( secADC, 1 );                                          //        |   
   	vADCarray->AddNode( secADC, 2, new TGeoRotation("rot", 0.  , 180., 0.) ); //   4    |   1
   	vADCarray->AddNode( secADC, 3, new TGeoRotation("rot", 180., 0.,   0.) ); // --------------->  x  
   	vADCarray->AddNode( secADC, 4, new TGeoRotation("rot", 180., 180., 0.) ); //   3    |   2
                                                                             //        |
 	// here I add ADC to AD volume

   	//const Float_t kPosADC = -1902.75;  // -1902.75 (with 4cm thick) puts the ADC just next to the YSAA3_CC_BLOCK
	const Float_t kPosADC = -1900.0;
	ad->AddNode(vADCarray,3,new TGeoTranslation(0,0,kPosADC));
	// add second array
	if (GetADCToInstalled())
	{
		const Float_t kPosADC2 = -1895.0;
		ad->AddNode(vADCarray,4,new TGeoTranslation(0,0,kPosADC2));
	}

	// at the end, I add "AD" volume into ALICE

	alice->AddNode(ad,1);
}

//_____________________________________________________________________________
void AliADv1::AddAlignableVolumes() const
{
   //
   // Create entries for alignable volumes associating the symbolic volume
   // name with the corresponding volume path. Needs to be syncronized with
   // eventual changes in the geometry.
   //
   // ADA and ADC 


   
   TString volpath1 = "/ALIC_1/AD_1/ADC_3";
   TString volpath2 = "/ALIC_1/AD_1/ADC_4";
   TString volpath3 = "/ALIC_1/AD_1/ADA_1";
   TString volpath4 = "/ALIC_1/AD_1/ADA_2";
 
   TString symname1 = "AD/ADC3";
   TString symname2 = "AD/ADC4"; 
   TString symname3 = "AD/ADA1";
   TString symname4 = "AD/ADA2"; 
   
   if ( !gGeoManager->SetAlignableEntry(symname1.Data(), volpath1.Data()) )
      AliFatal(Form( "Alignable entry %s not created. Volume path %s not valid", symname1.Data(), volpath1.Data()) );
   if ( GetADCToInstalled() && !gGeoManager->SetAlignableEntry(symname2.Data(), volpath2.Data()) )
      AliFatal(Form( "Alignable entry %s not created. Volume path %s not valid", symname2.Data(), volpath2.Data()) );
   if ( !gGeoManager->SetAlignableEntry(symname3.Data(), volpath3.Data()) )
      AliFatal(Form( "Alignable entry %s not created. Volume path %s not valid", symname3.Data(), volpath3.Data()) );
   if ( GetADAToInstalled() && !gGeoManager->SetAlignableEntry(symname4.Data(), volpath4.Data()) )
      AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", symname4.Data(), volpath4.Data()) );
   
}


//_____________________________________________________________________________
void AliADv1::StepManager()
{

   //
   // Routine called at every step in the AD
   //

   // ADA and ADC static Variables         //
   //static  Int_t   numStep_ad = 0;         //
//   static  Int_t   vol_ad[2];              //
  
   /////////////////////////////////////////////////////////////////////////
   // ADA and ADC
   /////////////////////////////////////////////////////////////////////////
      
      
   // Get sensitive volumes id (scintillator pads)
   static Int_t idADA = gMC->VolId( "ADApad" );
   static Int_t idADC = gMC->VolId( "ADCpad" );
   
   // We keep only charged tracks : 
   // if ( !gMC->TrackCharge() || !gMC->IsTrackAlive() ) return;   
   // We keep charged and non-charged tracks : 
   if ( !gMC->IsTrackAlive() ) return;   
   
   Int_t copy;
   Int_t current_volid = gMC->CurrentVolID( copy );

   // check is the track is in a sensitive volume
   if( current_volid != idADA && current_volid != idADC ) {
      return; // not in the sensitive volume 
   }
   
   // First read the position, otherwise weird reults! //ecv
   Double_t s[3];
   Float_t  x[3];
   gMC->TrackPosition( s[0], s[1], s[2] );
   for ( Int_t j=0; j<3; j++ ) x[j] = s[j];
   
   // Set detectro type: ADA or ADC
   Int_t detType = (current_volid == idADA ) ? 0 : 1;
   
   // Get sector copy (1,2,3,4) ( 1 level up from pad )
   Int_t sect;
   gMC->CurrentVolOffID( 1, sect );

   // Get Detector copy (1,2) ( 2 levels up from pad )
   Int_t detc;
   gMC->CurrentVolOffID( 2, detc );
   
   // Sector number 
   // ADA1 = 10-14
   // ADA2 = 20-24
   // ADC1 = 30-34
   // ADC2 = 40-44
   Int_t sectorNumber_AD = detType*20 + detc*10 + sect;
   
   Double_t lightYield_ad;
   Double_t photoCathodeEfficiency;
  
   if( detType == 1 )  {
      lightYield_ad          = fADCLightYield;
      photoCathodeEfficiency = fADCPhotoCathodeEfficiency;
   } else  {
      lightYield_ad          = fADALightYield;
      photoCathodeEfficiency = fADAPhotoCathodeEfficiency;
   }
      
   Float_t destep_ad = gMC->Edep();
   Float_t step_ad   = gMC->TrackStep();
   Int_t  nPhotonsInStep_ad = Int_t( destep_ad / (lightYield_ad * 1e-9) ); 
   nPhotonsInStep_ad = gRandom->Poisson( nPhotonsInStep_ad );
   
   static  Float_t eloss_ad    = 0.;
   static  Float_t tlength_ad  = 0.;   
   static  Int_t   nPhotons_ad = 0;      
   static  Float_t hits_ad[11];            
   static  Int_t   vol_ad[5];

   eloss_ad   += destep_ad;
   tlength_ad += step_ad;  
 
   if ( gMC->IsTrackEntering() ) { 
      nPhotons_ad = nPhotonsInStep_ad;
      Double_t p[4];
      gMC->TrackMomentum( p[0], p[1], p[2], p[3] );
      Float_t pt  = TMath::Sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ); 
      TParticle *par = gAlice->GetMCApp()->Particle(gAlice->GetMCApp()->GetCurrentTrackNumber());
      Int_t imo = par->GetFirstMother();
      Int_t pdgMo = 0;
      if ( imo > 0 ) {
         TParticle * pmot = gAlice->GetMCApp()->Particle(imo);
         pdgMo = pmot->GetPdgCode();
      }

      // Set integer values
      vol_ad[0]  = par->GetStatusCode();    // secondary flag //ecv
      vol_ad[1]  = par->GetPdgCode();       // PDG code
      vol_ad[2]  = pdgMo;                   // PDG of the mother
      // Set float values
      hits_ad[0]  = x[0];     // X
      hits_ad[1]  = x[1];     // Y 
      hits_ad[2]  = x[2];     // Z       
      hits_ad[3]  = p[3];     // kinetic energy of the entering particle
      hits_ad[4]  = pt;       // Pt
      hits_ad[5]  = p[0];     // Px
      hits_ad[6]  = p[1];     // Py
      hits_ad[7]  = p[2];     // Pz
      hits_ad[8]  = 1.0e09*gMC->TrackTime(); // in ns!
  
      tlength_ad = 0.0;
      eloss_ad   = 0.0; 
      
      return; // without return, we count 2 times nPhotonsInStep_ad !!!???
   }
   
   nPhotons_ad += nPhotonsInStep_ad;

   if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared() ) {

      // Set integer values
      vol_ad[3]  = nPhotons_ad;
	// brutal correction for ADA_1
      if (sectorNumber_AD==11) sectorNumber_AD=0;
      if (sectorNumber_AD==12) sectorNumber_AD=1;
      if (sectorNumber_AD==13) sectorNumber_AD=2;
      if (sectorNumber_AD==14) sectorNumber_AD=3;

	// same for ADA_2
      if (sectorNumber_AD==21) sectorNumber_AD=4;
      if (sectorNumber_AD==22) sectorNumber_AD=5;
      if (sectorNumber_AD==23) sectorNumber_AD=6;
      if (sectorNumber_AD==24) sectorNumber_AD=7;

	// brutal correction for ADC_3
      if (sectorNumber_AD==51) sectorNumber_AD=8;
      if (sectorNumber_AD==52) sectorNumber_AD=9;
      if (sectorNumber_AD==53) sectorNumber_AD=10;
      if (sectorNumber_AD==54) sectorNumber_AD=11;

	// same for ADC_4
      if (sectorNumber_AD==61) sectorNumber_AD=12;
      if (sectorNumber_AD==62) sectorNumber_AD=13;
      if (sectorNumber_AD==63) sectorNumber_AD=14;
      if (sectorNumber_AD==64) sectorNumber_AD=15;


      vol_ad[4]  = sectorNumber_AD;  // sector number (scintillator ID)
      // Set float values
      hits_ad[9]  = tlength_ad;    // track lenght inside ADC or ADA
      hits_ad[10] = eloss_ad;      // energy loss
      Int_t track = gAlice->GetMCApp()->GetCurrentTrackNumber();
      AddHit( track, vol_ad, hits_ad ); // <-- this is in AliAD.cxx
      tlength_ad        = 0.0;
      eloss_ad          = 0.0; 
      nPhotons_ad       = 0;
   }
       
   //   Do we need track reference ????
   // if( gMC->IsTrackEntering() || gMC->IsTrackExiting() ) {
   //    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), 49);
   // }
}
//_________________________________________________________
void AliADv1::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
	TClonesArray &lhits = *fHits;
	new(lhits[fNhits++]) AliADhit(fIshunt,track,vol,hits);
}
//_________________________________________________________
void AliADv1::AddDigits(Int_t* track, Int_t module, Float_t time)
{
	TClonesArray &ldigits = *fDigits;
	new(ldigits[fNdigits++]) AliADdigit(track,module,time);
}
//_________________________________________________________
void AliADv1::MakeBranch(Option_t *option)
{

	// Create branches in the current tree
	TString branchname(Form("%s",GetName()));
	AliDebug(2,Form("fBufferSize = %d",fBufferSize));
	const char *cH = strstr(option,"H");
	if (fHits && fLoader->TreeH() && cH)
	{
		fLoader->TreeH()->Branch(branchname.Data(),&fHits,fBufferSize);
		AliDebug(2,Form("Making Branch %s for hits",branchname.Data()));
	}
	const char *cD = strstr(option,"D");
  	if (fDigits   && fLoader->TreeD() && cD) 
	{
    		fLoader->TreeD()->Branch(branchname.Data(),&fDigits, fBufferSize);
    		AliDebug(2,Form("Making Branch %s for digits",branchname.Data()));
  	}  
}
