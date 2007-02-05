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

/* $Id$ */

/////////////////////////////////////////////////////////////////////
//                                                                 //
// T0 ( T-zero) detector  version 0                        //
//
//Begin Html       
/*
<img src="gif/AliT0v2Class.gif">
*/
//End Html
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#define RIGHT_ARRAY		1
#define LEFT_ARRAY		2

#include <Riostream.h>
#include <stdlib.h>

#include <TGeometry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TVirtualMC.h>

#include "AliMagF.h"
#include "AliRun.h"
#include "AliT0hit.h"
#include "AliT0v2.h"
#include "AliMC.h"
//#include "AliT0hitPhoton.h"
//#include "TGeant3.h"

ClassImp(AliT0v2)


//////////////////////////////////////////////////////////////////////
// Standart constructor for T0 Detector version 2
//////////////////////////////////////////////////////////////////////
AliT0v2::AliT0v2(const char *name, const char *title):
 AliT0(name,title)
{
  fIdSens1 = 0;
//  setBufferSize(128000);
}


//////////////////////////////////////////////////////////////////////
// Create the geometry of T0 Detector version 2
//////////////////////////////////////////////////////////////////////
void AliT0v2::CreateGeometry()
{
  /*
     <img src="gif/AliT0v2.gif">
  */

  Int_t 	*idtmed;
  Int_t 	idrotm[999];
  Int_t 	i, j;
  Float_t 	x, y, z;
  Float_t 	zRight;
  Float_t	zLeft;
  Float_t  	phi[3];
  Float_t 	theta;
  Double_t 	angel;
  
  Float_t 	pstart[3]	= {4.5,  10.7,   5.3};
  Float_t 	pinstart[3]	= {0.,    1.31,  5.25};
  Float_t 	ppmt[3]		= {0.,    1.31,  3.5};
  Float_t 	ptop[3]		= {0.,    1.3,   1.5};
//  Float_t 	preg[3]		= {0.,    1.3,   0.005};
  Float_t 	preg[3]		= {0.,    0.875, 0.005};
  Float_t 	pdes[3]		= {0.875, 1.3,   0.005};


  zRight = 75.;
  zLeft  = 350.;
  theta  = (180 / TMath::Pi()) * TMath::ATan(6.5 / zRight);
  angel  = 2 * TMath::Pi() / 12;
  idtmed = fIdtmed->GetArray();
    

  AliMatrix (idrotm[901], 90., 0., 90., 90., 180., 0.);
    
  gMC->Gsvolu ("0RST", "TUBE", idtmed[3], pstart, 3);
  gMC->Gsvolu ("0LST", "TUBE", idtmed[3], pstart, 3);
  gMC->Gspos ("0RST", 1, "ALIC", 0., 0., -zRight, 0, "ONLY");
  gMC->Gspos ("0LST", 1, "ALIC", 0., 0., zLeft, idrotm[901], "ONLY");

//  12 unit: PMT + divider
  gMC->Gsvolu("0INS", "TUBE", idtmed[3], pinstart, 3);   
  z = 0;
  for (i = 0; i < 12; i++)
  {
     x = 6.5 * TMath::Sin(i * angel);
     y = 6.5 * TMath::Cos(i * angel);
	
     phi[0] = -30 * i;
     phi[1] = 90 - i * 30;
     phi[2] = 90 - i * 30;
     for (j = 0; j < 3; j++)
	if (phi[j] < 0)  phi[j] += 360;
	
     AliMatrix (idrotm[902 + i], 90.,         phi[0],
		                 90. + theta, phi[1],
		                 theta,       phi[2]);  

     gMC->Gspos ("0INS", i + 1, "0RST", x, y, z, idrotm[902 + i], "ONLY");
     gMC->Gspos ("0INS", i + 13, "0LST", x, y, z, 0, "ONLY");
  }

  gMC->Gsvolu ("0PMT", "TUBE", idtmed[1], ppmt, 3);     
  x = y = 0;      
  z = -pinstart[2] + ppmt[2];
  gMC->Gspos ("0PMT", 1, "0INS", x, y, z, 0, "ONLY");
    
// PMT
  // Entry window (glass)
  gMC->Gsvolu ("0TOP", "TUBE", idtmed[6], ptop, 3);
  z = -ppmt[2] + ptop[2];
  gMC->Gspos ("0TOP", 1, "0PMT", 0, 0, z, 0, "ONLY");

  gMC->Gsvolu ("0REG", "TUBE", idtmed[6], preg, 3);
  z = -ppmt[2] + 2 * ptop[2] + preg[2];
  gMC->Gspos ("0REG", 1, "0PMT", 0, 0, z, 0, "ONLY");

  gMC->Gsvolu ("0DES", "TUBE", idtmed[6], pdes, 3);
  z = -ppmt[2] + 2 * ptop[2] + preg[2];
  gMC->Gspos ("0DES", 1, "0PMT", 0, 0, z, 0, "ONLY");
}   
 

//////////////////////////////////////////////////////////////////////
// Definition of avaible T0 materials
//////////////////////////////////////////////////////////////////////
void AliT0v2::CreateMaterials()
{
   Int_t isxfld   = gAlice->Field()->Integ();
   Float_t sxmgmx = gAlice->Field()->Max();
   
   Float_t a, z, d, radl, absl, buf[1];
   Int_t nbuf;

// Scintillator CH
   Float_t 	ascin[2] = {1.01, 12.01};
   Float_t 	zscin[2] = {1.,    6.};
   Float_t 	wscin[2] = {1.,    1.};
   Float_t 	denscin  = 1.03;
   
// PMT glass SiO2
   Float_t 	aglass[2] = {28.0855, 15.9994};
   Float_t 	zglass[2] = {14.,      8.};
   Float_t 	wglass[2] = { 1.,      2.};
   Float_t 	dglass    = 2.65;
   
// Ceramic   97.2% Al2O3 , 2.8% SiO2
   Float_t 	acer[2], zcer[2];
   Float_t 	wcer[2]	   = { 0.972,     0.028};
   Float_t 	aal2o3[2]  = {26.981539, 15.9994 };
   Float_t 	zal2o3[2]  = {13.,        8.};
   Float_t 	wal2o3[2]  = { 2.,        3.};
   Float_t 	denscer    = 3.6;

// Brass 80% Cu, 20% Zn
   Float_t abrass[2] = {63.546, 65.39};
   Float_t zbrass[2] = {29.,    30.};
   Float_t wbrass[2] = { 0.8,    0.2};
   Float_t denbrass  = 8.96;

//Ribber C6H12S
   Float_t aribber[3] = {12.,  1., 32.};
   Float_t zribber[3] = { 6.,  1., 16.};
   Float_t wribber[3] = { 6., 12.,  1.};
   Float_t denribber  = 0.8;
   
    
   AliMaterial (0, "T0 Steel$", 55.850, 26., 7.87, 1.76, 999);
   AliMaterial (1, "T0 Vacuum$", 1.e-16, 1.e-16, 1.e-16, 1.e16, 999);
   AliMaterial (2, "T0 Air$", 14.61, 7.3, .001205, 30423., 999); 
   
   AliMixture (3, "Al2O3   $", aal2o3, zal2o3, denscer, -2, wal2o3);
   AliMixture (4, "PMT glass   $", aglass, zglass, dglass, -2, wglass);
   char namate[21]="";
   gMC->Gfmate ((*fIdmate)[3], namate, a, z, d, radl, absl, buf, nbuf);
   acer[0] = a;
   zcer[0] = z;
   gMC->Gfmate ((*fIdmate)[4], namate, a, z, d, radl, absl, buf, nbuf);
   acer[1] = a;
   zcer[1] = z;
   
   AliMixture (5, "Scintillator$",ascin,zscin,denscin,-2,wscin);
   AliMixture (6, "Brass    $", abrass, zbrass, denbrass, 2, wbrass);
   AliMixture (7, "Ribber $",aribber,zribber,denribber,-3,wribber);
   AliMixture (9, "Ceramic    $", acer, zcer, denscer, 2, wcer);
   
   AliMedium (1, "T0 Air$", 2, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   AliMedium (2, "Scintillator$", 5, 1, isxfld, sxmgmx, 10., .01, 1., .003, .003);
   AliMedium (3, "Vacuum$", 1, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium (4, "Ceramic$", 9, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium (6, "Glass$", 4, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium (8, "Steel$", 0, 0, isxfld, sxmgmx, 1., .001, 1., .001, .001);
   AliMedium (9, "Ribber  $", 7, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(11, "Brass  $", 6, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
}


//////////////////////////////////////////////////////////////////////
// Draw a shaded view of the Forward multiplicity detector version 2
//////////////////////////////////////////////////////////////////////
void AliT0v2::DrawModule() const
{  
  //Set ALIC mother transparent
  gMC->Gsatt ("ALIC", "SEEN", 0);

  //Set volumes visible
  gMC->Gsatt ("0STA", "SEEN", 0);
  gMC->Gsatt ("0PMT", "SEEN", 1);
  gMC->Gsatt ("0DIV", "SEEN", 1);

  gMC->Gdopt ("hide", "on");
  gMC->Gdopt ("shad", "on");
  gMC->SetClipBox (".");
  gMC->SetClipBox ("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw ("alic", 40, 30, 0, 12, 9.5, 0.7, 0.7);
  gMC->Gdhead (1111, "T-Zero detector");
  gMC->Gdopt ("hide", "off");
}


//////////////////////////////////////////////////////////////////////
// Initialises version 2 of the Forward Multiplicity Detector
//////////////////////////////////////////////////////////////////////
void AliT0v2::Init()
{
//Int_t *idtmed  = gAlice->Idtmed();

  AliT0::Init();
  fIdSens1 = gMC->VolId ("0REG");
// Definition Cherenkov parameters
  const Int_t NUMENTRIES = 32;

  Float_t ppckov[NUMENTRIES] =
            { 2.034E-9, 2.068E-9, 2.103E-9, 2.139E-9,
              2.177E-9, 2.216E-9, 2.256E-9, 2.298E-9,
              2.341E-9, 2.386E-9, 2.433E-9, 2.481E-9,
              2.532E-9, 2.585E-9, 2.640E-9, 2.697E-9,
              2.757E-9, 2.820E-9, 2.885E-9, 2.954E-9,
              3.026E-9, 3.102E-9, 3.181E-9, 3.265E-9,
              3.353E-9, 3.446E-9, 3.545E-9, 3.649E-9,
              3.760E-9, 3.877E-9, 4.002E-9, 4.136E-9 };

  Float_t rindex_qwarz[NUMENTRIES] =
            { 1.458, 1.458, 1.458, 1.458, 1.458, 1.458, 1.458,
              1.458, 1.458, 1.458, 1.458, 1.458, 1.458, 1.458,
              1.458, 1.458, 1.458, 1.458, 1.458, 1.458, 1.458,
              1.458, 1.458, 1.458, 1.458, 1.458, 1.458, 1.458,
              1.458, 1.458, 1.458, 1.458 };

  Float_t rindex_air[NUMENTRIES] =
            { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1. };

  Float_t effic_all[NUMENTRIES] =
            { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1. };
              
  Float_t absor_qwarz[NUMENTRIES] =
  	    { 2000., 2000., 2000., 2000., 2000., 2000., 2000., 
  	      2000., 2000., 2000., 2000., 2000., 2000., 2000.,
  	      2000., 2000., 2000., 2000., 2000., 2000., 2000.,
  	      2000., 2000., 2000., 2000., 2000., 2000., 2000.,
  	      2000., 2000., 2000., 2000. };

  Float_t absor_air[NUMENTRIES] =
  	    { 200., 200., 200., 200., 200., 200., 200., 
  	      200., 200., 200., 200., 200., 200., 200.,
  	      200., 200., 200., 200., 200., 200., 200.,
  	      200., 200., 200., 200., 200., 200., 200.,
  	      200., 200., 200., 200. };
   
  Int_t *idtmed = fIdtmed->GetArray();

   gMC->SetCerenkov (idtmed[6], NUMENTRIES, ppckov, absor_qwarz, effic_all, rindex_qwarz);
   gMC->SetCerenkov (idtmed[1], NUMENTRIES, ppckov, absor_air, effic_all, rindex_air);

  printf ("*** T0 version 2 initialized ***\n");
}


//////////////////////////////////////////////////////////////////////
// Called for every step in the T0 Detector
//////////////////////////////////////////////////////////////////////
void AliT0v2::StepManager()
{
  Int_t 			id;
  Int_t				copy;
//  Int_t				copy1;
  Float_t			xyz[3];
  Float_t			XYZ[3];
  Float_t			hitPhoton[6];
  static Float_t 		hits[7];
  static Float_t 		edep;
  static Int_t 			vol[2];
  TLorentzVector 		pos;
  TLorentzVector		mom;
  
  
  if(!gMC->IsTrackAlive()) return; // particle has disappeared

//  TGeant3  *g3 = (TGeant3*) gMC;
//  g3->Gpcxyz();

  TClonesArray 	&lhits = *fHits;


///////////////////////////////////////////////
// If particles is photon then ...

  if (gMC->TrackPid() == 50000050)
  {
     id = gMC->CurrentVolID(copy);

     // Check the sensetive volume
     if (id == fIdSens1) 
     {
	if (gMC->IsTrackEntering()) 
	{
	   gMC->CurrentVolOffID(2,copy);
	   if (copy < 13) 
	     {
	       vol[0] = RIGHT_ARRAY;
	       vol[1] = copy;
	     }
	   else 
	     {
	       vol[0] = LEFT_ARRAY;
	       vol[1] = copy - 12;
	     }

	   gMC->TrackPosition(pos);
	   gMC->TrackMomentum(mom);
	   XYZ[0] = pos[0];
	   XYZ[1] = pos[1];
	   XYZ[2] = pos[2];
	   gMC->Gmtod (XYZ, xyz, 1);

	   hitPhoton[0] = sqrt (xyz[0] * xyz[0] + xyz[1] * xyz[1]);
	   hitPhoton[1] = mom[0];
	   hitPhoton[2] = mom[1];
	   hitPhoton[3] = mom[2];
	   hitPhoton[4] = 1e9 * gMC->TrackTime();
	   hitPhoton[5] = 1e9 * gMC->Etot();
	   
	   AddHitPhoton (gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hitPhoton);
	}
	gMC->StopTrack();
     }
     else  
	if (id == gMC->VolId ("0DES"))  gMC->StopTrack();
  }

// end photon particalary
///////////////////////////////////////////////


///////////////////////////////////////////////
// If particles is charge then ...
 
  // Float_t charge = gMC->TrackCharge();
  //  if(TMath::Abs(charge) <= 0.) return;
  id = gMC->CurrentVolID(copy);

  
  //  printf("gMC->ckine->ipart %d",gMC->ckine->ipart);
  // Check the sensetive volume
  if(id==gMC->VolId("0TOP") ) {
    if(gMC->IsTrackEntering()) {
      gMC->CurrentVolOffID(2,copy);
      if (copy < 13) 
	{
	  vol[0] = RIGHT_ARRAY;
	  vol[1] = copy;
	}
      else 
	{
	  vol[0] = LEFT_ARRAY;
	  vol[1] = copy - 12;
	}
      
      gMC->TrackPosition(pos);
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];
      Float_t etot = gMC->Etot();
      hits[4] = etot;
      Int_t part = gMC->TrackPid();
      hits[5] = part;
      Float_t ttime = gMC->TrackTime();
      hits[6] = ttime*1e9;
      edep = 0;
    }
    if(gMC->IsTrackInside()) 	{
      Float_t de = gMC->Edep(); 
      edep = edep + de;
      //       printf ("E deposition %f\n",edep);
      //    for (i=0; i<=6; i++){
      //    printf(" HITS on T0 inside %f\n",hits[i]); } 
    }
    if(gMC->IsTrackExiting())	{
      Float_t de = gMC->Edep(); 
      edep = edep + de;
      hits[3] = edep * 1e3;

      //       for (i=0; i<=6; i++){
      //	 printf(" HITS on T0 Exit %f\n",hits[i]); } 
      //for (i=0; i<=1; i++) { printf("T0 vol %d\n",vol[i]);}
     
      new(lhits[fNhits++]) AliT0hit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);      
    }
  }

// end charge particles particalary
///////////////////////////////////////////////

}
