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

#include "AliRICHv3.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliMagF.h"

#include "AliConst.h"
#include "AliPDG.h"

#include <iostream.h>
#include <TNode.h>
#include <TGeometry.h>
#include <TBRIK.h>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticle.h>


#include "AliRICHGeometry.h"
#include "AliRICHSegmentationV1.h"
#include "AliRICHResponseV0.h"
#include "AliRICHHit.h"

ClassImp(AliRICHv3)

//______________________________________________________________
//    Implementation of the RICH version 3 with azimuthal rotation


AliRICHv3::AliRICHv3(const char *sName, const char *sTitle)
      	  :AliRICH(sName,sTitle)
{
// The named ctor currently creates a single copy of 
// AliRICHGeometry AliRICHSegmentationV1 AliRICHResponseV0
// and initialises the corresponding models of all 7 chambers with these stuctures.
// Note: all chambers share the single copy of models. MUST be changed later (???).
   cout<<ClassName()<<"::named ctor(sName,sTitle)>\n"; // no way to control it as ctor is called before call to SetDebugXXXX()

   fCkovNumber=fFreonProd=fDebugLevel=0;
   
   AliRICHGeometry       *pRICHGeometry    =new AliRICHGeometry;           // ??? to be moved to AlRICHChamber::named ctor
   AliRICHSegmentationV1 *pRICHSegmentation=new AliRICHSegmentationV1;     // ??? to be moved to AlRICHChamber::named ctor
   AliRICHResponseV0     *pRICHResponse    =new AliRICHResponseV0;         // ??? to be moved to AlRICHChamber::named ctor
     
   fChambers = new TObjArray(kNCH);
   for (Int_t i=0; i<kNCH; i++){    
      fChambers->AddAt(new AliRICHChamber,i); // ??? to be changed to named ctor of AliRICHChamber
      SetGeometryModel(i,pRICHGeometry);
      SetSegmentationModel(i,pRICHSegmentation);
      SetResponseModel(i,pRICHResponse);
      ((AliRICHChamber*)fChambers->At(i))->Init(i); // ??? to be removed     
   }
}//AliRICHv3::ctor(const char *pcName, const char *pcTitle)

AliRICHv3::~AliRICHv3()
{
// Dtor deletes RICH models. In future (???) AliRICHChamber will be responsible for that.
   if(IsDebugStart()) cout<<ClassName()<<"::dtor()>\n";
      
   delete GetChamber(0)->GetGeometryModel();
   delete GetChamber(0)->GetResponseModel();
   delete GetChamber(0)->GetSegmentationModel();
}//AliRICHv3::dtor()

void AliRICHv3::CreateMaterials()
{
// Provides material definition for simulation (currently GEANT)    
   if(IsDebugStart()) cout<<ClassName()<<"::CreateMaterials()>\n";
   
   Int_t   isxfld = gAlice->Field()->Integ();
   Float_t sxmgmx = gAlice->Field()->Max();
   Int_t i;
    

    //Photons energy intervals
    Float_t ppckov[26];
    for (i=0;i<26;i++) 
    {
	ppckov[i] = (Float_t(i)*0.1+5.5)*1e-9;
	//printf ("Energy intervals: %e\n",ppckov[i]);
    }
    
    
    //Refraction index for quarz
    Float_t rIndexQuarz[26];
    Float_t  e1= 10.666;
    Float_t  e2= 18.125;
    Float_t  f1= 46.411;
    Float_t  f2= 228.71;
    for (i=0;i<26;i++)
    {
	Float_t ene=ppckov[i]*1e9;
	Float_t a=f1/(e1*e1 - ene*ene);
	Float_t b=f2/(e2*e2 - ene*ene);
	rIndexQuarz[i] = TMath::Sqrt(1. + a + b );
	//printf ("rIndexQuarz: %e\n",rIndexQuarz[i]);
    } 
    
    //Refraction index for opaque quarz, methane and grid
    Float_t rIndexOpaqueQuarz[26];
    Float_t rIndexMethane[26];
    Float_t rIndexGrid[26];
    for (i=0;i<26;i++)
    {
	rIndexOpaqueQuarz[i]=1;
	rIndexMethane[i]=1.000444;
	rIndexGrid[i]=1;
	//printf ("rIndexOpaqueQuarz , etc: %e, %e, %e\n",rIndexOpaqueQuarz[i], rIndexMethane[i], rIndexGrid[i]=1);
    } 
    
    //Absorption index for freon
    Float_t abscoFreon[26] = {179.0987, 179.0987, 179.0987, 179.0987, 179.0987,  179.0987, 179.0987, 179.0987, 
	 		       179.0987, 142.9206, 56.64957, 25.58622, 13.95293, 12.03905, 10.42953, 8.804196, 
			       7.069031, 4.461292, 2.028366, 1.293013, .577267,   .40746,  .334964, 0., 0., 0.};
    

    Float_t abscoQuarz [26] = {105.8, 65.52, 48.58, 42.85, 35.79, 31.262, 28.598, 27.527, 25.007, 22.815, 21.004,
				19.266, 17.525, 15.878, 14.177, 11.719, 9.282, 6.62, 4.0925, 2.601, 1.149, .667, .3627,
				.192, .1497, .10857};
    
    //Absorption index for methane
    Float_t abscoMethane[26];
    for (i=0;i<26;i++) 
    {
	abscoMethane[i]=AbsoCH4(ppckov[i]*1e9); 
	//printf("abscoMethane: %e for energy: %e\n", abscoMethane[i],ppckov[i]*1e9);
    }
    
    //Absorption index for opaque quarz, csi and grid, efficiency for all and grid
    Float_t abscoOpaqueQuarz[26];
    Float_t abscoCsI[26];
    Float_t abscoGrid[26];
    Float_t efficAll[26];
    Float_t efficGrid[26];
    for (i=0;i<26;i++)
    { 
	abscoOpaqueQuarz[i]=1e-5; 
	abscoCsI[i]=1e-4; 
	abscoGrid[i]=1e-4; 
	efficAll[i]=1; 
	efficGrid[i]=1;
    } 
    
    //Efficiency for csi 
    
    Float_t efficCsI[26] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983, 0.010125,
			     0.0242999997, 0.0405000001, 0.0688500032, 0.105299994, 0.121500008, 0.141749993, 0.157949999,
			     0.162, 0.166050002, 0.167669997, 0.174299985, 0.176789999, 0.179279998, 0.182599992, 0.18592,
			     0.187579989, 0.189239994, 0.190899998, 0.207499996, 0.215799987};
	
    

    //FRESNEL LOSS CORRECTION FOR PERPENDICULAR INCIDENCE AND
    //UNPOLARIZED PHOTONS

    for (i=0;i<26;i++)
    {
	efficCsI[i] = efficCsI[i]/(1.-Fresnel(ppckov[i]*1e9,1.,0)); 
    }
	
    
    Float_t afre[2], agri, amet[2], aqua[2], ahon, zfre[2], zgri, zhon, 
    zmet[2], zqua[2];
    Int_t nlmatfre;
    Float_t densquao;
    Int_t nlmatmet, nlmatqua;
    Float_t wmatquao[2], rIndexFreon[26];
    Float_t aquao[2], epsil, stmin, zquao[2];
    Int_t nlmatquao;
    Float_t radlal, densal, tmaxfd, deemax, stemax;
    Float_t aal, zal, radlgri, densfre, radlhon, densgri, denshon,densqua, densmet, wmatfre[2], wmatmet[2], wmatqua[2];
    
    Int_t *idtmed = fIdtmed->GetArray()-999;
    
    // --- Photon energy (GeV) 
    // --- Refraction indexes 
    for (i = 0; i < 26; ++i) {
      rIndexFreon[i] = ppckov[i] * .0172 * 1e9 + 1.177;
      //rIndexFreon[i] = 1;
	//printf ("rIndexFreon: %e \n efficCsI: %e for energy: %e\n",rIndexFreon[i], efficCsI[i], ppckov[i]);
    }
            
    // --- Detection efficiencies (quantum efficiency for CsI) 
    // --- Define parameters for honeycomb. 
    //     Used carbon of equivalent rad. lenght 
    
    ahon    = 12.01;
    zhon    = 6.;
    denshon = 0.1;
    radlhon = 18.8;
    
    // --- Parameters to include in GSMIXT, relative to Quarz (SiO2) 
    
    aqua[0]    = 28.09;
    aqua[1]    = 16.;
    zqua[0]    = 14.;
    zqua[1]    = 8.;
    densqua    = 2.64;
    nlmatqua   = -2;
    wmatqua[0] = 1.;
    wmatqua[1] = 2.;
    
    // --- Parameters to include in GSMIXT, relative to opaque Quarz (SiO2) 
    
    aquao[0]    = 28.09;
    aquao[1]    = 16.;
    zquao[0]    = 14.;
    zquao[1]    = 8.;
    densquao    = 2.64;
    nlmatquao   = -2;
    wmatquao[0] = 1.;
    wmatquao[1] = 2.;
    
    // --- Parameters to include in GSMIXT, relative to Freon (C6F14) 
    
    afre[0]    = 12.;
    afre[1]    = 19.;
    zfre[0]    = 6.;
    zfre[1]    = 9.;
    densfre    = 1.7;
    nlmatfre   = -2;
    wmatfre[0] = 6.;
    wmatfre[1] = 14.;
    
    // --- Parameters to include in GSMIXT, relative to methane (CH4) 
    
    amet[0]    = 12.01;
    amet[1]    = 1.;
    zmet[0]    = 6.;
    zmet[1]    = 1.;
    densmet    = 7.17e-4;
    nlmatmet   = -2;
    wmatmet[0] = 1.;
    wmatmet[1] = 4.;
    
    // --- Parameters to include in GSMIXT, relative to anode grid (Cu) 
  
    agri    = 63.54;
    zgri    = 29.;
    densgri = 8.96;
    radlgri = 1.43;
    
    // --- Parameters to include in GSMATE related to aluminium sheet 
    
    aal    = 26.98;
    zal    = 13.;
    densal = 2.7;
    radlal = 8.9;

    // --- Glass parameters

    Float_t aglass[5]={12.01, 28.09, 16.,   10.8,  23.};
    Float_t zglass[5]={ 6.,   14.,    8.,    5.,   11.};
    Float_t wglass[5]={ 0.5,  0.105, 0.355, 0.03,  0.01};
    Float_t dglass=1.74;

    
    AliMaterial(1, "Air     $", 14.61, 7.3, .001205, 30420., 67500);
    AliMaterial(6, "HON", ahon, zhon, denshon, radlhon, 0);
    AliMaterial(16, "CSI", ahon, zhon, denshon, radlhon, 0);
    AliMixture(20, "QUA", aqua, zqua, densqua, nlmatqua, wmatqua);
    AliMixture(21, "QUAO", aquao, zquao, densquao, nlmatquao, wmatquao);
    AliMixture(30, "FRE", afre, zfre, densfre, nlmatfre, wmatfre);
    AliMixture(40, "MET", amet, zmet, densmet, nlmatmet, wmatmet);
    AliMixture(41, "METG", amet, zmet, densmet, nlmatmet, wmatmet);
    AliMaterial(11, "GRI", agri, zgri, densgri, radlgri, 0);
    AliMaterial(50, "ALUM", aal, zal, densal, radlal, 0);
    AliMixture(32, "GLASS",aglass, zglass, dglass, 5, wglass);
    AliMaterial(31, "COPPER$",   63.54,    29.,   8.96,  1.4, 0.);
    
    tmaxfd = -10.;
    stemax = -.1;
    deemax = -.2;
    epsil  = .001;
    stmin  = -.001;
    
    AliMedium(1, "DEFAULT MEDIUM AIR$", 1, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(2, "HONEYCOMB$", 6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(3, "QUARZO$", 20, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(4, "FREON$", 30, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(5, "METANO$", 40, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(6, "CSI$", 16, 1, isxfld, sxmgmx,tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(7, "GRIGLIA$", 11, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(8, "QUARZOO$", 21, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(9, "GAP$", 41, 1, isxfld, sxmgmx,tmaxfd, .1, -deemax, epsil, -stmin);
    AliMedium(10, "ALUMINUM$", 50, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(11, "GLASS", 32, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(12, "PCB_COPPER", 31, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    

    gMC->SetCerenkov(idtmed[1000], 26, ppckov, abscoMethane, efficAll, rIndexMethane);
    gMC->SetCerenkov(idtmed[1001], 26, ppckov, abscoMethane, efficAll, rIndexMethane);
    gMC->SetCerenkov(idtmed[1002], 26, ppckov, abscoQuarz, efficAll,rIndexQuarz);
    gMC->SetCerenkov(idtmed[1003], 26, ppckov, abscoFreon, efficAll,rIndexFreon);
    gMC->SetCerenkov(idtmed[1004], 26, ppckov, abscoMethane, efficAll, rIndexMethane);
    gMC->SetCerenkov(idtmed[1005], 26, ppckov, abscoCsI, efficCsI, rIndexMethane);
    gMC->SetCerenkov(idtmed[1006], 26, ppckov, abscoGrid, efficGrid, rIndexGrid);
    gMC->SetCerenkov(idtmed[1007], 26, ppckov, abscoOpaqueQuarz, efficAll, rIndexOpaqueQuarz);
    gMC->SetCerenkov(idtmed[1008], 26, ppckov, abscoMethane, efficAll, rIndexMethane);
    gMC->SetCerenkov(idtmed[1009], 26, ppckov, abscoGrid, efficGrid, rIndexGrid);
    gMC->SetCerenkov(idtmed[1010], 26, ppckov, abscoOpaqueQuarz, efficAll, rIndexOpaqueQuarz);
}//AliRICHv3::CreateMaterials()


void AliRICHv3::CreateGeometry()
{
// Provides geometry structure for simulation (currently GEANT volumes tree)         
   if(IsDebugStart()) cout<<ClassName()<<"::CreateGeometry()>\n";

  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH"); 
  AliRICHSegmentationV0*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHChamber*       iChamber;

  iChamber = &(pRICH->Chamber(0));
  segmentation=(AliRICHSegmentationV0*) iChamber->GetSegmentationModel(0);
  geometry=iChamber->GetGeometryModel();

  Float_t distance;
  distance = geometry->GetFreonThickness()/2 + geometry->GetQuartzThickness() + geometry->GetGapThickness();
  geometry->SetRadiatorToPads(distance);
    
  //Opaque quartz thickness
  Float_t oqua_thickness = .5;
  //CsI dimensions


  Float_t csi_width = segmentation->Npx()*segmentation->Dpx() + segmentation->DeadZone();
  Float_t csi_length = segmentation->Npy()*segmentation->Dpy() + 2*segmentation->DeadZone();
  
  
  Int_t *idtmed = fIdtmed->GetArray()-999;
    
    Int_t i;
    Float_t zs;
    Int_t idrotm[1099];
    Float_t par[3];
    
    // --- Define the RICH detector 
    //     External aluminium box 
    par[0] = 68.8;
    par[1] = 13;                 //Original Settings
    par[2] = 70.86;
    gMC->Gsvolu("RICH", "BOX ", idtmed[1009], par, 3);
    
    //     Air 
    par[0] = 66.3;
    par[1] = 13;                 //Original Settings
    par[2] = 68.35;
    gMC->Gsvolu("SRIC", "BOX ", idtmed[1000], par, 3);
    
    //    Air 2 (cutting the lower part of the box)
    
    par[0] = 1.25;
    par[1] = 3;                 //Original Settings
    par[2] = 70.86;
    gMC->Gsvolu("AIR2", "BOX ", idtmed[1000], par, 3);

    //    Air 3 (cutting the lower part of the box)
    
    par[0] = 66.3;
    par[1] = 3;                 //Original Settings
    par[2] = 1.2505;
    gMC->Gsvolu("AIR3", "BOX ", idtmed[1000], par, 3);
    
    //     Honeycomb 
    par[0] = 66.3;
    par[1] = .188;                 //Original Settings
    par[2] = 68.35;
    gMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);
    
    //     Aluminium sheet 
    par[0] = 66.3;
    par[1] = .025;                 //Original Settings
    par[2] = 68.35;
    /*par[0] = 66.5;
    par[1] = .025;
    par[2] = 63.1;*/
    gMC->Gsvolu("ALUM", "BOX ", idtmed[1009], par, 3);
    
    //     Quartz 
    par[0] = geometry->GetQuartzWidth()/2;
    par[1] = geometry->GetQuartzThickness()/2;
    par[2] = geometry->GetQuartzLength()/2;
    gMC->Gsvolu("QUAR", "BOX ", idtmed[1002], par, 3);
    
    //     Spacers (cylinders) 
    par[0] = 0.;
    par[1] = .5;
    par[2] = geometry->GetFreonThickness()/2;
    gMC->Gsvolu("SPAC", "TUBE", idtmed[1002], par, 3);
    
    //     Feet (freon slabs supports)

    par[0] = .7;
    par[1] = .3;
    par[2] = 1.9;
    gMC->Gsvolu("FOOT", "BOX", idtmed[1009], par, 3);

    //     Opaque quartz 
    par[0] = geometry->GetQuartzWidth()/2;
    par[1] = .2;
    par[2] = geometry->GetQuartzLength()/2;
    gMC->Gsvolu("OQUA", "BOX ", idtmed[1007], par, 3);
  
    //     Frame of opaque quartz
    par[0] = geometry->GetOuterFreonWidth()/2;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetOuterFreonLength()/2; 
    gMC->Gsvolu("OQF1", "BOX ", idtmed[1007], par, 3);

    par[0] = geometry->GetInnerFreonWidth()/2;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetInnerFreonLength()/2; 
    gMC->Gsvolu("OQF2", "BOX ", idtmed[1007], par, 3);
    
    
    //     Freon 
    par[0] = geometry->GetOuterFreonWidth()/2 - oqua_thickness;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetOuterFreonLength()/2 - 2*oqua_thickness; 
    gMC->Gsvolu("FRE1", "BOX ", idtmed[1003], par, 3);

    par[0] = geometry->GetInnerFreonWidth()/2 - oqua_thickness;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetInnerFreonLength()/2 - 2*oqua_thickness; 
    gMC->Gsvolu("FRE2", "BOX ", idtmed[1003], par, 3);
    
    //     Methane 
    par[0] = csi_width/2;
    par[1] = geometry->GetGapThickness()/2;
    par[2] = csi_length/2;
    gMC->Gsvolu("META", "BOX ", idtmed[1004], par, 3);
    
    //     Methane gap 
    par[0] = csi_width/2;
    par[1] = geometry->GetProximityGapThickness()/2;
    par[2] = csi_length/2;
    gMC->Gsvolu("GAP ", "BOX ", idtmed[1008], par, 3);
    
    //     CsI photocathode 
    par[0] = csi_width/2;
    par[1] = .25;
    par[2] = csi_length/2;
    gMC->Gsvolu("CSI ", "BOX ", idtmed[1005], par, 3);
    
    //     Anode grid 
    par[0] = 0.;
    par[1] = .001;
    par[2] = 20.;
    gMC->Gsvolu("GRID", "TUBE", idtmed[1006], par, 3);

    // Wire supports
    // Bar of metal
    
    par[0] = csi_width/2;
    par[1] = 1.05;
    par[2] = 1.05;
    gMC->Gsvolu("WSMe", "BOX ", idtmed[1009], par, 3);

    // Ceramic pick up (base)
    
    par[0] =  csi_width/2;
    par[1] = .25;
    par[2] = 1.05;
    gMC->Gsvolu("WSG1", "BOX ", idtmed[1010], par, 3);

    // Ceramic pick up (head)

    par[0] = csi_width/2;
    par[1] = .1;
    par[2] = .1;
    gMC->Gsvolu("WSG2", "BOX ", idtmed[1010], par, 3);

    // Aluminium supports for methane and CsI
    // Short bar

    par[0] = csi_width/2;
    par[1] = geometry->GetGapThickness()/2 + .25;
    par[2] = (68.35 - csi_length/2)/2;
    gMC->Gsvolu("SMSH", "BOX", idtmed[1009], par, 3);
    
    // Long bar

    par[0] = (66.3 - csi_width/2)/2;
    par[1] = geometry->GetGapThickness()/2 + .25;
    par[2] = csi_length/2 + 68.35 - csi_length/2;
    gMC->Gsvolu("SMLG", "BOX", idtmed[1009], par, 3);
    
    // Aluminium supports for freon
    // Short bar

    par[0] = geometry->GetQuartzWidth()/2;
    par[1] = .3;
    par[2] = (68.35 - geometry->GetQuartzLength()/2)/2;
    gMC->Gsvolu("SFSH", "BOX", idtmed[1009], par, 3);
    
    // Long bar

    par[0] = (66.3 - geometry->GetQuartzWidth()/2)/2;
    par[1] = .3;
    par[2] = geometry->GetQuartzLength()/2 + 68.35 - geometry->GetQuartzLength()/2;
    gMC->Gsvolu("SFLG", "BOX", idtmed[1009], par, 3);
    
    // PCB backplane
    
    par[0] = csi_width/2;
    par[1] = .25;
    par[2] = csi_length/4 -.5025;
    gMC->Gsvolu("PCB ", "BOX", idtmed[1011], par, 3);

    
    // Backplane supports

    // Aluminium slab
    
    par[0] = 33.15;
    par[1] = 2;
    par[2] = 21.65;
    gMC->Gsvolu("BACK", "BOX", idtmed[1009], par, 3);
    
    // Big hole
    
    par[0] = 9.05;
    par[1] = 2;
    par[2] = 4.4625;
    gMC->Gsvolu("BKHL", "BOX", idtmed[1000], par, 3);

    // Small hole
    
    par[0] = 5.7;
    par[1] = 2;
    par[2] = 4.4625;
    gMC->Gsvolu("BKHS", "BOX", idtmed[1000], par, 3);

    // Place holes inside backplane support

    gMC->Gspos("BKHS", 1, "BACK", .8 + 5.7,0., .6 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 2, "BACK", -.8 - 5.7,0., .6 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 3, "BACK", .8 + 5.7,0., -.6 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 4, "BACK", -.8 - 5.7,0., -.6 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 5, "BACK", .8 + 5.7,0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 6, "BACK", -.8 - 5.7,0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 7, "BACK", .8 + 5.7,0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 8, "BACK", -.8 - 5.7,0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 1, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., .6 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 2, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., .6 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 3, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., -.6 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 4, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., -.6 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 5, "BACK", .8 + 11.4+ 1.6 + 9.05, 0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 6, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 7, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 8, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");

    
  
    // --- Places the detectors defined with GSVOLU 
    //     Place material inside RICH 
    gMC->Gspos("SRIC", 1, "RICH", 0.,0., 0., 0, "ONLY");
    gMC->Gspos("AIR2", 1, "RICH", 66.3 + 1.2505, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
    gMC->Gspos("AIR2", 2, "RICH", -66.3 - 1.2505, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
    gMC->Gspos("AIR3", 1, "RICH", 0.,  1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, -68.35 - 1.25, 0, "ONLY");
    gMC->Gspos("AIR3", 2, "RICH", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35,  68.35 + 1.25, 0, "ONLY");
    
      
    gMC->Gspos("ALUM", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.025, 0., 0, "ONLY");
    gMC->Gspos("HONE", 1, "SRIC", 0., 1.276- geometry->GetGapThickness()/2  - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .188, 0., 0, "ONLY");
    gMC->Gspos("ALUM", 2, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .025, 0., 0, "ONLY");
    gMC->Gspos("FOOT", 1, "SRIC", 64.95, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, 36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 2, "SRIC", 21.65, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3 , 36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 3, "SRIC", -21.65, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, 36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 4, "SRIC", -64.95, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, 36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 5, "SRIC", 64.95, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, -36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 6, "SRIC", 21.65, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, -36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 7, "SRIC", -21.65, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, -36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 8, "SRIC", -64.95, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, -36.9, 0, "ONLY");
    gMC->Gspos("OQUA", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .2, 0., 0, "ONLY");
    
    // Supports placing

    // Methane supports
    gMC->Gspos("SMLG", 1, "SRIC", csi_width/2 + (66.3 - csi_width/2)/2, 1.276 + .25, 0., 0, "ONLY");
    gMC->Gspos("SMLG", 2, "SRIC", - csi_width/2 - (66.3 - csi_width/2)/2, 1.276 + .25, 0., 0, "ONLY");
    gMC->Gspos("SMSH", 1, "SRIC", 0., 1.276 + .25, csi_length/2 + (68.35 - csi_length/2)/2, 0, "ONLY");
    gMC->Gspos("SMSH", 2, "SRIC", 0., 1.276 + .25, - csi_length/2 - (68.35 - csi_length/2)/2, 0, "ONLY");

    //Freon supports

    Float_t supp_y = 1.276 - geometry->GetGapThickness()/2- geometry->GetQuartzThickness() -geometry->GetFreonThickness() - .2 + .3; //y position of freon supports

    gMC->Gspos("SFLG", 1, "SRIC", geometry->GetQuartzWidth()/2 + (66.3 - geometry->GetQuartzWidth()/2)/2, supp_y, 0., 0, "ONLY");
    gMC->Gspos("SFLG", 2, "SRIC", - geometry->GetQuartzWidth()/2 - (66.3 - geometry->GetQuartzWidth()/2)/2, supp_y, 0., 0, "ONLY");
    gMC->Gspos("SFSH", 1, "SRIC", 0., supp_y, geometry->GetQuartzLength()/2 + (68.35 - geometry->GetQuartzLength()/2)/2, 0, "ONLY");
    gMC->Gspos("SFSH", 2, "SRIC", 0., supp_y, - geometry->GetQuartzLength()/2 - (68.35 - geometry->GetQuartzLength()/2)/2, 0, "ONLY");
    
    AliMatrix(idrotm[1019], 0., 0., 90., 0., 90., 90.);
    

    Int_t nspacers = 30;
    
    for (i = 0; i < nspacers/3; i++) {
	zs = -11.6/2 + (TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE1", 10.5, 0., zs, idrotm[1019], "ONLY");  //Original settings 
    }
    
    for (i = nspacers/3; i < (nspacers*2)/3; i++) {
	zs = -11.6/2 + (nspacers/3 + TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE1", 0, 0., zs, idrotm[1019], "ONLY");  //Original settings 
    }
    
    for (i = (nspacers*2)/3; i < nspacers; ++i) {
	zs = -11.6/2 + ((nspacers*2)/3 + TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE1", -10.5, 0., zs, idrotm[1019], "ONLY"); //Original settings  
    }

    for (i = 0; i < nspacers/3; i++) {
	zs = -11.6/2 + (TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE2", 10.5, 0., zs, idrotm[1019], "ONLY");  //Original settings 
    }
    
    for (i = nspacers/3; i < (nspacers*2)/3; i++) {
	zs = -11.6/2 + (nspacers/3 + TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE2", 0, 0., zs, idrotm[1019], "ONLY");  //Original settings 
    }
    
    for (i = (nspacers*2)/3; i < nspacers; ++i) {
	zs = -11.6/2 + ((nspacers*2)/3 + TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE2", -10.5, 0., zs, idrotm[1019], "ONLY"); //Original settings  
    }

    
    gMC->Gspos("FRE1", 1, "OQF1", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("FRE2", 1, "OQF2", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("OQF1", 1, "SRIC", geometry->GetOuterFreonWidth()/2 + geometry->GetInnerFreonWidth()/2 + 2, 1.276 - geometry->GetGapThickness()/2- geometry->GetQuartzThickness() -geometry->GetFreonThickness()/2, 0., 0, "ONLY"); //Original settings (31.3)
    gMC->Gspos("OQF2", 2, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()/2, 0., 0, "ONLY");          //Original settings 
    gMC->Gspos("OQF1", 3, "SRIC", - (geometry->GetOuterFreonWidth()/2 + geometry->GetInnerFreonWidth()/2) - 2, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()/2, 0., 0, "ONLY");       //Original settings (-31.3)
    gMC->Gspos("QUAR", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness()/2, 0., 0, "ONLY");
    gMC->Gspos("GAP ", 1, "META", 0., geometry->GetGapThickness()/2 - geometry->GetProximityGapThickness()/2 - 0.0001, 0., 0, "ONLY");
    gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
    gMC->Gspos("CSI ", 1, "SRIC", 0., 1.276 + geometry->GetGapThickness()/2 + .25, 0., 0, "ONLY");
    printf("CSI pos: %f\n",1.276 + geometry->GetGapThickness()/2 + .25);
   
    // Wire support placing

    gMC->Gspos("WSG2", 1, "GAP ", 0., geometry->GetProximityGapThickness()/2 - .1, 0., 0, "ONLY");
    gMC->Gspos("WSG1", 1, "CSI ", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("WSMe", 1, "SRIC ", 0., 1.276 + geometry->GetGapThickness()/2 + .5 + 1.05, 0., 0, "ONLY");

    // Backplane placing
    
    gMC->Gspos("BACK", 1, "SRIC ", -33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, 43.3, 0, "ONLY");
    gMC->Gspos("BACK", 2, "SRIC ", 33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2 , 43.3, 0, "ONLY");
    gMC->Gspos("BACK", 3, "SRIC ", -33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
    gMC->Gspos("BACK", 4, "SRIC ", 33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
    gMC->Gspos("BACK", 5, "SRIC ", 33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");
    gMC->Gspos("BACK", 6, "SRIC ", -33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");

    // PCB placing
    
    gMC->Gspos("PCB ", 1, "SRIC ", 0.,  1.276 + geometry->GetGapThickness()/2 + .5 + 1.05, csi_width/4 + .5025 + 2.5, 0, "ONLY");
    gMC->Gspos("PCB ", 2, "SRIC ", 0.,  1.276 + geometry->GetGapThickness()/2 + .5 + 1.05, -csi_width/4 - .5025 - 2.5, 0, "ONLY");

// Place chambers into mother volume ALIC
           
   Double_t dOffset        = geometry->GetOffset() - geometry->GetGapThickness()/2;  // distance from center of mother volume ALIC to methane
   
   Double_t dAlpha         = geometry->GetAlphaAngle(); // angle between centers of chambers - y-z plane
   Double_t dAlphaRad      = dAlpha*kDegrad;
   
   Double_t dBeta          = geometry->GetBetaAngle();   // angle between center of chambers - y-x plane
   Double_t dBetaRad       = dBeta*kDegrad;
   
   Double_t dRotAngle      = geometry->GetRotationAngle();     // the whole RICH is to be rotated in x-y plane + means clockwise rotation 
   Double_t dRotAngleRad   = dRotAngle*kDegrad;
   
    
   TRotMatrix *pRotMatrix; // tmp pointer
   
   TVector3 vector(0,dOffset,0); // Position of chamber 2 without rotation
    
// Chamber 0  standalone (no other chambers in this row) 
                   AliMatrix(idrotm[1000],      90, -dRotAngle           , 90-dAlpha     , 90-dRotAngle        , dAlpha        , -90           ); 
   pRotMatrix=new TRotMatrix("rot993","rot993", 90, -dRotAngle           , 90-dAlpha     , 90-dRotAngle        , dAlpha        , -90           );

   vector.SetXYZ(0,dOffset,0);  vector.RotateX(dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("RICH",1,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1000], "ONLY");           
   Chamber(0).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
// Chamber 1   
                   AliMatrix(idrotm[1001],      90, -dBeta-dRotAngle     , 90            , 90-dBeta-dRotAngle  , 0             ,   0           );  
   pRotMatrix=new TRotMatrix("rot994","rot994", 90, -dBeta-dRotAngle     , 90            , 90-dBeta-dRotAngle  , 0             ,   0           );  
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(-dBetaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("RICH",2,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1001], "ONLY");           
   Chamber(1).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
// Chamber 2   the top one with no Alpha-Beta rotation
                   AliMatrix(idrotm[1002],      90, -dRotAngle           , 90            , 90-dRotAngle        , 0             ,   0           );
   pRotMatrix=new TRotMatrix("rot995","rot995", 90, -dRotAngle           , 90            , 90-dRotAngle        , 0             ,   0           );
   
   vector.SetXYZ(0,dOffset,0);
   vector.RotateZ(-dRotAngleRad);
      
   gMC->Gspos("RICH",3,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1002], "ONLY");           
   Chamber(2).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
// Chamber 3
                   AliMatrix(idrotm[1003],      90,  dBeta-dRotAngle     , 90.           , 90+dBeta-dRotAngle  , 0             ,   0           );
   pRotMatrix=new TRotMatrix("rot996","rot996", 90,  dBeta-dRotAngle     , 90.           , 90+dBeta-dRotAngle  , 0             ,   0           );
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(dBetaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("RICH",4,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1003], "ONLY");           
   Chamber(3).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
// Chamber 4   
                   AliMatrix(idrotm[1004],      90,  360-dBeta-dRotAngle , 108.2         , 90-dBeta-dRotAngle  , 18.2          ,  90-dBeta     );
   pRotMatrix=new TRotMatrix("rot997","rot997", 90,  360-dBeta-dRotAngle , 108.2         , 90-dBeta-dRotAngle  , 18.2          ,  90-dBeta     );
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(-dBetaRad); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("RICH",5,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1004], "ONLY");
   Chamber(4).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
// Chamber 5   
                   AliMatrix(idrotm[1005],      90, -dRotAngle           , 90+dAlpha     , 90-dRotAngle        , dAlpha        ,  90           );     
   pRotMatrix=new TRotMatrix("rot998","rot998", 90, -dRotAngle           , 90+dAlpha     , 90-dRotAngle        , dAlpha        ,  90           );     
   
   vector.SetXYZ(0,dOffset,0); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
      
   gMC->Gspos("RICH",6,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1005], "ONLY");           
   Chamber(5).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
// Chamber 6           
                   AliMatrix(idrotm[1006],      90,  dBeta-dRotAngle     , 108.2         , 90+dBeta-dRotAngle  , 18.2          ,  90+dBeta     );    
   pRotMatrix=new TRotMatrix("rot999","rot999", 90,  dBeta-dRotAngle     , 108.2         , 90+dBeta-dRotAngle  , 18.2          ,  90+dBeta     );    
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(dBetaRad); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
      
   gMC->Gspos("RICH",7,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1006], "ONLY");
   Chamber(6).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
      
}//void AliRICHv3::CreateGeometry()




void AliRICHv3::Init()
{
// Makes nothing for a while   
   if(IsDebugStart()) cout<<ClassName()<<"::Init()>\n";
    
}


void AliRICHv3::BuildGeometry()    
{                                      
// Provides geometry structure for event display (ROOT TNode tree)
   
   if(IsDebugStart()) cout<<ClassName()<<"::BuildGeometry()>\n";
  
    TNode *node, *subnode, *top;
    
    const int kColorRICH = kRed;
    //
    top=gAlice->GetGeometry()->GetNode("alice");

    AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    AliRICHChamber*       iChamber;
    AliRICHGeometry*  geometry;
 
    iChamber = &(pRICH->Chamber(0));
    AliRICHSegmentationV1* segmentation=(AliRICHSegmentationV1*) iChamber->GetSegmentationModel(0);
    geometry=iChamber->GetGeometryModel();
    
    new TBRIK("S_RICH","S_RICH","void",71.09999,11.5,73.15);

    Float_t padplane_width = segmentation->GetPadPlaneWidth();
    Float_t padplane_length = segmentation->GetPadPlaneLength();


    new TBRIK("PHOTO","PHOTO","void", padplane_width/2,.1,padplane_length/2);

// Chamber 0             
    top->cd();
    node = new TNode("RICH1","RICH1","S_RICH",Chamber(0).GetX(),Chamber(0).GetY(),Chamber(0).GetZ(),"rot993");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);

// Chamber 1
    top->cd(); 
    node = new TNode("RICH2","RICH2","S_RICH",Chamber(1).GetX(),Chamber(1).GetY(),Chamber(1).GetZ(),"rot994");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);

// Chamber 2
    top->cd();
    node = new TNode("RICH3","RICH3","S_RICH",Chamber(2).GetX(),Chamber(2).GetY(),Chamber(2).GetZ(),"rot995");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);
    
// Chamber 3
    top->cd();
    node = new TNode("RICH4","RICH4","S_RICH",Chamber(3).GetX(),Chamber(3).GetY(),Chamber(3).GetZ(),"rot996");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);

// Chamber 4
    top->cd();
    node = new TNode("RICH5","RICH5","S_RICH",Chamber(4).GetX(),Chamber(4).GetY(),Chamber(4).GetZ(),"rot997");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);

// Chamber 5
    top->cd();
    node = new TNode("RICH6","RICH6","S_RICH",Chamber(5).GetX(),Chamber(5).GetY(),Chamber(5).GetZ(),"rot998");
    node->SetLineColor(kColorRICH);
    fNodes->Add(node);node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);

// Chamber 6
    top->cd();
    node = new TNode("RICH7","RICH7","S_RICH",Chamber(6).GetX(),Chamber(6).GetY(),Chamber(6).GetZ(),"rot999");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node); 
    
}//AliRICHv3::BuildGeometry()

void AliRICHv3::StepManager()
{
// The active Step Manager is realised currently in AliRICHv3 for debug
// leaving StepManager in AliRICH intact. To be removed in future. 

    Int_t          copy, id;
    static Int_t   idvol;
    static Int_t   vol[2];
    Int_t          ipart;
    static Float_t hits[22];
    static Float_t ckovData[19];
    TLorentzVector position;
    TLorentzVector momentum;
    Float_t        pos[3];
    Float_t        mom[4];
    Float_t        localPos[3];
    Float_t        localMom[4];
    Float_t        localTheta,localPhi;
    Float_t        theta,phi;
    Float_t        destep, step;
    Double_t       ranf[2];
    Int_t          nPads;
    Float_t        coscerenkov;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t kBig=1.e10;
       
    TClonesArray &lhits = *fHits;
    TParticle *current = (TParticle*)(*gAlice->Particles())[gAlice->CurrentTrack()];

 //if (current->Energy()>1)
   //{
        
    // Only gas gap inside chamber
    // Tag chambers and record hits when track enters 
    
    idvol=-1;
    id=gMC->CurrentVolID(copy);
    Float_t cherenkovLoss=0;
    //gAlice->KeepTrack(gAlice->CurrentTrack());
    
    gMC->TrackPosition(position);
    pos[0]=position(0);
    pos[1]=position(1);
    pos[2]=position(2);
    //bzero((char *)ckovData,sizeof(ckovData)*19);
    ckovData[1] = pos[0];                 // X-position for hit
    ckovData[2] = pos[1];                 // Y-position for hit
    ckovData[3] = pos[2];                 // Z-position for hit
    ckovData[6] = 0;                      // dummy track length
    //ckovData[11] = gAlice->CurrentTrack();
    
    //printf("\n+++++++++++\nTrack: %d\n++++++++++++\n",gAlice->CurrentTrack());

    //AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    
    /********************Store production parameters for Cerenkov photons************************/ 
//is it a Cerenkov photon? 
    if (gMC->TrackPid() == 50000050) { 

      //if (gMC->VolId("GAP ")==gMC->CurrentVolID(copy))
        //{                    
	  Float_t ckovEnergy = current->Energy();
	  //energy interval for tracking
	  if  (ckovEnergy > 5.6e-09 && ckovEnergy < 7.8e-09 )       
	    //if (ckovEnergy > 0)
	    {
	      if (gMC->IsTrackEntering()){        //is track entering?
		//printf("Track entered (1)\n");
		if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
		  {                                                          //is it in freo?
		    if (gMC->IsNewTrack()){                          //is it the first step?
		      //printf("I'm in!\n");
		      Int_t mother = current->GetFirstMother(); 
		      
		      //printf("Second Mother:%d\n",current->GetSecondMother());
		      
		      ckovData[10] = mother;
		      ckovData[11] = gAlice->CurrentTrack();
		      ckovData[12] = 1;             //Media where photon was produced 1->Freon, 2->Quarz
		      //printf("Produced in FREO\n");
		      fCkovNumber++;
		      fFreonProd=1;
		      //printf("Index: %d\n",fCkovNumber);
		    }    //first step question
		  }        //freo question
		
		if (gMC->IsNewTrack()){                                  //is it first step?
		  if (gMC->VolId("QUAR")==gMC->CurrentVolID(copy))             //is it in quarz?
		    {
		      ckovData[12] = 2;
		      //printf("Produced in QUAR\n");
		    }    //quarz question
		}        //first step question
		
		//printf("Before %d\n",fFreonProd);
	      }   //track entering question
	      
	      if (ckovData[12] == 1)                                        //was it produced in Freon?
		//if (fFreonProd == 1)
		{
		  if (gMC->IsTrackEntering()){                                     //is track entering?
		    //printf("Track entered (2)\n");
		    //printf("Current volume (should be META): %s\n",gMC->CurrentVolName());
		    //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("META"),gMC->CurrentVolID(copy));
		    if (gMC->VolId("META")==gMC->CurrentVolID(copy))                //is it in gap?      
		      {
			//printf("Got in META\n");
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);
			// Z-position for hit
			
			
			/**************** Photons lost in second grid have to be calculated by hand************/ 
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = (1. - .025 / cophi) * (1. - .05 /  cophi);
			//gMC->Rndm(ranf, 1);
			gMC->GetRandom()->RndmArray(1,ranf);
			//printf("grid calculation:%f\n",t);
			if (ranf[0] > t) {
			  gMC->StopTrack();
			  ckovData[13] = 5;
			  AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
			  //printf("Added One (1)!\n");
			  //printf("Lost one in grid\n");
			}
			/**********************************************************************************/
		      }    //gap
		    
		    //printf("Current volume (should be CSI) (1): %s\n",gMC->CurrentVolName());
		    //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("CSI "),gMC->CurrentVolID(copy));
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))             //is it in csi?      
		      {
			//printf("Got in CSI\n");
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);
			
			/********* Photons lost by Fresnel reflection have to be calculated by hand********/ 
			/***********************Cerenkov phtons (always polarised)*************************/
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = Fresnel(ckovEnergy*1e9,cophi,1);
			//gMC->Rndm(ranf, 1);
			gMC->GetRandom()->RndmArray(1,ranf);
			if (ranf[0] < t) {
			  gMC->StopTrack();
			  ckovData[13] = 6;
			  AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
			  //printf("Added One (2)!\n");
			  //printf("Lost by Fresnel\n");
			}
			/**********************************************************************************/
		      }
		  } //track entering?
		  
		  
		  /********************Evaluation of losses************************/
		  /******************still in the old fashion**********************/
		  
		  TArrayI procs;
		  Int_t i1 = gMC->StepProcesses(procs);            //number of physics mechanisms acting on the particle
		  for (Int_t i = 0; i < i1; ++i) {
		    //        Reflection loss 
		    if (procs[i] == kPLightReflection) {        //was it reflected
		      ckovData[13]=10;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			ckovData[13]=1;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			ckovData[13]=2;
		      //gMC->StopTrack();
		      //AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
		    } //reflection question
		     
		    //        Absorption loss 
		    else if (procs[i] == kPLightAbsorption) {              //was it absorbed?
		      //printf("Got in absorption\n");
		      ckovData[13]=20;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			ckovData[13]=11;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			ckovData[13]=12;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("META")) 
			ckovData[13]=13;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("GAP ")) 
			ckovData[13]=13;
		      
		      if (gMC->CurrentVolID(copy) == gMC->VolId("SRIC")) 
			ckovData[13]=15;
		      
		      //        CsI inefficiency 
		      if (gMC->CurrentVolID(copy) == gMC->VolId("CSI ")) {
			ckovData[13]=16;
		      }
		      gMC->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
		      //printf("Added One (3)!\n");
		      //printf("Added cerenkov %d\n",fCkovNumber);
		    } //absorption question 
		    
		    
		    //        Photon goes out of tracking scope 
		    else if (procs[i] == kPStop) {                 //is it below energy treshold?
		      ckovData[13]=21;
		      gMC->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
		      //printf("Added One (4)!\n");
		    }	// energy treshold question	    
		  }  //number of mechanisms cycle
		  /**********************End of evaluation************************/
		} //freon production question
	    } //energy interval question
	//}//inside the proximity gap question
    } //cerenkov photon question
      
    /**************************************End of Production Parameters Storing*********************/ 
    
    
    /*******************************Treat photons that hit the CsI (Ckovs and Feedbacks)************/ 
    
    if (gMC->TrackPid() == 50000050 || gMC->TrackPid() == 50000051) {
      //printf("Cerenkov\n");
      
      //if (gMC->TrackPid() == 50000051)
	//printf("Tracking a feedback\n");
      
      if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))
	{
	  //printf("Current volume (should be CSI) (2): %s\n",gMC->CurrentVolName());
	  //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("CSI "),gMC->CurrentVolID(copy));
	  //printf("Got in CSI\n");
	  //printf("Tracking a %d\n",gMC->TrackPid());
	  if (gMC->Edep() > 0.){
		gMC->TrackPosition(position);
		gMC->TrackMomentum(momentum);
		pos[0]=position(0);
		pos[1]=position(1);
		pos[2]=position(2);
		mom[0]=momentum(0);
		mom[1]=momentum(1);
		mom[2]=momentum(2);
		mom[3]=momentum(3);
		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		gMC->Gmtod(pos,localPos,1);                                                                    
		gMC->Gmtod(mom,localMom,2);
		
		gMC->CurrentVolOffID(2,copy);
		vol[0]=copy;
		idvol=vol[0]-1;

		//Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(localPos[0], localPos[2]);
		//printf("Sector:%d\n",sector);

		/*if (gMC->TrackPid() == 50000051){
		  fFeedbacks++;
		  printf("Feedbacks:%d\n",fFeedbacks);
		}*/	
		
        //PH		((AliRICHChamber*) (*fChambers)[idvol])
		((AliRICHChamber*)fChambers->At(idvol))
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if(idvol<kNCH) {	
		    ckovData[0] = gMC->TrackPid();        // particle type
		    ckovData[1] = pos[0];                 // X-position for hit
		    ckovData[2] = pos[1];                 // Y-position for hit
		    ckovData[3] = pos[2];                 // Z-position for hit
		    ckovData[4] = theta;                      // theta angle of incidence
		    ckovData[5] = phi;                      // phi angle of incidence 
		    ckovData[8] = (Float_t) fNSDigits;      // first sdigit
		    ckovData[9] = -1;                       // last pad hit
		    ckovData[13] = 4;                       // photon was detected
		    ckovData[14] = mom[0];
		    ckovData[15] = mom[1];
		    ckovData[16] = mom[2];
		    
		    destep = gMC->Edep();
		    gMC->SetMaxStep(kBig);
		    cherenkovLoss  += destep;
		    ckovData[7]=cherenkovLoss;
		    
		    nPads = Hits2SDigits(localPos[0],localPos[2],cherenkovLoss,idvol,kCerenkov);
		    		    
		    if (fNSDigits > (Int_t)ckovData[8]) {
			ckovData[8]= ckovData[8]+1;
			ckovData[9]= (Float_t) fNSDigits;
		    }

		    //printf("Cerenkov loss: %f\n", cherenkovLoss);

		    ckovData[17] = nPads;
		    //printf("nPads:%d",nPads);
		    
		    //TClonesArray *Hits = RICH->Hits();
		    AliRICHHit *mipHit =  (AliRICHHit*) (fHits->UncheckedAt(0));
		    if (mipHit)
		      {
			mom[0] = current->Px();
			mom[1] = current->Py();
			mom[2] = current->Pz();
			Float_t mipPx = mipHit->MomX();
			Float_t mipPy = mipHit->MomY();
			Float_t mipPz = mipHit->MomZ();
			
			Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
			Float_t rt = TMath::Sqrt(r);
			Float_t mipR = mipPx*mipPx + mipPy*mipPy + mipPz*mipPz;	
			Float_t mipRt = TMath::Sqrt(mipR);
			if ((rt*mipRt) > 0)
			  {
			    coscerenkov = (mom[0]*mipPx + mom[1]*mipPy + mom[2]*mipPz)/(rt*mipRt);
			  }
			else
			  {
			    coscerenkov = 0;
			  }
			Float_t cherenkov = TMath::ACos(coscerenkov);
			ckovData[18]=cherenkov;
		      }
		    //if (sector != -1)
		    //{
		    AddHit(gAlice->CurrentTrack(),vol,ckovData);
		    AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
		    //printf("Added One (5)!\n");
		    //}
		}
	    }
	}
    }
    
    /***********************************************End of photon hits*********************************************/
    

    /**********************************************Charged particles treatment*************************************/

    else if (gMC->TrackCharge())
    //else if (1 == 1)
      {
//If MIP
	/*if (gMC->IsTrackEntering())
	  {                
	    hits[13]=20;//is track entering?
	  }*/
	if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
	  {
	    gMC->TrackMomentum(momentum);
	    mom[0]=momentum(0);
	    mom[1]=momentum(1);
	    mom[2]=momentum(2);
	    mom[3]=momentum(3);
	    hits [19] = mom[0];
	    hits [20] = mom[1];
	    hits [21] = mom[2];
	    fFreonProd=1;
	  }

	if (gMC->VolId("GAP ")== gMC->CurrentVolID(copy)) {
// Get current particle id (ipart), track position (pos)  and momentum (mom)
	    
	    gMC->CurrentVolOffID(3,copy);
	    vol[0]=copy;
	    idvol=vol[0]-1;

	    //Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(localPos[0], localPos[2]);
	    //printf("Sector:%d\n",sector);
	    
	    gMC->TrackPosition(position);
	    gMC->TrackMomentum(momentum);
	    pos[0]=position(0);
	    pos[1]=position(1);
	    pos[2]=position(2);
	    mom[0]=momentum(0);
	    mom[1]=momentum(1);
	    mom[2]=momentum(2);
	    mom[3]=momentum(3);
	    gMC->Gmtod(pos,localPos,1);                                                                    
	    gMC->Gmtod(mom,localMom,2);
	    
	    ipart  = gMC->TrackPid();
	    //
	    // momentum loss and steplength in last step
	    destep = gMC->Edep();
	    step   = gMC->TrackStep();
  
	    //
	    // record hits when track enters ...
	    if( gMC->IsTrackEntering()) {
//		gMC->SetMaxStep(fMaxStepGas);
		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		

		Double_t localTc = localMom[0]*localMom[0]+localMom[2]*localMom[2];
		Double_t localRt = TMath::Sqrt(localTc);
		localTheta   = Float_t(TMath::ATan2(localRt,Double_t(localMom[1])))*kRaddeg;                       
		localPhi     = Float_t(TMath::ATan2(Double_t(localMom[2]),Double_t(localMom[0])))*kRaddeg;    
		
		hits[0] = Float_t(ipart);         // particle type
		hits[1] = localPos[0];                 // X-position for hit
		hits[2] = localPos[1];                 // Y-position for hit
		hits[3] = localPos[2];                 // Z-position for hit
		hits[4] = localTheta;                  // theta angle of incidence
		hits[5] = localPhi;                    // phi angle of incidence 
		hits[8] = (Float_t) fNSDigits;    // first sdigit
		hits[9] = -1;                     // last pad hit
		hits[13] = fFreonProd;           // did id hit the freon?
		hits[14] = mom[0];
		hits[15] = mom[1];
		hits[16] = mom[2];
		hits[18] = 0;               // dummy cerenkov angle

		tlength = 0;
		eloss   = 0;
		fFreonProd = 0;
	
		Chamber(idvol).LocaltoGlobal(localPos,hits+1);
	   
		
		//To make chamber coordinates x-y had to pass localPos[0], localPos[2]
		xhit    = localPos[0];
		yhit    = localPos[2];
		// Only if not trigger chamber
		if(idvol<kNCH) {
		    //
		    //  Initialize hit position (cursor) in the segmentation model 
          //PH		    ((AliRICHChamber*) (*fChambers)[idvol])
		    ((AliRICHChamber*)fChambers->At(idvol))
			->SigGenInit(localPos[0], localPos[2], localPos[1]);
		}
	    }
	    
	    // 
	    // Calculate the charge induced on a pad (disintegration) in case 
	    //
	    // Mip left chamber ...
	    if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
		gMC->SetMaxStep(kBig);
		eloss   += destep;
		tlength += step;
		
				
		// Only if not trigger chamber
		if(idvol<kNCH) {
		  if (eloss > 0) 
		    {
		      if(gMC->TrackPid() == kNeutron)
			printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		      nPads = Hits2SDigits(xhit,yhit,eloss,idvol,kMip);
		      hits[17] = nPads;
		      //printf("nPads:%d",nPads);
		    }
		}
		
		hits[6]=tlength;
		hits[7]=eloss;
		if (fNSDigits > (Int_t)hits[8]) {
		    hits[8]= hits[8]+1;
		    hits[9]= (Float_t) fNSDigits;
		}
		
		//if(sector !=-1)
		new(lhits[fNhits++]) AliRICHHit(fIshunt,gAlice->CurrentTrack(),vol,hits);
		eloss = 0; 
		//
		// Check additional signal generation conditions 
		// defined by the segmentation
		// model (boundary crossing conditions) 
	    } else if 
          //PH		(((AliRICHChamber*) (*fChambers)[idvol])
		(((AliRICHChamber*)fChambers->At(idvol))
		 ->SigGenCond(localPos[0], localPos[2], localPos[1]))
	    {
          //PH		((AliRICHChamber*) (*fChambers)[idvol])
		((AliRICHChamber*)fChambers->At(idvol))
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if (eloss > 0) 
		  {
		    if(gMC->TrackPid() == kNeutron)
		      printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		    nPads = Hits2SDigits(xhit,yhit,eloss,idvol,kMip);
		    hits[17] = nPads;
		    //printf("Npads:%d",NPads);
		  }
		xhit     = localPos[0];
		yhit     = localPos[2]; 
		eloss    = destep;
		tlength += step ;
		//
		// nothing special  happened, add up energy loss
	    } else {        
		eloss   += destep;
		tlength += step ;
	    }
	}
      }
    /*************************************************End of MIP treatment**************************************/
   //}
}// void AliRICHv3::StepManager()
