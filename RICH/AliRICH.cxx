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

#include <TArrayF.h>
#include <TBRIK.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TNode.h> 
#include <TObjArray.h>
#include <TObject.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TRandom.h> 
#include <TTUBE.h>
#include <TTree.h>
#include <TVector.h>
#include <AliMagF.h>
#include <AliPoints.h>
#include "AliRICH.h"
#include "AliRICHParam.h"
#include "AliRICHClusterFinder.h"
#include "AliRICHRawCluster.h"
#include "AliRICHRecHit1D.h"
#include "AliRICHRecHit3D.h"
#include "AliRICHSegmentationV0.h"
#include "AliRICHTransientDigit.h"
#include <AliRun.h>
#include <AliRunDigitizer.h>
#include "AliRICHSegmentationV1.h"
#include "AliRICHResponseV0.h"
 
ClassImp(AliRICHhit)
//__________________________________________________________________________________________________
void AliRICHhit::Print(Option_t*)const
{
  Info("","chamber=%2i, PID=%9i, TID=%6i, eloss=%9.3f eV",fChamber,fPid,fTrack,fEloss*1e9);
}//void AliRICHdigit::Print(Option_t *option)const
//__________________________________________________________________________________________________
ClassImp(AliRICHdigit)
//__________________________________________________________________________________________________
void AliRICHdigit::Print(Option_t*)const
{
  Info("","ID=%6i, chamber=%2i, PadX=%3i, PadY=%3i, Qdc=%4i, TID1=%5i, TID2=%5i, TID3=%5i",
         Id(),fChamber,fPadX,fPadY,fQdc,fTracks[0],fTracks[1],fTracks[2]);
}//void AliRICHdigit::Print(Option_t *option)const
//__________________________________________________________________________________________________

ClassImp(AliRICH)    
//__________________________________________________________________________________________________
// RICH manager class   
//BEGIN_HTML
/*
  <img src="gif/alirich.gif">
*/
//END_HTML
//__________________________________________________________________________________________________
AliRICH::AliRICH()
        :AliDetector() 
{//Default ctor should not contain any new operators
  fpParam     =0;
  fChambers   =0;
  
  
        //AliDetector ctor deals with Hits and Digits
  
  fSdigits    =0; fNsdigits   =0;
  fDigitsNew  =0; for(int i=0;i<kNCH;i++) fNdigitsNew[i]  =0;
  fClusters   =0; for(int i=0;i<kNCH;i++) fNclusters[i]=0;
  
  fCerenkovs  =0; fNcerenkovs =0;
  fSpecials   =0; fNspecials  =0;  
  fDchambers  =0; for(int i=0;i<kNCH;i++) fNdch[i]=0;
  fRecHits1D  =0; for(int i=0;i<kNCH;i++) fNrechits1D[i]=0;
  fRecHits3D  =0; for(int i=0;i<kNCH;i++) fNrechits3D[i]=0;
  fRawClusters=0; for(int i=0;i<kNCH;i++) fNrawch[i]=0;  
}//AliRICH::AliRICH()
//__________________________________________________________________________________________________
AliRICH::AliRICH(const char *name, const char *title)
        :AliDetector(name,title)
{//Named ctor
  if(GetDebug())Info("named ctor","Start.");
  fpParam     =   new AliRICHParam;
  fChambers = 0;  CreateChambers();
  
        //AliDetector ctor deals with Hits and Digits (reset them to 0, does not create them)
  fHits=       0;     CreateHits();          gAlice->AddHitList(fHits);
  fSdigits=    0;
  fDigitsNew=  0;
  fClusters=   0;
  
  fCerenkovs=  0;     CreateCerenkovsOld();  gAlice->AddHitList(fCerenkovs);
  fSpecials=   0;     CreateSpecialsOld();   
  fDchambers=  0;   //CreateDigitsOld();
  fRawClusters=0;   //CreateRawClustersOld();
  fRecHits1D=  0;   //CreateRecos1Old();
  fRecHits3D=  0;   //CreateRecos3Old();
  
  fCkovNumber=fFreonProd=0;  
  if(GetDebug())Info("named ctor","Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//__________________________________________________________________________________________________
AliRICH::~AliRICH()
{//dtor
  if(GetDebug()) Info("dtor","Start.");

  if(fpParam)    delete fpParam;
  if(fChambers)  delete fChambers;
  
  if(fHits)      delete fHits;
  if(fSdigits)   delete fSdigits;
  if(fDigits)    delete fDigits;
  if(fDigitsNew) {fDigitsNew->Delete();   delete fDigitsNew;}
  if(fClusters)  {fClusters->Delete();    delete fClusters;}
  
  if(fCerenkovs) delete fCerenkovs;
  if(fSpecials)  delete fSpecials;
  if(fDchambers)   {fDchambers->Delete();     delete fDchambers;}
  if(fRawClusters) {fRawClusters->Delete();   delete fRawClusters;}          
  if(fRecHits1D) {fRecHits1D->Delete();       delete fRecHits1D;}
  if(fRecHits3D) {fRecHits3D->Delete();       delete fRecHits3D;}                     
  if(GetDebug()) Info("dtor","Stop.");    
}//AliRICH::~AliRICH()
//__________________________________________________________________________________________________
void AliRICH::Hits2SDigits()
{//Create a list of sdigits corresponding to list of hits. Every hit generates one or more sdigits.
  if(GetDebug()) Info("Hit2SDigits","Start.");
  
  if(GetDebug()) Info("Hit2SDigits","Stop.");
}//void AliRICH::Hits2SDigits()
//__________________________________________________________________________________________________
void AliRICH::SDigits2Digits()
{//Generate digits from sdigits.
  if(GetDebug()) Info("SDigits2Digits","Start.");
    
  if(GetDebug()) Info("SDigits2Digits","Stop.");
}//void AliRICH::SDigits2Digits()
//__________________________________________________________________________________________________
void AliRICH::Digits2Reco()
{//Generate clusters from digits then generate recos from clusters or digits
  if(GetDebug()) Info("Digits2reco","Start.");

}//void AliRICH::Digits2Reco()  
//__________________________________________________________________________________________________
void AliRICH::AddRawCluster(Int_t id, const AliRICHRawCluster& c)
{// Add a RICH digit to the list
   
    TClonesArray &lrawcl = *((TClonesArray*)fRawClusters->At(id));
    new(lrawcl[fNrawch[id]++]) AliRICHRawCluster(c);
}
//_____________________________________________________________________________
void AliRICH::AddRecHit1D(Int_t id, Float_t *rechit, Float_t *photons, Int_t *padsx, Int_t* padsy)
{// Add a RICH reconstructed hit to the list

    TClonesArray &lrec1D = *((TClonesArray*)fRecHits1D->At(id));
    new(lrec1D[fNrechits1D[id]++]) AliRICHRecHit1D(id,rechit,photons,padsx,padsy);
}
//_____________________________________________________________________________
void AliRICH::AddRecHit3D(Int_t id, Float_t *rechit, Float_t omega, Float_t theta, Float_t phi)
{// Add a RICH reconstructed hit to the list

    TClonesArray &lrec3D = *((TClonesArray*)fRecHits3D->At(id));
    new(lrec3D[fNrechits3D[id]++]) AliRICHRecHit3D(id,rechit,omega,theta,phi);
}
//______________________________________________________________________________
void AliRICH::BuildGeometry() 
{//Builds a TNode geometry for event display
  if(GetDebug())Info("BuildGeometry","Start.");
  
  TNode *node, *subnode, *top;
  top=gAlice->GetGeometry()->GetNode("alice");
  
  new TBRIK("S_RICH","S_RICH","void",71.09999,11.5,73.15);

  Float_t wid=Param()->SectorSizeX();
  Float_t len=Param()->SectorSizeY();
  new TBRIK("PHOTO","PHOTO","void",wid/2,0.1,len/2);
  
  for(int i=1;i<=kNCH;i++){
    top->cd();
    node = new TNode(Form("RICH%i",i),Form("RICH%i",i),"S_RICH",C(i)->X(),C(i)->Y(),C(i)->Z(),C(i)->RotMatrixName());
    node->SetLineColor(kRed);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",wid+Param()->DeadZone(),5,len/2+Param()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,len/2+Param()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-wid-Param()->DeadZone(),5,len/2+Param()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",wid+Param()->DeadZone(),5,-len/2-Param()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-len/2 -Param()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-wid-Param()->DeadZone(),5,-len/2 - Param()->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);
  }  
  if(GetDebug())Info("BuildGeometry","Stop.");    
}//void AliRICH::BuildGeometry()

static Int_t kCSI=6;
static Int_t kGAP=9;
//______________________________________________________________________________
void AliRICH::CreateMaterials()
{
    //
    // *** DEFINITION OF AVAILABLE RICH MATERIALS *** 
    // ORIGIN    : NICK VAN EIJNDHOVEN 
    // Modified by:  N. Colonna (INFN - BARI, Nicola.Colonna@ba.infn.it) 
    //               R.A. Fini  (INFN - BARI, Rosanna.Fini@ba.infn.it) 
    //               R.A. Loconsole (Bari University, loco@riscom.ba.infn.it) 
    //
    Int_t i;
    
    //Photons energy intervals
    Float_t ppckov[26];
    for (i=0;i<26;i++) 
    {
	ppckov[i] = (Float_t(i)*0.1+5.5)*1e-9;
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
	
    /*******************************************End of rich_media.f***************************************/
    
  Float_t rIndexFreon[26];
    
    
    // --- Photon energy (GeV) 
    // --- Refraction indexes 
    for (i = 0; i < 26; ++i) {
      rIndexFreon[i] = ppckov[i] * .0172 * 1e9 + 1.177;
      //rIndexFreon[i] = 1;
    }
            
      
  Float_t a=0,z=0,den=0,radl=0,absl=0;
  Float_t tmaxfd=-10.0, deemax=-0.2, stemax=-0.1,epsil=0.001, stmin=-0.001; 
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
    
  AliMaterial( 1, "Air     $",a=14.61,z=7.3, den=0.001205,radl=30420.0,absl=67500);//(Air)
  AliMedium(1, "DEFAULT MEDIUM AIR$", 1, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial( 6, "HON",      a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);    //(C)-equivalent radl
  AliMedium(2, "HONEYCOMB$", 6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(16, "CSI",      a=12.01,z=6.0, den=0.1,     radl=18.8,   absl=0);    //CsI-radl equivalent
  AliMedium(kCSI, "CSI$", 16, 1, isxfld, sxmgmx,tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(11, "GRI",      a=63.54,z=29.0,den=8.96,    radl=1.43,   absl=0);    //anode grid (Cu) 
  AliMedium(7, "GRIGLIA$", 11, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(50, "ALUM",     a=26.98,z=13.0,den=2.7,     radl=8.9,    absl=0);    //aluminium sheet (Al)
  AliMedium(10, "ALUMINUM$", 50, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(31, "COPPER$",  a=63.54,z=29.0,den=8.96,    radl=1.4,    absl=0);    //(Cu)
  AliMedium(12, "PCB_COPPER", 31, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aQuartz[2]={28.09,16.0};  Float_t  zQuartz[2]={14.00, 8.0};  Float_t  wmatQuartz[2]={1,2};
  AliMixture (20, "QUA",aQuartz,zQuartz,den=2.64,-2, wmatQuartz);//Quarz (SiO2) - trasnparent 
  AliMedium(3, "QUARZO$", 20, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture (21, "QUAO",aQuartz, zQuartz, den=2.64, -2, wmatQuartz);//Quarz (SiO2) - opaque
  AliMedium(8, "QUARZOO$", 21, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t  aFreon[2]={12,19};  Float_t  zFreon[2]={6,9};  Float_t wmatFreon[2]={6,14};
  AliMixture (30, "FRE",aFreon,zFreon,den=1.7,-2,wmatFreon);//Freon (C6F14) 
  AliMedium(4, "FREON$", 30, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  Float_t aMethane[2]={12.01,1}; Float_t zMethane[2]={6,1}; Float_t wmatMethane[2]={1,4};
  AliMixture (40, "MET", aMethane, zMethane, den=7.17e-4,-2, wmatMethane);//methane (CH4)     
  AliMedium(5, "METANO$", 40, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMixture (41, "METG", aMethane, zMethane, den=7.17e-4, -2, wmatMethane);
  AliMedium(kGAP, "GAP$", 41, 1, isxfld, sxmgmx,tmaxfd, 0.1, -deemax, epsil, -stmin);
  
  Float_t aGlass[5]={12.01, 28.09, 16.,   10.8,  23.};
  Float_t zGlass[5]={ 6.,   14.,    8.,    5.,   11.};
  Float_t wGlass[5]={ 0.5,  0.105, 0.355, 0.03,  0.01};
  AliMixture (32, "GLASS",aGlass, zGlass, den=1.74, 5, wGlass);//Glass 50%C+10.5%Si+35.5%O+3% + 1%
  AliMedium(11, "GLASS", 32, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
            
  Int_t *idtmed = fIdtmed->GetArray()-999;
  gMC->SetCerenkov(idtmed[1000], 26, ppckov, abscoMethane,     efficAll,  rIndexMethane);
  gMC->SetCerenkov(idtmed[1001], 26, ppckov, abscoMethane,     efficAll,  rIndexMethane);
  gMC->SetCerenkov(idtmed[1002], 26, ppckov, abscoQuarz,       efficAll,  rIndexQuarz);
  gMC->SetCerenkov(idtmed[1003], 26, ppckov, abscoFreon,       efficAll,  rIndexFreon);
  gMC->SetCerenkov(idtmed[1004], 26, ppckov, abscoMethane,     efficAll,  rIndexMethane);
  gMC->SetCerenkov(idtmed[1005], 26, ppckov, abscoCsI,         efficCsI,  rIndexMethane);
  gMC->SetCerenkov(idtmed[1006], 26, ppckov, abscoGrid,        efficGrid, rIndexGrid);
  gMC->SetCerenkov(idtmed[1007], 26, ppckov, abscoOpaqueQuarz, efficAll,  rIndexOpaqueQuarz);
  gMC->SetCerenkov(idtmed[1008], 26, ppckov, abscoMethane,     efficAll,  rIndexMethane);
  gMC->SetCerenkov(idtmed[1009], 26, ppckov, abscoGrid,        efficGrid, rIndexGrid);
  gMC->SetCerenkov(idtmed[1010], 26, ppckov, abscoOpaqueQuarz, efficAll,  rIndexOpaqueQuarz);
}//void AliRICH::CreateMaterials()
//__________________________________________________________________________________________________
Float_t AliRICH::Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)
{

    //ENE(EV), PDOTI=COS(INC.ANG.), PDOTR=COS(POL.PLANE ROT.ANG.)
    
    Float_t en[36] = {5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,
		      6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,
		      7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5};
    Float_t csin[36] = {2.14,2.21,2.33,2.48,2.76,2.97,2.99,2.59,2.81,3.05,
			2.86,2.53,2.55,2.66,2.79,2.96,3.18,3.05,2.84,2.81,2.38,2.11,
			2.01,2.13,2.39,2.73,3.08,3.15,2.95,2.73,2.56,2.41,2.12,1.95,
			1.72,1.53};
    Float_t csik[36] = {0.,0.,0.,0.,0.,0.196,0.408,0.208,0.118,0.49,0.784,0.543,
	 		0.424,0.404,0.371,0.514,0.922,1.102,1.139,1.376,1.461,1.253,0.878,
			0.69,0.612,0.649,0.824,1.347,1.571,1.678,1.763,1.857,1.824,1.824,
			1.714,1.498};
    Float_t xe=ene;
    Int_t  j=Int_t(xe*10)-49;
    Float_t cn=csin[j]+((csin[j+1]-csin[j])/0.1)*(xe-en[j]);
    Float_t ck=csik[j]+((csik[j+1]-csik[j])/0.1)*(xe-en[j]);

    //FORMULAE FROM HANDBOOK OF OPTICS, 33.23 OR
    //W.R. HUNTER, J.O.S.A. 54 (1964),15 , J.O.S.A. 55(1965),1197

    Float_t sinin=TMath::Sqrt(1-pdoti*pdoti);
    Float_t tanin=sinin/pdoti;

    Float_t c1=cn*cn-ck*ck-sinin*sinin;
    Float_t c2=4*cn*cn*ck*ck;
    Float_t aO=TMath::Sqrt(0.5*(TMath::Sqrt(c1*c1+c2)+c1));
    Float_t b2=0.5*(TMath::Sqrt(c1*c1+c2)-c1);
    
    Float_t rs=((aO-pdoti)*(aO-pdoti)+b2)/((aO+pdoti)*(aO+pdoti)+b2);
    Float_t rp=rs*((aO-sinin*tanin)*(aO-sinin*tanin)+b2)/((aO+sinin*tanin)*(aO+sinin*tanin)+b2);
    

    //CORRECTION FACTOR FOR SURFACE ROUGHNESS
    //B.J. STAGG  APPLIED OPTICS, 30(1991),4113

    Float_t sigraf=18.;
    Float_t lamb=1240/ene;
    Float_t fresn;
 
    Float_t  rO=TMath::Exp(-(4*TMath::Pi()*pdoti*sigraf/lamb)*(4*TMath::Pi()*pdoti*sigraf/lamb));

    if(pola)
    {
	Float_t pdotr=0.8;                                 //DEGREE OF POLARIZATION : 1->P , -1->S
	fresn=0.5*(rp*(1+pdotr)+rs*(1-pdotr));
    }
    else
	fresn=0.5*(rp+rs);
      
    fresn = fresn*rO;
    return(fresn);
}//Float_t AliRICH::Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)
//__________________________________________________________________________________________________
Float_t AliRICH::AbsoCH4(Float_t x)
{

    //KLOSCH,SCH4(9),WL(9),EM(9),ALENGTH(31)
    Float_t sch4[9] = {.12,.16,.23,.38,.86,2.8,7.9,28.,80.};              //MB X 10^22
    //Float_t wl[9] = {153.,152.,151.,150.,149.,148.,147.,146.,145};
    Float_t em[9] = {8.1,8.158,8.212,8.267,8.322,8.378,8.435,8.493,8.55};
    const Float_t kLosch=2.686763E19;                                      // LOSCHMIDT NUMBER IN CM-3
    const Float_t kIgas1=100, kIgas2=0, kOxy=10., kWater=5., kPressure=750.,kTemperature=283.;                                      
    Float_t pn=kPressure/760.;
    Float_t tn=kTemperature/273.16;
    
	
// ------- METHANE CROSS SECTION -----------------
// ASTROPH. J. 214, L47 (1978)
	
    Float_t sm=0;
    if (x<7.75) 
	sm=.06e-22;
    
    if(x>=7.75 && x<=8.1)
    {
	Float_t c0=-1.655279e-1;
	Float_t c1=6.307392e-2;
	Float_t c2=-8.011441e-3;
	Float_t c3=3.392126e-4;
	sm=(c0+c1*x+c2*x*x+c3*x*x*x)*1.e-18;
    }
    
    if (x> 8.1)
    {
	Int_t j=0;
	while (x<=em[j] && x>=em[j+1])
	{
	    j++;
	    Float_t a=(sch4[j+1]-sch4[j])/(em[j+1]-em[j]);
	    sm=(sch4[j]+a*(x-em[j]))*1e-22;
	}
    }
    
    Float_t dm=(kIgas1/100.)*(1.-((kOxy+kWater)/1.e6))*kLosch*pn/tn;
    Float_t abslm=1./sm/dm;
    
//    ------- ISOBUTHANE CROSS SECTION --------------
//     i-C4H10 (ai) abs. length from curves in
//     Lu-McDonald paper for BARI RICH workshop .
//     -----------------------------------------------------------
    
    Float_t ai;
    Float_t absli;
    if (kIgas2 != 0) 
    {
	if (x<7.25)
	    ai=100000000.;
	
	if(x>=7.25 && x<7.375)
	    ai=24.3;
	
	if(x>=7.375)
	    ai=.0000000001;
	
	Float_t si = 1./(ai*kLosch*273.16/293.);                    // ISOB. CRO.SEC.IN CM2
	Float_t di=(kIgas2/100.)*(1.-((kOxy+kWater)/1.e6))*kLosch*pn/tn;
	absli =1./si/di;
    }
    else
	absli=1.e18;
//    ---------------------------------------------------------
//
//       transmission of O2
//
//       y= path in cm, x=energy in eV
//       so= cross section for UV absorption in cm2
//       do= O2 molecular density in cm-3
//    ---------------------------------------------------------
    
    Float_t abslo;
    Float_t so=0;
    if(x>=6.0)
    {
	if(x>=6.0 && x<6.5)
	{
	    so=3.392709e-13 * TMath::Exp(2.864104 *x);
	    so=so*1e-18;
	}
	
	if(x>=6.5 && x<7.0) 
	{
	    so=2.910039e-34 * TMath::Exp(10.3337*x);
	    so=so*1e-18;
	}
	    

	if (x>=7.0) 
	{
	    Float_t a0=-73770.76;
	    Float_t a1=46190.69;
	    Float_t a2=-11475.44;
	    Float_t a3=1412.611;
	    Float_t a4=-86.07027;
	    Float_t a5=2.074234;
	    so= a0+(a1*x)+(a2*x*x)+(a3*x*x*x)+(a4*x*x*x*x)+(a5*x*x*x*x*x);
	    so=so*1e-18;
	}
	
	Float_t dox=(kOxy/1e6)*kLosch*pn/tn;
	abslo=1./so/dox;
    }
    else
	abslo=1.e18;
//     ---------------------------------------------------------
//
//       transmission of H2O
//
//       y= path in cm, x=energy in eV
//       sw= cross section for UV absorption in cm2
//       dw= H2O molecular density in cm-3
//     ---------------------------------------------------------
    
    Float_t abslw;
    
    Float_t b0=29231.65;
    Float_t b1=-15807.74;
    Float_t b2=3192.926;
    Float_t b3=-285.4809;
    Float_t b4=9.533944;
    
    if(x>6.75)
    {    
	Float_t sw= b0+(b1*x)+(b2*x*x)+(b3*x*x*x)+(b4*x*x*x*x);
	sw=sw*1e-18;
	Float_t dw=(kWater/1e6)*kLosch*pn/tn;
	abslw=1./sw/dw;
    }
    else
    	abslw=1.e18;
	    
//    ---------------------------------------------------------
    
    Float_t alength=1./(1./abslm+1./absli+1./abslo+1./abslw);
    return (alength);
}
//__________________________________________________________________________________________________
void AliRICH::MakeBranch(Option_t* option)
{//Create Tree branches for the RICH.
  if(GetDebug())Info("MakeBranch","Start with option= %s.",option);
    
  const Int_t kBufferSize = 4000;
      
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");

  if(cH&&TreeH()){//H
    CreateHits();      //branch will be created in AliDetector::MakeBranch
    CreateCerenkovsOld(); MakeBranchInTree(TreeH(),"RICHCerenkov", &fCerenkovs, kBufferSize,0);
    CreateSpecialsOld();  MakeBranchInTree(TreeH(),"RICHSpecials", &fSpecials,kBufferSize,0); 
  }//H     
  AliDetector::MakeBranch(option);//this is after cH because we need to guarantee that fHits array is created
      
  if(cS&&fLoader->TreeS()){//S  
    CreateSdigits();   MakeBranchInTree(fLoader->TreeS(),"RICH",&fSdigits,kBufferSize,0) ;
  }//S
   
  if(cD&&fLoader->TreeD()){//D
    CreateDigitsOld();  
    for(int i=0;i<kNCH;i++) 
      MakeBranchInTree(fLoader->TreeD(),Form("%sDigits%d",GetName(),i+1),&((*fDchambers)[i]),kBufferSize,0);
  }//D
  
  if(cR&&fLoader->TreeR()){//R
    CreateRawClustersOld(); 
    for(int i=0; i<kNCH ;i++)
      MakeBranchInTree(fLoader->TreeR(),Form("%sRawClusters%d",GetName(),i+1), &((*fRawClusters)[i]), kBufferSize, 0);

    CreateRecos1Old();   
    for(int i=0; i<kNCH ;i++) 
      MakeBranchInTree(fLoader->TreeR(),Form("%sRecHits1D%d",GetName(),i+1),&((*fRecHits1D)[i]),kBufferSize,0);
    
    CreateRecos3Old();   
    for(int i=0; i<kNCH ;i++)
      MakeBranchInTree(fLoader->TreeR(),Form("%sRecHits3D%d",GetName(),i+1), &((*fRecHits3D)[i]), kBufferSize, 0);
   }//R
  if(GetDebug())Info("MakeBranch","Stop.");   
}//void AliRICH::MakeBranch(Option_t* option)
//__________________________________________________________________________________________________
void AliRICH::SetTreeAddress()
{//Set branch address for the Hits and Digits Tree.
  if(GetDebug())Info("SetTreeAddress","Start.");
      
  TBranch *branch;
    
  if(fLoader->TreeH()){//H
    if(GetDebug())Info("SetTreeAddress","tree H is requested.");
    CreateHits();//branch map will be in AliDetector::SetTreeAddress    
    branch=fLoader->TreeH()->GetBranch("RICHCerenkov");   if(branch){CreateCerenkovsOld(); branch->SetAddress(&fCerenkovs);}       
    branch=fLoader->TreeH()->GetBranch("RICHSpecials");   if(branch){CreateSpecialsOld();  branch->SetAddress(&fSpecials);}
  }//H
  AliDetector::SetTreeAddress();//this is after TreeH because we need to guarantee that fHits array is created

  if(fLoader->TreeS()){//S
    if(GetDebug())Info("SetTreeAddress","tree S is requested.");
    branch=fLoader->TreeS()->GetBranch(GetName());        if(branch){CreateSdigits();   branch->SetAddress(&fSdigits);}
  }//S
    
  if(fLoader->TreeD()){//D    
    if(GetDebug())Info("SetTreeAddress","tree D is requested.");
    for(int i=0;i<kNCH;i++){      
      branch=fLoader->TreeD()->GetBranch(Form("%s%d",GetName(),i+1)); 
      if(branch){CreateDigits(); branch->SetAddress(&((*fDigitsNew)[i]));}
      
      branch=fLoader->TreeD()->GetBranch(Form("%sDigits%d",GetName(),i+1)); 
      if(branch){CreateDigitsOld(); branch->SetAddress(&((*fDchambers)[i]));}
    }//for
  }//D
    
  if(fLoader->TreeR()){//R
    if(GetDebug())Info("SetTreeAddress","tree R is requested.");

    for(int i=0;i<kNCH;i++){         
      branch=fLoader->TreeR()->GetBranch(Form("%sClusters%d" ,GetName(),i+1));
      if(branch){CreateClusters(); branch->SetAddress(&((*fRawClusters)[i]));}
    }
    
    for(int i=0;i<kNCH;i++) {
      branch=fLoader->TreeR()->GetBranch(Form("%sRawClusters%d" ,GetName(),i+1));
      if(branch){CreateRawClustersOld(); branch->SetAddress(&((*fRawClusters)[i]));}
      
      branch=fLoader->TreeR()->GetBranch(Form("%sRecHits1D%d",GetName(),i+1));
      if(branch){CreateRecos1Old(); branch->SetAddress(&((*fRecHits1D)[i]));}
      
      branch=fLoader->TreeR()->GetBranch(Form("%sRecHits3D%d",GetName(),i+1));
      if(branch){CreateRecos3Old();branch->SetAddress(&((*fRecHits3D)[i]));}
    }
  }//R
  if(GetDebug())Info("SetTreeAddress","Stop.");
}//void AliRICH::SetTreeAddress()
//__________________________________________________________________________________________________
void AliRICH::Print(Option_t *option)const
{
  TObject::Print(option);
  Param()->Dump();
  fChambers->Print(option);  
}//void AliRICH::Print(Option_t *option)const
//__________________________________________________________________________________________________
void AliRICH::CreateGeometry()
{//Creates detailed geometry simulation (currently GEANT volumes tree)         
  if(GetDebug())Info("CreateGeometry","Start.");
//???????? to be removed to AliRICHParam?
  Param()->RadiatorToPads(Param()->FreonThickness()/2+Param()->QuartzThickness()+Param()->GapThickness());
    
//Opaque quartz thickness
  Float_t oqua_thickness = .5;
//CsI dimensions
  Float_t pcX=Param()->PcSizeX();
  Float_t pcY=Param()->PcSizeY();
  
  Int_t *idtmed = fIdtmed->GetArray()-999;
    
  Int_t i;
  Float_t zs;
  Int_t idrotm[1099];
  Float_t par[3];
    
//External aluminium box 
  par[0]=68.8*kcm;par[1]=13*kcm;par[2]=70.86*kcm;  gMC->Gsvolu("RICH", "BOX ", idtmed[1009], par, 3);
//Air 
  par[0]=66.3;   par[1] = 13; par[2] = 68.35;      gMC->Gsvolu("SRIC", "BOX ", idtmed[1000], par, 3); 
//Air 2 (cutting the lower part of the box)
  par[0]=1.25;    par[1] = 3;    par[2] = 70.86;   gMC->Gsvolu("AIR2", "BOX ", idtmed[1000], par, 3);
//Air 3 (cutting the lower part of the box)
  par[0]=66.3;    par[1] = 3;  par[2] = 1.2505;    gMC->Gsvolu("AIR3", "BOX ", idtmed[1000], par, 3);
//Honeycomb 
  par[0]=66.3;par[1]=0.188;  par[2] = 68.35;       gMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);
//Aluminium sheet 
  par[0]=66.3;par[1]=0.025;par[2]=68.35;           gMC->Gsvolu("ALUM", "BOX ", idtmed[1009], par, 3);
  //par[0] = 66.5; par[1] = .025; par[2] = 63.1;
//Quartz 
  par[0]=Param()->QuartzWidth()/2;par[1]=Param()->QuartzThickness()/2;par[2]=Param()->QuartzLength()/2;
  gMC->Gsvolu("QUAR", "BOX ", idtmed[1002], par, 3);
//Spacers (cylinders) 
  par[0]=0.;par[1]=.5;par[2]=Param()->FreonThickness()/2;  gMC->Gsvolu("SPAC", "TUBE", idtmed[1002], par, 3);    
//Feet (freon slabs supports)
  par[0] = .7;  par[1] = .3;  par[2] = 1.9;        gMC->Gsvolu("FOOT", "BOX", idtmed[1009], par, 3);
//Opaque quartz 
  par[0]=Param()->QuartzWidth()/2;par[1]= .2;par[2]=Param()->QuartzLength()/2;
  gMC->Gsvolu("OQUA", "BOX ", idtmed[1007], par, 3);
//Frame of opaque quartz
  par[0]=Param()->OuterFreonWidth()/2;par[1]=Param()->FreonThickness()/2;par[2]=Param()->OuterFreonLength()/2; 
  gMC->Gsvolu("OQF1", "BOX ", idtmed[1007], par, 3);
  par[0]=Param()->InnerFreonWidth()/2;par[1]=Param()->FreonThickness()/2;par[2]=Param()->InnerFreonLength()/2; 
  gMC->Gsvolu("OQF2", "BOX ", idtmed[1007], par, 3);
//Freon 
  par[0]=Param()->OuterFreonWidth()/2 - oqua_thickness;
  par[1]=Param()->FreonThickness()/2;
  par[2]=Param()->OuterFreonLength()/2 - 2*oqua_thickness; 
  gMC->Gsvolu("FRE1", "BOX ", idtmed[1003], par, 3);

  par[0]=Param()->InnerFreonWidth()/2 - oqua_thickness;
  par[1]=Param()->FreonThickness()/2;
  par[2]=Param()->InnerFreonLength()/2 - 2*oqua_thickness; 
  gMC->Gsvolu("FRE2", "BOX ", idtmed[1003], par, 3);    
//Methane 
  par[0]=pcX/2;par[1]=Param()->GapThickness()/2;par[2]=pcY/2;         gMC->Gsvolu("META","BOX ",idtmed[1004], par, 3);
//Methane gap 
  par[0]=pcX/2;par[1]=Param()->ProximityGapThickness()/2;par[2]=pcY/2;gMC->Gsvolu("GAP ","BOX ",(*fIdtmed)[kGAP],par,3);
//CsI PC
  par[0]=pcX/2;par[1]=.25;par[2]=pcY/2;  gMC->Gsvolu("CSI ", "BOX ", (*fIdtmed)[kCSI], par, 3);
//Anode grid 
  par[0] = 0.;par[1] = .001;par[2] = 20.;  gMC->Gsvolu("GRID", "TUBE", idtmed[1006], par, 3);

//Wire supports
//Bar of metal
  par[0]=pcX/2;par[1]=1.05;par[2]=1.05;  gMC->Gsvolu("WSMe", "BOX ", idtmed[1009], par, 3);
//Ceramic pick up (base)
  par[0]=pcX/2;par[1]= .25;par[2]=1.05;  gMC->Gsvolu("WSG1", "BOX ", idtmed[1010], par, 3);
//Ceramic pick up (head)
  par[0] = pcX/2;par[1] = .1;par[2] = .1;  gMC->Gsvolu("WSG2", "BOX ", idtmed[1010], par, 3);

//Aluminium supports for methane and CsI
//Short bar
  par[0]=pcX/2;par[1]=Param()->GapThickness()/2 + .25; par[2] = (68.35 - pcY/2)/2;
  gMC->Gsvolu("SMSH", "BOX", idtmed[1009], par, 3);
//Long bar
  par[0]=(66.3 - pcX/2)/2;par[1]=Param()->GapThickness()/2+.25;par[2]=pcY/2+68.35-pcY/2;
  gMC->Gsvolu("SMLG", "BOX", idtmed[1009], par, 3);
    
//Aluminium supports for freon
//Short bar
  par[0] = Param()->QuartzWidth()/2; par[1] = .3; par[2] = (68.35 - Param()->QuartzLength()/2)/2;
  gMC->Gsvolu("SFSH", "BOX", idtmed[1009], par, 3);    
//Long bar
  par[0] = (66.3 - Param()->QuartzWidth()/2)/2; par[1] = .3;
  par[2] = Param()->QuartzLength()/2 + 68.35 - Param()->QuartzLength()/2;
  gMC->Gsvolu("SFLG", "BOX", idtmed[1009], par, 3);    
//PCB backplane
  par[0] = pcX/2;par[1] = .25; par[2] = pcY/4 -.5025;  gMC->Gsvolu("PCB ", "BOX", idtmed[1011], par, 3);

//Backplane supports
//Aluminium slab
  par[0] = 33.15;par[1] = 2;par[2] = 21.65;  gMC->Gsvolu("BACK", "BOX", idtmed[1009], par, 3);    
//Big hole
  par[0] = 9.05; par[1] = 2; par[2] = 4.4625;  gMC->Gsvolu("BKHL", "BOX", idtmed[1000], par, 3);
//Small hole
  par[0] = 5.7;par[1] = 2;par[2] = 4.4625;  gMC->Gsvolu("BKHS", "BOX", idtmed[1000], par, 3);
//Place holes inside backplane support
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
//Place material inside RICH 
  gMC->Gspos("SRIC", 1, "RICH", 0.,0., 0., 0, "ONLY");
  gMC->Gspos("AIR2", 1, "RICH", 66.3 + 1.2505, 1.276-Param()->GapThickness()/2-Param()->QuartzThickness()-Param()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
  gMC->Gspos("AIR2", 2, "RICH", -66.3 - 1.2505,1.276-Param()->GapThickness()/2-Param()->QuartzThickness()-Param()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
  gMC->Gspos("AIR3", 1, "RICH", 0., 1.276-Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, -68.35 - 1.25, 0, "ONLY");
  gMC->Gspos("AIR3", 2, "RICH", 0., 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35,  68.35 + 1.25, 0, "ONLY");
  gMC->Gspos("ALUM", 1, "SRIC", 0., 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .6 - .05 - .376 -.025, 0., 0, "ONLY");
  gMC->Gspos("HONE", 1, "SRIC", 0., 1.276- Param()->GapThickness()/2  - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .6 - .05 - .188, 0., 0, "ONLY");
  gMC->Gspos("ALUM", 2, "SRIC", 0., 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .6 - .025, 0., 0, "ONLY");
  gMC->Gspos("FOOT", 1, "SRIC", 64.95, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 2, "SRIC", 21.65, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .3 , 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 3, "SRIC", -21.65, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 4, "SRIC", -64.95, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 5, "SRIC", 64.95, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 6, "SRIC", 21.65, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 7, "SRIC", -21.65, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 8, "SRIC", -64.95, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("OQUA", 1, "SRIC", 0., 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()- .2, 0., 0, "ONLY");
// Methane supports
  gMC->Gspos("SMLG", 1, "SRIC", pcX/2 + (66.3 - pcX/2)/2, 1.276 + .25, 0., 0, "ONLY");
  gMC->Gspos("SMLG", 2, "SRIC", - pcX/2 - (66.3 - pcX/2)/2, 1.276 + .25, 0., 0, "ONLY");
  gMC->Gspos("SMSH", 1, "SRIC", 0., 1.276 + .25, pcY/2 + (68.35 - pcY/2)/2, 0, "ONLY");
  gMC->Gspos("SMSH", 2, "SRIC", 0., 1.276 + .25, - pcY/2 - (68.35 - pcY/2)/2, 0, "ONLY");
//Freon supports
  Float_t supp_y = 1.276 - Param()->GapThickness()/2- Param()->QuartzThickness() -Param()->FreonThickness() - .2 + .3; //y position of freon supports
  gMC->Gspos("SFLG", 1, "SRIC", Param()->QuartzWidth()/2 + (66.3 - Param()->QuartzWidth()/2)/2, supp_y, 0., 0, "ONLY");
  gMC->Gspos("SFLG", 2, "SRIC", - Param()->QuartzWidth()/2 - (66.3 - Param()->QuartzWidth()/2)/2, supp_y, 0., 0, "ONLY");
  gMC->Gspos("SFSH", 1, "SRIC", 0., supp_y, Param()->QuartzLength()/2 + (68.35 - Param()->QuartzLength()/2)/2, 0, "ONLY");
  gMC->Gspos("SFSH", 2, "SRIC", 0., supp_y, - Param()->QuartzLength()/2 - (68.35 - Param()->QuartzLength()/2)/2, 0, "ONLY");
  AliMatrix(idrotm[1019], 0., 0., 90., 0., 90., 90.);
//Place spacers
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
  gMC->Gspos("OQF1", 1, "SRIC", Param()->OuterFreonWidth()/2 + Param()->InnerFreonWidth()/2 + 2, 1.276 - Param()->GapThickness()/2- Param()->QuartzThickness() -Param()->FreonThickness()/2, 0., 0, "ONLY"); //Original settings (31.3)
  gMC->Gspos("OQF2", 2, "SRIC", 0., 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()/2, 0., 0, "ONLY");          //Original settings 
  gMC->Gspos("OQF1", 3, "SRIC", - (Param()->OuterFreonWidth()/2 + Param()->InnerFreonWidth()/2) - 2, 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness() - Param()->FreonThickness()/2, 0., 0, "ONLY");       //Original settings (-31.3)
  gMC->Gspos("QUAR", 1, "SRIC", 0., 1.276 - Param()->GapThickness()/2 - Param()->QuartzThickness()/2, 0., 0, "ONLY");
  gMC->Gspos("GAP ", 1, "META", 0., Param()->GapThickness()/2 - Param()->ProximityGapThickness()/2 - 0.0001, 0., 0, "ONLY");
  gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
  gMC->Gspos("CSI ", 1, "SRIC", 0., 1.276 + Param()->GapThickness()/2 + .25, 0., 0, "ONLY");
//Wire support placing
  gMC->Gspos("WSG2", 1, "GAP ", 0., Param()->ProximityGapThickness()/2 - .1, 0., 0, "ONLY");
  gMC->Gspos("WSG1", 1, "CSI ", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("WSMe", 1, "SRIC ", 0., 1.276 + Param()->GapThickness()/2 + .5 + 1.05, 0., 0, "ONLY");
//Backplane placing
  gMC->Gspos("BACK", 1, "SRIC ", -33.15, 1.276 + Param()->GapThickness()/2 + .5 + 2.1 + 2, 43.3, 0, "ONLY");
  gMC->Gspos("BACK", 2, "SRIC ", 33.15, 1.276 + Param()->GapThickness()/2 + .5 + 2.1 + 2 , 43.3, 0, "ONLY");
  gMC->Gspos("BACK", 3, "SRIC ", -33.15, 1.276 + Param()->GapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
  gMC->Gspos("BACK", 4, "SRIC ", 33.15, 1.276 + Param()->GapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
  gMC->Gspos("BACK", 5, "SRIC ", 33.15, 1.276 + Param()->GapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");
  gMC->Gspos("BACK", 6, "SRIC ", -33.15, 1.276 + Param()->GapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");
//PCB placing
  gMC->Gspos("PCB ", 1, "SRIC ", 0.,  1.276 + Param()->GapThickness()/2 + .5 + 1.05, pcX/4 + .5025 + 2.5, 0, "ONLY");
  gMC->Gspos("PCB ", 2, "SRIC ", 0.,  1.276 + Param()->GapThickness()/2 + .5 + 1.05, -pcX/4 - .5025 - 2.5, 0, "ONLY");

//place chambers into mother volume ALIC
  for(int i=1;i<=kNCH;i++){
    AliMatrix(idrotm[1000+i],C(i)->ThetaXd(),C(i)->PhiXd(),
                             C(i)->ThetaYd(),C(i)->PhiYd(),
                             C(i)->ThetaZd(),C(i)->PhiZd());
    gMC->Gspos("RICH",i,"ALIC",C(i)->X(),C(i)->Y(),C(i)->Z(),idrotm[1000+i], "ONLY");
  }

  if(GetDebug())Info("CreateGeometry","Stop.");  
}//void AliRICH::CreateGeometry()
//__________________________________________________________________________________________________
void AliRICH::CreateChambers()
{//create all RICH Chambers, do nothing if already created.
  if(fChambers) return;//
  if(GetDebug())Info("CreateChambers","Creating RICH chambers.");
  fChambers=new TObjArray(kNCH);
  fChambers->SetOwner();
  for(int i=0;i<kNCH;i++)  fChambers->AddAt(new AliRICHChamber(i+1,Param()),i);  
}//void AliRICH::CreateChambers()
//__________________________________________________________________________________________________
void AliRICH::GenerateFeedbacks(Int_t iChamber,Float_t eloss)
{// Generate FeedBack photons
  Int_t j;
  Float_t cthf, phif, enfp = 0, sthf;
  Float_t e1[3], e2[3], e3[3];
  Float_t vmod, uswop;
  Float_t dir[3], phi;
  Float_t pol[3], mom[4];
//Determine number of feedback photons
  TLorentzVector x4;
  gMC->TrackPosition(x4);  
  Float_t charge=Param()->TotalCharge(gMC->TrackPid(),eloss,C(iChamber)->G2Ly(x4));//Total Charge
  Int_t iNphotons=gMC->GetRandom()->Poisson(Param()->AlphaFeedback()*charge);    
  Info("GenerateFeedbacks","N photons=%i",iNphotons);
//Generate photons
  for(Int_t i=0;i<iNphotons;i++){
    Double_t ranf[2];
    gMC->GetRandom()->RndmArray(2,ranf);    //Sample direction
    cthf=ranf[0]*2-1.0;
    if(cthf<0) continue;
    sthf = TMath::Sqrt((1 - cthf) * (1 + cthf));
    phif = ranf[1] * 2 * TMath::Pi();
    
    if(Double_t randomNumber=gMC->GetRandom()->Rndm()<=0.57)
      enfp = 7.5e-9;
    else if(randomNumber<=0.7)
      enfp = 6.4e-9;
    else
      enfp = 7.9e-9;
    

    dir[0] = sthf * TMath::Sin(phif);    dir[1] = cthf;    dir[2] = sthf * TMath::Cos(phif);
    gMC->Gdtom(dir, mom, 2);
    mom[0]*=enfp;    mom[1]*=enfp;    mom[2]*=enfp;
    mom[3] = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
    
    // Polarisation
    e1[0]=      0;    e1[1]=-dir[2];    e1[2]= dir[1];
    e2[0]=-dir[1];    e2[1]= dir[0];    e2[2]=      0;
    e3[0]= dir[1];    e3[1]=      0;    e3[2]=-dir[0];
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e1[j];
      e1[j]=e3[j];
      e3[j]=uswop;
    }
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e2[j];
      e2[j]=e3[j];
      e3[j]=uswop;
    }
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    vmod=TMath::Sqrt(1/vmod);
    for(j=0;j<3;j++) e1[j]*=vmod;
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    vmod=TMath::Sqrt(1/vmod);
    for(j=0;j<3;j++) e2[j]*=vmod;
    
    phi = gMC->GetRandom()->Rndm()* 2 * TMath::Pi();
    for(j=0;j<3;j++) pol[j]=e1[j]*TMath::Sin(phi)+e2[j]*TMath::Cos(phi);
    gMC->Gdtom(pol, pol, 2);
    Int_t outputNtracksStored;    
    gAlice->PushTrack(1,                             //do not transport
                     gAlice->GetCurrentTrackNumber(),//parent track 
                     kFeedback,                      //PID
		     mom[0],mom[1],mom[2],mom[3],    //track momentum  
                     x4.X(),x4.Y(),x4.Z(),x4.T(),    //track origin 
                     pol[0],pol[1],pol[2],           //polarization
		     kPFeedBackPhoton,outputNtracksStored,1.0);
    
  }
}//Int_t AliRICH::FeedBackPhotons()
//__________________________________________________________________________________________________
