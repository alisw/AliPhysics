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


#include <strings.h>

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
#include "AliMagF.h"
#include "AliPoints.h"
#include "AliRICH.h"
#include "AliRICHParam.h"
#include "AliRICHClusterFinder.h"
#include "AliRICHDigitizer.h"
#include "AliRICHHitMapA1.h"
#include "AliRICHMerger.h"
#include "AliRICHRawCluster.h"
#include "AliRICHRecHit1D.h"
#include "AliRICHRecHit3D.h"
#include "AliRICHSDigit.h"
#include "AliRICHSegmentationV0.h"
#include "AliRICHTransientDigit.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRICHSegmentationV1.h"
#include "AliRICHResponseV0.h"
 
ClassImp(AliRICHhit)
ClassImp(AliRICHdigit)
ClassImp(AliRICH)
    
//___________________________________________
// RICH manager class   
//Begin_Html
/*
  <img src="gif/alirich.gif">
*/
//End_Html

AliRICH::AliRICH()
{//Default ctor should not contain any new operators
  fIshunt     = 0;
  fHits       = 0;
  fSDigits    = 0;
  fNsdigits   = 0;
  fNcerenkovs = 0;
  fDchambers  = 0;
  fRecHits1D = 0;
  fRecHits3D = 0;
  fRawClusters = 0;
  fChambers = 0;
  fCerenkovs  = 0;
  for (Int_t i=0; i<kNCH; i++){
      fNdch[i]       = 0;
      fNrawch[i]     = 0;
      fNrechits1D[i] = 0;
      fNrechits3D[i] = 0;
  }
  fpParam=0;
}//AliRICH::AliRICH()
//__________________________________________________________________________________________________
AliRICH::AliRICH(const char *name, const char *title)
        :AliDetector(name,title)
{//Named ctor
  if(GetDebug())Info("named ctor","Start.");
  fpParam     =new AliRICHParam;
  fNcerenkovs =fNsdigits   =0;//fNhits and fNdigits reset in AliDetector ctor
  fHits       =new TClonesArray("AliRICHhit",1000  ); gAlice->AddHitList(fHits);//hits
  fCerenkovs  =new TClonesArray("AliRICHCerenkov",1000);gAlice->AddHitList(fCerenkovs);//cerenkovs ??? to be removed    
  fSDigits    =new TClonesArray("AliRICHdigit",100000);//sdigits
  fDigits     =new TClonesArray("AliRICHdigit",100000);//digits
  
  fIshunt     =0;
  fDchambers  =new TObjArray(kNCH);//digits             ??? to be removed
  fRawClusters=new TObjArray(kNCH);//clusters
  fRecHits1D  =new TObjArray(kNCH);//recos Bari
  fRecHits3D  =new TObjArray(kNCH);//recos Lisbon
  for(int i=0;i<kNCH;i++) {
    fNdch[i]=fNrawch[i]=0;
    fDchambers  ->AddAt(new TClonesArray("AliRICHDigit",10000), i); //??? to be removed
    fRawClusters->AddAt(new TClonesArray("AliRICHRawCluster",10000), i); 
    fRecHits1D  ->AddAt(new TClonesArray("AliRICHRecHit1D",1000), i);
    fRecHits3D  ->AddAt(new TClonesArray("AliRICHRecHit3D",1000), i);
  }
  SetMarkerColor(kRed);
  fCkovNumber=fFreonProd=0;
  CreateChambers();
  if(GetDebug())Info("named ctor","Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//__________________________________________________________________________________________________
AliRICH::~AliRICH()
{//dtor
  if(GetDebug()) Info("dtor","Start.");

    fIshunt  = 0;
    delete fHits;
    delete fSDigits;
    delete fCerenkovs;
    
    //PH Delete TObjArrays
    if (fChambers) {
      fChambers->Delete();
      delete fChambers;
    }
    if (fDchambers) {
      fDchambers->Delete();
      delete fDchambers;
    }
    if (fRawClusters) {
      fRawClusters->Delete();
      delete fRawClusters;
    }
    if (fRecHits1D) {
      fRecHits1D->Delete();
      delete fRecHits1D;
    }
    if (fRecHits3D) {
      fRecHits3D->Delete();
      delete fRecHits3D;
    }                     
  if(GetDebug()) Info("dtor","Stop.");    
}//AliRICH::~AliRICH()
//__________________________________________________________________________________________________
void AliRICH::Hits2SDigits()
{//Create a list of sdigits corresponding to list of hits. Every hit generates sdigit.
  if(GetDebug()) Info("Hit2SDigits","Start.");
  
  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//events loop
    fLoader->GetRunLoader()->GetEvent(iEventN);
  
    if(!fLoader->TreeH()) fLoader->LoadHits();
    if(!fLoader->TreeS()) fLoader->MakeTree("S");
    MakeBranch("S");
    
    AliRICHSegmentationV1 *pSeg=new AliRICHSegmentationV1;
    AliRICHResponseV0     *pRes=new AliRICHResponseV0;
    
    Float_t dx=Param()->SigmaIntegration()*Param()->ChargeSpreadX();
    Float_t dy=Param()->SigmaIntegration()*Param()->ChargeSpreadY();
    Float_t charge;
    
    for(int iPrimN=0;iPrimN<TreeH()->GetEntries();iPrimN++){//prims loop
      fLoader->TreeH()->GetEntry(iPrimN); 
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){//hits loop
        AliRICHhit *pHit=(AliRICHhit*)Hits()->At(iHitN);//get current hit
        
        if(pHit->Y()>0)
          charge=Param()->TotalCharge(pHit->Particle(),pHit->Eloss(),pHit->Y()-Param()->SectorSizeY());
        else
          charge=Param()->TotalCharge(pHit->Particle(),pHit->Eloss(),pHit->Y()+Param()->SectorSizeY());
                  
        for(pSeg->FirstPad(pHit->X(),pHit->Y(),0,dx,dy);pSeg->MorePads();pSeg->NextPad()){//pads loop
            
          AddSDigit(pHit->Chamber(),pSeg->Ix(),pSeg->Iy(),
                                    Int_t(charge*TMath::Abs(pRes->IntXY(pSeg))),
                iPrimN);//chamber-xpad-ypad-qdc-track1-2-3
        }//pads loop
      }//hits loop
    }//prims loop
    
    delete pSeg;
    
    fLoader->TreeS()->Fill();
    fLoader->WriteSDigits("OVERWRITE");
  }//events loop
  
  fLoader->UnloadHits();
  fLoader->UnloadSDigits();  
  if(GetDebug()) Info("Hit2SDigits","Stop.");
}//void AliRICH::Hits2SDigits()
//__________________________________________________________________________________________________
void AliRICH::SDigits2Digits()
{//Generate digits from sdigits.
  if(GetDebug()) Info("SDigits2Digits","Start.");
    
  AliRunDigitizer *pManager = new AliRunDigitizer(1,1);
  pManager->SetInputStream(0,"galice.root");
  pManager->Exec("deb");
  if(GetDebug()) Info("SDigits2Digits","Stop.");
}//void AliRICH::SDigits2Digits()
//__________________________________________________________________________________________________
void AliRICH::Digits2Reco()
{
// Generate clusters
// Called from alirun, single event only.     
  if(GetDebug()) Info("Digits2Reco","Start.");

  int nparticles = gAlice->GetNtrack();
  cout << "Particles (RICH):" <<nparticles<<endl;
  if (nparticles > 0) FindClusters(0);

}//void AliRICH::Digits2Reco()  
//__________________________________________________________________________________________________


void AliRICH::AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{// Add a RICH digit to the list   

   TClonesArray &ldigits = *((TClonesArray*)fDchambers->At(id));
   new(ldigits[fNdch[id]++]) AliRICHDigit(tracks,charges,digits);
}

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

  Float_t wid=fpParam->SectorSizeX();
  Float_t len=fpParam->SectorSizeY();
  new TBRIK("PHOTO","PHOTO","void",wid/2,0.1,len/2);
  
  for(int i=1;i<=kNCH;i++){
    top->cd();
    node = new TNode(Form("RICH%i",i),Form("RICH%i",i),"S_RICH",C(i)->X(),C(i)->Y(),C(i)->Z(),C(i)->RotMatrixName());
    node->SetLineColor(kRed);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",wid+fpParam->DeadZone(),5,len/2+fpParam->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,len/2+fpParam->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-wid-fpParam->DeadZone(),5,len/2+fpParam->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",wid+fpParam->DeadZone(),5,-len/2-fpParam->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-len/2 -fpParam->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-wid-fpParam->DeadZone(),5,-len/2 - fpParam->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);
  }  
  if(GetDebug())Info("BuildGeometry","Stop.");    
}//void AliRICH::BuildGeometry()
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
    Int_t   isxfld = gAlice->Field()->Integ();
    Float_t sxmgmx = gAlice->Field()->Max();
    Int_t i;

    /************************************Antonnelo's Values (14-vectors)*****************************************/
    /*
    Float_t ppckov[14] = { 5.63e-9,5.77e-9,5.9e-9,6.05e-9,6.2e-9,6.36e-9,6.52e-9,
			   6.7e-9,6.88e-9,7.08e-9,7.3e-9,7.51e-9,7.74e-9,8e-9 };
    Float_t rIndexQuarz[14] = { 1.528309,1.533333,
				 1.538243,1.544223,1.550568,1.55777,
				 1.565463,1.574765,1.584831,1.597027,
			       1.611858,1.6277,1.6472,1.6724 };
    Float_t rIndexOpaqueQuarz[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t rIndexMethane[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t rIndexGrid[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t abscoFreon[14] = { 179.0987,179.0987,
				179.0987,179.0987,179.0987,142.92,56.65,13.95,10.43,7.07,2.03,.5773,.33496,0. };
    //Float_t abscoFreon[14] = { 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,
	//			 1e-5,1e-5,1e-5,1e-5,1e-5 };
    Float_t abscoQuarz[14] = { 64.035,39.98,35.665,31.262,27.527,22.815,21.04,17.52,
				14.177,9.282,4.0925,1.149,.3627,.10857 };
    Float_t abscoOpaqueQuarz[14] = { 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,
				 1e-5,1e-5,1e-5,1e-5,1e-5 };
    Float_t abscoCsI[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			      1e-4,1e-4,1e-4,1e-4 };
    Float_t abscoMethane[14] = { 1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,
				  1e6,1e6,1e6 };
    Float_t abscoGrid[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			      1e-4,1e-4,1e-4,1e-4 };
    Float_t efficAll[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t efficCsI[14] = { 6e-4,.005,.0075,.01125,.045,.117,.135,.16575,
			      .17425,.1785,.1836,.1904,.1938,.221 };
    Float_t efficGrid[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    */
   
    
    /**********************************End of Antonnelo's Values**********************************/
    
    /**********************************Values from rich_media.f (31-vectors)**********************************/
    

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
    
    //Absorption index for quarz
    /*Float_t Qzt [21] = {.0,.0,.005,.04,.35,.647,.769,.808,.829,.844,.853,.858,.869,.887,.903,.902,.902,
	 		.906,.907,.907,.907};
    Float_t Wavl2[] = {150.,155.,160.0,165.0,170.0,175.0,180.0,185.0,190.0,195.0,200.0,205.0,210.0,
	 	       215.0,220.0,225.0,230.0,235.0,240.0,245.0,250.0};		 		 
    Float_t abscoQuarz[31];	     
    for (Int_t i=0;i<31;i++)
    {
	Float_t Xlam = 1237.79 / (ppckov[i]*1e9);
	if (Xlam <= 160) abscoQuarz[i] = 0;
	if (Xlam > 250) abscoQuarz[i] = 1;
	else 
	{
	    for (Int_t j=0;j<21;j++)
	    {
		if (Xlam > Wavl2[j] && Xlam < Wavl2[j+1])
		{
		    Float_t Dabs = (Qzt[j+1] - Qzt[j])/(Wavl2[j+1] - Wavl2[j]);
		    Float_t Abso = Qzt[j] + Dabs*(Xlam - Wavl2[j]);
		    abscoQuarz[i] = -5.0/(TMath::Log(Abso));
		} 
	    }
	}
    }*/

    /*Float_t abscoQuarz[31] = {49.64211, 48.41296, 47.46989, 46.50492, 45.13682, 44.47883, 43.1929 , 41.30922, 40.5943 ,
			       39.82956, 38.98623, 38.6247 , 38.43448, 37.41084, 36.22575, 33.74852, 30.73901, 24.25086, 
			       17.94531, 11.88753, 5.99128,  3.83503,  2.36661,  1.53155, 1.30582, 1.08574, .8779708, 
			       .675275, 0., 0., 0.};
    
    for (Int_t i=0;i<31;i++)
    {
	abscoQuarz[i] = abscoQuarz[i]/10;
    }*/

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
}
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
void AliRICH::ResetDigits()
{//Reset number of digits and the digits array for this detector
  for ( int i=0;i<kNCH;i++ ) {
    if (fDchambers && fDchambers->At(i))   fDchambers->At(i)->Clear();
    if (fNdch)  fNdch[i]=0;
  }
}
//__________________________________________________________________________________________________
void AliRICH::ResetRawClusters()
{//Reset number of raw clusters and the raw clust array for this detector
  for ( int i=0;i<kNCH;i++ ) {
    if (fRawClusters->At(i))    ((TClonesArray*)fRawClusters->At(i))->Clear();
    if (fNrawch)  fNrawch[i]=0;
  }
}
//__________________________________________________________________________________________________
void AliRICH::ResetRecHits1D()
{//Reset number of raw clusters and the raw clust array for this detector
  for ( int i=0;i<kNCH;i++ ) {
    if (fRecHits1D->At(i))    ((TClonesArray*)fRecHits1D->At(i))->Clear();
    if (fNrechits1D)  fNrechits1D[i]=0;
  }
}

//__________________________________________________________________________________________________
void AliRICH::ResetRecHits3D()
{// Reset number of raw clusters and the raw clust array for this detector
  for ( int i=0;i<kNCH;i++ ) {
    if (fRecHits3D->At(i))    ((TClonesArray*)fRecHits3D->At(i))->Clear();
    if (fNrechits3D)  fNrechits3D[i]=0;
  }
}
//__________________________________________________________________________________________________
void AliRICH::FindClusters(Int_t nev /*kir,Int_t lastEntry*/)
{// Loop on chambers and on cathode planes
    for (Int_t icat=1;icat<2;icat++) {
	gAlice->ResetDigits();
	gAlice->TreeD()->GetEvent(0);
	for (Int_t ich=0;ich<kNCH;ich++) {
	  AliRICHChamber* iChamber=(AliRICHChamber*)fChambers->At(ich);
	  TClonesArray *pRICHdigits  = this->DigitsAddress(ich);
	  if (pRICHdigits == 0)	      
	      continue;
	  //
	  // Get ready the current chamber stuff
	  //
	  AliRICHResponse* response = iChamber->GetResponseModel();
	  AliSegmentation*  seg = iChamber->GetSegmentationModel();
	  AliRICHClusterFinder* rec = iChamber->GetReconstructionModel();
	  if (seg) {	  
	      rec->SetSegmentation(seg);
	      rec->SetResponse(response);
	      rec->SetDigits(pRICHdigits);
	      rec->SetChamber(ich);
	      if (nev==0) rec->CalibrateCOG(); 
	      rec->FindRawClusters();
	  }  
	  TClonesArray *fRch;
	  fRch=RawClustAddress(ich);
	  fRch->Sort();
	} // for ich

	gAlice->TreeR()->Fill();
	TClonesArray *fRch;
	for (int i=0;i<kNCH;i++) {
	    fRch=RawClustAddress(i);
	    fRch->GetEntriesFast();
	}
	
	ResetRawClusters();
	
    } // for icat
    
    char hname[30];
    sprintf(hname,"TreeR%d",nev);
    gAlice->TreeR()->Write(hname,kOverwrite,0);
    gAlice->TreeR()->Reset();    
}//void AliRICH::FindClusters(Int_t nev)
//__________________________________________________________________________________________________
void AliRICH::MakeBranchInTreeD(TTree *treeD, const char *file)
{// Create TreeD branches for the RICH.
  if(GetDebug())Info("MakeBranchInTreeD","Start.");

  const Int_t kBufferSize = 4000;
  char branchname[30];
    
  //
  // one branch for digits per chamber
  // 
  for (Int_t i=0; i<kNCH ;i++) {
    sprintf(branchname,"%sDigits%d",GetName(),i+1);	
    if (fDchambers && treeD) {
      MakeBranchInTree(treeD,branchname, &((*fDchambers)[i]), kBufferSize, file);
//      printf("Making Branch %s for digits in chamber %d\n",branchname,i+1);
    }
  }
}
//__________________________________________________________________________________________________
void AliRICH::MakeBranch(Option_t* option)
{//Create Tree branches for the RICH.
  if(GetDebug())Info("MakeBranch","Start with option= %s.",option);
    
  const Int_t kBufferSize = 4000;
  char branchname[20];
      
   
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");

  if(cH&&TreeH()){
    if(!fHits) fHits=new TClonesArray("AliRICHhit",1000  );
    if(!fCerenkovs) fCerenkovs  = new TClonesArray("AliRICHCerenkov",1000);
    MakeBranchInTree(TreeH(),"RICHCerenkov", &fCerenkovs, kBufferSize, 0) ;

    //kir if(!fSDigits) fSDigits    = new TClonesArray("AliRICHdigit",100000);
    //kir MakeBranchInTree(TreeH(),"RICHSDigits", &fSDigits, kBufferSize, 0) ;
  }     
  AliDetector::MakeBranch(option);//this is after cH because we need to guarantee that fHits array is created
      
  if(cS&&fLoader->TreeS()){  
    if(!fSDigits) fSDigits=new TClonesArray("AliRICHdigit",100000);
    MakeBranchInTree(fLoader->TreeS(),"RICH",&fSDigits,kBufferSize,0) ;
  }
   
  int i;
  if (cD&&fLoader->TreeD()){
    if(!fDchambers){
      fDchambers=new TObjArray(kNCH);    // one branch for digits per chamber
      for(i=0;i<kNCH;i++){ 
        fDchambers->AddAt(new TClonesArray("AliRICHDigit",10000), i); 
      }       
    }
    for (i=0; i<kNCH ;i++) 
      {
        sprintf(branchname,"%sDigits%d",GetName(),i+1);	
        MakeBranchInTree(fLoader->TreeD(),branchname, &((*fDchambers)[i]), kBufferSize, 0);
      }
   }

  if (cR&&gAlice->TreeR()){//one branch for raw clusters per chamber
    Int_t i;
    if (fRawClusters == 0x0 ) 
     {
       fRawClusters = new TObjArray(kNCH);
       for (i=0; i<kNCH ;i++) 
         {
           fRawClusters->AddAt(new TClonesArray("AliRICHRawCluster",10000), i); 
         }
     }
     
    if (fRecHits1D == 0x0) 
     {
        fRecHits1D = new TObjArray(kNCH);
        for (i=0; i<kNCH ;i++) 
         {
          fRecHits1D->AddAt(new TClonesArray("AliRICHRecHit1D",1000), i);
         }
     }

    if (fRecHits3D == 0x0) 
     {
        fRecHits3D = new TObjArray(kNCH);
        for (i=0; i<kNCH ;i++) 
         {
          fRecHits3D->AddAt(new TClonesArray("AliRICHRecHit3D",1000), i);
         }
     }
       
    for (i=0; i<kNCH ;i++){
       sprintf(branchname,"%sRawClusters%d",GetName(),i+1);      
       MakeBranchInTree(gAlice->TreeR(),branchname, &((*fRawClusters)[i]), kBufferSize, 0);
       sprintf(branchname,"%sRecHits1D%d",GetName(),i+1);
       MakeBranchInTree(fLoader->TreeR(),branchname, &((*fRecHits1D)[i]), kBufferSize, 0);
       sprintf(branchname,"%sRecHits3D%d",GetName(),i+1);  
       MakeBranchInTree(fLoader->TreeR(),branchname, &((*fRecHits3D)[i]), kBufferSize, 0);
     }
   }//if (cR && gAlice->TreeR())
  if(GetDebug())Info("MakeBranch","Stop.");   
}//void AliRICH::MakeBranch(Option_t* option)
//______________________________________________________________________________
void AliRICH::SetTreeAddress()
{//Set branch address for the Hits and Digits Tree.
  if(GetDebug())Info("SetTreeAddress","Start.");
  
  char branchname[20];
  Int_t i;

    
  TBranch *branch;
  TTree *treeH = fLoader->TreeH();
  TTree *treeD = fLoader->TreeD();
  TTree *treeR = fLoader->TreeR();
    
  if(treeH){
    if(GetDebug())Info("SetTreeAddress","tree H is requested.");
    if(fHits==0x0) fHits=new TClonesArray("AliRICHhit",1000); 
    
    branch = treeH->GetBranch("RICHCerenkov");
    if(branch){
      if (fCerenkovs == 0x0) fCerenkovs  = new TClonesArray("AliRICHCerenkov",1000); 
        branch->SetAddress(&fCerenkovs);
    }
       
//kir      branch = treeH->GetBranch("RICHSDigits");
//kir      if (branch) 
//kir       {
//kir         if (fSDigits == 0x0) fSDigits    = new TClonesArray("AliRICHdigit",100000);
//kir         branch->SetAddress(&fSDigits);
//kir       }
  }//if(treeH)
 
  AliDetector::SetTreeAddress();//this is after TreeH because we need to guarantee that fHits array is created

    
  if(fLoader->TreeS()){
    if(GetDebug())Info("SetTreeAddress","tree S is requested.");
    if(!fSDigits) fSDigits=new TClonesArray("AliRICHdigit",100000);
    fLoader->TreeS()->GetBranch("RICH")->SetAddress(&fSDigits);
  }
    
    
  if(treeD){
    if(GetDebug())Info("SetTreeAddress","tree D is requested.");

      if (fDchambers == 0x0) 
        {
           fDchambers = new TObjArray(kNCH);
           for (i=0; i<kNCH ;i++) 
             {
               fDchambers->AddAt(new TClonesArray("AliRICHDigit",10000), i); 
             }
        }
      
      for (i=0; i<kNCH; i++) {
        sprintf(branchname,"%sDigits%d",GetName(),i+1);
        if (fDchambers) {
           branch = treeD->GetBranch(branchname);
           if (branch) branch->SetAddress(&((*fDchambers)[i]));
        }
      }
    }
    
  if(treeR){
    if(GetDebug())Info("SetTreeAddress","tree R is requested.");

    if (fRawClusters == 0x0 ) 
     {
       fRawClusters = new TObjArray(kNCH);
       for (i=0; i<kNCH ;i++) 
         {
           fRawClusters->AddAt(new TClonesArray("AliRICHRawCluster",10000), i); 
         }
     }
     
    if (fRecHits1D == 0x0) 
     {
        fRecHits1D = new TObjArray(kNCH);
        for (i=0; i<kNCH ;i++) 
         {
          fRecHits1D->AddAt(new TClonesArray("AliRICHRecHit1D",1000), i);
         }
     }

    if (fRecHits3D == 0x0) 
     {
        fRecHits3D = new TObjArray(kNCH);
        for (i=0; i<kNCH ;i++) 
         {
          fRecHits3D->AddAt(new TClonesArray("AliRICHRecHit3D",1000), i);
         }
     }
    
    for (i=0; i<kNCH; i++) {
	  sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	  if (fRawClusters) {
	      branch = treeR->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	  }
    }
      
    for (i=0; i<kNCH; i++) {
	sprintf(branchname,"%sRecHits1D%d",GetName(),i+1);
	if (fRecHits1D) {
	  branch = treeR->GetBranch(branchname);
	  if (branch) branch->SetAddress(&((*fRecHits1D)[i]));
	  }
     }
      
     for (i=0; i<kNCH; i++) {
	sprintf(branchname,"%sRecHits3D%d",GetName(),i+1);
	if (fRecHits3D) {
	  branch = treeR->GetBranch(branchname);
	  if (branch) branch->SetAddress(&((*fRecHits3D)[i]));
	  }
      } 
      
  }//if(treeR)
  if(GetDebug())Info("SetTreeAddress","Stop.");
}//void AliRICH::SetTreeAddress()
//__________________________________________________________________________________________________
void AliRICH::Print(Option_t *option)const
{
  TObject::Print(option);
  fpParam->Dump();
  fChambers->Print(option);  
}//void AliRICH::Print(Option_t *option)const
//__________________________________________________________________________________________________
void AliRICH::CreateGeometry()
{//Creates detailed geometry simulation (currently GEANT volumes tree)         
  if(GetDebug())Info("CreateGeometry","Start.");
//???????? to be removed to AliRICHParam?
  fpParam->RadiatorToPads(fpParam->FreonThickness()/2+fpParam->QuartzThickness()+fpParam->GapThickness());
    
//Opaque quartz thickness
  Float_t oqua_thickness = .5;
//CsI dimensions
  Float_t csi_width =fpParam->Nx()*fpParam->PadSizeX()+fpParam->DeadZone();
  Float_t csi_length=fpParam->Ny()*fpParam->PadSizeY()+2*fpParam->DeadZone();
  
  Int_t *idtmed = fIdtmed->GetArray()-999;
    
  Int_t i;
  Float_t zs;
  Int_t idrotm[1099];
  Float_t par[3];
    
//External aluminium box 
  par[0]=68.8*kcm;par[1]=13*kcm;par[2]=70.86*kcm;//Original Settings
  gMC->Gsvolu("RICH", "BOX ", idtmed[1009], par, 3);
//Air 
  par[0]=66.3;   par[1] = 13; par[2] = 68.35; //Original Settings
  gMC->Gsvolu("SRIC", "BOX ", idtmed[1000], par, 3); 
//Air 2 (cutting the lower part of the box)
  par[0]=1.25;    par[1] = 3;    par[2] = 70.86; //Original Settings
  gMC->Gsvolu("AIR2", "BOX ", idtmed[1000], par, 3);
//Air 3 (cutting the lower part of the box)
  par[0]=66.3;    par[1] = 3;  par[2] = 1.2505; //Original Settings
  gMC->Gsvolu("AIR3", "BOX ", idtmed[1000], par, 3);
//Honeycomb 
  par[0]=66.3;par[1]=0.188;  par[2] = 68.35;  //Original Settings
  gMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);
//Aluminium sheet 
  par[0]=66.3;par[1]=0.025;par[2]=68.35; //Original Settings
  //par[0] = 66.5; par[1] = .025; par[2] = 63.1;
  gMC->Gsvolu("ALUM", "BOX ", idtmed[1009], par, 3);
//Quartz 
  par[0]=fpParam->QuartzWidth()/2;par[1]=fpParam->QuartzThickness()/2;par[2]=fpParam->QuartzLength()/2;
  gMC->Gsvolu("QUAR", "BOX ", idtmed[1002], par, 3);
//Spacers (cylinders) 
  par[0]=0.;par[1]=.5;par[2]=fpParam->FreonThickness()/2;
  gMC->Gsvolu("SPAC", "TUBE", idtmed[1002], par, 3);    
//Feet (freon slabs supports)
  par[0] = .7;  par[1] = .3;  par[2] = 1.9;
  gMC->Gsvolu("FOOT", "BOX", idtmed[1009], par, 3);
//Opaque quartz 
  par[0]=fpParam->QuartzWidth()/2;par[1]= .2;par[2]=fpParam->QuartzLength()/2;
  gMC->Gsvolu("OQUA", "BOX ", idtmed[1007], par, 3);
//Frame of opaque quartz
  par[0]=fpParam->OuterFreonWidth()/2;par[1]=fpParam->FreonThickness()/2;par[2]=fpParam->OuterFreonLength()/2; 
  gMC->Gsvolu("OQF1", "BOX ", idtmed[1007], par, 3);
  par[0]=fpParam->InnerFreonWidth()/2;par[1]=fpParam->FreonThickness()/2;par[2]=fpParam->InnerFreonLength()/2; 
  gMC->Gsvolu("OQF2", "BOX ", idtmed[1007], par, 3);
//Freon 
  par[0]=fpParam->OuterFreonWidth()/2 - oqua_thickness;
  par[1]=fpParam->FreonThickness()/2;
  par[2]=fpParam->OuterFreonLength()/2 - 2*oqua_thickness; 
  gMC->Gsvolu("FRE1", "BOX ", idtmed[1003], par, 3);

  par[0]=fpParam->InnerFreonWidth()/2 - oqua_thickness;
  par[1]=fpParam->FreonThickness()/2;
  par[2]=fpParam->InnerFreonLength()/2 - 2*oqua_thickness; 
  gMC->Gsvolu("FRE2", "BOX ", idtmed[1003], par, 3);    
//Methane 
  par[0]=csi_width/2;par[1]=fpParam->GapThickness()/2;par[2]=csi_length/2;
  gMC->Gsvolu("META", "BOX ", idtmed[1004], par, 3);
//Methane gap 
  par[0]=csi_width/2;par[1]=fpParam->ProximityGapThickness()/2;par[2] = csi_length/2;
  gMC->Gsvolu("GAP ", "BOX ", idtmed[1008], par, 3);
//CsI photocathode 
  par[0]=csi_width/2;par[1]=.25;par[2]=csi_length/2;
  gMC->Gsvolu("CSI ", "BOX ", idtmed[1005], par, 3);
//Anode grid 
  par[0] = 0.;par[1] = .001;par[2] = 20.;
  gMC->Gsvolu("GRID", "TUBE", idtmed[1006], par, 3);

//Wire supports
//Bar of metal
  par[0]=csi_width/2;par[1]=1.05;par[2]=1.05;
  gMC->Gsvolu("WSMe", "BOX ", idtmed[1009], par, 3);
//Ceramic pick up (base)
  par[0]=csi_width/2;par[1]= .25;par[2]=1.05;
  gMC->Gsvolu("WSG1", "BOX ", idtmed[1010], par, 3);
//Ceramic pick up (head)
  par[0] = csi_width/2;par[1] = .1;par[2] = .1;
  gMC->Gsvolu("WSG2", "BOX ", idtmed[1010], par, 3);

//Aluminium supports for methane and CsI
//Short bar
  par[0]=csi_width/2;par[1]=fpParam->GapThickness()/2 + .25; par[2] = (68.35 - csi_length/2)/2;
  gMC->Gsvolu("SMSH", "BOX", idtmed[1009], par, 3);
//Long bar
  par[0]=(66.3 - csi_width/2)/2;par[1]=fpParam->GapThickness()/2+.25;par[2]=csi_length/2+68.35-csi_length/2;
  gMC->Gsvolu("SMLG", "BOX", idtmed[1009], par, 3);
    
//Aluminium supports for freon
//Short bar
  par[0] = fpParam->QuartzWidth()/2; par[1] = .3; par[2] = (68.35 - fpParam->QuartzLength()/2)/2;
  gMC->Gsvolu("SFSH", "BOX", idtmed[1009], par, 3);    
//Long bar
  par[0] = (66.3 - fpParam->QuartzWidth()/2)/2; par[1] = .3;
  par[2] = fpParam->QuartzLength()/2 + 68.35 - fpParam->QuartzLength()/2;
  gMC->Gsvolu("SFLG", "BOX", idtmed[1009], par, 3);    
//PCB backplane
  par[0] = csi_width/2;par[1] = .25; par[2] = csi_length/4 -.5025;
  gMC->Gsvolu("PCB ", "BOX", idtmed[1011], par, 3);

//Backplane supports
//Aluminium slab
  par[0] = 33.15;par[1] = 2;par[2] = 21.65;
  gMC->Gsvolu("BACK", "BOX", idtmed[1009], par, 3);    
//Big hole
  par[0] = 9.05; par[1] = 2; par[2] = 4.4625;
  gMC->Gsvolu("BKHL", "BOX", idtmed[1000], par, 3);
//Small hole
  par[0] = 5.7;par[1] = 2;par[2] = 4.4625;
  gMC->Gsvolu("BKHS", "BOX", idtmed[1000], par, 3);
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
  gMC->Gspos("AIR2", 1, "RICH", 66.3 + 1.2505, 1.276-fpParam->GapThickness()/2-fpParam->QuartzThickness()-fpParam->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
  gMC->Gspos("AIR2", 2, "RICH", -66.3 - 1.2505,1.276-fpParam->GapThickness()/2-fpParam->QuartzThickness()-fpParam->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
  gMC->Gspos("AIR3", 1, "RICH", 0., 1.276-fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, -68.35 - 1.25, 0, "ONLY");
  gMC->Gspos("AIR3", 2, "RICH", 0., 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35,  68.35 + 1.25, 0, "ONLY");
  gMC->Gspos("ALUM", 1, "SRIC", 0., 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .6 - .05 - .376 -.025, 0., 0, "ONLY");
  gMC->Gspos("HONE", 1, "SRIC", 0., 1.276- fpParam->GapThickness()/2  - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .6 - .05 - .188, 0., 0, "ONLY");
  gMC->Gspos("ALUM", 2, "SRIC", 0., 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .6 - .025, 0., 0, "ONLY");
  gMC->Gspos("FOOT", 1, "SRIC", 64.95, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 2, "SRIC", 21.65, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .3 , 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 3, "SRIC", -21.65, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 4, "SRIC", -64.95, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .3, 36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 5, "SRIC", 64.95, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 6, "SRIC", 21.65, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 7, "SRIC", -21.65, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("FOOT", 8, "SRIC", -64.95, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .4 - .3, -36.9, 0, "ONLY");
  gMC->Gspos("OQUA", 1, "SRIC", 0., 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()- .2, 0., 0, "ONLY");
// Methane supports
  gMC->Gspos("SMLG", 1, "SRIC", csi_width/2 + (66.3 - csi_width/2)/2, 1.276 + .25, 0., 0, "ONLY");
  gMC->Gspos("SMLG", 2, "SRIC", - csi_width/2 - (66.3 - csi_width/2)/2, 1.276 + .25, 0., 0, "ONLY");
  gMC->Gspos("SMSH", 1, "SRIC", 0., 1.276 + .25, csi_length/2 + (68.35 - csi_length/2)/2, 0, "ONLY");
  gMC->Gspos("SMSH", 2, "SRIC", 0., 1.276 + .25, - csi_length/2 - (68.35 - csi_length/2)/2, 0, "ONLY");
//Freon supports
  Float_t supp_y = 1.276 - fpParam->GapThickness()/2- fpParam->QuartzThickness() -fpParam->FreonThickness() - .2 + .3; //y position of freon supports
  gMC->Gspos("SFLG", 1, "SRIC", fpParam->QuartzWidth()/2 + (66.3 - fpParam->QuartzWidth()/2)/2, supp_y, 0., 0, "ONLY");
  gMC->Gspos("SFLG", 2, "SRIC", - fpParam->QuartzWidth()/2 - (66.3 - fpParam->QuartzWidth()/2)/2, supp_y, 0., 0, "ONLY");
  gMC->Gspos("SFSH", 1, "SRIC", 0., supp_y, fpParam->QuartzLength()/2 + (68.35 - fpParam->QuartzLength()/2)/2, 0, "ONLY");
  gMC->Gspos("SFSH", 2, "SRIC", 0., supp_y, - fpParam->QuartzLength()/2 - (68.35 - fpParam->QuartzLength()/2)/2, 0, "ONLY");
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
  gMC->Gspos("OQF1", 1, "SRIC", fpParam->OuterFreonWidth()/2 + fpParam->InnerFreonWidth()/2 + 2, 1.276 - fpParam->GapThickness()/2- fpParam->QuartzThickness() -fpParam->FreonThickness()/2, 0., 0, "ONLY"); //Original settings (31.3)
  gMC->Gspos("OQF2", 2, "SRIC", 0., 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()/2, 0., 0, "ONLY");          //Original settings 
  gMC->Gspos("OQF1", 3, "SRIC", - (fpParam->OuterFreonWidth()/2 + fpParam->InnerFreonWidth()/2) - 2, 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness() - fpParam->FreonThickness()/2, 0., 0, "ONLY");       //Original settings (-31.3)
  gMC->Gspos("QUAR", 1, "SRIC", 0., 1.276 - fpParam->GapThickness()/2 - fpParam->QuartzThickness()/2, 0., 0, "ONLY");
  gMC->Gspos("GAP ", 1, "META", 0., fpParam->GapThickness()/2 - fpParam->ProximityGapThickness()/2 - 0.0001, 0., 0, "ONLY");
  gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
  gMC->Gspos("CSI ", 1, "SRIC", 0., 1.276 + fpParam->GapThickness()/2 + .25, 0., 0, "ONLY");
//Wire support placing
  gMC->Gspos("WSG2", 1, "GAP ", 0., fpParam->ProximityGapThickness()/2 - .1, 0., 0, "ONLY");
  gMC->Gspos("WSG1", 1, "CSI ", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("WSMe", 1, "SRIC ", 0., 1.276 + fpParam->GapThickness()/2 + .5 + 1.05, 0., 0, "ONLY");
//Backplane placing
  gMC->Gspos("BACK", 1, "SRIC ", -33.15, 1.276 + fpParam->GapThickness()/2 + .5 + 2.1 + 2, 43.3, 0, "ONLY");
  gMC->Gspos("BACK", 2, "SRIC ", 33.15, 1.276 + fpParam->GapThickness()/2 + .5 + 2.1 + 2 , 43.3, 0, "ONLY");
  gMC->Gspos("BACK", 3, "SRIC ", -33.15, 1.276 + fpParam->GapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
  gMC->Gspos("BACK", 4, "SRIC ", 33.15, 1.276 + fpParam->GapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
  gMC->Gspos("BACK", 5, "SRIC ", 33.15, 1.276 + fpParam->GapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");
  gMC->Gspos("BACK", 6, "SRIC ", -33.15, 1.276 + fpParam->GapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");
//PCB placing
  gMC->Gspos("PCB ", 1, "SRIC ", 0.,  1.276 + fpParam->GapThickness()/2 + .5 + 1.05, csi_width/4 + .5025 + 2.5, 0, "ONLY");
  gMC->Gspos("PCB ", 2, "SRIC ", 0.,  1.276 + fpParam->GapThickness()/2 + .5 + 1.05, -csi_width/4 - .5025 - 2.5, 0, "ONLY");

//place chambers into mother volume ALIC
  for(int i=1;i<=kNCH;i++){
    AliMatrix(idrotm[1000+i],C(i)->ThetaXd(),C(i)->PhiXd(),
                             C(i)->ThetaYd(),C(i)->PhiYd(),
                             C(i)->ThetaZd(),C(i)->PhiZd());
    gMC->Gspos("RICH",i,"ALIC",C(i)->X(),C(i)->Y(),C(i)->Z(),idrotm[1000+i], "ONLY");
  }

  if(GetDebug())Info("CreateGeometry","Stop.");  
}//void AliRICH::CreateGeometry()
//______________________________________________________________________________
void AliRICH::CreateChambers()
{//(re)create all RICH Chambers
  if(GetDebug())Info("CreateChambers","Start.");

  if(fChambers) delete fChambers;//recreate chambers
  fChambers=new TObjArray(kNCH);
  fChambers->SetOwner();
  for(int i=0;i<kNCH;i++){
    fChambers->AddAt(new AliRICHChamber(i+1,fpParam),i);
  }

  if(GetDebug())Info("CreateChambers","Stop.");
}//void AliRICH::CreateChambers()
//__________________________________________________________________________________________________
void AliRICH::GenerateFeedbacks(Float_t eloss)
{// Generate FeedBack photons
  Int_t j;
  Float_t cthf, phif, enfp = 0, sthf;
  Float_t e1[3], e2[3], e3[3];
  Float_t vmod, uswop;
  Float_t dir[3], phi;
  Float_t pol[3], mom[4];
//Determine number of feedback photons
  TLorentzVector x4;
  gMC->TrackPosition(x4);//This sould return the current track position
  Float_t charge=Param()->TotalCharge(gMC->TrackPid(),eloss,5*kcm);//??? Totsl Charge
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
static Int_t sMaxIterPad=0;    // Static variables for the pad-hit iterator routines
static Int_t sCurIterPad=0;

//__________________________________________________________________________________________________
AliRICHSDigit* AliRICH::FirstPad(AliRICHhit*  hit,TClonesArray *clusters ) 
{// Initialise the pad iterator Return the address of the first sdigit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->PHlast() > 0) {
	sMaxIterPad=Int_t(hit->PHlast());
	sCurIterPad=Int_t(hit->PHfirst());
	return (AliRICHSDigit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
    
}
//__________________________________________________________________________________________________
AliRICHSDigit* AliRICH::NextPad(TClonesArray *clusters) 
{// Iterates over pads
  
    sCurIterPad++;
    if (sCurIterPad <= sMaxIterPad) {
	return (AliRICHSDigit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}
//__________________________________________________________________________________________________
