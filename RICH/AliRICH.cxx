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

////////////////////////////////////////////////
//  Manager and hits classes for set:RICH     //
////////////////////////////////////////////////

#include <Riostream.h>
#include <strings.h>

#include <TArrayF.h>
#include <TBRIK.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TH1.h>
#include <TH2.h>
#include <TNode.h> 
#include <TObjArray.h>
#include <TObject.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TRandom.h> 
#include <TStyle.h>
#include <TTUBE.h>
#include <TTree.h>
#include <TVector.h>
#include "AliConst.h"
#include "AliMagF.h"
#include "AliPoints.h"
#include "AliRICH.h"
#include "AliRICHClusterFinder.h"
#include "AliRICHDigit.h"
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
#include "AliSegmentation.h"
#include "AliRICHParam.h"

static Int_t sMaxIterPad=0;    // Static variables for the pad-hit iterator routines
static Int_t sCurIterPad=0;
 
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
//kir  fFileName = 0;
//kir  fMerger = 0;
}//AliRICH::AliRICH()
//______________________________________________________________________________
AliRICH::AliRICH(const char *name, const char *title)
        :AliDetector(name,title)
{//Named ctor
  if(GetDebug())Info("named ctor","Start.");
  fpParam     =new AliRICHParam;
  fHits       =new TClonesArray("AliRICHhit",1000  );
  fCerenkovs  =new TClonesArray("AliRICHCerenkov",1000);
  fSDigits    =new TClonesArray("AliRICHdigit",100000);
  gAlice->AddHitList(fHits);
  gAlice->AddHitList(fCerenkovs);
  fNsdigits   =0;
  fNcerenkovs =0;
  fIshunt     =0;
  fDchambers  =new TObjArray(kNCH);
  fRawClusters=new TObjArray(kNCH);
  fRecHits1D  =new TObjArray(kNCH);
  fRecHits3D  =new TObjArray(kNCH);
  for(int i=0;i<kNCH;i++) {
    fDchambers->AddAt(new TClonesArray("AliRICHDigit",10000), i); 
    fRawClusters->AddAt(new TClonesArray("AliRICHRawCluster",10000), i); 
    fRecHits1D->AddAt(new TClonesArray("AliRICHRecHit1D",1000), i);
    fRecHits3D->AddAt(new TClonesArray("AliRICHRecHit3D",1000), i);
    fNdch[i]=0;
    fNrawch[i]=0;
  }
  SetMarkerColor(kRed);
  fCkovNumber=fFreonProd=0;
//kir  fFileName = 0;
//kir  fMerger = 0;
  if(GetDebug())Info("named ctor","Stop.");
}//AliRICH::AliRICH(const char *name, const char *title)
//______________________________________________________________________________
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
//______________________________________________________________________________
Int_t AliRICH::Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res)
{//calls the charge disintegration method of the current chamber and adds all generated sdigits to the list of digits
   
   Int_t iChamber,iPadX,iPadY,iAdc,iTrack;
   Float_t list[4][500];
   Int_t iNdigits;
        

  ((AliRICHChamber*)fChambers->At(idvol))->DisIntegration(eloss, xhit, yhit, iNdigits, list, res);
    Int_t ic=0;
    
  for(Int_t i=0; i<iNdigits; i++) {
    if(Int_t(list[0][i]) > 0) {
	    ic++;
	    iAdc = Int_t(list[0][i]);
	    iPadX = Int_t(list[1][i]);
	    iPadY = Int_t(list[2][i]);
	    iChamber = Int_t(list[3][i]);

	    
	    AddSDigit(iChamber,iPadX,iPadY,iAdc,iTrack);
	}
    }
    
   if(gAlice->TreeS()){
	gAlice->TreeS()->Fill();
	gAlice->TreeS()->Write(0,TObject::kOverwrite);
   }
   return iNdigits;
}//Int_t AliRICH::Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res)
//______________________________________________________________________________
void AliRICH::Hits2SDigits()
{//Create a list of sdigits corresponding to list of hits. Every hit generates sdigit.
  if(GetDebug()) Info("Hit2SDigits","Start.");
  
  for(Int_t iEventN=0;iEventN<gAlice->GetEventsPerRun();iEventN++){//loop on events
    fLoader->GetRunLoader()->GetEvent(iEventN);
  
    if(!fLoader->TreeH()) fLoader->LoadHits();
    if(!fLoader->TreeS()) fLoader->MakeTree("S");
    MakeBranch("S");
  
    for(int iPrimN=0;iPrimN<TreeH()->GetEntries();iPrimN++){//loop on primary tracks
      fLoader->TreeH()->GetEntry(iPrimN); 
      for(Int_t iHitN=0;iHitN<Hits()->GetEntries();iHitN++){//loop on hits for given primary track  
        AddSDigit(4,13,24,55,4);//chamber-xpad-ypad-qdc-track1-2-3
      }//loop on hits for given primary track
    }//loop on primary tracks
  
    fLoader->TreeS()->Fill();
    fLoader->WriteSDigits("OVERWRITE");
  }//loop on events
  
  if(GetDebug()) Info("Hit2SDigits","Stop.");
}
//______________________________________________________________________________
void AliRICH::SDigits2Digits()
{//Generate digits from sdigits.
  if(GetDebug()) Info("SDigits2Digits","Start.");
   //AliRICHChamber*       iChamber;
  
  
   //for(Int_t i=0;i<7;i++) {
   //iChamber = &(Chamber(i));
   //iChamber->GenerateTresholds();
   //}
  
   //int nparticles = gAlice->GetNtrack();
   //cout << "Particles (RICH):" <<nparticles<<endl;
   //if (nparticles <= 0) return;
   //if (!fMerger) {
   //fMerger = new AliRICHMerger();
   //}


   //fMerger->Init();
   //fMerger->Digitise(nev,flag);

   AliRunDigitizer * manager = new AliRunDigitizer(1,1);
   manager->SetInputStream(0,"galice.root");
   //AliRICHDigitizer *dRICH  = new AliRICHDigitizer(manager);
   manager->Exec("deb");
  if(GetDebug()) Info("SDigits2Digits","Stop.");
}//void AliRICH::SDigits2Digits()
//______________________________________________________________________________
void AliRICH::Digits2Reco()
{
// Generate clusters
// Called from alirun, single event only.     
  if(GetDebug()) Info("Digits2Reco","Start.");

  int nparticles = gAlice->GetNtrack();
  cout << "Particles (RICH):" <<nparticles<<endl;
  if (nparticles > 0) FindClusters(0);

}//void AliRICH::Digits2Reco()  


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

  Float_t wid=fpParam->PadPlaneWidth();
  Float_t len=fpParam->PadPlaneLength();
  new TBRIK("PHOTO","PHOTO","void",wid/2,0.1,len/2);
  
  for(int i=0;i<kNCH;i++){
    top->cd();
    node = new TNode(Form("RICH%i",i+1),Form("RICH%i",i+1),"S_RICH",C(i)->X(),C(i)->Y(),C(i)->Z(),C(i)->RotMatrixName());
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
//______________________________________________________________________________
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
}

//__________________________________________
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



//___________________________________________
//____________________________________________
void AliRICH::ResetDigits()
{//Reset number of digits and the digits array for this detector
  for ( int i=0;i<kNCH;i++ ) {
    if (fDchambers && fDchambers->At(i))   fDchambers->At(i)->Clear();
    if (fNdch)  fNdch[i]=0;
  }
}
//____________________________________________
void AliRICH::ResetRawClusters()
{//Reset number of raw clusters and the raw clust array for this detector
  for ( int i=0;i<kNCH;i++ ) {
    if (fRawClusters->At(i))    ((TClonesArray*)fRawClusters->At(i))->Clear();
    if (fNrawch)  fNrawch[i]=0;
  }
}
//____________________________________________
void AliRICH::ResetRecHits1D()
{//Reset number of raw clusters and the raw clust array for this detector
  for ( int i=0;i<kNCH;i++ ) {
    if (fRecHits1D->At(i))    ((TClonesArray*)fRecHits1D->At(i))->Clear();
    if (fNrechits1D)  fNrechits1D[i]=0;
  }
}

//____________________________________________
void AliRICH::ResetRecHits3D()
{// Reset number of raw clusters and the raw clust array for this detector
  for ( int i=0;i<kNCH;i++ ) {
    if (fRecHits3D->At(i))    ((TClonesArray*)fRecHits3D->At(i))->Clear();
    if (fNrechits3D)  fNrechits3D[i]=0;
  }
}
//______________________________________________________________________________
void AliRICH::FindClusters(Int_t nev /*kir,Int_t lastEntry*/)
{// Loop on chambers and on cathode planes
    for (Int_t icat=1;icat<2;icat++) {
	gAlice->ResetDigits();
	gAlice->TreeD()->GetEvent(0);
	for (Int_t ich=0;ich<kNCH;ich++) {
      //PH	  AliRICHChamber* iChamber=(AliRICHChamber*) (*fChambers)[ich];
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
//______________________________________________________________________________
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
//______________________________________________________________________________
AliRICHSDigit* AliRICH::NextPad(TClonesArray *clusters) 
{// Iterates over pads
  
    sCurIterPad++;
    if (sCurIterPad <= sMaxIterPad) {
	return (AliRICHSDigit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}


void AliRICH::DiagnosticsFE(Int_t evNumber1,Int_t evNumber2)
{
  
  Int_t NpadX = 162;                 // number of pads on X
  Int_t NpadY = 162;                 // number of pads on Y
  
  Int_t Pad[162][162];
  for (Int_t i=0;i<NpadX;i++) {
    for (Int_t j=0;j<NpadY;j++) {
      Pad[i][j]=0;
    }
  }
  
  //  Create some histograms

  TH1F *pionspectra1 = new TH1F("pionspectra1","Pion Spectra",200,-4,2);
  TH1F *pionspectra2 = new TH1F("pionspectra2","Pion Spectra",200,-4,2);
  TH1F *pionspectra3 = new TH1F("pionspectra3","Pion Spectra",200,-4,2);
  TH1F *protonspectra1 = new TH1F("protonspectra1","Proton Spectra",200,-4,2);
  TH1F *protonspectra2 = new TH1F("protonspectra2","Proton Spectra",200,-4,2);
  TH1F *protonspectra3 = new TH1F("protonspectra3","Proton Spectra",200,-4,2);
  TH1F *kaonspectra1 = new TH1F("kaonspectra1","Kaon Spectra",100,-4,2);
  TH1F *kaonspectra2 = new TH1F("kaonspectra2","Kaon Spectra",100,-4,2);
  TH1F *kaonspectra3 = new TH1F("kaonspectra3","Kaon Spectra",100,-4,2);
  TH1F *electronspectra1 = new TH1F("electronspectra1","Electron Spectra",100,-4,2);
  TH1F *electronspectra2 = new TH1F("electronspectra2","Electron Spectra",100,-4,2);
  TH1F *electronspectra3 = new TH1F("electronspectra3","Electron Spectra",100,-4,2);
  TH1F *muonspectra1 = new TH1F("muonspectra1","Muon Spectra",100,-4,2);
  TH1F *muonspectra2 = new TH1F("muonspectra2","Muon Spectra",100,-4,2);
  TH1F *muonspectra3 = new TH1F("muonspectra3","Muon Spectra",100,-4,2);
  TH1F *neutronspectra1 = new TH1F("neutronspectra1","Neutron Spectra",100,-4,2);
  TH1F *neutronspectra2 = new TH1F("neutronspectra2","Neutron Spectra",100,-4,2);
  TH1F *neutronspectra3 = new TH1F("neutronspectra2","Neutron Spectra",100,-4,2);
  TH1F *chargedspectra1 = new TH1F("chargedspectra1","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *chargedspectra2 = new TH1F("chargedspectra2","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *chargedspectra3 = new TH1F("chargedspectra2","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *pionptspectrafinal = new TH1F("pionptspectrafinal","Primary Pions Transverse Momenta at HMPID",20,0,5);
  TH1F *pionptspectravertex = new TH1F("pionptspectravertex","Primary Pions Transverse Momenta at vertex",20,0,5);
  TH1F *kaonptspectrafinal = new TH1F("kaonptspectrafinal","Primary Kaons Transverse Momenta at HMPID",20,0,5);
  TH1F *kaonptspectravertex = new TH1F("kaonptspectravertex","Primary Kaons Transverse Momenta at vertex",20,0,5);
  //TH1F *hitsPhi = new TH1F("hitsPhi","Distribution of phi angle of incidence",100,-180,180);
  TH1F *hitsTheta = new TH1F("hitsTheta","Distribution of Theta angle of incidence, all tracks",100,0,50);
  TH1F *hitsTheta500MeV = new TH1F("hitsTheta500MeV","Distribution of Theta angle of incidence, 0.5-1 GeV primary tracks",100,0,50);
  TH1F *hitsTheta1GeV = new TH1F("hitsTheta1GeV","Distribution of Theta angle of incidence, 1-2 GeV primary tracks",100,0,50);
  TH1F *hitsTheta2GeV = new TH1F("hitsTheta2GeV","Distribution of Theta angle of incidence, 2-3 GeV primary tracks",100,0,50);
  TH1F *hitsTheta3GeV = new TH1F("hitsTheta3GeV","Distribution of Theta angle of incidence, >3 GeV primary tracks",100,0,50);
  TH2F *production = new TH2F("production","Mother production vertices",100,-300,300,100,0,600);
   
   
   

//   Start loop over events 

  Int_t pion=0, kaon=0, proton=0, electron=0, positron=0, neutron=0, highneutrons=0, muon=0;
  Int_t chargedpions=0,primarypions=0,highprimarypions=0,chargedkaons=0,primarykaons=0,highprimarykaons=0;
  Int_t photons=0, primaryphotons=0, highprimaryphotons=0;
  TRandom* random=0;

   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       

       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
// Get pointers to RICH detector and Hits containers
       
       AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
     
       TTree *treeH = TreeH();
       Int_t ntracks =(Int_t) treeH->GetEntries();
            
// Start loop on tracks in the hits containers
       
       for (Int_t track=0; track<ntracks;track++) {
	   printf ("Processing Track: %d\n",track);
	   gAlice->ResetHits();
	   treeH->GetEvent(track);
	   	   	   
	   for(AliRICHhit* mHit=(AliRICHhit*)pRICH->FirstHit(-1); 
	       mHit;
	       mHit=(AliRICHhit*)pRICH->NextHit()) 
	     {
	       //Int_t nch  = mHit->fChamber;              // chamber number
	       //Float_t x  = mHit->X();                    // x-pos of hit
	       //Float_t y  = mHit->Z();                    // y-pos
	       //Float_t z  = mHit->Y();
	       //Float_t phi = mHit->Phi();                 //Phi angle of incidence
	       Float_t theta = mHit->Theta();             //Theta angle of incidence
	       Float_t px = mHit->MomX();
	       Float_t py = mHit->MomY();
	       Int_t index = mHit->Track();
	       Int_t particle = (Int_t)(mHit->Particle());    
	       Float_t R;
	       Float_t PTfinal;
	       Float_t PTvertex;

	      TParticle *current = gAlice->Particle(index);
	      
	      //Float_t energy=current->Energy(); 

	      R=TMath::Sqrt(current->Vx()*current->Vx() + current->Vy()*current->Vy());
	      PTfinal=TMath::Sqrt(px*px + py*py);
	      PTvertex=TMath::Sqrt(current->Px()*current->Px() + current->Py()*current->Py());
	      
	      

	      if (TMath::Abs(particle) < 10000000)
		{
		  hitsTheta->Fill(theta,(float) 1);
		  if (R<5)
		    {
		      if (PTvertex>.5 && PTvertex<=1)
			{
			  hitsTheta500MeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>1 && PTvertex<=2)
			{
			  hitsTheta1GeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>2 && PTvertex<=3)
			{
			  hitsTheta2GeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>3)
			{
			  hitsTheta3GeV->Fill(theta,(float) 1);
			}
		    }
		  
		}

	      //if (nch == 3)
		//{
	      
	      if (TMath::Abs(particle) < 50000051)
		{
		  //if (TMath::Abs(particle) == 50000050 || TMath::Abs(particle) == 2112)
		  if (TMath::Abs(particle) == 2112 || TMath::Abs(particle) == 50000050)
		    {
		      //gMC->Rndm(&random, 1);
		      if (random->Rndm() < .1)
			production->Fill(current->Vz(),R,(float) 1);
		      if (TMath::Abs(particle) == 50000050)
			//if (TMath::Abs(particle) > 50000000)
			{
			  photons +=1;
			  if (R<5)
			    {
			      primaryphotons +=1;
			      if (current->Energy()>0.001)
				highprimaryphotons +=1;
			    }
			}	
		      if (TMath::Abs(particle) == 2112)
			{
			  neutron +=1;
			  if (current->Energy()>0.0001)
			    highneutrons +=1;
			}
		    }
		  if (TMath::Abs(particle) < 50000000)
		    {
		      production->Fill(current->Vz(),R,(float) 1);
		    }
		  //mip->Fill(x,y,(float) 1);
		}
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  if (R<5)
		    {
		      pionptspectravertex->Fill(PTvertex,(float) 1);
		      pionptspectrafinal->Fill(PTfinal,(float) 1);
		    }
		}
	      
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  if (R<5)
		    {
		      kaonptspectravertex->Fill(PTvertex,(float) 1);
		      kaonptspectrafinal->Fill(PTfinal,(float) 1);
		    }
		}
	      
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  pionspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    pionspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    {
		      pionspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		  pion +=1;
		  if (TMath::Abs(particle)==211)
		    {
		      chargedpions +=1;
		      if (R<5)
			{
			  primarypions +=1;
			  if (current->Energy()>1)
			    highprimarypions +=1;
			}
		    }	
		}
	      if (TMath::Abs(particle)==2212)
		{
		  protonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    protonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    protonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  proton +=1;
		}
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  kaonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    kaonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    kaonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  kaon +=1;
		  if (TMath::Abs(particle)==321)
		    {
		      chargedkaons +=1;
		      if (R<5)
			{
			  primarykaons +=1;
			  if (current->Energy()>1)
			    highprimarykaons +=1;
			}
		    }
		}
	      if (TMath::Abs(particle)==11)
		{
		  electronspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    electronspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    electronspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (particle == 11)
		    electron +=1;
		  if (particle == -11)
		    positron +=1;
		}
	      if (TMath::Abs(particle)==13)
		{
		  muonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    muonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    muonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  muon +=1;
		}
	      if (TMath::Abs(particle)==2112)
		{
		  neutronspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    neutronspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    {
		      neutronspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		  neutron +=1;
		}
	      if(TMath::Abs(particle)==211 || TMath::Abs(particle)==2212 || TMath::Abs(particle)==321)
		{
		  if (current->Energy()-current->GetCalcMass()>1)
		    {
		      chargedspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
			chargedspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      if (R>250 && R<450)
			chargedspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		}
	      // Fill the histograms
	      //Nh1+=nhits;
	      //h->Fill(x,y,(float) 1);
	      //}
	      //}
	   }          
	   
       }
       
   }
   //   }

   TStyle *mystyle=new TStyle("Plain","mystyle");
   mystyle->SetPalette(1,0);
   mystyle->cd();
   
   //Create canvases, set the view range, show histograms

    TCanvas *c2 = new TCanvas("c2","Angles of incidence",150,150,100,150);
    c2->Divide(2,2);
    //c2->SetFillColor(42);
    
    c2->cd(1);
    hitsTheta500MeV->SetFillColor(5);
    hitsTheta500MeV->Draw();
    c2->cd(2);
    hitsTheta1GeV->SetFillColor(5);
    hitsTheta1GeV->Draw();
    c2->cd(3);
    hitsTheta2GeV->SetFillColor(5);
    hitsTheta2GeV->Draw();
    c2->cd(4);
    hitsTheta3GeV->SetFillColor(5);
    hitsTheta3GeV->Draw();
    
            
   
    TCanvas *c15 = new TCanvas("c15","Mothers Production Vertices",50,50,600,600);
    c15->cd();
    production->SetFillColor(42);
    production->SetXTitle("z (m)");
    production->SetYTitle("R (m)");
    production->Draw();

    TCanvas *c10 = new TCanvas("c10","Pt Spectra",50,50,600,700);
    c10->Divide(2,2);
    c10->cd(1);
    pionptspectravertex->SetFillColor(5);
    pionptspectravertex->SetXTitle("Pt (GeV)");
    pionptspectravertex->Draw();
    c10->cd(2);
    pionptspectrafinal->SetFillColor(5);
    pionptspectrafinal->SetXTitle("Pt (GeV)");
    pionptspectrafinal->Draw();
    c10->cd(3);
    kaonptspectravertex->SetFillColor(5);
    kaonptspectravertex->SetXTitle("Pt (GeV)");
    kaonptspectravertex->Draw();
    c10->cd(4);
    kaonptspectrafinal->SetFillColor(5);
    kaonptspectrafinal->SetXTitle("Pt (GeV)");
    kaonptspectrafinal->Draw();
   
  
   TCanvas *c16 = new TCanvas("c16","Particles Spectra II",150,150,600,350);
   c16->Divide(2,1);
   
   c16->cd(1);
   //TCanvas *c13 = new TCanvas("c13","Electron Spectra",400,10,600,700);
   electronspectra1->SetFillColor(5);
   electronspectra1->SetXTitle("log(GeV)");
   electronspectra2->SetFillColor(46);
   electronspectra2->SetXTitle("log(GeV)");
   electronspectra3->SetFillColor(10);
   electronspectra3->SetXTitle("log(GeV)");
   //c13->SetLogx();
   electronspectra1->Draw();
   electronspectra2->Draw("same");
   electronspectra3->Draw("same");
   
   c16->cd(2);
   //TCanvas *c14 = new TCanvas("c14","Muon Spectra",400,10,600,700);
   muonspectra1->SetFillColor(5);
   muonspectra1->SetXTitle("log(GeV)");
   muonspectra2->SetFillColor(46);
   muonspectra2->SetXTitle("log(GeV)");
   muonspectra3->SetFillColor(10);
   muonspectra3->SetXTitle("log(GeV)");
   //c14->SetLogx();
   muonspectra1->Draw();
   muonspectra2->Draw("same");
   muonspectra3->Draw("same");
   
   //c16->cd(3);
   //TCanvas *c16 = new TCanvas("c16","Neutron Spectra",400,10,600,700);
   //neutronspectra1->SetFillColor(42);
   //neutronspectra1->SetXTitle("log(GeV)");
   //neutronspectra2->SetFillColor(46);
   //neutronspectra2->SetXTitle("log(GeV)");
   //neutronspectra3->SetFillColor(10);
   //neutronspectra3->SetXTitle("log(GeV)");
   //c16->SetLogx();
   //neutronspectra1->Draw();
   //neutronspectra2->Draw("same");
   //neutronspectra3->Draw("same");

   TCanvas *c9 = new TCanvas("c9","Particles Spectra",150,150,600,700);
   //TCanvas *c9 = new TCanvas("c9","Pion Spectra",400,10,600,700);
   c9->Divide(2,2);
   
   c9->cd(1);
   pionspectra1->SetFillColor(5);
   pionspectra1->SetXTitle("log(GeV)");
   pionspectra2->SetFillColor(46);
   pionspectra2->SetXTitle("log(GeV)");
   pionspectra3->SetFillColor(10);
   pionspectra3->SetXTitle("log(GeV)");
   //c9->SetLogx();
   pionspectra1->Draw();
   pionspectra2->Draw("same");
   pionspectra3->Draw("same");
   
   c9->cd(2);
   //TCanvas *c10 = new TCanvas("c10","Proton Spectra",400,10,600,700);
   protonspectra1->SetFillColor(5);
   protonspectra1->SetXTitle("log(GeV)");
   protonspectra2->SetFillColor(46);
   protonspectra2->SetXTitle("log(GeV)");
   protonspectra3->SetFillColor(10);
   protonspectra3->SetXTitle("log(GeV)");
   //c10->SetLogx();
   protonspectra1->Draw();
   protonspectra2->Draw("same");
   protonspectra3->Draw("same");
   
   c9->cd(3);
   //TCanvas *c11 = new TCanvas("c11","Kaon Spectra",400,10,600,700); 
   kaonspectra1->SetFillColor(5);
   kaonspectra1->SetXTitle("log(GeV)");
   kaonspectra2->SetFillColor(46);
   kaonspectra2->SetXTitle("log(GeV)");
   kaonspectra3->SetFillColor(10);
   kaonspectra3->SetXTitle("log(GeV)");
   //c11->SetLogx();
   kaonspectra1->Draw();
   kaonspectra2->Draw("same");
   kaonspectra3->Draw("same");
   
   c9->cd(4);
   //TCanvas *c12 = new TCanvas("c12","Charged Particles Spectra",400,10,600,700);
   chargedspectra1->SetFillColor(5);
   chargedspectra1->SetXTitle("log(GeV)");
   chargedspectra2->SetFillColor(46);
   chargedspectra2->SetXTitle("log(GeV)");
   chargedspectra3->SetFillColor(10);
   chargedspectra3->SetXTitle("log(GeV)");
   //c12->SetLogx();
   chargedspectra1->Draw();
   chargedspectra2->Draw("same");
   chargedspectra3->Draw("same");
   


   printf("*****************************************\n");
   printf("* Particle                   *  Counts  *\n");
   printf("*****************************************\n");

   printf("* Pions:                     *   %4d   *\n",pion);
   printf("* Charged Pions:             *   %4d   *\n",chargedpions);
   printf("* Primary Pions:             *   %4d   *\n",primarypions);
   printf("* Primary Pions (p>1GeV/c):  *   %4d   *\n",highprimarypions);
   printf("* Kaons:                     *   %4d   *\n",kaon);
   printf("* Charged Kaons:             *   %4d   *\n",chargedkaons);
   printf("* Primary Kaons:             *   %4d   *\n",primarykaons);
   printf("* Primary Kaons (p>1GeV/c):  *   %4d   *\n",highprimarykaons);
   printf("* Muons:                     *   %4d   *\n",muon);
   printf("* Electrons:                 *   %4d   *\n",electron);
   printf("* Positrons:                 *   %4d   *\n",positron);
   printf("* Protons:                   *   %4d   *\n",proton);
   printf("* All Charged:               *   %4d   *\n",(chargedpions+chargedkaons+muon+electron+positron+proton));
   printf("*****************************************\n");
   //printf("* Photons:                   *   %3.1f   *\n",photons); 
   //printf("* Primary Photons:           *   %3.1f   *\n",primaryphotons);
   //printf("* Primary Photons (p>1MeV/c):*   %3.1f   *\n",highprimaryphotons);
   //printf("*****************************************\n");
   //printf("* Neutrons:                  *   %3.1f   *\n",neutron);
   //printf("* Neutrons (p>100keV/c):     *   %3.1f   *\n",highneutrons);
   //printf("*****************************************\n");

   if (gAlice->TreeD())
     {
       gAlice->TreeD()->GetEvent(0);
   
       Float_t occ[7]; 
       Float_t sum=0;
       Float_t mean=0; 
       printf("\n*****************************************\n");
       printf("* Chamber   * Digits      * Occupancy   *\n");
       printf("*****************************************\n");
       
       for (Int_t ich=0;ich<7;ich++)
	 {
	   TClonesArray *Digits = DigitsAddress(ich);    //  Raw clusters branch
	   Int_t ndigits = Digits->GetEntriesFast();
	   occ[ich] = Float_t(ndigits)/(160*144);
	   sum += Float_t(ndigits)/(160*144);
	   printf("*   %d      *    %d      *   %3.1f%%     *\n",ich,ndigits,occ[ich]*100);
	 }
       mean = sum/7;
       printf("*****************************************\n");
       printf("* Mean occupancy          *   %3.1f%%     *\n",mean*100);
       printf("*****************************************\n");
     }
 
  printf("\nEnd of analysis\n");
   
}//void AliRICH::DiagnosticsFE(Int_t evNumber1,Int_t evNumber2)
//______________________________________________________________________________
void AliRICH::DiagnosticsSE(Int_t diaglevel,Int_t evNumber1,Int_t evNumber2)
{

AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
   AliRICHSegmentationV0*  segmentation;
   AliRICHChamber*       chamber;
   
   chamber = &(pRICH->Chamber(0));
   segmentation=(AliRICHSegmentationV0*) chamber->GetSegmentationModel();

   Int_t NpadX = segmentation->Npx();                 // number of pads on X
   Int_t NpadY = segmentation->Npy();                 // number of pads on Y
    
   //Int_t Pad[144][160];
   /*for (Int_t i=0;i<NpadX;i++) {
     for (Int_t j=0;j<NpadY;j++) {
       Pad[i][j]=0;
     }
   } */


   Int_t xmin= -NpadX/2;  
   Int_t xmax=  NpadX/2;
   Int_t ymin= -NpadY/2;
   Int_t ymax=  NpadY/2;

   Float_t PTfinal = 0;
   Int_t pionCount = 0;
   Int_t kaonCount = 0;
   Int_t protonCount = 0;
   
   TH2F *feedback = 0;
   TH2F *mip = 0;
   TH2F *cerenkov = 0;
   TH2F *h = 0;
   TH1F *hitsX = 0;
   TH1F *hitsY = 0;

   TH2F *hc0 = new TH2F("hc0","Zoom on center of central chamber",150,-25,25,150,-45,5);

   if (diaglevel == 1)
     {
       printf("Single Ring Hits\n");
       feedback = new TH2F("feedback","Feedback hit distribution",150,-20,20,150,-35,5);
       mip = new TH2F("mip","Mip hit distribution",150,-20,20,150,-35,5);
       cerenkov = new TH2F("cerenkov","Cerenkov hit distribution",150,-20,20,150,-35,5);
       h = new TH2F("h","Detector hit distribution",150,-20,20,150,-35,5);
       hitsX = new TH1F("hitsX","Distribution of hits along x-axis",150,-50,50);
       hitsY = new TH1F("hitsY","Distribution of hits along z-axis",150,-50,50);
     }       
   else
     {
       printf("Full Event Hits\n");
       
       feedback = new TH2F("feedback","Feedback hit distribution",150,-300,300,150,-300,300);
       mip = new TH2F("mip","Mip hit distribution",150,-300,300,150,-300,300);
       cerenkov = new TH2F("cerenkov","Cerenkov hit distribution",150,-300,300,150,-300,300);
       h = new TH2F("h","Detector hit distribution",150,-300,300,150,-300,300); 
       hitsX = new TH1F("digitsX","Distribution of hits along x-axis",200,-300,300);
       hitsY = new TH1F("digitsY","Distribution of hits along z-axis",200,-300,300);
     }
   


   TH2F *hc1 = new TH2F("hc1","Chamber 1 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc2 = new TH2F("hc2","Chamber 2 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc3 = new TH2F("hc3","Chamber 3 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc4 = new TH2F("hc4","Chamber 4 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc5 = new TH2F("hc5","Chamber 5 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc6 = new TH2F("hc6","Chamber 6 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc7 = new TH2F("hc7","Chamber 7 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
      
   TH1F *Clcharge = new TH1F("Clcharge","Cluster Charge Distribution",500,0.,500.);
   TH1F *ckovangle = new TH1F("ckovangle","Cerenkov angle per photon",100,.35,.8);
   TH1F *hckphi = new TH1F("hckphi","Cerenkov phi angle per photon",620,-3.1,3.1);
   TH1F *mother = new TH1F("mother","Cerenkovs per Mip",75,0.,75.);
   TH1F *radius = new TH1F("radius","Mean distance to Mip",100,0.,20.);
   TH1F *phspectra1 = new TH1F("phspectra1","Detected Photon Spectra",200,5.,10.);
   TH1F *phspectra2 = new TH1F("phspectra2","Produced Photon Spectra",200,5.,10.);
   TH1F *totalphotonstrack = new TH1F("totalphotonstrack","Produced Photons per Mip",100,200,700.);
   TH1F *totalphotonsevent = new TH1F("totalphotonsevent","Produced Photons per Mip",100,200,700.);
   //TH1F *feedbacks = new TH1F("feedbacks","Produced Feedbacks per Mip",50,0.5,50.);
   TH1F *padnumber = new TH1F("padnumber","Number of pads per cluster",50,-0.5,50.);
   TH1F *padsev = new TH1F("padsev","Number of pads hit per MIP",50,0.5,100.);
   TH1F *clusev = new TH1F("clusev","Number of clusters per MIP",50,0.5,50.);
   TH1F *photev = new TH1F("photev","Number of detected photons per MIP",50,0.5,50.);
   TH1F *feedev = new TH1F("feedev","Number of feedbacks per MIP",50,0.5,50.);
   TH1F *padsmip = new TH1F("padsmip","Number of pads per event inside MIP region",50,0.5,50.);
   TH1F *padscl = new TH1F("padscl","Number of pads per event from cluster count",50,0.5,100.);
   TH1F *pionspectra = new TH1F("pionspectra","Pion Spectra",200,.5,10.);
   TH1F *protonspectra = new TH1F("protonspectra","Proton Spectra",200,.5,10.);
   TH1F *kaonspectra = new TH1F("kaonspectra","Kaon Spectra",100,.5,10.);
   TH1F *chargedspectra = new TH1F("chargedspectra","Charged particles above 1 GeV Spectra",100,.5,10.);
   TH1F *hitsPhi = new TH1F("hitsPhi","Distribution of phi angle of incidence",50,0,360);
   TH1F *hitsTheta = new TH1F("hitsTheta","Distribution of theta angle of incidence",50,0,15);
   TH1F *Omega1D = new TH1F("omega","Reconstructed Cerenkov angle per track",50,.5,1);
   TH1F *Theta = new TH1F("theta","Reconstructed theta incidence angle per track",100,0,15);
   TH1F *Phi = new TH1F("phi","Reconstructed phi incidence per track",100,0,360);
   TH1F *Omega3D = new TH1F("omega","Reconstructed Cerenkov angle per track",100,.35,.8);
   TH1F *PhotonCer = new TH1F("photoncer","Reconstructed Cerenkov angle per photon",100,.35,.8);
   TH2F *PadsUsed = new TH2F("padsused","Pads Used for Reconstruction",100,-30,30,100,-30,30);
   TH1F *MeanRadius = new TH1F("radius","Mean Radius for reconstructed track",100,0.,20.);
   TH2F *identification = new TH2F("identification","Particle Identification",100,1,5,100,0,.8);
   TH1F *OriginalOmega = new TH1F("Original Omega","Cerenkov angle per track",100,.35,.8);
   TH1F *OriginalPhi = new TH1F("Original Phi","Distribution of phi angle of incidence per track",100,0,360);
   TH1F *OriginalTheta = new TH1F("Original Theta","Distribution of theta angle per track",100,0,15);
   TH1F *OmegaError = new TH1F("Omega Error","Difference between original an reconstructed cerenkov angle",100,0,.2);
   TH1F *PhiError = new TH1F("Phi Error","Difference between original an reconstructed phi angle",100,0,360);
   TH1F *ThetaError = new TH1F("Theta Error","Difference between original an reconstructed phi angle",100,0,15);


//   Start loop over events 

   Int_t Nh=0;
   Int_t pads=0;
   Int_t Nh1=0;
   Int_t mothers[80000];
   Int_t mothers2[80000];
   Float_t mom[3];
   Int_t nraw=0;
   Int_t phot=0;
   Int_t feed=0;
   Int_t padmip=0;
   Float_t x=0,y=0;

   Float_t chiSquareOmega = 0;
   Float_t chiSquareTheta = 0;
   Float_t chiSquarePhi = 0;

   Float_t recEffEvent = 0;
   Float_t recEffTotal = 0;

   Float_t trackglob[3];
   Float_t trackloc[3];

   
   for (Int_t i=0;i<100;i++) mothers[i]=0;

   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       

       //cout<<"nev  "<<nev<<endl;
       printf ("\n**********************************\nProcessing Event: %d\n",nev);
       //cout<<"nparticles  "<<nparticles<<endl;
       printf ("Particles       : %d\n\n",nparticles);
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
// Get pointers to RICH detector and Hits containers
       

       TTree *TH = TreeH(); 
       Stat_t ntracks = TH->GetEntries();

       // Start loop on tracks in the hits containers
       //Int_t Nc=0;
       for (Int_t track=0; track<ntracks;track++) {
	   
	 printf ("\nProcessing Track: %d\n",track);
	 gAlice->ResetHits();
	 TH->GetEvent(track);
	 Int_t nhits = pRICH->Hits()->GetEntriesFast();
	 if (nhits) Nh+=nhits;
	 printf("Hits            : %d\n",nhits);
	 for(AliRICHhit* mHit=(AliRICHhit*)pRICH->FirstHit(-1); 
	     mHit;
	     mHit=(AliRICHhit*)pRICH->NextHit()) 
	   {
	     Int_t nch  = mHit->Chamber();              // chamber number
	     trackglob[0] = mHit->X();                 // x-pos of hit
	     trackglob[1] = mHit->Y();
	     trackglob[2] = mHit->Z();                 // y-pos of hit
	     //x  = mHit->X();                           // x-pos of hit
	     //y  = mHit->Z();                           // y-pos
	     Float_t phi = mHit->Phi();                 //Phi angle of incidence
	     Float_t theta = mHit->Theta();             //Theta angle of incidence
	     Int_t index = mHit->Track();
	     Int_t particle = (Int_t)(mHit->Particle());        
	     //Int_t freon = (Int_t)(mHit->fLoss);    
	     Float_t px = mHit->MomX();
	     Float_t py = mHit->MomY();
	     
	     if (TMath::Abs(particle) < 10000000)
	       {
		 PTfinal=TMath::Sqrt(px*px + py*py);
	       }
	
	     chamber = &(pRICH->Chamber(nch-1));
	     
	     
	     chamber->GlobaltoLocal(trackglob,trackloc);
	     
	     chamber->LocaltoGlobal(trackloc,trackglob);
	     
       
	     x=trackloc[0];
	     y=trackloc[2];
	     
	     hitsX->Fill(x,(float) 1);
	     hitsY->Fill(y,(float) 1);
	       
	      
	      TParticle *current = (TParticle*)gAlice->Particle(index);
	      //printf("Particle type: %d\n",sizeoff(Particles));

	      hitsTheta->Fill(theta,(float) 1);
	      //hitsPhi->Fill(phi,(float) 1);
	      //if (pRICH->GetDebugLevel() == -1)
	     
	      if (current->GetPdgCode() < 10000000)
		{
		  mip->Fill(x,y,(float) 1);
		  //printf("adding mip\n");
		  //if (current->Energy() - current->GetCalcMass()>1 && freon==1)
		  //{
		  hitsPhi->Fill(TMath::Abs(phi),(float) 1);
		  //hitsTheta->Fill(theta,(float) 1);
		  //printf("Theta:%f, Phi:%f\n",theta,phi);
		  //}
		}
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  pionspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if (TMath::Abs(particle)==2212)
		{
		  protonspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  kaonspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if(TMath::Abs(particle)==211 || TMath::Abs(particle)==2212 || TMath::Abs(particle)==321)
		{
		  if (current->Energy() - current->GetCalcMass()>1)
		    chargedspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      //printf("Hits:%d\n",hit);
	      //printf ("Chamber number:%d x:%f y:%f\n",nch,x,y);
              // Fill the histograms
	      Nh1+=nhits;
	      h->Fill(x,y,(float) 1);
		  //}
              //}
	   }
	   
	   Int_t ncerenkovs = pRICH->Cerenkovs()->GetEntriesFast();
	   //if (current->GetPdgCode() < 50000051 && current->GetPdgCode() > 50000040)
	   //totalphotonsevent->Fill(ncerenkovs,(float) 1);

 	   if (ncerenkovs) {
	     printf("Cerenkovs       : %d\n",ncerenkovs);
	     totalphotonsevent->Fill(ncerenkovs,(float) 1);
	     for (Int_t hit=0;hit<ncerenkovs;hit++) {
	       AliRICHCerenkov* cHit = (AliRICHCerenkov*) pRICH->Cerenkovs()->UncheckedAt(hit);
	       Int_t nchamber = cHit->fChamber;     // chamber number
	       Int_t index =    cHit->Track();
	       //Int_t pindex =   (Int_t)(cHit->fIndex);
	       trackglob[0] = cHit->X();                 // x-pos of hit
	       trackglob[1] = cHit->Y();
	       trackglob[2] = cHit->Z();                 // y-pos of hit
	       //Float_t cx  =      cHit->X();                // x-position
	       //Float_t cy  =      cHit->Z();                // y-position
	       Int_t cmother =  cHit->fCMother;      // Index of mother particle
	       Int_t closs =    (Int_t)(cHit->fLoss);           // How did the particle get lost? 
	       Float_t cherenkov = cHit->fCerenkovAngle;   //production cerenkov angle
	       
	       chamber = &(pRICH->Chamber(nchamber-1));
	     
	       //printf("Nch:%d\n",nch);
	       
	       chamber->GlobaltoLocal(trackglob,trackloc);
	     
	       chamber->LocaltoGlobal(trackloc,trackglob);
	     
       
	       Float_t cx=trackloc[0];
	       Float_t cy=trackloc[2];
	       
	       //printf ("Cerenkov hit number %d/%d, X:%f, Y:%f\n",hit,ncerenkovs,cx,cy); 


	       //printf("Particle:%9d\n",index);
		 		 
	       TParticle *current = (TParticle*)gAlice->Particle(index);
	       Float_t energyckov = current->Energy();
	       
	       if (current->GetPdgCode() == 50000051)
		 {
		   if (closs==4)
		     {
		       feedback->Fill(cx,cy,(float) 1);
		       feed++;
		     }
		 }
	       if (current->GetPdgCode() == 50000050)
		 {
		   
		   if (closs !=4)
		     {
		       phspectra2->Fill(energyckov*1e9,(float) 1);
		     }
		       
		   if (closs==4)
		     {
		       cerenkov->Fill(cx,cy,(float) 1); 
		       
		       //printf ("Cerenkov hit number %d/%d, X:%d, Y:%d\n",hit,ncerenkovs,cx,cy); 
		       
		       //TParticle *MIP = (TParticle*)gAlice->Particle(cmother);
		       AliRICHhit* mipHit = (AliRICHhit*) pRICH->Hits()->UncheckedAt(0);
		       mom[0] = current->Px();
		       mom[1] = current->Py();
		       mom[2] = current->Pz();
		       //mom[0] = cHit->fMomX;
		       // mom[1] = cHit->fMomZ;
		       //mom[2] = cHit->fMomY;
		       //Float_t energymip = MIP->Energy();
		       //Float_t Mip_px = mipHit->fMomFreoX;
		       //Float_t Mip_py = mipHit->fMomFreoY;
		       //Float_t Mip_pz = mipHit->fMomFreoZ;
		       //Float_t Mip_px = MIP->Px();
		       //Float_t Mip_py = MIP->Py();
		       //Float_t Mip_pz = MIP->Pz();
		       
		       
		       
		       //Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
		       //Float_t rt = TMath::Sqrt(r);
		       //Float_t Mip_r = Mip_px*Mip_px + Mip_py*Mip_py + Mip_pz*Mip_pz;	
		       //Float_t Mip_rt = TMath::Sqrt(Mip_r);
		       //Float_t coscerenkov = (mom[0]*Mip_px + mom[1]*Mip_py + mom[2]*Mip_pz)/(rt*Mip_rt+0.0000001);
		       //Float_t cherenkov = TMath::ACos(coscerenkov);
		       ckovangle->Fill(cherenkov,(float) 1);                           //Cerenkov angle calculus
		       //printf("Cherenkov: %f\n",cherenkov);
		       Float_t ckphi=TMath::ATan2(mom[0], mom[2]);
		       hckphi->Fill(ckphi,(float) 1);
		       
		       
		       //Float_t mix = MIP->Vx();
		       //Float_t miy = MIP->Vy();
		       Float_t mx = mipHit->X();
		       Float_t my = mipHit->Z();
		       //printf("FX %e, FY %e, VX %e, VY %e\n",cx,cy,mx,my);
		       Float_t dx = trackglob[0] - mx;
		       Float_t dy = trackglob[2] - my;
		       //printf("Dx:%f, Dy:%f\n",dx,dy);
		       Float_t final_radius = TMath::Sqrt(dx*dx+dy*dy);
		       //printf("Final radius:%f\n",final_radius);
		       radius->Fill(final_radius,(float) 1);
		       
		       phspectra1->Fill(energyckov*1e9,(float) 1);
		       phot++;
		     }
		   for (Int_t nmothers=0;nmothers<=ntracks;nmothers++){
		     if (cmother == nmothers){
		       if (closs == 4)
			 mothers2[cmother]++;
		       mothers[cmother]++;
		     }
		   } 
		 }
	     }
	   }
	   

	   if(gAlice->TreeR())
	     {
	       Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
	       gAlice->TreeR()->GetEvent(nent-1);
	       TClonesArray *Rawclusters = pRICH->RawClustAddress(2);    //  Raw clusters branch
	       //printf ("Rawclusters:%p",Rawclusters);
	       Int_t nrawclusters = Rawclusters->GetEntriesFast();
	       	       
	       if (nrawclusters) {
		 printf("Raw Clusters    : %d\n",nrawclusters);
		 for (Int_t hit=0;hit<nrawclusters;hit++) {
		   AliRICHRawCluster* rcHit = (AliRICHRawCluster*) pRICH->RawClustAddress(2)->UncheckedAt(hit);
		   //Int_t nchamber = rcHit->fChamber;     // chamber number
		   //Int_t nhit = cHit->fHitNumber;        // hit number
		   Int_t qtot = rcHit->fQ;                 // charge
		   Float_t fx  =  rcHit->fX;                 // x-position
		   Float_t fy  =  rcHit->fY;                 // y-position
		   //Int_t type = rcHit->fCtype;             // cluster type ?   
		   Int_t mult = rcHit->fMultiplicity;      // How many pads form the cluster
		   pads += mult;
		   if (qtot > 0) {
		     //printf ("fx: %d, fy: %d\n",fx,fy);
		     if (fx>(x-4) && fx<(x+4)  && fy>(y-4) && fy<(y+4)) {
		       //printf("There %d \n",mult);
		       padmip+=mult;
		     } else {
		       padnumber->Fill(mult,(float) 1);
		       nraw++;
		       if (mult<4) Clcharge->Fill(qtot,(float) 1);
		     }
		     
		   }
		 }
	       }
	       
	       
	       TClonesArray *RecHits1D = pRICH->RecHitsAddress1D(2);
	       Int_t nrechits1D = RecHits1D->GetEntriesFast();
	       //printf (" nrechits:%d\n",nrechits);
	       
	       if(nrechits1D)
		 {
		   for (Int_t hit=0;hit<nrechits1D;hit++) {
		     AliRICHRecHit1D* recHit1D = (AliRICHRecHit1D*) pRICH->RecHitsAddress1D(2)->UncheckedAt(hit);
		     Float_t r_omega = recHit1D->fOmega;                  // Cerenkov angle
		     Float_t *cer_pho = recHit1D->fCerPerPhoton;        // Cerenkov angle per photon
		     Int_t *padsx = recHit1D->fPadsUsedX;           // Pads Used fo reconstruction (x)
		     Int_t *padsy = recHit1D->fPadsUsedY;           // Pads Used fo reconstruction (y)
		     Int_t goodPhotons = recHit1D->fGoodPhotons;    // Number of pads used for reconstruction
		     
		     Omega1D->Fill(r_omega,(float) 1);
		     
		     for (Int_t i=0; i<goodPhotons; i++)
		       {
			 PhotonCer->Fill(cer_pho[i],(float) 1);
			 PadsUsed->Fill(padsx[i],padsy[i],1);
			 //printf("Angle:%f, pad: %d %d\n",cer_pho[i],padsx[i],padsy[i]);
		       }
		     
		     //printf("Omega: %f, Theta: %f, Phi: %f\n",r_omega,r_theta,r_phi);
		   }
		 }

	       
	       TClonesArray *RecHits3D = pRICH->RecHitsAddress3D(2);
	       Int_t nrechits3D = RecHits3D->GetEntriesFast();
	       //printf (" nrechits:%d\n",nrechits);
	       
	       if(nrechits3D)
		 {
		   recEffEvent = 0;
		   
		   //for (Int_t hit=0;hit<nrechits3D;hit++) {
		   AliRICHRecHit3D* recHit3D = (AliRICHRecHit3D*) pRICH->RecHitsAddress3D(2)->UncheckedAt(track);
		   Float_t r_omega    = recHit3D->fOmega;                  // Cerenkov angle
		   Float_t r_theta    = recHit3D->fTheta;                  // Theta angle of incidence
		   Float_t r_phi      = recHit3D->fPhi;                    // Phi angle if incidence
		   Float_t meanradius = recHit3D->fMeanRadius;              // Mean radius for reconstructed point
		   Float_t originalOmega = recHit3D->fOriginalOmega;       // Real Cerenkov angle
		   Float_t originalTheta = recHit3D->fOriginalTheta;       // Real incidence angle
		   Float_t originalPhi = recHit3D->fOriginalPhi;           // Real azimuthal angle
		   
		   
		   //correction to track cerenkov angle
		   originalOmega = (Float_t) ckovangle->GetMean();
		   
		   if(diaglevel == 4)
		     {
		       printf("\nMean cerenkov angle: %f\n", originalOmega);
		       printf("Reconstructed cerenkov angle: %f\n",r_omega);
		     }
		   
		   Float_t omegaError = TMath::Abs(originalOmega - r_omega);
		   Float_t thetaError = TMath::Abs(originalTheta - r_theta);
		   Float_t phiError   = TMath::Abs(originalPhi - r_phi);
		   
		   //chiSquareOmega += (omegaError/originalOmega)*(omegaError/originalOmega); 
		   //chiSquareTheta += (thetaError/originalTheta)*(thetaError/originalTheta); 
		   //chiSquarePhi += (phiError/originalPhi)*(phiError/originalPhi); 
		   
		   if(TMath::Abs(omegaError) < 0.015)
		     recEffEvent += 1;
		   
		   
		   
		   //printf("rechit %f %f %f %f %f\n",recHit3D->fOmega,recHit3D->fTheta,recHit3D->fPhi, recHit3D->fX,recHit3D->fY);  
		   
		   Omega3D->Fill(r_omega,(float) 1);
		   Theta->Fill(r_theta*180/TMath::Pi(),(float) 1);
		   Phi->Fill(r_phi*180/TMath::Pi()-180,(float) 1);
		   MeanRadius->Fill(meanradius,(float) 1);
		   identification->Fill(PTfinal, r_omega,1);
		   OriginalOmega->Fill(originalOmega, (float) 1);
		   OriginalTheta->Fill(originalTheta, (float) 1);
		   OriginalPhi->Fill(TMath::Abs(originalPhi), (float) 1);
		   OmegaError->Fill(omegaError, (float) 1);
		   ThetaError->Fill(thetaError, (float) 1);
		   PhiError->Fill(phiError, (float) 1);
		   
		   recEffEvent = recEffEvent;
		   recEffTotal += recEffEvent;
		   
		   Float_t pioncer = acos(sqrt((.139*.139+PTfinal*PTfinal)/(PTfinal*PTfinal*1.285*1.285)));
		   Float_t kaoncer = acos(sqrt((.439*.439+PTfinal*PTfinal)/(PTfinal*PTfinal*1.285*1.285)));
		   Float_t protoncer = acos(sqrt((.938*.938+PTfinal*PTfinal)/(PTfinal*PTfinal*1.285*1.285)));

		   Float_t piondist = TMath::Abs(r_omega - pioncer);
		   Float_t kaondist = TMath::Abs(r_omega - kaoncer);
		   Float_t protondist = TMath::Abs(r_omega - protoncer);

		   if(diaglevel == 4)
		     {
		       if(pioncer<r_omega)
			 {
			   printf("Identified as a PION!\n");
			   pionCount += 1;
			 }
		       if(kaoncer<r_omega && pioncer>r_omega)
			 {
			   if(kaondist>piondist)
			     {
			       printf("Identified as a PION!\n");
			       pionCount += 1;
			     }
			   else
			     {
			       printf("Identified as a KAON!\n");
			       kaonCount += 1;
			     }
			 }			 }
		       if(protoncer<r_omega && kaoncer>r_omega)
			 {
			   if(kaondist>protondist)
			     {
			       printf("Identified as a PROTON!\n");
			       protonCount += 1;
			     }
			   else
			     {
			       printf("Identified as a KAON!\n");
			       pionCount += 1;
			     }
			 }
		       if(protoncer>r_omega)
			 {
			   printf("Identified as a PROTON!\n");
			   protonCount += 1;
			 }

		       printf("\nReconstruction efficiency: %5.2f%%\n", recEffEvent*100);
		 }
	     }
       }
   
       
       for (Int_t nmothers=0;nmothers<ntracks;nmothers++){
	 totalphotonstrack->Fill(mothers[nmothers],(float) 1);
	 mother->Fill(mothers2[nmothers],(float) 1);
	 //printf ("Entries in %d : %d\n",nmothers, mothers[nmothers]);
       }
       
       clusev->Fill(nraw,(float) 1);
       photev->Fill(phot,(float) 1);
       feedev->Fill(feed,(float) 1);
       padsmip->Fill(padmip,(float) 1);
       padscl->Fill(pads,(float) 1);
       //printf("Photons:%d\n",phot);
       phot = 0;
       feed = 0;
       pads = 0;
       nraw=0;
       padmip=0;
       
       
       
       gAlice->ResetDigits();
       //Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
       gAlice->TreeD()->GetEvent(0);
       
       if (diaglevel < 4)
	 {
	   
	   
	   TClonesArray *Digits  = pRICH->DigitsAddress(2);
	   Int_t ndigits = Digits->GetEntriesFast();
	   printf("Digits          : %d\n",ndigits);
	   padsev->Fill(ndigits,(float) 1);
	   for (Int_t hit=0;hit<ndigits;hit++) {
	     AliRICHDigit* dHit = (AliRICHDigit*) Digits->UncheckedAt(hit);
	     Int_t qtot = dHit->Signal();                // charge
	     Int_t ipx  = dHit->PadX();               // pad number on X
	     Int_t ipy  = dHit->PadY();               // pad number on Y
	     //printf("%d, %d\n",ipx,ipy);
	     if( ipx<=100 && ipy <=100) hc0->Fill(ipx,ipy,(float) qtot);
	   }
	 }
       
       if (diaglevel == 5)
	 {
	   for (Int_t ich=0;ich<7;ich++)
	     {
	       TClonesArray *Digits = pRICH->DigitsAddress(ich);    //  Raw clusters branch
	       Int_t ndigits = Digits->GetEntriesFast();
	       //printf("Digits:%d\n",ndigits);
	       padsev->Fill(ndigits,(float) 1); 
	       if (ndigits) {
		 for (Int_t hit=0;hit<ndigits;hit++) {
		   AliRICHDigit* dHit = (AliRICHDigit*) Digits->UncheckedAt(hit);
		   //Int_t nchamber = dHit->GetChamber();     // chamber number
		   //Int_t nhit = dHit->fHitNumber;          // hit number
		   Int_t qtot = dHit->Signal();                // charge
		   Int_t ipx  = dHit->PadX();               // pad number on X
		   Int_t ipy  = dHit->PadY();               // pad number on Y
		   //Int_t iqpad  = dHit->fQpad;           // charge per pad
		   //Int_t rpad  = dHit->fRSec;            // R-position of pad
		   //printf ("Pad hit, PadX:%d, PadY:%d\n",ipx,ipy);
		   if( ipx<=100 && ipy <=100 && ich==2) hc0->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==0) hc1->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==1) hc2->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==2) hc3->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==3) hc4->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==4) hc5->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==5) hc6->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==6) hc7->Fill(ipx,ipy,(float) qtot);
		 }
	       }
	     }
	 }
   }
   
   if(diaglevel == 4)
     {

       Stat_t omegaE;
       Stat_t thetaE;
       Stat_t phiE;
       
       Stat_t omegaO;
       Stat_t thetaO;
       Stat_t phiO;
       
       for(Int_t i=0;i<99;i++)
	 {
	   omegaE = OriginalOmega->GetBinContent(i);
	   if(omegaE != 0)
	     {
	       omegaO = Omega3D->GetBinContent(i);
	       chiSquareOmega += (TMath::Power(omegaE,2) - TMath::Power(omegaO,2))/omegaO;
	     }

  	   thetaE = OriginalTheta->GetBinContent(i);
	   if(thetaE != 0)
	     {
	       thetaO = Theta->GetBinContent(i);
	       chiSquareTheta += (TMath::Power(thetaE,2) - TMath::Power(thetaO,2))/thetaO;
	     }

	   phiE = OriginalPhi->GetBinContent(i);
	   if(phiE != 0)
	     {
	       phiO = Phi->GetBinContent(i);
	       chiSquarePhi += (TMath::Power(phiE,2) - TMath::Power(phiO,2))/phiO;
	     }
	   
	   //printf(" o: %f  t: %f  p: %f\n", OriginalOmega->GetBinContent(i), OriginalTheta->GetBinContent(i),OriginalPhi->GetBinContent(i));

	 }

       

       printf("\nChi square test values:   Omega - %f\n", chiSquareOmega);
       printf("                          Theta - %f\n", chiSquareTheta);
       printf("                          Phi   - %f\n", chiSquarePhi);
       
       printf("\nKolmogorov test values:   Omega - %5.4f\n", Omega3D->KolmogorovTest(OriginalOmega));
       printf("                          Theta - %5.4f\n", Theta->KolmogorovTest(OriginalTheta));
       printf("                          Phi   - %5.4f\n", Phi->KolmogorovTest(OriginalPhi));

       recEffTotal = recEffTotal/evNumber2;
       printf("\nTotal reconstruction efficiency: %5.2f%%\n", recEffTotal*100);
       printf("\n Pions: %d\n Kaons: %d\n Protons:%d\n",pionCount, kaonCount, protonCount);

     }
   
   
   //Create canvases, set the view range, show histograms

   TCanvas *c1 = 0;
   TCanvas *c2 = 0;
   TCanvas *c3 = 0;
   TCanvas *c4 = 0;
   TCanvas *c5 = 0;
   TCanvas *c6 = 0;
   TCanvas *c7 = 0;
   TCanvas *c8 = 0;
   TCanvas *c9 = 0;
   TCanvas *c10 = 0;
   TCanvas *c11 = 0;
   TCanvas *c12 = 0;
   TCanvas *c13 = 0;

   //TF1* expo = 0;
   //TF1* gaus = 0;
   
   TStyle *mystyle=new TStyle("Plain","mystyle");
   mystyle->SetPalette(1,0);
   //mystyle->SetTitleYSize(0.2);
   //mystyle->SetStatW(0.19);
   //mystyle->SetStatH(0.1);
   //mystyle->SetStatFontSize(0.01);
   //mystyle->SetTitleYSize(0.3);
   mystyle->SetFuncColor(2);
   //mystyle->SetOptStat(0111);
   mystyle->SetDrawBorder(0);
   mystyle->SetTitleBorderSize(0);
   mystyle->SetOptFit(1111);
   mystyle->cd();

   
   TClonesArray *RecHits3D = pRICH->RecHitsAddress3D(2);
   Int_t nrechits3D = RecHits3D->GetEntriesFast();
   TClonesArray *RecHits1D = pRICH->RecHitsAddress1D(2);
   Int_t nrechits1D = RecHits1D->GetEntriesFast();

  switch(diaglevel)
     {
     case 1:
       
       c1 = new TCanvas("c1","Alice RICH digits",50,50,300,350);
       hc0->SetXTitle("ix (npads)");
       hc0->Draw("colz");
	
//
       c2 = new TCanvas("c2","Hits per type",100,100,600,700);
       c2->Divide(2,2);
       //c4->SetFillColor(42);

       c2->cd(1);
       feedback->SetXTitle("x (cm)");
       feedback->SetYTitle("y (cm)");
       feedback->Draw("colz");
       
       c2->cd(2);
       //mip->SetFillColor(5);
       mip->SetXTitle("x (cm)");
       mip->SetYTitle("y (cm)");
       mip->Draw("colz");
       
       c2->cd(3);
       //cerenkov->SetFillColor(5);
       cerenkov->SetXTitle("x (cm)");
       cerenkov->SetYTitle("y (cm)"); 
       cerenkov->Draw("colz");
       
       c2->cd(4);
       //h->SetFillColor(5);
       h->SetXTitle("x (cm)");
       h->SetYTitle("y (cm)");
       h->Draw("colz");

       c3 = new TCanvas("c3","Hits distribution",150,150,600,350);
       c3->Divide(2,1);
       //c10->SetFillColor(42);
       
       c3->cd(1);
       hitsX->SetFillColor(5);
       hitsX->SetXTitle("(cm)");
       hitsX->Draw();
       
       c3->cd(2);
       hitsY->SetFillColor(5);
       hitsY->SetXTitle("(cm)");
       hitsY->Draw();
       
      
       break;
//
     case 2:
       
       c4 = new TCanvas("c4","Photon Spectra",50,50,600,350);
       c4->Divide(2,1);
       //c6->SetFillColor(42);
       
       c4->cd(1);
       phspectra2->SetFillColor(5);
       phspectra2->SetXTitle("energy (eV)");
       phspectra2->Draw();
       c4->cd(2);
       phspectra1->SetFillColor(5);
       phspectra1->SetXTitle("energy (eV)");
       phspectra1->Draw();
       
       c5 = new TCanvas("c5","Particles Spectra",100,100,600,700);
       c5->Divide(2,2);
       //c9->SetFillColor(42);
       
       c5->cd(1);
       pionspectra->SetFillColor(5);
       pionspectra->SetXTitle("(GeV)");
       pionspectra->Draw();
       
       c5->cd(2);
       protonspectra->SetFillColor(5);
       protonspectra->SetXTitle("(GeV)");
       protonspectra->Draw();
       
       c5->cd(3);
       kaonspectra->SetFillColor(5);
       kaonspectra->SetXTitle("(GeV)");
       kaonspectra->Draw();
       
       c5->cd(4);
       chargedspectra->SetFillColor(5);
       chargedspectra->SetXTitle("(GeV)");
       chargedspectra->Draw();

       break;
       
     case 3:

       
       if(gAlice->TreeR())
	 {
	   c6=new TCanvas("c6","Clusters Statistics",50,50,600,700);
	   c6->Divide(2,2);
	   //c3->SetFillColor(42);
	   
	   c6->cd(1);
	   //TPad* c6_1;
	   //c6_1->SetLogy();
	   Clcharge->SetFillColor(5);
	   Clcharge->SetXTitle("ADC counts");
	   if (evNumber2>10)
	     {
	       Clcharge->Fit("expo");
	       //expo->SetLineColor(2);
	       //expo->SetLineWidth(3);
	     }
	   Clcharge->Draw();
	   
	   c6->cd(2);
	   padnumber->SetFillColor(5);
	   padnumber->SetXTitle("(counts)");
	   padnumber->Draw();
	   
	   c6->cd(3);
	   clusev->SetFillColor(5);
	   clusev->SetXTitle("(counts)");
	   if (evNumber2>10)
	     {
	       clusev->Fit("gaus");
	       //gaus->SetLineColor(2);
	       //gaus->SetLineWidth(3);
	     }
	   clusev->Draw();
	   
	   c6->cd(4);
	   padsmip->SetFillColor(5);
	   padsmip->SetXTitle("(counts)");
	   padsmip->Draw(); 
	 }
       
       if(evNumber2<1)
	 {
	   c11 = new TCanvas("c11","Cherenkov per Mip",400,10,600,700);
	   mother->SetFillColor(5);
	   mother->SetXTitle("counts");
	   mother->Draw();
	 }

       c7 = new TCanvas("c7","Production Statistics",100,100,600,700);
       c7->Divide(2,2);
       //c7->SetFillColor(42);
       
       c7->cd(1);
       totalphotonsevent->SetFillColor(5);
       totalphotonsevent->SetXTitle("Photons (counts)");
       if (evNumber2>10)
	   {
	     totalphotonsevent->Fit("gaus");
	     //gaus->SetLineColor(2);
	     //gaus->SetLineWidth(3);
	   }
       totalphotonsevent->Draw();
       
       c7->cd(2);
       photev->SetFillColor(5);
       photev->SetXTitle("(counts)");
       if (evNumber2>10)
	 {
	   photev->Fit("gaus");
	   //gaus->SetLineColor(2);
	   //gaus->SetLineWidth(3);
	 }
       photev->Draw();
       
       c7->cd(3);
       feedev->SetFillColor(5);
       feedev->SetXTitle("(counts)");
       if (evNumber2>10)
	 {
	   feedev->Fit("gaus");
	   //gaus->SetLineColor(2);
	   //gaus->SetLineWidth(3);
	 }
       feedev->Draw();

       c7->cd(4);
       padsev->SetFillColor(5);
       padsev->SetXTitle("(counts)");
       if (evNumber2>10)
	 {
	   padsev->Fit("gaus");
	   //gaus->SetLineColor(2);
	   //gaus->SetLineWidth(3);
	 }
       padsev->Draw();

       break;

     case 4:
       

       if(nrechits3D)
	 {
	   c8 = new TCanvas("c8","3D reconstruction of Phi angle",50,50,300,1050);
	   c8->Divide(1,3);
	   //c2->SetFillColor(42);
	   
	   
	   // data per hit
	   c8->cd(1);
	   hitsPhi->SetFillColor(5);
	   if (evNumber2>10)
	     hitsPhi->Fit("gaus");
	   hitsPhi->Draw();
	   
	    //data per track
	   c8->cd(2);
	   OriginalPhi->SetFillColor(5);
	   if (evNumber2>10)
	     OriginalPhi->Fit("gaus");
	   OriginalPhi->Draw();

	   //recontructed data
	   c8->cd(3);
	   Phi->SetFillColor(5);
	   if (evNumber2>10)
	     Phi->Fit("gaus");
	   Phi->Draw();

	   c9 = new TCanvas("c9","3D reconstruction of theta angle",75,75,300,1050);
	   c9->Divide(1,3);

	   // data per hit
	   c9->cd(1);
	   hitsTheta->SetFillColor(5);
	   if (evNumber2>10)
	     hitsTheta->Fit("gaus");
	   hitsTheta->Draw();
	   
	   //data per track
	   c9->cd(2);
	   OriginalTheta->SetFillColor(5);
	   if (evNumber2>10)
	     OriginalTheta->Fit("gaus");
	   OriginalTheta->Draw();

	   //recontructed data
	   c9->cd(3);
	   Theta->SetFillColor(5);
	   if (evNumber2>10)
	     Theta->Fit("gaus");
	   Theta->Draw();

	   c10 = new TCanvas("c10","3D reconstruction of cherenkov angle",100,100,300,1050);
	   c10->Divide(1,3);

	   // data per hit
	   c10->cd(1);
	   ckovangle->SetFillColor(5);
	   ckovangle->SetXTitle("angle (radians)");
	   if (evNumber2>10)
	     ckovangle->Fit("gaus");
	   ckovangle->Draw();
	   
	   //data per track
	   c10->cd(2);
	   OriginalOmega->SetFillColor(5);
	   OriginalOmega->SetXTitle("angle (radians)");
	   if (evNumber2>10)
	     OriginalOmega->Fit("gaus");
	   OriginalOmega->Draw();

	   //recontructed data
	   c10->cd(3);
	   Omega3D->SetFillColor(5);
	   Omega3D->SetXTitle("angle (radians)");
	   if (evNumber2>10)
	     Omega3D->Fit("gaus");
	   Omega3D->Draw(); 


	   c11 = new TCanvas("c11","3D reconstruction of mean radius",125,125,300,700);
	   c11->Divide(1,2);

	   // data per hit
	   c11->cd(1);
	   radius->SetFillColor(5);
	   radius->SetXTitle("radius (cm)");
	   radius->Draw();

	   //recontructed data
	   c11->cd(2);
	   MeanRadius->SetFillColor(5);
	   MeanRadius->SetXTitle("radius (cm)");
	   MeanRadius->Draw();

	   
	   c12 = new TCanvas("c12","Cerenkov angle vs. Momentum",150,150,550,350);

	   c12->cd(1);
	   identification->SetFillColor(5);
	   identification->SetXTitle("Momentum (GeV/c)");
	   identification->SetYTitle("Cherenkov angle (radians)");
	   
	   //Float_t pionmass=.139;
	   //Float_t kaonmass=.493;
	   //Float_t protonmass=.938;
	   //Float_t n=1.295;
	   
	   TF1 *pionplot = new TF1("pion","acos(sqrt((.139*.139+x*x)/(x*x*1.285*1.285)))",1,5);
	   TF1 *kaonplot = new TF1("kaon","acos(sqrt((.439*.439+x*x)/(x*x*1.285*1.285)))",1,5);
	   TF1 *protonplot = new TF1("proton","acos(sqrt((.938*.938+x*x)/(x*x*1.285*1.285)))",1,5);
	   
	   identification->Draw();

	   pionplot->SetLineColor(5);
	   pionplot->Draw("same");

	   kaonplot->SetLineColor(4);
	   kaonplot->Draw("same");

	   protonplot->SetLineColor(3);
	   protonplot->Draw("same");
	   //identification->Draw("same");



	   c13 = new TCanvas("c13","Reconstruction Errors",200,200,900,350);
	   c13->Divide(3,1);

	   c13->cd(1);
	   PhiError->SetFillColor(5);
	   if (evNumber2>10)
	     PhiError->Fit("gaus");
	   PhiError->Draw();
	   c13->cd(2);
	   ThetaError->SetFillColor(5);
	   if (evNumber2>10)
	     ThetaError->Fit("gaus");
	   ThetaError->Draw();
	   c13->cd(3);
	   OmegaError->SetFillColor(5);
	   OmegaError->SetXTitle("angle (radians)");
	   if (evNumber2>10)
	     OmegaError->Fit("gaus");
	   OmegaError->Draw();
	   
	 }
       
       if(nrechits1D)
	 {
	   c9 = new TCanvas("c9","1D Reconstruction",100,100,1100,700);
	   c9->Divide(3,2);
	   //c5->SetFillColor(42);
	   
	   c9->cd(1);
	   ckovangle->SetFillColor(5);
	   ckovangle->SetXTitle("angle (radians)");
	   ckovangle->Draw();
	   
	   c9->cd(2);
	   radius->SetFillColor(5);
	   radius->SetXTitle("radius (cm)");
	   radius->Draw();
	   
	   c9->cd(3);
	   hc0->SetXTitle("pads");
	   hc0->Draw("box"); 
	   
	   c9->cd(5);
	   Omega1D->SetFillColor(5);
	   Omega1D->SetXTitle("angle (radians)");
	   Omega1D->Draw();
	   
	   c9->cd(4);
	   PhotonCer->SetFillColor(5);
	   PhotonCer->SetXTitle("angle (radians)");
	   PhotonCer->Draw();
	   
	   c9->cd(6);
	   PadsUsed->SetXTitle("pads");
	   PadsUsed->Draw("box"); 
	 }
       
       break;
       
     case 5:
       
       printf("Drawing histograms.../n");

       //if (ndigits)
	 //{
       c10 = new TCanvas("c10","Alice RICH digits",50,50,1200,700);
       c1->Divide(4,2);
       //c1->SetFillColor(42);
       
       c10->cd(1);
       hc1->SetXTitle("ix (npads)");
       hc1->Draw("box");
       c10->cd(2);
       hc2->SetXTitle("ix (npads)");
       hc2->Draw("box");
       c10->cd(3);
       hc3->SetXTitle("ix (npads)");
       hc3->Draw("box");
       c10->cd(4);
       hc4->SetXTitle("ix (npads)");
       hc4->Draw("box");
       c10->cd(5);
       hc5->SetXTitle("ix (npads)");
       hc5->Draw("box");
       c10->cd(6);
       hc6->SetXTitle("ix (npads)");
       hc6->Draw("box");
       c10->cd(7);
       hc7->SetXTitle("ix (npads)");
       hc7->Draw("box");
       c10->cd(8);
       hc0->SetXTitle("ix (npads)");
       hc0->Draw("box");
	 //}
//
       c11 = new TCanvas("c11","Hits per type",100,100,600,700);
       c11->Divide(2,2);
       //c4->SetFillColor(42);
       
       c11->cd(1);
       feedback->SetXTitle("x (cm)");
       feedback->SetYTitle("y (cm)");
       feedback->Draw();
       
       c11->cd(2);
       //mip->SetFillColor(5);
       mip->SetXTitle("x (cm)");
       mip->SetYTitle("y (cm)");
       mip->Draw();
       
       c11->cd(3);
       //cerenkov->SetFillColor(5);
       cerenkov->SetXTitle("x (cm)");
       cerenkov->SetYTitle("y (cm)"); 
       cerenkov->Draw();
       
       c11->cd(4);
       //h->SetFillColor(5);
       h->SetXTitle("x (cm)");
       h->SetYTitle("y (cm)");
       h->Draw();

       c12 = new TCanvas("c12","Hits distribution",150,150,600,350);
       c12->Divide(2,1);
       //c10->SetFillColor(42);
       
       c12->cd(1);
       hitsX->SetFillColor(5);
       hitsX->SetXTitle("(cm)");
       hitsX->Draw();
       
       c12->cd(2);
       hitsY->SetFillColor(5);
       hitsY->SetXTitle("(cm)");
       hitsY->Draw();
       
       break;
       
     }
       

   // calculate the number of pads which give a signal


   //Int_t Np=0;
   /*for (Int_t i=0;i< NpadX;i++) {
       for (Int_t j=0;j< NpadY;j++) {
	   if (Pad[i][j]>=6){
	       Np+=1;
	   }
       }
   }*/
   //printf("The total number of pads which give a signal: %d %d\n",Nh,Nh1);
   printf("\nEnd of analysis\n");
   printf("**********************************\n");
}//void AliRICH::DiagnosticsSE(Int_t diaglevel,Int_t evNumber1,Int_t evNumber2)
//______________________________________________________________________________
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
//______________________________________________________________________________
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
}
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
  TTree *treeS = fLoader->TreeS();
    
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
 
   //this is after TreeH because we need to guarantee that fHits array is created
  AliDetector::SetTreeAddress();
    
  if(treeS){
    if(GetDebug())Info("SetTreeAddress","tree S is requested.");
    branch = treeS->GetBranch("RICH");
    if(branch){
      if(!fSDigits) fSDigits=new TClonesArray("AliRICHdigit",100000);
      branch->SetAddress(&fSDigits);
    }
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
//______________________________________________________________________________
void AliRICH::Print(Option_t *option)const
{
  TObject::Print(option);
  fpParam->Dump();
  Chambers()->Print(option);  
}//void AliRICH::Print(Option_t *option)const
//______________________________________________________________________________
void AliRICH::CreateGeometry()
{//Creates detailed geometry simulation (currently GEANT volumes tree)         
  if(GetDebug())Info("CreateGeometry","Start.");
//???????? to be removed to AliRICHParam?
  fpParam->RadiatorToPads(fpParam->FreonThickness()/2+fpParam->QuartzThickness()+fpParam->GapThickness());
    
//Opaque quartz thickness
  Float_t oqua_thickness = .5;
//CsI dimensions
  Float_t csi_width =fpParam->Nx()*fpParam->PadX()+fpParam->DeadZone();
  Float_t csi_length=fpParam->Ny()*fpParam->PadY()+2*fpParam->DeadZone();
  
  Int_t *idtmed = fIdtmed->GetArray()-999;
    
  Int_t i;
  Float_t zs;
  Int_t idrotm[1099];
  Float_t par[3];
    
//External aluminium box 
  par[0]=68.8;par[1]=13;par[2]=70.86;//Original Settings
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
  CreateChambers();

  for(int i=0;i<kNCH;i++){
    AliMatrix(idrotm[1000+i],C(i)->ThetaXd(),C(i)->PhiXd(),
                             C(i)->ThetaYd(),C(i)->PhiYd(),
                             C(i)->ThetaZd(),C(i)->PhiZd());
    gMC->Gspos("RICH",i+1,"ALIC",C(i)->X(),C(i)->Y(),C(i)->Z(),idrotm[1000+i], "ONLY");
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
