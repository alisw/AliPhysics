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

//_________________________________________________________________________
// This is a TTask that constructs ReconstParticles (reconstructed particles) 
// out of Digits
// 
//-- Authors: Evgeny Karpechev(INR) and Alla Maevsksia
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"
#include "TFolder.h"

// --- Standard library ---
#include <stdlib.h>
#include <Riostream.h>

// --- AliRoot header files ---

#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliFMDdigit.h"
#include "AliFMDhit.h"
#include "AliFMDReconstParticles.h"
#include "AliFMD.h"
#include "AliFMDv1.h"
#include "AliFMDReconstruction.h"
#include "AliRun.h"
#include "AliConfig.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"

ClassImp(AliFMDReconstruction)

        
//____________________________________________________________________________ 

AliFMDReconstruction::AliFMDReconstruction():TTask("AliFMDReconstruction","") 
{
  fNevents = 0 ;  // Number of events to rreconnstraction, 0 means all events in current file
  fRunLoader = 0x0;
  
}
//____________________________________________________________________________ 

AliFMDReconstruction::AliFMDReconstruction(AliRunLoader* rl):TTask("AliFMDReconstruction","")
{

  if (rl == 0x0)
   {
     Fatal("AliFMDReconstruction","Argument AliRunLoader* is null!");
     return;
   }
   
  fNevents = 0 ;    // Number of events to rreconnstraction, 0 means all events in current file

  fRunLoader = rl;
  AliLoader* gime = fRunLoader->GetLoader("FMDLoader");
  if (gime == 0x0)
   {
     Fatal("AliFMDReconstruction","Can not find FMD (loader) in specified event");
     return;//never reached
   }
  //add Task to //root/Tasks folder
  gime->PostReconstructioner(this);
}

//____________________________________________________________________________ 

AliFMDReconstruction::~AliFMDReconstruction()
{
}

//____________________________________________________________________________

void AliFMDReconstruction::Exec(Option_t *option) 
{ 
 //Collects all digits in the same active volume into number of particles
  /*
    Reconstruct number of particles 
    in given group of pads for given FMDvolume
    determine by numberOfVolume , 
    numberOfMinSector,numberOfMaxSector,
    numberOfMinRing, numberOgMaxRing
    Reconstruction method choose dependence on number of empty pads  
  */


  cout<<"\nStart AliFMDReconstruction::Exec(...)"<<endl;
  Int_t const knumVolumes=5;
  Int_t const knumSectors=40;
  Int_t const knumRings=768;
  Int_t padADC[10][50][800];
  Float_t eta, etain,etaout,rad,theta;
  Int_t ivol, iSector, iRing;
  Float_t rin[5]={4.2,15.4,4.2,15.4,4.2};
  Float_t rout[5]={17.4,28.4,17.4,28.4,17.4};
  Float_t z[5]={62.8, 75.2, -83.4, -75.2, -340.};
  Int_t numberOfRings[5]={768,384,768,384,768};
  Int_t numberOfSectors[5]=  {20,40,20,40,20};
  Int_t numberOfEtaIntervals[5];
  // number of ring for boundary 0.1 eta

  
  if (fRunLoader == 0x0)
   {
    Error("Exec","Run Loader loader is NULL - Session not opened");
    return;
   }
  AliLoader* gime = fRunLoader->GetLoader("FMDLoader");
  if (gime == 0x0)
   {
     Fatal("AliFMDReconstruction","Can not find FMD (loader) in specified event");
     return;//never reached
   }
   
  fRunLoader->LoadgAlice();
  Int_t retval;
  
  retval = gime->LoadHits("READ"); 
  if (retval)
   {
     Error("Exec","Error occured while loading hits. Exiting.");
     return;
   }

  retval = gime->LoadDigits("READ"); 
  if (retval)
   {
     Error("Exec","Error occured while loading digits. Exiting.");
     return;
   }
  
  AliFMD * fFMD = (AliFMD *) gAlice->GetDetector("FMD");
  TClonesArray *fReconParticles=fFMD->ReconParticles();
  
  TTree* treeD = gime->TreeD();
  if (treeD == 0x0)
   {
     Error("Exec","Can not get Tree with Digits. Nothing to reconstruct - Exiting");
     return;
   }
  if(fNevents == 0) fNevents=(Int_t)treeD->GetEntries();
  //PH Do we use TreeE (kinematics), or TreeD (digits) toaccess the number
  //PH of events?
//PH   if(fNevents == 0) fNevents=(Int_t)gAlice->TreeE()->GetEntries(); 
//PH   cout<<" fNevents "<<fNevents<<endl;

  TBranch *brDigits=0;

  for(Int_t ievent=0;ievent<fNevents;ievent++)
    { 
      fRunLoader->GetEvent(ievent) ;

      TTree* treeH = gime->TreeH();
      if (treeH == 0x0)
       {
         Error("Exec","Can not get TreeH");
         return;
       }
      cout<<" ievent "<<ievent<<endl;
      for (Int_t i=0; i<knumVolumes; i++)
	for(Int_t j=0; j<knumSectors; j++)
	  for(Int_t ij=0; ij<knumRings; ij++)
	    padADC[i][j][ij]=0;                    //zhachem ???

	brDigits=treeD->GetBranch("FMD");
	if (!brDigits){
	  cerr<<"EXEC Branch FMD digits not found"<<endl;
        return;
      } 
 
      if(gime->TreeR()==0) gime->MakeTree("R");

      //Make branches
      fFMD->MakeBranch("R");

      
      Int_t zeroADC=1;
 
      AliFMDdigit  *fmdDigit;
       if (fFMD)
	{
	  gime->TreeD()->GetEvent(0); 
	  TClonesArray *fFMDdigits=fFMD->Digits();
	  Int_t nDigits=fFMDdigits->GetEntries();
	  cout<<" nDigits "<<nDigits<<endl;
	   Int_t recParticles[6];
	   Int_t nRecPart=0 ;
	   Int_t zeroPads=0;
	   Int_t numberOfPads=0; //To avoid warning
	   Int_t pedestal;
	   Float_t channelWidth=(22400*50)/1024;
	   for (Int_t digit=0;digit<nDigits;digit++) 
	     {
	       fmdDigit=(AliFMDdigit*)fFMDdigits->UncheckedAt(digit);    
	       ivol=fmdDigit->Volume();
	       iSector=fmdDigit->NumberOfSector();
	       iRing=fmdDigit->NumberOfRing();
	       pedestal=Int_t(gRandom->Gaus(500,250));
	       padADC[ivol-1][iSector-1][iRing-1]=
		 fmdDigit->ADCsignal()
		 -Int_t(Float_t(pedestal)/channelWidth);
	       if (padADC[ivol-1][iSector-1][iRing-1]<0) 
		 padADC[ivol-1][iSector-1][iRing-1]=0;
	     } //digit loop
	   Int_t rmin=0; Int_t rmax=0; //To avoid warning
	   Int_t smin=0; Int_t smax=0; //To avoid warning
	   AliHeader *header = fRunLoader->GetHeader();
	   AliGenEventHeader* genHeader = header->GenEventHeader();
	   TArrayF *o = new TArrayF(3); 
	   genHeader->PrimaryVertex(*o);
	   Float_t zVertex=o->At(2);
 	   for (ivol=0; ivol<knumVolumes; ivol++)
	     {
	       Float_t realZ=zVertex+z[ivol];
	       theta=TMath::ATan(rin[ivol]/TMath::Abs(realZ));
	       etain = - TMath::Log( TMath::Tan(theta/2.));
	       theta=TMath::ATan(rout[ivol]/TMath::Abs(realZ));
	       etaout=- TMath::Log( TMath::Tan(theta/2.));
	       numberOfEtaIntervals[ivol]=Int_t ((etain-etaout)*10);
	       eta=etain;
	       for (Int_t e1=0;e1<=numberOfEtaIntervals[ivol];e1++) 
		 {
		   theta = 2.*TMath::ATan (TMath::Exp (-eta));
		   Float_t radmin = TMath::Abs(realZ) * (TMath::Tan (theta));
		   rmin= Int_t ( (radmin-rin[ivol])*numberOfRings[ivol]/(rout[ivol]-rin[ivol]));
		   eta=eta-0.1;
		   theta = 2.*TMath::ATan (TMath::Exp (-eta));
		   rad = TMath::Abs(realZ) * (TMath::Tan (theta));
		   rmax=Int_t( (rad-rin[ivol])*numberOfRings[ivol]/(rout[ivol]-rin[ivol]));
		   zeroPads=0;
		   smin=0;
		   smax=numberOfSectors[ivol]; 
		   for (Int_t iring=rmin; iring<rmax; iring++) 
		     {
		       numberOfPads=(rmax-rmin)*(smax-smin);
		       for 
			 (Int_t isector=1;
			  isector<=numberOfSectors[ivol]; 
			  isector++) 
			 if(padADC[ivol][isector-1][iring-1]<zeroADC)
			   zeroPads++;
		     } //ring
		   Float_t zerosRatio= 
		     (Float_t)zeroPads/(Float_t)numberOfPads;
		   recParticles[0]=ivol;
		   recParticles[1]=smin;
		   recParticles[2]=smax;
		   recParticles[3]=rmin;
		   recParticles[4]=rmax;
		   if (zerosRatio>0.1 )
		     recParticles[5]=
		       DeterminationByPoisson
		       (padADC,ivol+1,rmin,rmax,smin,smax);
		   else
		     recParticles[5]=
		       DeterminationByThresholds
		       (padADC,ivol+1,rmin,rmax,smin,smax);
		   new((*fReconParticles)[nRecPart++]) 
		 AliFMDReconstParticles(recParticles);   	       
		 } // eta
	     } // volume
	   
	}//if FMD
       gime->TreeR()->Reset();
       gime->TreeR()->Fill(); 
       gime->WriteRecPoints("OVERWRITE");
    } //event loop
  cout<<"\nAliFMDReconstruction::Exec finished"<<endl;
}

//------------------------------------------------------------
Int_t AliFMDReconstruction::DeterminationByThresholds
(Int_t padAdc[][50][800],Int_t volume, Int_t rmin, Int_t rmax, 
 Int_t smin, Int_t smax)
{
  /*
    reconstruct number of particles by 
    energy deposited threshold method 
    Using if number of empty pads less then 10%

*/
  cout<<"\nStart threshold method\n";
  
  Int_t thresholdInner[30]={
    0,     14,  28,    42,   57,     
    72,    89,  104,  124,  129, 
    164,  174,  179,  208,  228, 
    248,  268,   317,  337,  357, 
    392,  407,  416,  426,  436, 
    461,  468,  493,  506,  515}; 

  Int_t thresholdOuter[30]={0, 18, 48, 77, 105,
			    132, 165, 198, 231, 
			    264, 286, 308, 334, 
			    352, 374, 418, 440,
			    462, 484, 506, 528,
			    550, 572, 594, 616}; 
  Int_t threshold[30];
  for (Int_t it=0; it<30; it++) {
    if(volume==1||volume==3||volume==5) threshold[it]=thresholdInner[it];
    if(volume==2||volume==4) threshold[it]=thresholdOuter[it];
  }
  Int_t thresholdArraySize = 30;     
  Int_t numPart=0;
  //Inner Si
  for (Int_t iring=rmin; iring<rmax; iring++) 
    {
      for (Int_t isector=smin; isector<smax; isector++) 
	{
	  for (int i=0;i<thresholdArraySize-1;i++)
	    {	      
	      if(padAdc[volume-1][isector][iring]>threshold[i]
		 &&padAdc[volume-1][isector][iring]<=threshold[i+1])
		{
		  numPart+=i;
		  break;
		}; // if in threshol interval
	    }; //threshold_array_size
	}; //iring
    }; //sector
  cout<<"\nEnd threshol method"<<endl;
  return numPart;
}

//__________________________________________________________________

Int_t AliFMDReconstruction::DeterminationByPoisson 
(Int_t padAdc[][50][800],
 Int_t vol, Int_t rmin, Int_t rmax, 
 Int_t secmin, Int_t secmax)
{
  /* 
     reconstruct number of particles by Poisson statistic method 
     according number of empty pad in chosen group of pads 
     using if number of empty pads more then 10% 

  */

  //  cout<<"\nStart Poisson method";
  Int_t thresholdADC=1;   
  Int_t zeropads=0;
  for (Int_t i=rmin;i<rmax;i++)
    {
      for (Int_t j=secmin;j<secmax;j++)
	{
	  if (padAdc[vol-1][j][i]<thresholdADC) zeropads++;
	};
    };
  Float_t lambda=-TMath::Log(Float_t(zeropads)/
			     ( Float_t(secmax-secmin)*
			       Float_t(rmax-rmin))); //+1 zdes' ne nado
  Int_t fRecon=(Int_t)(lambda*(secmax-secmin)*(rmax-rmin)); //+1 zdes' ne nado
  //  cout<<"\nEnd Poisson method"<<endl;
  return fRecon;
}


