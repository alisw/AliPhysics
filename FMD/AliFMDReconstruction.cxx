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
#include <Riostream.h>

// --- Standard library ---
#include <stdlib.h>

// --- AliRoot header files ---

#include "AliFMDdigit.h"
#include "AliFMDReconstParticles.h"
#include "AliFMD.h"
#include "AliFMDv1.h"
#include "AliFMDReconstruction.h"
#include "AliRun.h"
ClassImp(AliFMDReconstruction)

        
//____________________________________________________________________________ 

AliFMDReconstruction::AliFMDReconstruction():TTask("AliFMDReconstruction","") 
{
  fNevents = 0 ;  // Number of events to rreconnstraction, 0 means all events in current file
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}
//____________________________________________________________________________ 

AliFMDReconstruction::AliFMDReconstruction(char* HeaderFile, char *ReconstParticlesFile):TTask("AliFMDReconstruction","")
{
  fNevents = 0 ;    // Number of events to rreconnstraction, 0 means all events in current file
  fReconstParticlesFile=ReconstParticlesFile ;
  fHeadersFile=HeaderFile ;
  //add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ;     
}

//____________________________________________________________________________ 

AliFMDReconstruction::~AliFMDReconstruction()
{
}
//----------------------------------------------------------------------------

void AliFMDReconstruction::Exec(Option_t *option) 
{ 
  printf (" AliFMDReconstruction starting \n");
 //Collects all digits in the same active volume into number of particles
  
  Int_t const NumVolums=5;
  Int_t const NumSectors=40;
  Int_t const NumRings=768;
  Int_t PadADC[10][50][800];
  Int_t ivol, iSector, iRing;
  Int_t Ne1; 
  //Int_t NumberOfRings[5]=  {256,128,256,128,256};
  Int_t NumberOfSectors[5]=  {20,40,20,40,20};
  // number of ring for boundary 0.1 eta
   Int_t EtaIntervalInner []=
      {0, 55, 110, 165, 221, 276, 331, 386, 442,
	497, 552, 607, 663, 718, 767 };
   /*
     {0,  18, 36,   55,  73,
      92, 110, 128, 147, 165,
      184, 202, 221, 239, 255};*/
   Int_t  EtaIntervalOuter []=  //{0, 21, 43, 65, 86, 108, 127};
     {0, 65, 130, 195, 260, 325, 383};
 Int_t EtaInterval[20];
  

  Int_t size_EtaIntervalInner=sizeof (EtaIntervalInner)/sizeof(EtaIntervalInner[0]);

  Int_t size_EtaIntervalOuter=sizeof (EtaIntervalOuter)/sizeof(EtaIntervalOuter[0]);

  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD");  
  TClonesArray *fReconParticles=FMD->ReconParticles();
  if(fNevents == 0) fNevents=(Int_t)gAlice->TreeE()->GetEntries(); 
  for(Int_t ievent=0;ievent<fNevents;ievent++)
    { 
      for (Int_t i=0; i<NumVolums; i++)
	for(Int_t j=0; j<NumSectors; j++)
	  for(Int_t ij=0; ij<NumRings; ij++)
	    PadADC[i][j][ij]=0;                    //zhachem ???
      gAlice->GetEvent(ievent) ;
      if(gAlice->TreeH()==0) return; 
      if(gAlice->TreeR()==0) gAlice->MakeTree("R");
      //Make branches
      AliFMDdigit  *fmdDigit;
      FMD->MakeBranch("R");
      
      Int_t zeroADC=1;
      // Int_t  threshold_array_size=30;

      // cout<<" AliFMDdigit "<<AliFMDdigit<<endl;
      if (FMD)
	{
	  gAlice->TreeD()->GetEvent(0); 
	  TClonesArray *FMDdigits=FMD->Digits();
	  Int_t nDigits=FMDdigits->GetEntries();
	   Int_t RecParticles[6];
	   Int_t nRecPart=0 ;
	   Int_t ZeroPads=0;
	   Int_t NumberOfPads=0; //To avoid warning
	   Int_t pedestal;
	   Float_t channelWidth=(22400*50)/1024;
	   for (Int_t digit=0;digit<nDigits;digit++) 
	     {
	       fmdDigit=(AliFMDdigit*)FMDdigits->UncheckedAt(digit);    
	       ivol=fmdDigit->Volume();
	       iSector=fmdDigit->NumberOfSector();
	       iRing=fmdDigit->NumberOfRing();
	       pedestal=Int_t(gRandom->Gaus(500,250));
	       PadADC[ivol-1][iSector-1][iRing-1]=
		 fmdDigit->ADCsignal()
		 -Int_t(Float_t(pedestal)/channelWidth);
	       if (PadADC[ivol-1][iSector-1][iRing-1]<0) 
		 PadADC[ivol-1][iSector-1][iRing-1]=0;
	     } //digit loop
	   Int_t Rmin=0; Int_t Rmax=0; //To avoid warning
	   Int_t Smin=0; Int_t Smax=0; //To avoid warning
		   
	   for (ivol=1; ivol<=NumVolums; ivol++)
	     {
	       if (ivol==1||ivol==3||ivol==5)
		 {
		   Ne1=size_EtaIntervalInner;
		   for(Int_t ieta=0; ieta<Ne1; ieta++)
		     EtaInterval[ieta]=EtaIntervalInner[ieta];
		 }
	       if (ivol==2||ivol==4)
		 {
		   Ne1=size_EtaIntervalOuter;
		   for( Int_t ieta=0; ieta<Ne1; ieta++)
		     EtaInterval[ieta]=EtaIntervalOuter[ieta];
		 }


		   for (Int_t e1=0;e1<Ne1-1;e1++) // vol<=NumVolums
		     {
		       Rmin=EtaInterval[e1];
		   Rmax=EtaInterval[e1+1];
		   ZeroPads=0;
		   Smin=0;
		   Smax=NumberOfSectors[ivol-1]; 
		   for (Int_t iring=Rmin; iring<Rmax; iring++) 
		     {
		       NumberOfPads=(Rmax-Rmin)*(Smax-Smin);
		       for 
			 (Int_t isector=1;
			  isector<=NumberOfSectors[ivol-1]; 
			  isector++) 
			 if(PadADC[ivol-1][isector-1][iring-1]<zeroADC)
			   ZeroPads++;
		     } //ring
		   /*
		   cout<<"\nRmin="<<Rmin;
		   cout<<"\nRmax="<<Rmax;
		   cout<<"\nSmin="<<Smin;
		   cout<<"\nSmax="<<Smax;
	       
		   cout<<"\nvolume "<<ivol<<" zero "<<ZeroPads<<
		     " NumberOfPads "<<NumberOfPads<<endl;
		   */
		   Float_t zerosRatio= 
		     (Float_t)ZeroPads/(Float_t)NumberOfPads;
		   //  cout<<"\nzerosRatio="<<zerosRatio;		      
		   RecParticles[0]=ivol;
		   RecParticles[1]=Smin;
		   RecParticles[2]=Smax;
		   RecParticles[3]=Rmin;
		   RecParticles[4]=Rmax;
		   if (zerosRatio>0.1 ||(ivol==2||ivol==4))
		     RecParticles[5]=
		       Determination_by_Poisson
		       (PadADC,ivol,Rmin,Rmax,Smin,Smax);
		   else
		     RecParticles[5]=
		       Determination_by_thresholds
		       (PadADC,ivol,Rmin,Rmax,Smin,Smax);
		   //  cout<<" nDeterm "<<RecParticles[5]<<endl;
		   new((*fReconParticles)[nRecPart++]) 
		 AliFMDReconstParticles(RecParticles);   	       
	} // eta
    } // volume
	   
	   // if(zerosRatio<0.01) Determination_by_counting (ZeroPads) ;	   
	}//if FMD
      gAlice->TreeR()->Reset();
      gAlice->TreeR()->Fill(); 
      gAlice->TreeR()->Write(0,TObject::kOverwrite);
    } //event loop
  // cout<<"\nAliFMDReconstruction::Exec finished"<<endl;
};

//------------------------------------------------------------
Int_t AliFMDReconstruction::Determination_by_thresholds
(Int_t PadADC[][50][800],Int_t volume, Int_t Rmin, Int_t Rmax, 
 Int_t Smin, Int_t Smax)
{
  
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
  Int_t threshold_array_size = 30;     
  Int_t NumPart=0;
  //Inner Si
  for (Int_t iring=Rmin; iring<Rmax; iring++) 
    {
      for (Int_t isector=Smin; isector<Smax; isector++) 
	{
	  for (int i=0;i<threshold_array_size-1;i++)
	    {	      
	      if(PadADC[volume-1][isector][iring]>threshold[i]
		 &&PadADC[volume-1][isector][iring]<=threshold[i+1])
		{
		  NumPart+=i;
		  break;
		}; // if in threshol interval
	    }; //threshold_array_size
	}; //iring
    }; //sector
  // cout<<"\nEnd threshol method"<<endl;
  return NumPart;
}


 //__________________________________________________________________

Int_t AliFMDReconstruction::Determination_by_Poisson (Int_t PadADC[][50][800],
					   Int_t vol, Int_t rmin, Int_t rmax, 
					   Int_t secmin, Int_t secmax)
{
  Int_t threshold_adc=1;   
  Int_t zeropads=0;
  for (Int_t i=rmin;i<rmax;i++)
    {
      for (Int_t j=secmin;j<secmax;j++)
	{
	  if (PadADC[vol-1][j][i]<threshold_adc) zeropads++;
	};
    };
  Float_t lambda=-TMath::Log(Float_t(zeropads)/
			     ( Float_t(secmax-secmin)*
			       Float_t(rmax-rmin))); //+1 zdes' ne nado
  Int_t Recon=(Int_t)(lambda*(secmax-secmin)*(rmax-rmin)); //+1 zdes' ne nado
  //  cout<<"\nEnd Poisson method"<<endl;
  return Recon;
};

//------------------------------------------------------------------




























