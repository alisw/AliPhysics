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
#include <iostream.h>

// --- AliRoot header files ---

#include "AliFMDdigit.h"
#include "AliFMDhit.h"
#include "AliFMDReconstParticles.h"
#include "AliFMD.h"
#include "AliFMDv1.h"
#include "AliFMDReconstruction.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
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
 //Collects all digits in the same active volume into number of particles
  cout<<"\nStart AliFMDReconstruction::Exec(...)"<<endl;
  Int_t const NumVolums=5;
  Int_t const NumSectors=40;
  Int_t const NumRings=768;
  Int_t PadADC[10][50][800];
  Float_t eta, etain,etaout,rad,theta;
  Int_t ivol, iSector, iRing;
  Float_t rin[5]={4.2,15.4,4.2,15.4,4.2};
  Float_t rout[5]={17.4,28.4,17.4,28.4,17.4};
  Float_t z[5]={62.8, 75.2, -83.4, -75.2, -340.};
  Int_t NumberOfRings[5]={768,384,768,384,768};
  Int_t NumberOfSectors[5]=  {20,40,20,40,20};
  Int_t NumberOfEtaIntervals[5];
  // number of ring for boundary 0.1 eta
  TBranch *brDigits=0;

  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD");
  TClonesArray *fReconParticles=FMD->ReconParticles();
  if(fNevents == 0) fNevents=(Int_t)gAlice->TreeE()->GetEntries(); 
  cout<<" fNevents "<<fNevents<<endl;
  for(Int_t ievent=0;ievent<fNevents;ievent++)
    { 
      cout<<" ievent "<<ievent<<endl;
      for (Int_t i=0; i<NumVolums; i++)
	for(Int_t j=0; j<NumSectors; j++)
	  for(Int_t ij=0; ij<NumRings; ij++)
	    PadADC[i][j][ij]=0;                    //zhachem ???
      gAlice->GetEvent(ievent) ;
      if(gAlice->TreeH()==0) return; 
      if(gAlice->TreeD()==0) return;
	brDigits=gAlice->TreeD()->GetBranch("FMD");
	if (!brDigits){
	  cerr<<"EXEC Branch FMD digits not found"<<endl;
        return;
      } 
 
      if(gAlice->TreeR()==0) gAlice->MakeTree("R");
      //Make branches
      FMD->MakeBranch("R");

      
      Int_t zeroADC=1;
 
      AliFMDdigit  *fmdDigit;
       if (FMD)
	{
	  gAlice->TreeD()->GetEvent(0); 
	  TClonesArray *FMDdigits=FMD->Digits();
	  Int_t nDigits=FMDdigits->GetEntries();
	  cout<<" nDigits "<<nDigits<<endl;
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
	   AliHeader *header = gAlice->GetHeader();
	   AliGenEventHeader* genHeader = header->GenEventHeader();
	   TArrayF *o = new TArrayF(3); 
	   genHeader->PrimaryVertex(*o);
	   Float_t zVertex=o->At(2);
 	   for (ivol=0; ivol<NumVolums; ivol++)
	     {
	       Float_t realZ=zVertex+z[ivol];
	       theta=TMath::ATan(rin[ivol]/TMath::Abs(realZ));
	       etain = - TMath::Log( TMath::Tan(theta/2.));
	       theta=TMath::ATan(rout[ivol]/TMath::Abs(realZ));
	       etaout=- TMath::Log( TMath::Tan(theta/2.));
	       NumberOfEtaIntervals[ivol]=Int_t ((etain-etaout)*10);
	       eta=etain;
	       for (Int_t e1=0;e1<=NumberOfEtaIntervals[ivol];e1++) 
		 {
		   theta = 2.*TMath::ATan (TMath::Exp (-eta));
		   Float_t radmin = TMath::Abs(realZ) * (TMath::Tan (theta));
		   Rmin= Int_t ( (radmin-rin[ivol])*NumberOfRings[ivol]/(rout[ivol]-rin[ivol]));
		   eta=eta-0.1;
		   theta = 2.*TMath::ATan (TMath::Exp (-eta));
		   rad = TMath::Abs(realZ) * (TMath::Tan (theta));
		   Rmax=Int_t( (rad-rin[ivol])*NumberOfRings[ivol]/(rout[ivol]-rin[ivol]));
		   ZeroPads=0;
		   Smin=0;
		   Smax=NumberOfSectors[ivol]; 
		   for (Int_t iring=Rmin; iring<Rmax; iring++) 
		     {
		       NumberOfPads=(Rmax-Rmin)*(Smax-Smin);
		       for 
			 (Int_t isector=1;
			  isector<=NumberOfSectors[ivol]; 
			  isector++) 
			 if(PadADC[ivol][isector-1][iring-1]<zeroADC)
			   ZeroPads++;
		     } //ring
		   Float_t zerosRatio= 
		     (Float_t)ZeroPads/(Float_t)NumberOfPads;
		   RecParticles[0]=ivol;
		   RecParticles[1]=Smin;
		   RecParticles[2]=Smax;
		   RecParticles[3]=Rmin;
		   RecParticles[4]=Rmax;
		   if (zerosRatio>0.1 ||(ivol==2||ivol==4))
		     RecParticles[5]=
		       Determination_by_Poisson
		       (PadADC,ivol+1,Rmin,Rmax,Smin,Smax);
		   else
		     RecParticles[5]=
		       Determination_by_thresholds
		       (PadADC,ivol+1,Rmin,Rmax,Smin,Smax);
		   new((*fReconParticles)[nRecPart++]) 
		 AliFMDReconstParticles(RecParticles);   	       
		 } // eta
	     } // volume
	   
	}//if FMD
       gAlice->TreeR()->Reset();
       gAlice->TreeR()->Fill(); 
       gAlice->TreeR()->Write(0,TObject::kOverwrite);
    } //event loop
  cout<<"\nAliFMDReconstruction::Exec finished"<<endl;
};

//------------------------------------------------------------
Int_t AliFMDReconstruction::Determination_by_thresholds
(Int_t PadADC[][50][800],Int_t volume, Int_t Rmin, Int_t Rmax, 
 Int_t Smin, Int_t Smax)
{
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
  cout<<"\nEnd threshol method"<<endl;
  return NumPart;
}


 //__________________________________________________________________

Int_t AliFMDReconstruction::Determination_by_Poisson (Int_t PadADC[][50][800],
					   Int_t vol, Int_t rmin, Int_t rmax, 
					   Int_t secmin, Int_t secmax)
{
  //  cout<<"\nStart Poisson method";
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





























