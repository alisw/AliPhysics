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
// This is a class that constructs ReconstParticles (reconstructed particles) 
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
#include "TH2F.h"
#include "TRandom.h"

// --- Standard library ---
#include <stdlib.h>
#include <Riostream.h>

// --- AliRoot header files ---

#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliFMDdigit.h"
#include "AliFMDReconstParticles.h"
#include "AliFMDReconstructor.h"
#include "AliRun.h"
#include "AliConfig.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"

ClassImp(AliFMDReconstructor)

        
//____________________________________________________________________________

void AliFMDReconstructor::Reconstruct(AliRunLoader* runLoader) const
{ 
 //Collects all digits in the same active volume into number of particles
  /*
    Reconstruct number of particles 
    in given group of pads for given FMDvolume
    determine by numberOfVolume , 
    numberOfMinSector,numberOfMaxSector,
    numberOfMinRing, numberOgMaxRing
    Reconstructor method choose dependence on number of empty pads  
  */


#ifdef DEBUG
  Info("Exec ","Start");
#endif

  Int_t const knumVolumes=5;
  Int_t padADC;
  Float_t eta, etain,etaout,rad,theta;
  Int_t ivol, iSector, iRing;
  Float_t rin[5]={4.2,15.4,4.2,15.4,4.2};
  Float_t rout[5]={17.4,28.4,17.4,28.4,17.4};
  Float_t z[5]={-62.8, -75.2, 83.4, 75.2, 340.};
  Int_t numberOfRings[5]={512,256,512,256,512};
  Int_t numberOfSectors[5]=  {20,40,20,40,20};
  Int_t numberOfEtaIntervals[5];
  // number of ring for boundary 0.1 eta

  
  if (runLoader == 0x0)
   {
    Error("Exec","Run Loader loader is NULL - Session not opened");
    return;
   }

  AliLoader* plFMD = runLoader->GetLoader("FMDLoader");
  if (plFMD == 0x0)
   {
     Fatal("AliFMDReconstructor","Can not find FMD (loader) in specified event");
     return;//never reached
   }
   
  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  if (!runLoader->TreeE()) runLoader->LoadHeader();

  TDirectory* cwd = gDirectory;
  gDirectory = 0x0;
  Text_t buf1[20];
  TH2F* hTotal[10];
  for (Int_t j=1; j<=5; j++){
    sprintf(buf1,"hTotal%d",j);
    
    hTotal[j] = new TH2F(buf1," Number of primary particles ",
			 numberOfSectors[j-1],1,numberOfSectors[j-1],
			 numberOfRings[j-1],1,numberOfRings[j-1]);
  }   
  gDirectory = cwd;
 
  plFMD->LoadRecPoints("RECREATE");
  TClonesArray* reconParticles = new TClonesArray("AliFMDReconstParticles"); 
  TClonesArray* digits = new TClonesArray("AliFMDDigit");
 
  Int_t retval=0;     
  Int_t nevents=Int_t (runLoader->TreeE()->GetEntries()); 
#ifdef DEBUG
  cout<<" nevents "<<nevents<<endl;
#endif
   for(Int_t ievent=0;ievent<nevents;ievent++)
    { 
#ifdef DEBUG
      cout<<" *** event "<<ievent<<endl; 
#endif
      runLoader->GetEvent(ievent) ;
      //event z-vertex for correction eta-rad dependence      
      AliHeader *header = runLoader->GetHeader();
      AliGenEventHeader* genHeader = header->GenEventHeader();
      TArrayF *o = new TArrayF(3); 
      if (genHeader) genHeader->PrimaryVertex(*o);
      Float_t zVertex=o->At(2);

      for (Int_t i=0; i<5; i++)
        hTotal[i+1]->Reset();
      
      retval = plFMD->LoadDigits("READ"); 
      if (retval)
	{
	  Error("Exec","Error occured while loading digits. Exiting.");
	  return;
	}
      
      TTree* treeD = plFMD->TreeD();
      if (treeD == 0x0)
	{
	  Error("Exec","Can not get Tree with Digits. Nothing to reconstruct - Exiting");
	  return;
	}
      
      TBranch *brDigits=0;
	    
      brDigits=treeD->GetBranch("FMD");

      if (brDigits) {
	brDigits->SetAddress(&digits);
      }else{
	cerr<<"EXEC Branch FMD digits not found"<<endl;
	return;
      } 
	  
      if(plFMD->TreeR()==0) plFMD->MakeTree("R");

      //Make branches
      const Int_t kBufferSize = 16000;
      plFMD->TreeR()->Branch("FMD", &reconParticles, kBufferSize);

      
      Int_t zeroADC=6;
      // read Digits 
      AliFMDdigit  *fmdDigit;
       if (plFMD->TreeD()->GetEvent(0))
	{	  	  
	  Int_t nDigits=digits->GetEntries();
	  Int_t recParticles[6];
	  Int_t nRecPart=0 ;
	  Int_t zeroPads=0;
	  Int_t numberOfPads=0;
	  Int_t pedestal;
	  Float_t channelWidth=(22400*50)/1024;
	  for (Int_t digit=0;digit<nDigits;digit++) 
	    {
	      fmdDigit=(AliFMDdigit*)digits->UncheckedAt(digit);    
	      ivol=fmdDigit->Volume();
	      iSector=fmdDigit->NumberOfSector();
	      iRing=fmdDigit->NumberOfRing();
	      pedestal=Int_t(gRandom->Gaus(500,250));
	      padADC= fmdDigit->ADCsignal()
		-Int_t(Float_t(pedestal)/channelWidth);
	      if (padADC<0) padADC=0;
              hTotal[ivol]->Fill(iSector,iRing,padADC);
	    } //digit loop

	  //reconstruct multiplicity in 0.1 eta according Poisson distribution 

	  Int_t rmin=0; Int_t rmax=0; 
	  Int_t smin=0; Int_t smax=0; 
	  for (ivol=0; ivol<knumVolumes; ivol++)
	    {
	      Float_t ring2number=Float_t (numberOfRings[ivol])/
		(rout[ivol]-rin[ivol]);
	      Float_t realZ=zVertex+z[ivol];
	      theta=TMath::ATan(rout[ivol]/TMath::Abs(realZ));
	      etain = - TMath::Log( TMath::Tan(theta/2.));
	      theta=TMath::ATan(rin[ivol]/TMath::Abs(realZ));
	      etaout=- TMath::Log( TMath::Tan(theta/2.));
	      numberOfEtaIntervals[ivol]=Int_t((etaout-etain)*10)-1;
	      eta=etain;
	      for (Int_t e1=0;e1<=numberOfEtaIntervals[ivol];e1++) 
		{
		  theta = 2.*TMath::ATan (TMath::Exp (-eta));
		  Float_t radmin = TMath::Abs(realZ) * (TMath::Tan (theta));
		  rmax= Int_t ( (radmin-rin[ivol])*ring2number+0.5);
		  eta=eta+0.1;
		  theta = 2.*TMath::ATan (TMath::Exp (-eta));
		  rad = TMath::Abs(realZ) * (TMath::Tan (theta));
		  rmin=Int_t( (rad-rin[ivol])*ring2number+0.5);
		  
		  zeroPads=0;
		  smin=0;
		  smax=numberOfSectors[ivol]; 
		  numberOfPads=(rmax-rmin)*(smax-smin);
		  for (Int_t iring=rmin; iring<rmax; iring++) 
		    {
		      for 
			(Int_t isector=0;
			 isector<numberOfSectors[ivol]; 
			 isector++) 
			{
			  if(hTotal[ivol+1]->GetBinContent(isector+1,iring+1)
			     <zeroADC) zeroPads++;}
		      
		    } //ring
		  Float_t numberOfPads=Float_t(smax-smin)*Float_t(rmax-rmin);

		  Double_t lambda=-TMath::Log(Double_t(zeroPads)/numberOfPads);
		  Int_t fRecon=Int_t (lambda*numberOfPads+0.5);
		  recParticles[0]=ivol+1;
		  recParticles[1]=smin;
		  recParticles[2]=smax;
		  recParticles[3]=rmin;
		  recParticles[4]=rmax;
		  recParticles[5]= fRecon;
		  new((*reconParticles)[nRecPart++])   AliFMDReconstParticles(recParticles);             
		  

		 } // eta
	     } // volume
	   
	}//if (plFMD->TreeD()->GetEvent(0))
       plFMD->TreeR()->Reset();
       plFMD->TreeR()->Fill(); 
       plFMD->WriteRecPoints("OVERWRITE");
       plFMD->UnloadDigits();
    } //event loop
  plFMD->UnloadRecPoints();
#ifdef DEBUG
  Info(" Exec"," finished");
#endif
  //  delete hTotal[10]; 

}


//_____________________________________________________________________________
void AliFMDReconstructor::FillESD(AliRunLoader* /*runLoader*/, 
				  AliESD* /*esd*/) const
{
// nothing to be done

}

