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

//____________________________________________________________________________

void AliFMDReconstruction::Exec(Option_t *option) 
{ 
 //Collects all digits in the same active volume into number of particles

  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD");
  TClonesArray *fReconParticles=FMD->ReconParticles();
  if(fNevents == 0) fNevents=(Int_t)gAlice->TreeD()->GetEntries(); 
  for(Int_t ievent=0;ievent<fNevents;ievent++)
    { 
      gAlice->GetEvent(ievent) ;
      if(gAlice->TreeH()==0) return; 
      if(gAlice->TreeR()==0) gAlice->MakeTree("R");
      //Make branches
      FMD->MakeBranch("R");
      
      Int_t threshold[]={ 0,     14,  28,    42,   57, 
			  72,    89,  104,  124,  129, 
			  164,  174,  179,  208,  228, 
			  248,  268,   317,  337,  357, 
			  392,  407,  416,  426,  436, 
			  461,  468,  493,  506,  515, 
			  541,  566,  592,  617,  642, 
			  668,  693,  719,  744,  770, 
			  795,  821,  846,  871,  897, 
			  922,  948,  973,  999, 1024};
      

      /*
	Int_t threshold[]={    0,  18,  37,  56,   76, 
	                      96, 119, 138, 165,  172, 
               	             218, 231, 238, 277,  304, 
	                     330, 357, 423, 449,  476, 
	                     522, 542, 555, 568,  581, 
                             614, 624, 657, 674,  687, 
	                     700, 713, 720, 727,  733, 
	                     740, 759, 778, 797,  816, 
	                     834, 853, 872, 891,  910, 
	                     929, 948, 967, 986,  1024};
      */
      
      int threshold_array_size=sizeof(threshold)/sizeof(threshold[0]);
      AliFMDdigit  *fmdDigit;
      //    cout<<" AliFMDdigit "<<AliFMDdigit<<endl;
      if (FMD)
	{
	  gAlice->TreeD()->GetEvent(0); 
	  TClonesArray *FMDdigits=FMD->Digits();
	  Int_t nDigits=FMDdigits->GetEntries();
	   Int_t RecParticles[4];
	   Int_t nRecPart=0 ;
	   for (Int_t digit=0;digit<nDigits;digit++) 
	     {
	      fmdDigit=(AliFMDdigit*)FMDdigits->UncheckedAt(digit);    
	      RecParticles[0] = fmdDigit->Volume();
	      RecParticles[1] = fmdDigit->NumberOfSector();
	      RecParticles[2] = fmdDigit->NumberOfRing();
	      Int_t ADC=fmdDigit->ADCsignal(); 
	      RecParticles[3]=0; //case when fmdDigit->ADCsignal()==0
	      for (int i=0;i<threshold_array_size-1;i++)
		    {
		      if(ADC>threshold[i]&&ADC<=threshold[i+1])
			 RecParticles[3]=i;		 
		    }
	      new((*fReconParticles)[nRecPart++]) AliFMDReconstParticles(RecParticles); 
	     } //digit loop
	 }//if FMD
       gAlice->TreeR()->Reset();
       gAlice->TreeR()->Fill(); 
       gAlice->TreeR()->Write(0,TObject::kOverwrite);
    } //event loop
  cout<<"\nAliFMDReconstruction::Exec finished"<<endl;
}
//__________________________________________________________________

void AliFMDReconstruction::SetReconstParticlesFile(char * file )
{
   if (!fReconstParticlesFile.IsNull()) 
    cout<<"\nChanging reconstructed particles file from "<<
      (char *) fReconstParticlesFile.Data()<< " to "<<file<<endl;
    fReconstParticlesFile=file;
}
//__________________________________________________________________
void AliFMDReconstruction::Print(Option_t* option)const
{
  cout<<"------------------- "<<GetName() <<" -------------"<< endl ;
  if(fReconstParticlesFile.IsNull())
    cout<<"\nWriting reconstructed particles to file"<<endl ;
  else
    cout<<"\nWriting reconstructed particles to file  "<<
      (char*) fReconstParticlesFile.Data() << endl ;
}























