// --- ROOT system ---
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
// --- Standard library ---

// --- AliRoot header files ---

#include "AliFMDdigit.h"
#include "AliFMDReconstParticles.h"
#include "AliFMD.h"
#include "AliFMDv1.h"
#include "AliFMDReconstruction.h"
#include "AliRun.h"
#include "AliDetector.h"

#include "TROOT.h"
#include "TFolder.h"
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

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

void AliFMDReconstruction::Exec(TClonesArray *fReconParticles,Option_t *option) 
{ 
 //Collects all digits in the same active volume into number of particles

  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD");
  if(fNevents == 0) fNevents=(Int_t)gAlice->TreeD()->GetEntries(); 
  for(Int_t ievent=0;ievent<fNevents;ievent++)
    { 
      gAlice->GetEvent(ievent) ;
      if(gAlice->TreeH()==0) return; 
      if(gAlice->TreeR()==0) gAlice->MakeTree("R");
      //Make branches
      FMD->MakeBranch("R");
      Int_t threshold[]={   0,   18,  37,  56,   76, 
			    96, 119, 138, 165,  172, 
			    218, 231, 238, 277,  304, 
			    330, 357, 423, 449,  476, 
			    522, 542, 555, 568,  581, 
			    614, 624, 657, 674,  687, 
			    700, 713, 720, 727,  733, 
			    740, 759, 778, 797,  816, 
			    834, 853, 872, 891,  910, 
			    929, 948, 967, 986,  1024};
                         
      
      int threshold_array_size=sizeof(threshold)/sizeof(threshold[0]);
      AliFMDdigit  *fmdDigit;
       
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
  cout<<"\nAliFMDReconstruction::Exec(TClonesArray *fReconParticles,Option_t *option) finished"<<endl;
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





