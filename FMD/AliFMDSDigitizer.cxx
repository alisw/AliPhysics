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
// This is a TTask that constructs SDigits out of Hits
// A Summable Digits is the sum of all hits in a cell
// A threshold is applied 
//
//-- Author: Alla Maevskaia(INR)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
// --- Standard library ---

// --- AliRoot header files ---

#include "AliFMDdigit.h"
#include "AliFMDhit.h"
#include "AliFMD.h"
#include "AliFMDv1.h"
#include "AliFMDSDigitizer.h"
#include "AliRun.h"
#include "AliDetector.h"

#include "TFile.h"
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include <stdlib.h>
#include <Riostream.h>

ClassImp(AliFMDSDigitizer)

//____________________________________________________________________________ 
  AliFMDSDigitizer::AliFMDSDigitizer():TTask("AliFMDSDigitizer","") 
{
  // ctor
  fNevents = 0 ;     
  fSDigits = 0 ;
  fHits = 0 ;

}
           
//____________________________________________________________________________ 
  AliFMDSDigitizer::AliFMDSDigitizer(char* HeaderFile,char *SdigitsFile ):TTask("AliFMDSDigitizer","") 
{
  fNevents = 0 ;     // Number of events to digitize, 0 means all evens in current file
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}

//____________________________________________________________________________ 
  AliFMDSDigitizer::~AliFMDSDigitizer()
{
  // dtor
}

//---------------------------------------------------------------------
void AliFMDSDigitizer::SetRingsSi1(Int_t ringsSi1)
{
  fRingsSi1=ringsSi1;
}
void AliFMDSDigitizer::SetSectorsSi1(Int_t sectorsSi1)
{
  fSectorsSi1=sectorsSi1;
}
void AliFMDSDigitizer::SetRingsSi2(Int_t ringsSi2)
{
  fRingsSi2=ringsSi2;
}
void AliFMDSDigitizer::SetSectorsSi2(Int_t sectorsSi2)
{
  fSectorsSi2=sectorsSi2;
}

//____________________________________________________________________________
void AliFMDSDigitizer::Exec(Option_t *option) { 




  Int_t NumberOfRings[5]=
  {fRingsSi1,fRingsSi2,fRingsSi1,fRingsSi2,fRingsSi1};
  Int_t NumberOfSectors[5]=
  {fSectorsSi1,fSectorsSi2,fSectorsSi1,fSectorsSi2,fSectorsSi1};

  // Initialise Hit array
  fHits = new TClonesArray ("AliFMDhit", 1000);
  fSDigits = new TClonesArray ("AliFMDdigit", 1000);

  AliFMD *FMD = (AliFMD *) gAlice->GetDetector ("FMD");

  if (fNevents == 0)
    fNevents = (Int_t) gAlice->TreeE ()->GetEntries ();

  for (Int_t ievent = 0; ievent < fNevents; ievent++)
    {
      gAlice->GetEvent (ievent);
      if (gAlice->TreeH () == 0)
	return;
      if (gAlice->TreeS () == 0)
	gAlice->MakeTree ("S");



      
            //Make branches
      char branchname[20];
       sprintf (branchname, "%s", FMD->GetName ());
      //Make branch for digits
        FMD->MakeBranch ("S");
    
       //Now made SDigits from hits, for PHOS it is the same
      Int_t volume, sector, ring, charge;
      Float_t e;
      Float_t de[10][50][150];
      Int_t ivol, isec, iring;
      Int_t hit, nbytes;
      TParticle *particle;
      AliFMDhit *fmdHit;
      TClonesArray *FMDhits = FMD->Hits ();

      // Event ------------------------- LOOP  

      for (ivol = 0; ivol < 10; ivol++)
	for (isec = 0; isec < 50; isec++)
	  for (iring = 0; iring < 150; iring++)
	    de[ivol][isec][iring] = 0;

      if (FMD)
	{
	  FMDhits = FMD->Hits ();
	  TTree *TH = gAlice->TreeH ();
	  Stat_t ntracks = TH->GetEntries ();
	  for (Int_t track = 0; track < ntracks; track++)
	    {
	      gAlice->ResetHits ();
	      nbytes += TH->GetEvent (track);
	      particle = gAlice->Particle (track);
	      Int_t nhits = FMDhits->GetEntriesFast ();

	      for (hit = 0; hit < nhits; hit++)
		{
		  fmdHit = (AliFMDhit *) FMDhits->UncheckedAt (hit);

		  volume = fmdHit->Volume ();
		  sector = fmdHit->NumberOfSector ();
		  ring = fmdHit->NumberOfRing ();
		  e = fmdHit->Edep ();
		  de[volume][sector][ring] += e;
		}		//hit loop
	    }			//track loop
	}			//if FMD


      Int_t digit[5];
      Float_t I = 1.664 * 0.04 * 2.33 / 22400;	// = 0.69e-6;
      for (ivol = 1; ivol < 6; ivol++)
	{ 
	  for (isec = 1; isec <= NumberOfSectors[ivol-1]; isec++)
	    { 
	      for (iring = 1; iring <= NumberOfRings[ivol-1]; iring++)
		{
		      digit[0] = ivol;
		      digit[1] = isec;
		      digit[2] = iring;
		      charge = Int_t (de[ivol][isec][iring] / I);

		      digit[3] = charge;
		      //dinamic diapason from MIP(0.155MeV) to 30MIP(4.65MeV)
		      //1024 ADC channels 
		      Float_t channelWidth = (22400 * 30) / 1024;

		      digit[4] = Int_t (digit[3] / channelWidth);
		      FMD->AddSDigit(digit);

		}		// iring loop
	    }			//sector loop
	}			// volume loop
      
      gAlice->TreeS()->Reset();
      gAlice->TreeS()->Fill();
      gAlice->TreeS()->Write(0,TObject::kOverwrite) ;
    }				//event loop


}
 
//__________________________________________________________________
void AliFMDSDigitizer::SetSDigitsFile(char * file ){
  if(!fSDigitsFile.IsNull())
    cout << "Changing SDigits file from " <<(char *)fSDigitsFile.Data() << " to " << file << endl ;
  fSDigitsFile=file ;
}
//__________________________________________________________________
void AliFMDSDigitizer::Print(Option_t* option)const
{
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  if(fSDigitsFile.IsNull())
    cout << " Writing SDigitis to file galice.root "<< endl ;
  else
    cout << "    Writing SDigitis to file  " << (char*) fSDigitsFile.Data() << endl ;

}
