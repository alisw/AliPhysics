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

//_________________________________________________________________________
// This is a TTask that constructs SDigits out of Hits
// A Summable Digits is the sum of all hits in a cell
// A threshold is applied 
//
//-- Author: Alla Maevskaia(INR)
//////////////////////////////////////////////////////////////////////////////

// --- Standard library ---
#include <Riostream.h>
#include <stdlib.h>

// --- ROOT system ---
#include <TFile.h>
#include <TFolder.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTask.h>
#include <TTree.h>
#include <TVirtualMC.h>

// --- AliRoot header files ---
#include "AliConfig.h"
#include "AliDetector.h"
#include "AliFMD.h"
#include "AliFMDSDigitizer.h"
#include "AliFMDdigit.h"
#include "AliFMDhit.h"
#include "AliFMDv1.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"

ClassImp(AliFMDSDigitizer)

//____________________________________________________________________________ 
  AliFMDSDigitizer::AliFMDSDigitizer():TTask("AliFMDSDigitizer","") 
{
  // ctor
  fNevents = 0 ;     
  fSDigits = 0 ;
  fHits = 0 ;
  fRunLoader = 0;
}
           
//____________________________________________________________________________ 
AliFMDSDigitizer::AliFMDSDigitizer(const char* HeaderFile,char *SdigitsFile ):TTask("AliFMDSDigitizer","") 
{
  fNevents = 0 ;     // Number of events to digitize, 0 means all evens in current file
  // add Task to //root/Tasks folder
  fRunLoader = AliRunLoader::Open(HeaderFile);//Load event in default folder
  if (fRunLoader == 0x0)
   {
     Fatal("AliFMDSDigitizer","Can not open session. Header File is %s ",HeaderFile);
     return;//never reached
   }
  AliLoader* gime = fRunLoader->GetLoader("FMDLoader");
  if (gime == 0x0)
   {
     Fatal("AliFMDSDigitizer","Can not find FMD (loader) in specified event");
     return;//never reached
   }
  //add Task to //root/Tasks folder
  gime->PostSDigitizer(this);
}

//____________________________________________________________________________ 
  AliFMDSDigitizer::~AliFMDSDigitizer()
{
  // dtor
  AliLoader* gime = fRunLoader->GetLoader("FMDLoader");
  gime->CleanSDigitizer();
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
void AliFMDSDigitizer::Exec(Option_t *option) 
 { 
  Int_t NumberOfRings[5]=
  {fRingsSi1,fRingsSi2,fRingsSi1,fRingsSi2,fRingsSi1};
  Int_t NumberOfSectors[5]=
  {fSectorsSi1,fSectorsSi2,fSectorsSi1,fSectorsSi2,fSectorsSi1};

  if (fRunLoader)
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
  fRunLoader->LoadHeader();
  fRunLoader->LoadKinematics("READ");
  
  Int_t retval;

  retval = gime->LoadHits("READ");
  if (retval)
   {
     Error("Exec","Error occured while loading hits. Exiting.");
     return;
   }


  // Initialise Hit array
  fHits = new TClonesArray ("AliFMDhit", 1000);
  fSDigits = new TClonesArray ("AliFMDdigit", 1000);

  AliFMD *FMD = (AliFMD *) gAlice->GetDetector ("FMD");

  if (fNevents == 0)
    fNevents = (Int_t) fRunLoader->TreeE ()->GetEntries ();

  for (Int_t ievent = 0; ievent < fNevents; ievent++)
    {
      fRunLoader->GetEvent (ievent);

      TTree* TH = gime->TreeH();
      if (TH == 0x0)
       {
         Error("Exec","Can not get TreeH");
         return;
       }
      if (gime->TreeS () == 0)  gime->MakeTree("S");

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

         
         Stat_t ntracks = TH->GetEntries ();
         for (Int_t track = 0; track < ntracks; track++)
           {
             gAlice->ResetHits ();
             nbytes += TH->GetEvent(track);
             particle = fRunLoader->Stack()->Particle (track);
             Int_t nhits = FMDhits->GetEntriesFast ();

             for (hit = 0; hit < nhits; hit++)
              {
                fmdHit = (AliFMDhit *) FMDhits->UncheckedAt (hit);

                volume = fmdHit->Volume ();
                sector = fmdHit->NumberOfSector ();
                ring = fmdHit->NumberOfRing ();
                e = fmdHit->Edep ();
                de[volume][sector][ring] += e;
              }              //hit loop
           }                     //track loop
       }                     //if FMD


      Int_t digit[5];
      Float_t I = 1.664 * 0.04 * 2.33 / 22400;       // = 0.69e-6;
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

              }              // iring loop
           }                     //sector loop
       }                     // volume loop
      
      gime->TreeS()->Reset();
      gime->TreeS()->Fill();
      gime->WriteSDigits("OVERWRITE");
    }  //event loop


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
