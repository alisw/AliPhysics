/**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
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
//Piotr.Skowronski@cern.ch
//Fast fixes to be able to compile with new Folder and I/O structure
//To be implemented correctly by the responsible person

//PH 20/05/2003 It seems this class is obsolete and not used anymore

/* $Id$ */

#include "AliFMDMerger.h"


#include <TTree.h> 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>

#include "AliDetector.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliFMD.h"
#include "AliFMDSDigitizer.h"
#include "AliFMDhit.h"
#include "AliFMDdigit.h"

#include "AliRun.h"
#include "AliPDG.h"

#include <stdlib.h>
#include <Riostream.h>
#include <Riostream.h>

ClassImp(AliFMDMerger)

//___________________________________________
  AliFMDMerger::AliFMDMerger() 
{
// Default constructor    
    fEvNrSig = 0;
    fEvNrBgr = 0;
    fBgrLoader = 0x0;
    fSigLoader = 0x0;
    fMerge   = kDigitize;
    fDigits  = 0;
    fSDigits = 0;
    fFnBgr   = 0;
    fFnSig   = 0;
    fBgrFile = 0;
}

//------------------------------------------------------------------------
AliFMDMerger::~AliFMDMerger()
{
// Destructor
  if(fSDigits)  {
    fSDigits->Delete();
    delete fSDigits ;
    fSDigits = 0;
  }
}


//---------------------------------------------------------------------
void AliFMDMerger::SetRingsSi1(Int_t ringsSi1)
{
  fRingsSi1=ringsSi1;
}
void AliFMDMerger::SetSectorsSi1(Int_t sectorsSi1)
{
  fSectorsSi1=sectorsSi1;
}
void AliFMDMerger::SetRingsSi2(Int_t ringsSi2)
{
  fRingsSi2=ringsSi2;
}
void AliFMDMerger::SetSectorsSi2(Int_t sectorsSi2)
{
  fSectorsSi2=sectorsSi2;
}
//---------------------------------------------------------------------

//------------------------------------------------------------------------
void AliFMDMerger::Init()
{
// Initialisation
    if (fMerge) fBgrFile = InitBgr();
}



//------------------------------------------------------------------------
TFile* AliFMDMerger::InitBgr()
{
// Initialise background event
    fBgrLoader= AliRunLoader::Open(fFnBgr);
    TFile *file = TFile::Open(fFnBgr);
// add error checking later
    printf("\n AliFMDMerger has opened %s file with background event \n", fFnBgr);
    return file;
}

//------------------------------------------------------------------------
void AliFMDMerger::Digitise()
{

    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !



#ifdef DEBUG
  cout<<"ALiFMDMerger::>SDigits2Digits start...\n";
#endif
  if (fBgrLoader == 0x0)
   {
     cerr<<"AliFMDMerger::Digitise : Background Run Loader is NULL"<<endl;
     return;
   }
   
  fBgrLoader->LoadgAlice();
  fBgrLoader->LoadHeader();
  fBgrLoader->LoadKinematics();
  
  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD") ;

  Int_t chargeSum[10][30][150];
  Int_t digit[5];
  Int_t ivol, iSector, iRing;

  Int_t NumberOfRings[5]=
  {fRingsSi1,fRingsSi2,fRingsSi1,fRingsSi2,fRingsSi1};
  Int_t NumberOfSectors[5]=
  {fSectorsSi1,fSectorsSi2,fSectorsSi1,fSectorsSi2,fSectorsSi1};

  TFile *f1 =0;
  TTree *TK = fBgrLoader->TreeK();
  if (TK) f1 = TK->GetCurrentFile();

  //just patches to be able to compile
  

  fBgrLoader->GetEvent(fEvNrSig) ;
  AliLoader* loader = fBgrLoader->GetLoader("FMDLoader");
  if (loader == 0x0)
   {
     cerr<<"AliFMDMerger::Digitise : Can not find loader for FMD. Exiting"<<endl;
     return;
   }
  
  Int_t retval;
  retval = loader->LoadDigits("UPDATE");
  if (retval == 0x0)
   {
     cerr<<"AliFMDMerger::Digitise : Error occured while loading digits. Exiting"<<endl;
     return;
   }
  
  if(loader->TreeD() == 0)           
    loader->MakeTree("D") ;
  
  loader->TreeD()->Reset();

  //Make branches 
   ReadDigit( chargeSum, fEvNrSig);

   if(fMerge){ 
//    fBgrFile->cd();
    // gAlice->TreeS()->Reset();
    gAlice = fBgrLoader->GetAliRun();
    Int_t chargeBgr[10][30][150];
    ReadDigit( chargeBgr,fEvNrBgr);
    for ( ivol=1; ivol<=5; ivol++)
      for ( iSector=1; iSector<=NumberOfSectors[ivol-1]; iSector++)
       for (  iRing=1; iRing<=NumberOfRings[ivol-1]; iRing++)
         chargeSum[ivol][iSector][iRing]=
           chargeBgr[ivol][iSector][iRing]+
           chargeSum[ivol][iSector][iRing];
    
   } //if merge


  // Put noise and make ADC signal
  for ( ivol=1; ivol<=5; ivol++){
    for ( iSector=1; iSector<=NumberOfSectors[ivol-1]; iSector++){
      for ( iRing=1; iRing<=NumberOfRings[ivol-1]; iRing++){
       digit[0]=ivol;
       digit[1]=iSector;
       digit[2]=iRing;
       digit[3]=PutNoise(chargeSum[ivol][iSector][iRing]);
       if(chargeSum[ivol][iSector][iRing] <= 500) digit[3]=500; 

    //dinamic diapason from MIP(0.155MeV) to 30MIP(4.65MeV)
    //1024 ADC channels 
       Float_t channelWidth=(22400*50)/1024;
       digit[4]=Int_t(digit[3]/channelWidth);
       if (digit[4]>1024) digit[4]=1024; 

       FMD->AddDigit(digit);

      } //ivol
    } //iSector
  } //iRing

  f1->cd();
   
  //Make branch for digits
  FMD->MakeBranch("D");

  loader->TreeD()->Reset();
  loader->TreeD()->Fill();
  
  fDigits   = FMD->Digits();//should be moved to specialized loader (AliFMDLoader)
  
  loader->WriteDigits("OVERWRITE");
  
  gAlice->ResetDigits();

}

//---------------------------------------------------------------------

void AliFMDMerger::ReadDigit(Int_t chargeSum[][30][150], Int_t iEvNum)
{
  AliFMDdigit *fmddigit;
  
  for (Int_t i=0; i<10; i++)
    for(Int_t j=0; j<30; j++)
      for(Int_t ij=0; ij<150; ij++)
       chargeSum[i][j][ij]=0;

  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD") ;
  
  gAlice->GetEvent(iEvNum) ;
  if(gAlice->TreeS()==0) {
    cout<<" TreeS==0 -> return"<<gAlice->TreeS()<<endl; 
    return ;}
  
  Int_t ndig, k;
  gAlice->ResetDigits();
  gAlice->TreeS()->GetEvent(iEvNum);
  TClonesArray * FMDdigits   = FMD->SDigits();
  
  ndig=FMDdigits->GetEntries();

  for (k=0; k<ndig; k++) {
    fmddigit= (AliFMDdigit*) FMDdigits->UncheckedAt(k);

    Int_t iVolume=fmddigit->Volume();
    Int_t iNumberOfSector =fmddigit->NumberOfSector();
    Int_t iNumberOfRing=fmddigit->NumberOfRing();
    chargeSum[iVolume][iNumberOfSector][iNumberOfRing]=fmddigit->Charge();
  }
}


