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

/*
$Log$
Revision 1.3  2001/10/21 18:36:31  hristov
Several pointers were set to zero in the default constructors to avoid memory management problems

Revision 1.2  2001/09/24 16:41:15  hristov
New version of FMD code (A.Maevskaia)

Revision 1.1  2001/05/29 12:01:06  hristov
Last minute changes and new code for event mixing and reconstruction (A.Maevskaia)

Revision 1.3  2001/03/05 23:57:44  morsch
Writing of digit tree moved to macro.

Revision 1.2  2001/03/05 08:40:25  morsch
Method SortTracks(..) imported from AliMUON.

Revision 1.1  2001/02/02 14:11:53  morsch
AliMUONMerger prototype to be called by the merge manager.

*/

#include <TTree.h> 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>


#include "AliFMDMerger.h"
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
    TFile *file = new TFile(fFnBgr);
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

  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD") ;

  Int_t chargeSum[10][30][150];
  Int_t digit[5];
  Int_t ivol, iSector, iRing;

  Int_t NumberOfRings[5]=
  {fRingsSi1,fRingsSi2,fRingsSi1,fRingsSi2,fRingsSi1};
  Int_t NumberOfSectors[5]=
  {fSectorsSi1,fSectorsSi2,fSectorsSi1,fSectorsSi2,fSectorsSi1};

  TFile *f1 =0;
  TTree *TK = gAlice->TreeK();
  if (TK) f1 = TK->GetCurrentFile();

  gAlice->GetEvent(fEvNrSig) ;
  
  if(gAlice->TreeD() == 0)    	
    gAlice->MakeTree("D") ;
  gAlice->TreeD()->Reset();

  //Make branches 
   ReadDigit( chargeSum, fEvNrSig);

   if(fMerge){ 
    fBgrFile->cd();
    // gAlice->TreeS()->Reset();
    gAlice = (AliRun*)fBgrFile->Get("gAlice");
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

  gAlice->TreeD()->Reset();
  gAlice->TreeD()->Fill();
  
  fDigits   = FMD->Digits();
  
  gAlice->TreeD()->Write(0,TObject::kOverwrite) ;
  
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


