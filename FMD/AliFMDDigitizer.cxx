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


#include <TTree.h> 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TRandom.h>


#include "AliFMDDigitizer.h"
#include "AliFMD.h"
#include "AliFMDSDigitizer.h"
#include "AliFMDhit.h"
#include "AliFMDdigit.h"
#include "AliRunDigitizer.h"

#include "AliRun.h"
#include "AliPDG.h"

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

ClassImp(AliFMDDigitizer)

//___________________________________________
  AliFMDDigitizer::AliFMDDigitizer()  :AliDigitizer()
{
// Default ctor - don't use it
  ;
}

//___________________________________________
AliFMDDigitizer::AliFMDDigitizer(AliRunDigitizer* manager) 
    :AliDigitizer(manager) 
{
	cout<<"AliFMDDigitizer::AliFMDDigitizer"<<endl;
// ctor which should be used
//  fDebug =0;
 // if (GetDebug()>2)
  //  cerr<<"AliFMDDigitizer::AliFMDDigitizer"
   //     <<"(AliRunDigitizer* manager) was processed"<<endl;

}

//------------------------------------------------------------------------
AliFMDDigitizer::~AliFMDDigitizer()
{
// Destructor
  if(fSDigits)  {
    fSDigits->Delete();
    delete fSDigits ;
    fSDigits = 0;
  }
}

 //------------------------------------------------------------------------
Bool_t AliFMDDigitizer::Init()
{
// Initialization
 cout<<"AliFMDDigitizer::Init"<<endl;
 return kTRUE;
}
 

//---------------------------------------------------------------------

void AliFMDDigitizer::Exec(Option_t* option)
{



#ifdef DEBUG
  cout<<"AliFMDDigitizer::>SDigits2Digits start...\n";
#endif

  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD") ;
  cout<<" FMD "<<FMD<<endl;
  if (gAlice->TreeD () == 0)
    gAlice->MakeTree ("D");
  fDigits   = FMD->Digits();
  

  Int_t chargeSum[10][50][300];
  Int_t digit[5];
  Int_t ivol, iSector, iRing;
  for (Int_t i=0; i<10; i++)
    for(Int_t j=0; j<50; j++)
      for(Int_t ij=0; ij<300; ij++)
	chargeSum[i][j][ij]=0;
  Int_t NumberOfRings[5]=
  {256,128,256,128,256};
  Int_t NumberOfSectors[5]=
  {20,40,20,40,20}; 
  

// Loop over files to digitize

  for (Int_t inputFile=0; inputFile<fManager->GetNinputs();
       inputFile++) {

    ReadDigit( chargeSum, inputFile);

  }
  // Put noise and make ADC signal
  for ( ivol=1; ivol<=5; ivol++){
    for ( iSector=1; iSector<=NumberOfSectors[ivol-1]; iSector++){
      for ( iRing=1; iRing<=NumberOfRings[ivol-1]; iRing++){
	digit[0]=ivol;
	digit[1]=iSector;
	digit[2]=iRing;
	digit[3]=PutNoise(chargeSum[ivol][iSector][iRing]);
	if(chargeSum[ivol][iSector][iRing] <= 500) digit[3]=500; 
    //dynamic range from MIP(0.155MeV) to 30MIP(4.65MeV)
    //1024 ADC channels 
	Float_t channelWidth=(22400*50)/1024;
	digit[4]=Int_t(digit[3]/channelWidth);
	if (digit[4]>1024) digit[4]=1024; 
	FMD->AddDigit(digit);

      } //ivol
    } //iSector
  } //iRing
  TTree* treeD = fManager->GetTreeD();
  treeD->Reset();
  FMD->MakeBranch("D");
  treeD->Fill();
 
  fManager->GetTreeD()->Write(0,TObject::kOverwrite);
  
  gAlice->ResetDigits();

}

//---------------------------------------------------------------------

void AliFMDDigitizer::ReadDigit(Int_t chargeSum[][50][300], Int_t inputFile )
{
  cout<<" AliFMDDigitizer::ReadDigit "<<endl;
  AliFMDdigit *fmddigit;
  gAlice->GetEvent(0) ;
   AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD") ;
  Int_t ndig, k;
  gAlice->ResetDigits();
  treeS->GetEvent(0);

  treeS->GetEvent(0);
  TClonesArray * FMDSdigits   = FMD->SDigits();
  
  ndig=FMDSdigits->GetEntries();


  for (k=0; k<ndig; k++) {
    fmddigit= (AliFMDdigit*) FMDSdigits->UncheckedAt(k);
    Int_t iVolume=fmddigit->Volume();
    Int_t iNumberOfSector =fmddigit->NumberOfSector();
    Int_t iNumberOfRing=fmddigit->NumberOfRing();
    chargeSum[iVolume][iNumberOfSector][iNumberOfRing]+=fmddigit->Charge();
  }
}


