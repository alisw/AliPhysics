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
#include <Riostream.h>
#include <Riostream.h>

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
}

 //------------------------------------------------------------------------
/*
Bool_t AliFMDDigitizer::Init() 
{
 cout<<"AliFMDDigitizer::Init"<<endl;
 return kTRUE;
}
 
*/
//---------------------------------------------------------------------

void AliFMDDigitizer::Exec(Option_t* option)
{
  /*
   Conver hits to digits:
   - number of detector;
   - number of ring;
   - number of sector;
   - ADC signal in this channel
  */


#ifdef DEBUG
  cout<<"AliFMDDigitizer::>SDigits2Digits start...\n";
#endif

  //  cout<<" FMD "<<FMD<<endl;

   Int_t volume, sector, ring, charge;
  Float_t e;
  Float_t de[10][50][800];
  Int_t hit;
  Int_t digit[5];
  Int_t ivol, iSector, iRing;
  for (Int_t i=0; i<10; i++)
    for(Int_t j=0; j<50; j++)
      for(Int_t ij=0; ij<800; ij++)
	de[i][j][ij]=0;
  Int_t numberOfRings[5]=
  {768,384,768,384,768};
  Int_t numberOfSector[5]=
  {20,40,20,40,20}; 
  
  AliFMDhit *fmdHit=0;
  TTree *tH=0;
  TBranch *brHits=0;

  AliFMD * fFMD = (AliFMD *) gAlice->GetDetector("FMD") ;

// Loop over files to digitize

  Int_t nFiles=GetManager()->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles;
       inputFile++) {

    cout<<" event "<<fManager->GetOutputEventNr()<<endl;
    if (fFMD)
    {
      TClonesArray *FMDhits = fFMD->Hits ();
      tH = fManager->GetInputTreeH(inputFile);
      brHits = tH->GetBranch("FMD");
      if (brHits) {
	fFMD->SetHitsAddressBranch(brHits);
      }else{
	cerr<<"EXEC Branch FMD hit not found"<<endl;
	exit(111);
      } 
      Int_t ntracks    = (Int_t) tH->GetEntries();

        for (Int_t track = 0; track < ntracks; track++)
	{
	  brHits->GetEntry(track);
	  Int_t nhits = FMDhits->GetEntries ();

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

 
  // Put noise and make ADC signal
   Float_t iP = 1.664 * 0.04 * 2.33 / 22400;	// = 0.69e-6;
 for ( ivol=1; ivol<=5; ivol++){
    for ( iSector=1; iSector<=numberOfSector[ivol-1]; iSector++){
      for ( iRing=1; iRing<=numberOfRings[ivol-1]; iRing++){
	digit[0]=ivol;
	digit[1]=iSector;
	digit[2]=iRing;
	charge = Int_t (de[ivol][iSector][iRing] / iP);
	digit[3]=PutNoise(charge);
	if(digit[3]<= 500) digit[3]=500; 
    //dynamic range from MIP(0.155MeV) to 30MIP(4.65MeV)
    //1024 ADC channels 
	Float_t channelWidth=(22400*50)/1024;
	digit[4]=Int_t(digit[3]/channelWidth);
	if (digit[4]>1024) digit[4]=1024; 
	fFMD->AddDigit(digit);
      } //ivol
    } //iSector
  } //iRing

  TTree* treeD = fManager->GetTreeD();
  treeD->Clear();
  treeD->Reset();
  fFMD->MakeBranchInTreeD(treeD);
  treeD->Fill();
 
  fManager->GetTreeD()->Write(0,TObject::kOverwrite);
  
  gAlice->ResetDigits();
  }
}

