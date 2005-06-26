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

#include <Riostream.h>

#include <TDirectory.h>

#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliSTARTLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARTReconstructor.h"
#include <AliESD.h>
#include <TClonesArray.h>
#include "AliSTARTRecPoint.h"
#include "AliRawReader.h"
#include "AliSTARTRawReader.h"
#include "AliLog.h"
ClassImp(AliSTARTReconstructor)

  void  AliSTARTReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
  //START raw data-> digits conversion
 // reconstruct time information from raw data
  AliSTARTRawReader myrawreader(rawReader,digitsTree);
  myrawreader.NextThing();
}
  void AliSTARTReconstructor::Reconstruct(TTree*digitsTree, TTree*clustersTree) const
{
// START digits reconstruction
// STARTRecPoint writing 

  AliDebug(1,Form("Start DIGITS reconstruction "));
   Int_t channelWigth=25; //ps
  TArrayI* fSumMult = new TArrayI(6);
  Float_t ph2mV = 150./500.;
  Float_t mV2channel=200000/(25*25);  //5V -> 200ns

  TBranch *brDigits=digitsTree->GetBranch("START");
  AliSTARTdigit *fDigits = new AliSTARTdigit();
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    cerr<<"EXEC Branch START digits not found"<<endl;
    return;
  }
  brDigits->GetEntry(0);
  Int_t besttimeright=channelWigth * (fDigits->BestTimeRight());
  Int_t besttimeleft=channelWigth * (fDigits->BestTimeLeft());

  //folding with experimental time distribution
  //  Float_t c = 29.9792; // cm/ns
  Float_t c = 0.0299792; // cm/ps
  Float_t lenr=TMath::Sqrt(350*350 + 6.5*6.5);
  Float_t lenl=TMath::Sqrt(69.7*69.7 + 6.5*6.5);
  Float_t timeDiff=channelWigth * (fDigits->TimeDiff());
  Int_t meanTime=channelWigth * (fDigits->MeanTime());
  Float_t ds=(c*(timeDiff)-(lenr-lenl))/2;
  AliDebug(2,Form(" timediff in ns %f  real point%f",timeDiff,ds));
  

  fDigits->GetSumMult(*fSumMult);
  Int_t multipl[6]; 
 for (Int_t i=0; i<6; i++)
    {
      Float_t  mult=Float_t (fSumMult->At(i));
      Float_t   realMultmV=TMath::Exp(mult/mV2channel);
      multipl[i]=Int_t ((realMultmV/ph2mV)/500+0.5);
    }
  AliDebug(2,Form(" multiplicity Abs side %i  multiplicity non-Abs side %i",multipl[1],multipl[2]));

  AliSTARTRecPoint* frecpoints= new AliSTARTRecPoint ();
  clustersTree->Branch( "START", "AliSTARTRecPoint" ,&frecpoints, 405,1);
  frecpoints->SetTimeBestRight(besttimeright);
  frecpoints->SetTimeBestLeft(besttimeleft);
  frecpoints->SetVertex(ds);
  frecpoints->SetMeanTime(meanTime);
  frecpoints->SetMult(multipl[0]);
  frecpoints->SetMultA(multipl[2]);
  frecpoints->SetMultC(multipl[1]);

  clustersTree->Fill();
}


void AliSTARTReconstructor::FillESD(AliRunLoader* runLoader, AliESD *pESD) const
{

  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/
  
  //  Float_t c = 0.3;  //speed of light mm/ps
  // Float_t Zposition=0;
  
  if (!runLoader) {
    AliError("Reconstruct >> No run loader");
    return;
  }
  
  AliDebug(1,Form("Start FillESD START"));

  AliSTARTLoader* pStartLoader = (AliSTARTLoader*) runLoader->GetLoader("STARTLoader");
 
  pStartLoader->LoadRecPoints("READ");
  
    TTree *treeR = pStartLoader->TreeR();
  
   AliSTARTRecPoint* frecpoints= new AliSTARTRecPoint ();
    if (!frecpoints) {
    AliError("Reconstruct Fill ESD >> no recpoints found");
    return;
  }
  
  AliDebug(1,Form("Start FillESD START"));
   TBranch *brRec = treeR->GetBranch("START");
    if (brRec) {
      brRec->SetAddress(&frecpoints);
    }else{
      cerr<<"EXEC Branch START rec not found"<<endl;
      exit(111);
    } 
 
    brRec->GetEntry(0);
    Float_t Zposition=frecpoints->GetVertex();
    
    pESD->SetT0zVertex(Zposition);
    pStartLoader->UnloadRecPoints();
    
} // vertex in 3 sigma






