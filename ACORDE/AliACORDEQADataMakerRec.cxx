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
//---
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.


//  Authors:
//
//  Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//  Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP)
//  Arturo Fernandez Tellez <afernan@mail.cern.ch (FCFM-BUAP)
//
//  Created: June 13th 2008
//---
// Last Update: Aug. 27th 2008 --> Implementation to declare QA expert histogram


// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TDirectory.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliACORDEdigit.h" 
#include "AliACORDEhit.h"
#include "AliACORDEQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliACORDERawReader.h"
#include "AliACORDERawStream.h"
ClassImp(AliACORDEQADataMakerRec)
           
//____________________________________________________________________________ 
AliACORDEQADataMakerRec::AliACORDEQADataMakerRec():AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kACORDE), "ACORDE Quality Assurance Data Maker")
{

}
//____________________________________________________________________________ 
AliACORDEQADataMakerRec::AliACORDEQADataMakerRec(const AliACORDEQADataMakerRec& qadm):AliQADataMakerRec() 
{
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}
//__________________________________________________________________
AliACORDEQADataMakerRec& AliACORDEQADataMakerRec::operator = (const AliACORDEQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliACORDEQADataMakerRec();
  new(this) AliACORDEQADataMakerRec(qadm);
  return *this;
}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kACORDE, task, list) ;
}

//____________________________________________________________________________
void AliACORDEQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

}
 
//____________________________________________________________________________ 
void AliACORDEQADataMakerRec::InitRaws()
{
  // create Raw histograms in Raw subdir

  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  char *acoModule[60]={"0_0","0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9",
                        "1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9",
                        "2_0","2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9",
                        "3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9",
                        "4_0","4_1","4_2","4_3","4_4","4_5","4_6","4_7","4_8","4_9",
                        "5_0","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9"};


  TH1F *fhACORDEBitPattern[4];
  fhACORDEBitPattern[0] = new TH1F("ACORDEBitPatternfromRAWSingle","Distribution of ACORDE fired modules from RAW-Single;Modules;Counts",60,1,60);//AcordeSingleMuon BitPattern
  fhACORDEBitPattern[1] = new TH1F("ACORDEBitPatternfromRAWMulti","Distribution of ACORDE fired modules from RAW-Multi;Modules;Counts",60,1,60);//AcordeMultiMuon BitPattern
  fhACORDEBitPattern[2] = new TH1F("ACORDEMultiplicityfromRAWSingle","Number of fired ACORDE modules;No. of fired ACORDE modules;No. of events in ACORDE",60,1,60);//AcordeSingleMuon Multiplicity
  fhACORDEBitPattern[3] = new TH1F("ACORDEMultiplicityfromRAWMulti","Number of fired ACORDE modules; No. of fired ACORDE modules;No. of events in ACORDE",60,1,60);//AcordeMultiMuon Multiplicity
  for(Int_t i=0;i<4;i++) 
    Add2RawsList(fhACORDEBitPattern[i],i,!expert, image, !saveCorr);
  
  for (Int_t iModule = 0; iModule<60; iModule++)
  {
	fhACORDEBitPattern[0]->GetXaxis()->SetBinLabel(iModule+1,acoModule[iModule]);
	fhACORDEBitPattern[1]->GetXaxis()->SetBinLabel(iModule+1,acoModule[iModule]);
  }

}
//____________________________________________________________________________ 
void AliACORDEQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  TH1F *    fhDigitsModule;
   char *acoModule[60]={"0_0","0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9",
                        "1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9",
                        "2_0","2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9",
                        "3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9",
                        "4_0","4_1","4_2","4_3","4_4","4_5","4_6","4_7","4_8","4_9",
                        "5_0","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9"};


  fhDigitsModule = new TH1F("ACORDEBitPatternfromDigits","Distribution of ACORDE from DIGITS;Modules;Counts",60,1,60);
  Add2DigitsList(fhDigitsModule,0,!expert,image);
  for (Int_t i=0;i<60;i++) fhDigitsModule->GetXaxis()->SetBinLabel(i+1,acoModule[i]); 
}

//____________________________________________________________________________ 

void AliACORDEQADataMakerRec::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
  // Not needed for ACORDE by now !!!
}

//____________________________________________________________________________
void AliACORDEQADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *    fhESDsSingle;
  TH1F *    fhESDsMulti;
   char *acoModule[60]={"0_0","0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9",
                        "1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9",
                        "2_0","2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9",
                        "3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9",
                        "4_0","4_1","4_2","4_3","4_4","4_5","4_6","4_7","4_8","4_9",
                        "5_0","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9"};


   fhESDsSingle = new TH1F("ACORDEBitPatternfromESDsSingle","Distribution of ACORDE fired modules from ESDs-Single;Modules;Counts",60,1,60);
   Add2ESDsList(fhESDsSingle,0,!expert,image);

   fhESDsMulti = new TH1F("ACORDEBitPatternfromESDsMulti","Distribution of ACORDE fired modules from ESDs-Multi;Modules;Counts",60,1,60);
   Add2ESDsList(fhESDsMulti,1,!expert,image);
	
   for (Int_t i=0;i<60;i++)
   {
	fhESDsSingle->GetXaxis()->SetBinLabel(i+1,acoModule[i]);
	fhESDsMulti->GetXaxis()->SetBinLabel(i+1,acoModule[i]);
   }


}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //fills QA histos for RAW

  rawReader->Reset();
  AliACORDERawStream rawStream(rawReader);
  size_t contSingle=0;
  size_t contMulti=0;
  UInt_t dy[4];

  bool kroSingle[60],kroMulti[60];
  UInt_t tmpDy;

  for(Int_t m=0;m<60;m++) {kroSingle[m]=0;kroMulti[m]=0;}

if(rawStream.Next())
{
        dy[0]=rawStream.GetWord(0);
        dy[1]=rawStream.GetWord(1);
        dy[2]=rawStream.GetWord(2);
        dy[3]=rawStream.GetWord(3);
        tmpDy=dy[0];
        for(Int_t r=0;r<30;++r)
        {
                kroSingle[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[1];
        for(Int_t r=30;r<60;++r)
        {
                kroSingle[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[2];
        for(Int_t r=0;r<30;++r)
        {
                kroMulti[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[3];
        for(Int_t r=30;r<60;++r)
        {
                kroMulti[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        contSingle=0;
	contMulti=0;
        for(Int_t r=0;r<60;++r)
        {
                if(kroSingle[r]==1)
                {
                        GetRawsData(0)->Fill(r+1);
                        contSingle=contSingle+1;
                }
		if(kroMulti[r]==1)
		{
			GetRawsData(1)->Fill(r+1);
			contMulti++;
		}

        }GetRawsData(2)->Fill(contSingle);GetRawsData(3)->Fill(contMulti);
}
}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::MakeDigits( TTree *digitsTree)
{
  //fills QA histos for Digits
  if (fDigitsArray) 
    fDigitsArray->Clear() ; 
  else 
    fDigitsArray = new TClonesArray("AliACORDEdigit",1000);
  TBranch * branch = digitsTree->GetBranch("ACORDEdigit");
  if (!branch) {
    AliWarning("ACORDE branch in Digits Tree not found");
  } else {
    branch->SetAddress(&fDigitsArray);
    for(Int_t track = 0 ; track < branch->GetEntries() ; track++) {
      branch->GetEntry(track);
      for(Int_t idigit = 0 ; idigit < fDigitsArray->GetEntriesFast() ; idigit++) {
        AliACORDEdigit *AcoDigit = (AliACORDEdigit*) fDigitsArray->UncheckedAt(idigit);
        if (!AcoDigit) {
          AliError("The unchecked digit doesn't exist");
          continue ;
        }
        GetDigitsData(0)->Fill(AcoDigit->GetModule()-1);
      }
    }
  }
}

//____________________________________________________________________________
void AliACORDEQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD

  AliESDACORDE * fESDACORDE= esd->GetACORDEData();
       
  for(int i=0;i<60;i++)
  {
  	if(fESDACORDE->GetHitChannel(i)) 
	{
		GetESDsData(0)->Fill(i);
		GetESDsData(1)->Fill(i);
	}
  }

}
