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
// Last update: Sept. 29th. 2011 (by Mario RC: mrodrigu@mail.cern.ch) 
//	--> ACOMultiSL0_DQM_Shifter filling histogram fixed
//	--> Expert histogram updated: 2 histograms (Checks the hits for SL0 and AMU mode)
//	--> To be include in the next update: threshold settings from AliACORDEQAThreshold class (not yet)
// Last update: May 5th. 2011 (by Mario RC: mrodrigu@mail.cern.ch) -->Creates QA expert histograms 
// and QA-shifter histograms also with threshold lines and visual alarm
// Last Update: Aug. 27th 2008 --> Implementation to declare QA expert histogram
// Last update: Nov. 14t 2009 --> MRC <mrodrigu@mail.cern.ch> (FCFM-BUAP) 
//			|--> Change in Multiplicity histogram for AMORE (to detect empty triggers events of ACORDE)



// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TDirectory.h>
#include <TPaveText.h>
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
AliACORDEQADataMakerRec::AliACORDEQADataMakerRec(const AliACORDEQADataMakerRec& qadm):
  AliQADataMakerRec()
{
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliACORDEQADataMakerRec::~AliACORDEQADataMakerRec()
{
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
  ResetEventTrigClasses(); // reset triggers list to select all histos
  // Update for DQM GUI
  //

  for (Int_t specie = 0; specie < AliRecoParam::kNSpecies ; specie++) {
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) continue ;
    // 
    // RS Set event specie
    SetEventSpecie(AliRecoParam::ConvertIndex(specie));
    //
    for (int itc=-1;itc<GetNTrigClasses();itc++) { // RS Loop over the trigger classes
      //
      	TObjArray * parr = GetRawsDataOfTrigClass(itc);
      	if (!parr) continue;
      	TObjArray &harr = *parr;
      	TH1* h0 = (TH1*)harr[0];
      	TH1* h2 = (TH1*)harr[2];
      	if (!h0 || !h2) continue;

      	Double_t integralSL0 = 0;
      	Double_t integralAMU = 0;

	if ((itc==-1) && ( (integralAMU=h2->Integral()==0) || (integralSL0=h0->Integral()==0) ) ) continue;

	Float_t maxSL0 = 1.1*h0->GetMaximum();
	Float_t scaleSL0 = 0.;
	if (maxSL0!=0) scaleSL0 = 1./maxSL0;
	else scaleSL0 = 1.;
	h0->Scale(scaleSL0);


	Float_t maxAMU = 1.1*h2->GetMaximum();
	Float_t scaleAMU = 0.;
	if (maxAMU!=0) scaleAMU = 1./maxAMU;
	else scaleAMU = 1.;
	h2->Scale(scaleAMU);
      
	Int_t indexActiveModuleSL0 = 0;
	Int_t indexActiveModuleAMU = 0;

	for(Int_t iModule = 0; iModule < 60; iModule++){
		if (h0->GetBinContent(iModule) > 0) indexActiveModuleSL0++;
		if (h2->GetBinContent(iModule) > 0) indexActiveModuleAMU++;
	}
	
	Float_t meanHitsSL0 = 0.;
	Float_t meanHitsAMU = 0.;
	if ((indexActiveModuleSL0==0) || (indexActiveModuleAMU == 0)) continue;

	meanHitsSL0 = h0->Integral()/indexActiveModuleSL0;
	meanHitsAMU = h2->Integral()/indexActiveModuleAMU;

        TH1* h4 = (TH1*)harr[4];
        TH1* h5 = (TH1*)harr[5];

	if (h4 && h5 && meanHitsAMU!=0 && meanHitsSL0!=0){
		for (Int_t iModule = 0; iModule < 60; iModule++){
			h4->Fill(h0->GetBinContent(iModule)/meanHitsSL0-1);
			h5->Fill(h2->GetBinContent(iModule)/meanHitsAMU-1);
		}
	}

    } // end of trigger classes loop
  } // end specie loop
  // QA Checker standar (to be updated)
  //
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
	TH1F * fhACORDEBitPatternDQM = new TH1F("ACOHitsSL0_DQM_Shifter","Distribution of ACORDE fired modules for DQM shifter; No. of module; Counts",60,-0.5,59.5);// Hits histogram for QA-shifter ACO-SL0 trigger mode
        TH1F * fhACORDEMultiplicitySL0DQM = new TH1F("ACOMultiSL0_DQM_Shifter","Multiplicity of ACORDE fired modules for DQM shifter; No. of fired modules; No. of events",62,-1,60); // Multiplicity histo. for QA-shifter ACO-SL0 trigger mode
        TH1F * fhACORDEBitPatternAMUDQM = new TH1F("ACOHitsAMU_DQM_Shifter","Distribution of ACORDE fired modules for DQM shifter; No. of module; Counts",60,-0.5,59.5);// Hits histogram for QA-shifter ACO-SL0 trigger mode
        TH1F * fhACORDEMultiplicityAMUDQM = new TH1F("ACOMultiAMU_DQM_Shifter","Multiplicity of ACORDE fired modules for DQM shifter; No. of fired modules; No. of events",62,-1,60); // Multiplicity histo. for QA-shifter ACO-SL0 trigger mode
        TH1F * fhACORDEBitPatternCheckDQMSL0 = new TH1F("ACOHitsTriggerCheck_DQMExpertSL0","Check the activity of ACORDE's modules; Hits per module/mean of Hits; Counts",100,-3,5); // Check the trigger status of ACORDE (SL0 vs AMU)
        TH1F * fhACORDEBitPatternCheckDQMAMU = new TH1F("ACOHitsTriggerCheck_DQMExpertAMU","Check the activity of ACORDE's modules; Hits per module/mean of Hits; Counts",100,-3,5); // Check the trigger status of ACORDE (SL0 vs AMU)
         // Expert histograms
         // Check the hits multiplicity from trigger configuration
         Add2RawsList(fhACORDEBitPatternCheckDQMSL0,4,expert,image,!saveCorr);
         fhACORDEBitPatternCheckDQMSL0->SetFillColor(kOrange);
	 Add2RawsList(fhACORDEBitPatternCheckDQMAMU,5,expert,image,!saveCorr);
         fhACORDEBitPatternCheckDQMAMU->SetFillColor(kRed+2);

	
	// AMORE diplay settings for shifter on GUI
 
        // For SL0 ACO trigger mode
 
         Add2RawsList(fhACORDEBitPatternDQM,0,!expert,image,!saveCorr);
	 ForbidCloning(fhACORDEBitPatternDQM);
         Add2RawsList(fhACORDEMultiplicitySL0DQM,1,!expert,image,!saveCorr);
 	 ForbidCloning(fhACORDEMultiplicitySL0DQM);

         // For Hits distribution on ACORDE
 
         fhACORDEBitPatternDQM->SetFillColor(kMagenta+2);
 
         // For ACORDE Multiplicity distribution of fired modules
 
         fhACORDEMultiplicitySL0DQM->SetFillColor(kMagenta);
 
         // For AMU ACO trigger mode
 
         Add2RawsList(fhACORDEBitPatternAMUDQM,2,!expert,image,!saveCorr);
         Add2RawsList(fhACORDEMultiplicityAMUDQM,3,!expert,image,!saveCorr);
 
         // For Hits distribution on ACORDE
 
         fhACORDEBitPatternAMUDQM->SetFillColor(kViolet+7);
 
         // For ACORDE Multiplicity distribution of fired modules
 
         fhACORDEMultiplicityAMUDQM->SetFillColor(kViolet+6);
 
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}

//____________________________________________________________________________ 
void AliACORDEQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  TH1F *    fhDigitsModule;
  const char *acoModule[60]={"0_0","0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9",
                        "1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9",
                        "2_0","2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9",
                        "3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9",
                        "4_0","4_1","4_2","4_3","4_4","4_5","4_6","4_7","4_8","4_9",
                        "5_0","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9"};


  fhDigitsModule = new TH1F("ACORDEBitPatternfromDigits","Distribution of ACORDE from DIGITS;Modules;Counts",60,1,60);
  Add2DigitsList(fhDigitsModule,0,!expert,image);
  for (Int_t i=0;i<60;i++) fhDigitsModule->GetXaxis()->SetBinLabel(i+1,acoModule[i]); 
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 

void AliACORDEQADataMakerRec::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
  // Not needed for ACORDE by now !!!
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
}

//____________________________________________________________________________
void AliACORDEQADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *    fhESDsSingle;
  TH1F *    fhESDsMulti;
  TH1F *	fhESDsMultiplicity;
  const char *acoModule[60]={"0_0","0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9",
                        "1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9",
                        "2_0","2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9",
                        "3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9",
                        "4_0","4_1","4_2","4_3","4_4","4_5","4_6","4_7","4_8","4_9",
                        "5_0","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9"};


   fhESDsSingle = new TH1F("ACORDEBitPatternfromESDsSingle","Distribution of ACORDE fired modules from ESDs-Single;Modules;Counts",60,1,60);
   Add2ESDsList(fhESDsSingle,0,!expert,image);

   fhESDsMulti = new TH1F("ACORDEBitPatternfromESDsMulti","Distribution of ACORDE fired modules from ESDs-Multi;Modules;Counts",60,1,60);
   Add2ESDsList(fhESDsMulti,1,!expert,image);
   
   fhESDsMultiplicity = new TH1F("ACORDEMultiplicityfromESD","Number of fired ACORDE modules; No. of fired ACORDE modules;No. of events in ACORDE",60,-0.5,60);
   Add2ESDsList(fhESDsMultiplicity,2,!expert,image);	
   for (Int_t i=0;i<60;i++)
   {
	fhESDsSingle->GetXaxis()->SetBinLabel(i+1,acoModule[i]);
	fhESDsMulti->GetXaxis()->SetBinLabel(i+1,acoModule[i]);
   }
   //
   ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line
}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  //fills QA histos for RAW
  IncEvCountCycleRaws();
  IncEvCountTotalRaws();
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
        for(Int_t r=0;r<30;r++)
        {
                kroSingle[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[1];
        for(Int_t r=30;r<60;r++)
        {
                kroSingle[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[2];
        for(Int_t r=0;r<30;r++)
        {
                kroMulti[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        tmpDy=dy[3];
        for(Int_t r=30;r<60;r++)
        {
                kroMulti[r] = tmpDy & 1;
                tmpDy>>=1;
        }
        contSingle=0;
	contMulti=0;
        for(Int_t r=0;r<60;r++)
        {
			if(kroSingle[r]==1)
			{
			  FillRawsData(0,r);
			  //FillRawsData(4,r);
			  contSingle++;
			}
			if(kroMulti[r]==1)
			{
			  FillRawsData(2,r);
			  //FillRawsData(6,r);
			  contMulti++;
			}
			
        } 
	FillRawsData(1,contSingle); 
	FillRawsData(3,contMulti); 
	//	FillRawsData(7,contMulti);
}
}
//____________________________________________________________________________
void AliACORDEQADataMakerRec::MakeDigits( TTree *digitsTree)
{
  //fills QA histos for Digits
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();

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
        FillDigitsData(0,AcoDigit->GetModule()-1);
      }
    }
  }
}

//____________________________________________________________________________
void AliACORDEQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  //
  AliESDACORDE * fESDACORDE= esd->GetACORDEData();
  Int_t acoMulti=0;
  for(int i=0;i<60;i++)
  {
    if(fESDACORDE->GetHitChannel(i)) 
	  {
	  FillESDsData(0,i+1);
	  FillESDsData(1,i+1);
	  acoMulti++;
	}
  } FillESDsData(2,acoMulti);

}
