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
 AliACORDEQADataMakerRec::AliACORDEQADataMakerRec():AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kACORDE), "ACORDE Quality Assurance Data Maker"),
  fhACOMean(new TLine(0.,4.,60.,4.)),
  fhACOMin(new TLine(0.,4.,60.,4.)),
  fhACOMax(new TLine(0.,4.,60.,4.)),
  fhACOMulti(new TLine(0.,4.,60.,4.)),
  fhACOMeanAMU(new TLine(0.,4.,60.,4.)),
  fhACOMinAMU(new TLine(0.,4.,60.,4.)),
  fhACOMaxAMU(new TLine(0.,4.,60.,4.)),
  fhACOMultiAMU(new TLine(0.,4.,60.,4.))
{

}
//____________________________________________________________________________ 
AliACORDEQADataMakerRec::AliACORDEQADataMakerRec(const AliACORDEQADataMakerRec& qadm):
  AliQADataMakerRec(),
  fhACOMean(qadm.fhACOMean),
  fhACOMin(qadm.fhACOMin),
  fhACOMax(qadm.fhACOMax),
  fhACOMulti(qadm.fhACOMulti),
  fhACOMeanAMU(qadm.fhACOMeanAMU),
  fhACOMinAMU(qadm.fhACOMinAMU),
  fhACOMaxAMU(qadm.fhACOMaxAMU),
  fhACOMultiAMU(qadm.fhACOMultiAMU)
{
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliACORDEQADataMakerRec::~AliACORDEQADataMakerRec()
{
  delete fhACOMean;
  delete fhACOMin;
  delete fhACOMax;
  delete fhACOMulti;
  delete fhACOMeanAMU;
  delete fhACOMinAMU;
  delete fhACOMaxAMU;
  delete fhACOMultiAMU;
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

  // Thresholds and alarms definitions
  /*******************************************************************************************

	We check the performance of the ACORDE's modules respect to the mean
	value over the integral for hits for all the modules (in SL0 and AMU), then
	for the threshold lines we don't use a fixed threshold value.

	For the alarms, we need two fixed values (one for each configuration, SL0 and AMU):

	   -) SL0_ThresholdAlarm: minimum number of ACORDE's modules woriking properly
           -) AMU_ThresholdAlarm: minimum number of ACORDE's modules woriking properly

	This should work for p-p*, HI* and cosmic** data taking.

	* only as readout detector (if ACORDE is triggering, in some rare trigger cluster, 
				    an improvement should be done)
	** as readout and trigger detector.

  *********************************************************************************************/

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
      	TH1* h1 = (TH1*)harr[1];
      	if (!h0 || !h1) continue;
      	TH1* h2 = (TH1*)harr[2];
      	TH1* h3 = (TH1*)harr[3];
      	if (!h2 || !h3) continue;

      	Double_t integralSL0 = 0;
      	Double_t integralAMU = 0;

      	// maximimum and minimum for Pads

      	Double_t maxPad = h0->GetMaximum();
	
      	if ((itc==-1 && !(integralAMU=h2->Integral()) || h2->GetMaximum()==0)||(itc==-1 && !(integralSL0=h0->Integral()) || h0->GetMaximum()==0)) { // default clone
		printf("No entries in ACORDE Hits histograms for trigger class %d, fatal error, please check !!!\n",itc);
		Float_t maxPadAMU = h3->GetMaximum();
		TPaveText *acoBoxFatalAMU = new TPaveText(35,0,55,1,"b");
		acoBoxFatalAMU->SetFillColor(kRed);
		acoBoxFatalAMU->SetLineColor(kRed);
		acoBoxFatalAMU->SetLineWidth(2);
		acoBoxFatalAMU->AddText("ACO: Not O.K.");
		acoBoxFatalAMU->AddText("Call the experts");

		TPaveText *acoBoxFatalMultiAMU = new TPaveText(35,0,55,maxPadAMU,"b");
		acoBoxFatalMultiAMU->SetFillColor(kRed);
		acoBoxFatalMultiAMU->SetLineColor(kRed);
		acoBoxFatalMultiAMU->SetLineWidth(2);
		acoBoxFatalMultiAMU->AddText("ACO: Not O.K.");
		acoBoxFatalMultiAMU->AddText("Call the experts");

		h2->GetListOfFunctions()->Add(acoBoxFatalAMU);
		h3->GetListOfFunctions()->Add(acoBoxFatalMultiAMU);

		Float_t maxPadSL0 = h1->GetMaximum();

        	printf("No entries in ACORDE Hits histograms for trigger class %d, fatal error, please check !!!\n",itc);
        	TPaveText *acoBoxFatal = new TPaveText(35,0,55,1,"b");
        	acoBoxFatal->SetFillColor(kRed);
        	acoBoxFatal->SetLineColor(kRed);
        	acoBoxFatal->SetLineWidth(2);
        	acoBoxFatal->AddText("ACO: Not O.K.");
		acoBoxFatal->AddText("Call the experts");

        	TPaveText *acoBoxFatalMulti = new TPaveText(35,0,55,maxPadSL0,"b");
        	acoBoxFatalMulti->SetFillColor(kRed);
        	acoBoxFatalMulti->SetLineColor(kRed);
        	acoBoxFatalMulti->SetLineWidth(2);
        	acoBoxFatalMulti->AddText("ACO: Not O.K.");
		acoBoxFatalMulti->AddText("Call the experts");

        	h0->GetListOfFunctions()->Add(acoBoxFatal);
        	h1->GetListOfFunctions()->Add(acoBoxFatalMulti);
		continue;
	}
      
	// Check DQM

	// setting the thresholds

	Int_t SL0_ThresholdAlarm = 45; // default value until the creation of AliACORDEThreshold class
	Int_t AMU_ThresholdAlarm = 45; // default value until the creation of AliACORDEThreshold class

	// SL0 - histograms

	Int_t indexActiveModuleSL0 = 0;
        for(Int_t iModule=0;iModule<60;iModule++){	
		if (h0->GetBinContent(iModule)>0) indexActiveModuleSL0++;
	}
	
	Float_t meanHitsSL0 = 0.;
	if (indexActiveModuleSL0!=0) meanHitsSL0=h0->Integral()/indexActiveModuleSL0;

	Int_t indexSL0 = 0;
	
	// set the threshold lines: minimum, maximum and mean

	fhACOMean->SetX1(0);
	fhACOMean->SetY1(meanHitsSL0);
	fhACOMean->SetX2(59);
	fhACOMean->SetY2(meanHitsSL0);
	
	fhACOMin->SetX1(0);
	fhACOMin->SetX2(59);
	fhACOMin->SetY1(meanHitsSL0-0.80*meanHitsSL0);
	fhACOMin->SetY2(meanHitsSL0-0.80*meanHitsSL0);

	fhACOMax->SetX1(0);
	fhACOMax->SetX2(59);
	fhACOMax->SetY1(meanHitsSL0+0.80*meanHitsSL0);
	fhACOMax->SetY2(meanHitsSL0+0.80*meanHitsSL0);

	fhACOMulti->SetX1(0);
	fhACOMulti->SetY1(0);
	fhACOMulti->SetX2(0);
	Float_t maxMulti = 0;
	if (h1->GetMaximum()>0) maxMulti = h1->GetMaximum();
	fhACOMulti->SetY2(maxMulti);

	TPaveText *acoBoxOkHitsSL0 = new TPaveText(35,meanHitsSL0+0.5*meanHitsSL0,55,maxPad,"b");
	acoBoxOkHitsSL0->SetFillColor(kGreen);
	acoBoxOkHitsSL0->SetLineColor(kGreen);
	acoBoxOkHitsSL0->SetLineWidth(2);
	acoBoxOkHitsSL0->AddText("ACO: O.K.");

	TPaveText *acoBoxErrorHitsSL0 = new TPaveText(35,meanHitsSL0+0.5*meanHitsSL0,55,maxPad,"b");
	acoBoxErrorHitsSL0->SetFillColor(kRed);
	acoBoxErrorHitsSL0->SetLineColor(kRed);
	acoBoxErrorHitsSL0->SetLineWidth(2);
	acoBoxErrorHitsSL0->AddText("ACO: Not O.K.");

	Float_t maxPadMulti = h1->GetMaximum();

	TPaveText *acoBoxOkMultiSL0 = new TPaveText(35,maxPadMulti-0.3*maxPadMulti,55,maxPadMulti,"b");
	acoBoxOkMultiSL0->SetFillColor(kGreen);
	acoBoxOkMultiSL0->SetLineColor(kGreen);
	acoBoxOkMultiSL0->SetLineWidth(2);
	acoBoxOkMultiSL0->AddText("ACO: O.K.");

	TPaveText *acoBoxErrorMultiSL0 = new TPaveText(35,maxPadMulti-0.3*maxPadMulti,55,maxPadMulti,"b");
	acoBoxErrorMultiSL0->SetFillColor(kRed);
	acoBoxErrorMultiSL0->SetLineColor(kRed);
	acoBoxErrorMultiSL0->SetLineWidth(2);
	acoBoxErrorMultiSL0->AddText("ACO: Not O.K.");
        TH1* h4 = (TH1*)harr[4];

	for (Int_t iModule = 0; iModule < 60; iModule++){
		if (meanHitsSL0!=0){
			if (TMath::Abs(h0->GetBinContent(iModule)/meanHitsSL0-1) < 1) indexSL0++;
			if (h4){
				h4->Fill(h0->GetBinContent(iModule)/meanHitsSL0-1);
			}
		}
	}

	if (indexSL0>=SL0_ThresholdAlarm){
		h0->GetListOfFunctions()->Add(acoBoxOkHitsSL0);
		h1->GetListOfFunctions()->Add(acoBoxOkMultiSL0);
	}
	else{
		h0->GetListOfFunctions()->Add(acoBoxErrorHitsSL0);
		h1->GetListOfFunctions()->Add(acoBoxErrorMultiSL0);
	}


	// AMU - histograms	

	Int_t indexActiveModuleAMU = 0;
        for(Int_t iModule=0;iModule<60;iModule++){
		if (h2->GetBinContent(iModule)>0) indexActiveModuleAMU++;
	}

	Float_t meanHitsAMU = 0.;
	if (indexActiveModuleAMU!=0) meanHitsAMU=h2->Integral()/indexActiveModuleAMU;

	// setting the line's thresholds: min, max and mean of hits

	fhACOMeanAMU->SetX1(0);
	fhACOMeanAMU->SetY1(meanHitsAMU);
	fhACOMeanAMU->SetX2(59);
	fhACOMeanAMU->SetY2(meanHitsAMU);

	fhACOMinAMU->SetX1(0);
	fhACOMinAMU->SetX2(59);
	fhACOMinAMU->SetY1(meanHitsAMU-0.80*meanHitsAMU);
	fhACOMinAMU->SetY2(meanHitsAMU-0.80*meanHitsAMU);


	fhACOMaxAMU->SetX1(0);
	fhACOMaxAMU->SetX2(59);
	fhACOMaxAMU->SetY1(meanHitsAMU+0.80*meanHitsAMU);
	fhACOMaxAMU->SetY2(meanHitsAMU+0.80*meanHitsAMU);


	fhACOMultiAMU->SetX1(0);
	fhACOMultiAMU->SetY1(0);
	fhACOMultiAMU->SetX2(0);
	Float_t maxMultiAMU = 0;
	if (h3->GetMaximum()>0) maxMultiAMU = h3->GetMaximum();
	fhACOMultiAMU->SetY2(maxMultiAMU);
	
	// setting the alarms

	TPaveText *acoBoxOkHitsAMU = new TPaveText(35,meanHitsAMU+0.5*meanHitsAMU,55,maxPad,"b");
	acoBoxOkHitsAMU->SetFillColor(kGreen);
	acoBoxOkHitsAMU->SetLineColor(kGreen);
	acoBoxOkHitsAMU->SetLineWidth(2);
	acoBoxOkHitsAMU->AddText("ACO: O.K.");

	TPaveText *acoBoxErrorHitsAMU = new TPaveText(35,meanHitsAMU+0.5*meanHitsAMU,55,maxPad,"b");
	acoBoxErrorHitsAMU->SetFillColor(kRed);
	acoBoxErrorHitsAMU->SetLineColor(kRed);
	acoBoxErrorHitsAMU->SetLineWidth(2);
	acoBoxErrorHitsAMU->AddText("ACO: Not O.K.");

	Float_t maxPadMultiAMU = h3->GetMaximum();

	TPaveText *acoBoxOkMultiAMU = new TPaveText(35,maxPadMultiAMU-0.3*maxPadMultiAMU,55,maxPadMultiAMU,"b");
	acoBoxOkMultiAMU->SetFillColor(kGreen);
	acoBoxOkMultiAMU->SetLineColor(kGreen);
	acoBoxOkMultiAMU->SetLineWidth(2);
	acoBoxOkMultiAMU->AddText("ACO: O.K.");

	TPaveText *acoBoxErrorMultiAMU = new TPaveText(35,maxPadMultiAMU-0.3*maxPadMultiAMU,55,maxPadMultiAMU,"b");
	acoBoxErrorMultiAMU->SetFillColor(kRed);
	acoBoxErrorMultiAMU->SetLineColor(kRed);
	acoBoxErrorMultiAMU->SetLineWidth(2);
	acoBoxErrorMultiAMU->AddText("ACO: Not O.K.");
        TH1* h5 = (TH1*)harr[5];

	Int_t indexAMU=0;

	for (Int_t iModule = 0; iModule < 60; iModule++){
		if (meanHitsAMU!=0){
			if (TMath::Abs(h2->GetBinContent(iModule)/meanHitsAMU-1) < 1) indexAMU++;
			if (h5){
				h5->Fill(h2->GetBinContent(iModule)/meanHitsAMU-1);
			}
		}
	}

	if (indexAMU>=AMU_ThresholdAlarm) {
		h2->GetListOfFunctions()->Add(acoBoxOkHitsAMU);
		h3->GetListOfFunctions()->Add(acoBoxOkMultiAMU);
	}
	else{
		h2->GetListOfFunctions()->Add(acoBoxErrorHitsAMU);
		h3->GetListOfFunctions()->Add(acoBoxErrorMultiAMU);
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
        TH1F * fhACORDEBitPatternCheckDQMSL0 = new TH1F("ACOHitsTriggerCheck_DQMExpertSL0","Check the activity of ACORDE's modules; |Hits per module/mean of Hits - 1|; Counts",100,-3,5); // Check the trigger status of ACORDE (SL0 vs AMU)
        TH1F * fhACORDEBitPatternCheckDQMAMU = new TH1F("ACOHitsTriggerCheck_DQMExpertAMU","Check the activity of ACORDE's modules; |Hits per module/mean of Hits - 1|; Counts",100,-3,5); // Check the trigger status of ACORDE (SL0 vs AMU)
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
         fhACOMean->SetLineColor(kBlue);
         fhACOMean->SetLineStyle(2);
         fhACOMean->SetLineWidth(4);
         fhACORDEBitPatternDQM->GetListOfFunctions()->Add(fhACOMean);
         fhACOMin->SetLineColor(kGreen);
         fhACOMin->SetLineStyle(2);
         fhACOMin->SetLineWidth(4);
         fhACORDEBitPatternDQM->GetListOfFunctions()->Add(fhACOMin);
         fhACOMax->SetLineColor(kGreen);
         fhACOMax->SetLineStyle(2);
         fhACOMax->SetLineWidth(4);
         fhACORDEBitPatternDQM->GetListOfFunctions()->Add(fhACOMax);
 
         // For ACORDE Multiplicity distribution of fired modules
 
         fhACORDEMultiplicitySL0DQM->SetFillColor(kMagenta);
         fhACOMulti->SetLineColor(kMagenta);
         fhACOMulti->SetLineStyle(2);
         fhACOMulti->SetLineWidth(4);
         fhACORDEMultiplicitySL0DQM->GetListOfFunctions()->Add(fhACOMulti);
 
         // For AMU ACO trigger mode
 
         Add2RawsList(fhACORDEBitPatternAMUDQM,2,!expert,image,!saveCorr);
         Add2RawsList(fhACORDEMultiplicityAMUDQM,3,!expert,image,!saveCorr);
 
         // For Hits distribution on ACORDE
 
         fhACORDEBitPatternAMUDQM->SetFillColor(kViolet+7);
         fhACOMeanAMU->SetLineColor(kBlue);
         fhACOMeanAMU->SetLineStyle(2);
         fhACOMeanAMU->SetLineWidth(4);
         fhACORDEBitPatternAMUDQM->GetListOfFunctions()->Add(fhACOMeanAMU);
         fhACOMinAMU->SetLineColor(kGreen);
         fhACOMinAMU->SetLineStyle(2);
         fhACOMinAMU->SetLineWidth(4);
         fhACORDEBitPatternAMUDQM->GetListOfFunctions()->Add(fhACOMinAMU);
         fhACOMaxAMU->SetLineColor(kGreen);
         fhACOMaxAMU->SetLineStyle(2);
         fhACOMaxAMU->SetLineWidth(4);
         fhACORDEBitPatternAMUDQM->GetListOfFunctions()->Add(fhACOMaxAMU);
 
         // For ACORDE Multiplicity distribution of fired modules
 
         fhACORDEMultiplicityAMUDQM->SetFillColor(kViolet+6);
         fhACOMultiAMU->SetLineColor(kAzure-2);
         fhACOMultiAMU->SetLineStyle(2);
         fhACOMultiAMU->SetLineWidth(4);
         fhACORDEMultiplicityAMUDQM->GetListOfFunctions()->Add(fhACOMultiAMU);
 
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
