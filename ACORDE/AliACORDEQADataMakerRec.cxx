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
  fhACOMultiAMU(new TLine(0.,4.,60.,4.)),
  fhACOTriggerCheck(new TLine(0.,4.,60.,4.))
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
  fhACOMultiAMU(qadm.fhACOMultiAMU),
  fhACOTriggerCheck(qadm.fhACOTriggerCheck)
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
  delete fhACOTriggerCheck;
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
      TObjArray &harr = *GetRawsDataOfTrigClass(itc);
      //
      TH1* h0 = (TH1*)harr[0];
      TH1* h1 = (TH1*)harr[1];
      if (!h0 || !h1) continue;
      double integral = 0;
      if (itc==-1 && !(integral=h0->Integral())) { // default clone
	printf("No entries in ACORDE Hits histograms for trigger class %d, fatal error, please check !!!\n",itc);
	TPaveText *acoBoxFatal=new TPaveText(35,0.5,55,1,"b");
	acoBoxFatal->SetFillColor(kRed);
	acoBoxFatal->SetLineColor(kRed);
	acoBoxFatal->SetLineWidth(2);
	//acoBox2->SetTextSize(3);
	//acoBox2->SetTextColor(kBlack);
	acoBoxFatal->AddText("FLAG MESSAGE: ACO. Not Ok, Call the expert !!!");
	acoBoxFatal->AddText("Blue line: mean of hits");
	acoBoxFatal->AddText("Between GREEN lines: ACO. O.K.");
	h0->GetListOfFunctions()->Add(acoBoxFatal);
	
	TPaveText *acoMultiBoxFatal = new TPaveText(20,0.5,40,1,"b");
	acoMultiBoxFatal->SetFillColor(kRed);
	acoMultiBoxFatal->SetLineColor(kRed);
	acoMultiBoxFatal->SetLineWidth(2);
	acoMultiBoxFatal->AddText("ACO. Not O.K., Call the experts");
	h1->GetListOfFunctions()->Add(acoMultiBoxFatal);
	continue;
      }
      Double_t mean = integral/60;
      fhACOMean->SetY1(mean);
      fhACOMean->SetY2(mean);
      fhACOMin->SetY1(0.05*mean);
      fhACOMin->SetY2(0.05*mean);
      fhACOMax->SetY1(2.25*mean);
      fhACOMax->SetY2(2.25*mean);
      
      // alarms
      Double_t max = h0->GetMaximum();
      if (max == 0) {
	printf("Maximum of hits equal to ZERO, please check the status of ACORDE !!\n");
	continue;
      }
      // Multiplicity histogram threshold
      Double_t maxMulti = h1->GetMaximum();
      if (maxMulti == 0) {
	printf("Maximum of entries equal to ZERO, please check the status of ACORDE !!\n");
	continue;
      }
      fhACOMulti->SetX1(1);
      fhACOMulti->SetY1(1);
      fhACOMulti->SetX2(1);
      fhACOMulti->SetY2(maxMulti);
      TPaveText *acoBox=new TPaveText(35,max-0.20*max,55,max,"b");
      //acoBox->SetFillStyle(0);
      TPaveText *acoBox1=new TPaveText(35,max-0.20*max,55,max,"b");
      //acoBox1->SetFillStyle(0);
      TPaveText *acoBox2=new TPaveText(35,max-0.20*max,55,max,"b");
      //acoBox2->SetFillStyle(0);
      Int_t flagACO_DQM = 0;
      Bool_t isACOOk = kTRUE;
      Bool_t isACOWarning = kFALSE;
      for(Int_t iModule=0;iModule<60;iModule++)	{
	if ((h0->GetBinContent(iModule))/max < 0.5) flagACO_DQM++;
      }
      if (flagACO_DQM < 15) {
	acoBox->SetFillColor(kGreen);
	acoBox->SetLineColor(kGreen);
	acoBox->SetLineWidth(2);
	//acoBox->SetTextSize(3);
	//acoBox->SetTextColor(kBlack);
	acoBox->AddText("FLAG MESSAGE: O.K. !!!");
	acoBox->AddText("Blue line: mean of hits");
	acoBox->AddText("Between GREEN lines: ACO. O.K.");
	h0->GetListOfFunctions()->Add(acoBox);	
	//
      } 
      else if (flagACO_DQM > 15 && flagACO_DQM<30) {
	acoBox1->SetFillColor(kYellow);
	acoBox1->SetLineColor(kYellow);
	acoBox1->SetLineWidth(2);
	//acoBox1->SetTextSize(3);
	//acoBox1->SetTextColor(kBlack);
	acoBox1->AddText("FLAG MESSAGE: Warning, some modules are not working properly !!!");
	acoBox1->AddText("Blue line: mean of hits");
	acoBox1->AddText("Between GREEN lines: ACO. O.K.");
	h0->GetListOfFunctions()->Add(acoBox1);
	isACOWarning=kTRUE;	
      }
      else if (flagACO_DQM > 30) {
	acoBox2->SetFillColor(kRed);
	acoBox2->SetLineColor(kRed);
	acoBox2->SetLineWidth(2);
	//acoBox2->SetTextSize(3);
	//acoBox2->SetTextColor(kBlack);
	acoBox2->AddText("FLAG MESSAGE: ACO. Not Ok, Call the expert !!!");
	acoBox2->AddText("Blue line: mean of hits");
	acoBox2->AddText("Between GREEN lines: ACO. O.K.");
	h0->GetListOfFunctions()->Add(acoBox2);
	isACOOk=kFALSE;	
      }
      //
      
      TPaveText *acoMultiBox = new TPaveText(20,maxMulti-0.20*maxMulti,40,maxMulti,"b");
      if (h1->Integral()==0 || isACOOk==kFALSE) {
	acoMultiBox->SetFillColor(kRed);
	acoMultiBox->SetLineColor(kRed);
	acoMultiBox->SetLineWidth(2);
	acoMultiBox->AddText("ACO. Not O.K., Call the experts");
	h1->GetListOfFunctions()->Add(acoMultiBox);
      }
      /*    if (GetRawsData(5)->GetBinContent(1) > 0 || isACOOk && GetRawsData(5)->Integral()!=0 && isACOOk==kTRUE){
	    acoMultiBox->SetFillColor(kYellow);
	    acoMultiBox->SetLineColor(kYellow);
	    acoMultiBox->SetLineWidth(2);
	    acoMultiBox->AddText("Warning: possible empy events only IF ACORDE is triggering, else: O.K.");
	    GetRawsData(5)->GetListOfFunctions()->Add(acoMultiBox);
	    }
      */
      if (isACOOk==kTRUE) {
	acoMultiBox->SetFillColor(kGreen);
	acoMultiBox->SetLineColor(kGreen);
	acoMultiBox->SetLineWidth(2);
	acoMultiBox->AddText("FLAG MESSAGE: ACO. O.K.");
	//acoMultiBox->AddText("NOTE: if entries below the pink line and ACO is triggering, then call the expert (possible empty events)");
	h1->GetListOfFunctions()->Add(acoMultiBox);
      }
      if (isACOWarning==kTRUE) {
	acoMultiBox->SetFillColor(kYellow);
	acoMultiBox->SetLineColor(kYellow);
	acoMultiBox->SetLineWidth(2);
	acoMultiBox->AddText("FLAG MESSAGE: ACO. O.K., warning, some modules are not working properly");
	//acoMultiBox->AddText("NOTE: if entries below the pink line and ACO is triggering, then call the expert (possible empty events)");
	h1->GetListOfFunctions()->Add(acoMultiBox);
      }
      
      // for AMU ACORDE trigger option
      TH1* h2 = (TH1*)harr[2];
      TH1* h3 = (TH1*)harr[3];
      if (!h2 || !h3) continue;
      Double_t integral1 = h2->Integral();
      if (integral1==0) {
	printf("No entries in ACORDE Hits histograms for trigger class %d --> fatal error, please check !!!\n",itc);
	TPaveText *acoBoxFatalAMU=new TPaveText(35,0.5,55,1,"b");
	acoBoxFatalAMU->SetFillColor(kRed);
	acoBoxFatalAMU->SetLineColor(kRed);
	acoBoxFatalAMU->SetLineWidth(2);
	//acoBox2->SetTextSize(3);
	//acoBox2->SetTextColor(kBlack);
	acoBoxFatalAMU->AddText("FLAG MESSAGE: ACO. Not Ok, Call the expert !!!");
	acoBoxFatalAMU->AddText("Blue line: mean of hits");
	acoBoxFatalAMU->AddText("Between GREEN lines: ACO. O.K.");
	h2->GetListOfFunctions()->Add(acoBoxFatalAMU);
	
	TPaveText *acoMultiBoxFatalAMU = new TPaveText(20,0.5,40,1,"b");
	acoMultiBoxFatalAMU->SetFillColor(kRed);
	acoMultiBoxFatalAMU->SetLineColor(kRed);
	acoMultiBoxFatalAMU->SetLineWidth(2);
	acoMultiBoxFatalAMU->AddText("ACO. Not O.K., Call the experts");
	h3->GetListOfFunctions()->Add(acoMultiBoxFatalAMU);
	
	continue;
      }
      Double_t mean1 = integral1/60;
      fhACOMeanAMU->SetY1(mean1);
      fhACOMeanAMU->SetY2(mean1);
      fhACOMinAMU->SetY1(0.05*mean1);
      fhACOMinAMU->SetY2(0.05*mean1);
      fhACOMaxAMU->SetY1(2.25*mean1);
      fhACOMaxAMU->SetY2(2.25*mean1);
      
      // alarms
      Double_t max1 = h2->GetMaximum();
      if (max1 == 0) {
	printf("Maximum of hits equal to ZERO, please check the status of ACORDE !!\n");
	continue;
      }
      // Multiplicity histogram threshold
      Double_t maxMulti1 = h3->GetMaximum();
      if (maxMulti1 == 0) {
	printf("Maximum of entries equal to ZERO, please check the status of ACORDE !!\n");
	continue;
      }
      fhACOMultiAMU->SetX1(1);
      fhACOMultiAMU->SetY1(1);
      fhACOMultiAMU->SetX2(1);
      fhACOMultiAMU->SetY2(maxMulti1);
      TPaveText *acoBoxAMU=new TPaveText(35,max1-0.20*max1,55,max1,"b");
      //acoBox->SetFillStyle(0);
      TPaveText *acoBox1AMU=new TPaveText(35,max1-0.20*max1,55,max1,"b");
      //acoBox1->SetFillStyle(0);
      TPaveText *acoBox2AMU=new TPaveText(35,max1-0.20*max1,55,max1,"b");
      //acoBox2->SetFillStyle(0);
      Int_t flagACO_DQMAMU = 0;
      Bool_t isACOOkAMU = kTRUE;
      Bool_t isACOWarningAMU = kFALSE;
      for(Int_t iModule=0;iModule<60;iModule++) {
	if ((h2->GetBinContent(iModule))/max1 < 0.5) flagACO_DQMAMU++;
      }
      if (flagACO_DQMAMU < 15) {
	acoBoxAMU->SetFillColor(kGreen);
	acoBoxAMU->SetLineColor(kGreen);
	acoBoxAMU->SetLineWidth(2);
	//acoBox->SetTextSize(3);
	//acoBox->SetTextColor(kBlack);
	acoBoxAMU->AddText("FLAG MESSAGE: O.K. !!!");
	acoBoxAMU->AddText("Blue line: mean of hits");
	acoBoxAMU->AddText("Between GREEN lines: ACO. O.K.");
	h2->GetListOfFunctions()->Add(acoBoxAMU);	
	//
      }
      else if (flagACO_DQMAMU > 15 && flagACO_DQMAMU<30) {
	acoBox1AMU->SetFillColor(kYellow);
	acoBox1AMU->SetLineColor(kYellow);
	acoBox1AMU->SetLineWidth(2);
	//acoBox1->SetTextSize(3);
	//acoBox1->SetTextColor(kBlack);
	acoBox1AMU->AddText("FLAG MESSAGE: Warning, some modules are not working properly !!!");
	acoBox1AMU->AddText("Blue line: mean of hits");
	acoBox1AMU->AddText("Between GREEN lines: ACO. O.K.");
	h2->GetListOfFunctions()->Add(acoBox1AMU);
	isACOWarningAMU=kTRUE;
	//
      } 
      else if (flagACO_DQMAMU > 30) {
	acoBox2AMU->SetFillColor(kRed);
	acoBox2AMU->SetLineColor(kRed);
	acoBox2AMU->SetLineWidth(2);
	//acoBox2->SetTextSize(3);
	//acoBox2->SetTextColor(kBlack);
	acoBox2AMU->AddText("FLAG MESSAGE: ACO. Not Ok, Call the expert !!!");
	acoBox2AMU->AddText("Blue line: mean of hits");
	acoBox2AMU->AddText("Between GREEN lines: ACO. O.K.");
	h2->GetListOfFunctions()->Add(acoBox2AMU);
	isACOOkAMU=kFALSE;
      }
      //
      TPaveText *acoMultiBoxAMU = new TPaveText(20,maxMulti1-0.20*maxMulti1,40,maxMulti1,"b");
      if (h3->Integral()==0 || isACOOkAMU==kFALSE) {
	acoMultiBoxAMU->SetFillColor(kRed);
	acoMultiBoxAMU->SetLineColor(kRed);
	acoMultiBoxAMU->SetLineWidth(2);
	acoMultiBoxAMU->AddText("ACO. Not O.K., Call the experts");
	h3->GetListOfFunctions()->Add(acoMultiBoxAMU);
      }
      /*              if (GetRawsData(5)->GetBinContent(1) > 0 || isACOOk && GetRawsData(5)->Integral()!=0 && isACOOk==kTRUE){
		      acoMultiBox->SetFillColor(kYellow);
		      acoMultiBox->SetLineColor(kYellow);
		      acoMultiBox->SetLineWidth(2);
		      acoMultiBox->AddText("Warning: possible empy events only IF ACORDE is triggering, else: O.K.");
		      GetRawsData(5)->GetListOfFunctions()->Add(acoMultiBox);
		      }
      */
      if (isACOOkAMU==kTRUE) {
	acoMultiBoxAMU->SetFillColor(kGreen);
	acoMultiBoxAMU->SetLineColor(kGreen);
	acoMultiBoxAMU->SetLineWidth(2);
	acoMultiBoxAMU->AddText("FLAG MESSAGE: ACO. O.K.");
	//acoMultiBox->AddText("NOTE: if entries below the pink line and ACO is triggering, then call the expert (possible empty events)");
	h3->GetListOfFunctions()->Add(acoMultiBoxAMU);
      }
      if (isACOWarningAMU==kTRUE) {
	acoMultiBoxAMU->SetFillColor(kYellow);
	acoMultiBoxAMU->SetLineColor(kYellow);
	acoMultiBoxAMU->SetLineWidth(2);
	acoMultiBoxAMU->AddText("FLAG MESSAGE: ACO. O.K., warning, some modules are not working properly");
	//acoMultiBox->AddText("NOTE: if entries below the pink line and ACO is triggering, then call the expert (possible empty events)");
	h3->GetListOfFunctions()->Add(acoMultiBoxAMU);
      }
      
      // Checks if hits distribution from SL0 and AMU are equal
      Float_t eff = 0.0;
      Int_t effFlag = 0;
      //
      TH1* h4 = (TH1*)harr[4];
      if (h4) {
	for (Int_t iModule = 0; iModule < 60; iModule++) {
	  if (h2->GetBinContent(iModule)==0) {
	    eff = 0.0;
	    continue;
	  }
	  else {
	    eff = h0->GetBinContent(iModule)/h2->GetBinContent(iModule);
	    h4->Fill(iModule,eff);
	    if (eff!=1) effFlag++;
	  }
	}
	
	if (effFlag == 0)	{
	  TPaveText *checkTriggerBox = new TPaveText(20,0.6,40,0.8,"b");
	  checkTriggerBox->SetFillColor(kGreen);
	  checkTriggerBox->SetLineColor(kGreen);
	  checkTriggerBox->SetLineWidth(2);
	  checkTriggerBox->AddText("FLAG MESSAGE: ACO. trigger O.K.");
	  h4->GetListOfFunctions()->Add(checkTriggerBox);
	}
	else {
	  TPaveText *checkTriggerBox1 = new TPaveText(20,0.6,40,0.8,"b");
	  checkTriggerBox1->SetFillColor(kYellow);
	  checkTriggerBox1->SetLineColor(kYellow);
	  checkTriggerBox1->SetLineWidth(2);
	  checkTriggerBox1->AddText("FLAG MESSAGE: Warning, please check the ACO trigger configuration");
	  h4->GetListOfFunctions()->Add(checkTriggerBox1);
	}
      } // h4
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
  /*
  const char *acoModule[60]={"0_0","0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9",
                        "1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9",
                        "2_0","2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9",
                        "3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9",
                        "4_0","4_1","4_2","4_3","4_4","4_5","4_6","4_7","4_8","4_9",
                        "5_0","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9"};
  */
  // TH1F *fhACORDEBitPattern[4];
         //TH1F *fhACORDEBitPatternDQM;
 //  fhACORDEBitPattern[0] = new TH1F("ACORDEBitPatternfromRAWSingle","Distribution of ACORDE fired modules from RAW-Single;Modules;Counts",60,-0.5,59.5);//AcordeSingleMuon BitPattern
 //  fhACORDEBitPattern[1] = new TH1F("ACORDEBitPatternfromRAWMulti","Distribution of ACORDE fired modules from RAW-Multi;Modules;Counts",60,-0.5,59.5);//AcordeMultiMuon BitPattern
 //  fhACORDEBitPattern[2] = new TH1F("ACORDEMultiplicityfromRAWSingle","Number of fired ACORDE modules;No. of fired ACORDE modules;No. of events in ACORDE",61,-1,60);//AcordeSingleMuon Multiplicity
 //  fhACORDEBitPattern[3] = new TH1F("ACORDEMultiplicityfromRAWMulti","Number of fired ACORDE modules; No. of fired ACORDE modules;No. of events in ACORDE",61,-1,60);//AcordeMultiMuon Multiplicity
         TH1F * fhACORDEBitPatternDQM = new TH1F("ACOHitsSL0_DQM_Shifter","Distribution of ACORDE fired modules for DQM shifter; No. of module; Counts",60,-0.5,59.5);// Hits histogram for QA-shifter ACO-SL0 trigger mode
         TH1F * fhACORDEMultiplicitySL0DQM = new TH1F("ACOMultiSL0_DQM_Shifter","Multiplicity of ACORDE fired modules for DQM shifter; No. of fired modules; No. of events",62,-1,60); // Multiplicity histo. for QA-shifter ACO-SL0 trigger mode
         TH1F * fhACORDEBitPatternAMUDQM = new TH1F("ACOHitsAMU_DQM_Shifter","Distribution of ACORDE fired modules for DQM shifter; No. of module; Counts",60,-0.5,59.5);// Hits histogram for QA-shifter ACO-SL0 trigger mode
         TH1F * fhACORDEMultiplicityAMUDQM = new TH1F("ACOMultiAMU_DQM_Shifter","Multiplicity of ACORDE fired modules for DQM shifter; No. of fired modules; No. of events",62,-1,60); // Multiplicity histo. for QA-shifter ACO-SL0 trigger mode
         TH1F * fhACORDEBitPatternCheckDQM = new TH1F("ACOHitsTriggerCheck_DQMExpert","Check distribution for ACORDE trigger configuration; No. of module; Trigger hits difference",60,-0.5,59.5); // Check the trigger status of ACORDE (SL0 vs AMU)
         // Expert histograms
 //      for(Int_t i=0;i<4;i++)
 //    Add2RawsList(fhACORDEBitPattern[i],i,expert, !image, !saveCorr);
         // Check the hits multiplicity from trigger configuration
         Add2RawsList(fhACORDEBitPatternCheckDQM,4,expert,image,!saveCorr);
         fhACORDEBitPatternCheckDQM->SetFillColor(kOrange);
         // AMORE diplay settings for shifter on GUI
 
         // For SL0 ACO trigger mode
 
         Add2RawsList(fhACORDEBitPatternDQM,0,!expert,image,!saveCorr);
	 ForbidCloning(fhACORDEBitPatternDQM);
         Add2RawsList(fhACORDEMultiplicitySL0DQM,1,!expert,image,!saveCorr);
 	 ForbidCloning(fhACORDEMultiplicitySL0DQM);
         // For Hits distribution on ACORDE
 
         fhACORDEBitPatternDQM->SetFillColor(kCyan-7);
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
 
         fhACORDEMultiplicitySL0DQM->SetFillColor(kBlue+1);
         fhACOMulti->SetLineColor(kMagenta);
         fhACOMulti->SetLineStyle(2);
         fhACOMulti->SetLineWidth(4);
         fhACORDEMultiplicitySL0DQM->GetListOfFunctions()->Add(fhACOMulti);
 
         // For AMU ACO trigger mode
 
         Add2RawsList(fhACORDEBitPatternAMUDQM,2,!expert,image,!saveCorr);
         Add2RawsList(fhACORDEMultiplicityAMUDQM,3,!expert,image,!saveCorr);
 
         // For Hits distribution on ACORDE
 
         fhACORDEBitPatternAMUDQM->SetFillColor(kCyan-7);
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
 
         fhACORDEMultiplicityAMUDQM->SetFillColor(kBlue+1);
         fhACOMultiAMU->SetLineColor(kMagenta);
         fhACOMultiAMU->SetLineStyle(2);
         fhACOMultiAMU->SetLineWidth(4);
         fhACORDEMultiplicityAMUDQM->GetListOfFunctions()->Add(fhACOMultiAMU);
 
	 /*
  for (Int_t iModule = 0; iModule<60; iModule++)
  {
    fhACORDEBitPattern[0]->GetXaxis()->SetBinLabel(iModule+1,acoModule[iModule]);
    fhACORDEBitPattern[1]->GetXaxis()->SetBinLabel(iModule+1,acoModule[iModule]);
  }
	 */
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
	FillRawsData(3,contSingle); 
	FillRawsData(7,contMulti);
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
